// Выполнение расчёта по конфигурации: общий код CLI и GUI.
#include "runner.hpp"

#include "output.hpp"

#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/heom.hpp"
#include "ssp/excitons/plexciton.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/excitons/spectra.hpp"
#include "ssp/excitons/vibronic.hpp"
#include "ssp/units.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace ssp;
using namespace ssp::excitons;

namespace ssp::app {

namespace {

// Поток журнала: отдаёт в Log завершённые строки.
class LogStream {
 public:
  explicit LogStream(const Log& log) : log_(log) {}
  ~LogStream() {
    if (!buf_.empty()) log_(buf_);
  }
  template <class T>
  LogStream& operator<<(const T& v) {
    std::ostringstream os;
    os << v;
    buf_ += os.str();
    std::size_t p;
    while ((p = buf_.find('\n')) != std::string::npos) {
      log_(buf_.substr(0, p));
      buf_.erase(0, p + 1);
    }
    return *this;
  }

 private:
  const Log& log_;
  std::string buf_;
};


std::string fmt(double v, int prec = 6) {
  std::ostringstream os;
  os.precision(prec);
  os << v;
  return os.str();
}

int pick_state(const ExcitonStates& st, const std::string& how) {
  if (how == "highest") return static_cast<int>(st.size()) - 1;
  Eigen::Index ib = 0;
  st.mu2.maxCoeff(&ib);
  return static_cast<int>(ib);
}

void write_geometry(const std::string& dir, const Aggregate& a) {
  app::CsvWriter w(app::join_path(dir, "geometry.csv"), {"molecule positions (nm) and unit transition dipoles"},
                   {"x_nm", "y_nm", "z_nm", "mu_x", "mu_y", "mu_z"});
  for (std::size_t i = 0; i < a.size(); ++i)
    w.row({a.positions[i].x(), a.positions[i].y(), a.positions[i].z(), a.dipoles[i].x(), a.dipoles[i].y(),
           a.dipoles[i].z()});
}

}  // namespace

RunResult run(const RunConfig& rc, const std::string& config_json, const Log& log, const std::atomic<bool>* cancel,
              const Progress& progress) {
  LogStream out(log);
  RunResult result;
  const auto t_start = std::chrono::steady_clock::now();
  app::make_dirs(rc.output_dir);
  {
    std::ofstream f(app::join_path(rc.output_dir, "config_used.json"));
    f << config_json;
  }

  std::unique_ptr<plasmonics::Material> material;
  if (rc.nanoparticle.enabled) {
    const auto& np = rc.nanoparticle;
    material = std::make_unique<plasmonics::Material>(
        np.use_drude ? plasmonics::Material::drude(np.drude)
                     : plasmonics::Material::from_refractiveindex_yaml(np.material_file));
    const double r_eff = std::cbrt(np.spec.a_xy_nm * np.spec.a_xy_nm * np.spec.a_z_nm);
    material->set_size_damping(np.kreibig_A, np.fermi_velocity_m_s, r_eff, np.kreibig_omega_p_eV,
                               np.kreibig_gamma_bulk_eV);
    out << "NP material: " << material->name() << "\n";
  }

  // Сводная таблица по формам
  std::vector<std::string> obs_names;
  for (const auto& [k, d] : observable_descriptions()) obs_names.push_back(k);
  std::vector<std::string> cols = {"name", "geometry_type", "N", "realizations", "converged"};
  for (const auto& k : obs_names) {
    cols.push_back(k);
    cols.push_back(k + "_se");
  }
  if (rc.nanoparticle.enabled) {
    cols.push_back("np_max_field_enhancement");
    cols.push_back("np_min_gap_nm");
    cols.push_back("np_multipole_truncation");
    cols.push_back("qy_hybrid");
    cols.push_back("qy_free");
  }
  app::CsvWriter summary(app::join_path(rc.output_dir, "summary.csv"),
                         {"aggregate shape summary: mean and standard error over the disorder ensemble",
                          "observable descriptions — observables.csv in each shape folder"},
                         cols);

  for (const auto& cfg : rc.configurations) {
    if (cancel && cancel->load()) {
      out << "\nStopped by user.\n";
      break;
    }
    if (progress) progress(static_cast<int>(&cfg - rc.configurations.data()), static_cast<int>(rc.configurations.size()), cfg.name);
    out << "\n=== " << cfg.name << " (" << cfg.type << ") ===\n";
    Aggregate agg = build_geometry(cfg.type, cfg.params);
    if (cfg.placement.center) transform(agg, -centroid(agg), 0, 0, 0);
    transform(agg, cfg.placement.shift_nm, cfg.placement.euler_deg[0], cfg.placement.euler_deg[1],
              cfg.placement.euler_deg[2]);
    const double dmin = min_pair_distance(agg);
    out << "molecules: " << agg.size() << ", min. distance " << fmt(dmin, 4) << " nm\n";
    if (rc.model.coupling.model == CouplingModel::PointDipole && dmin < 1.5)
      out << "  ! warning: the point-dipole approximation is unreliable for r < ~1.5 nm; "
                   "consider extended_dipole\n";

    const std::string dir = app::join_path(rc.output_dir, cfg.name);
    app::make_dirs(dir);
    write_geometry(dir, agg);

    if (rc.model.dispersion.enabled) {
      // Зависящий от окружения дисперсионный сдвиг энергий сайтов: сводка и файл по молекулам
      const auto D = dispersion_site_shifts(agg, rc.model.dispersion);
      const Eigen::VectorXd& s = D.shift_cm;
      Eigen::VectorXd srt = s;
      std::sort(srt.data(), srt.data() + srt.size());
      out << "site dispersion shift: continuum " << fmt(D.continuum_cm, 5) << " cm⁻¹ (n = "
          << fmt(rc.model.dispersion.n_aggregate, 3) << "), back to gas phase +" << fmt(D.to_gas_cm, 5)
          << " cm⁻¹; site shifts: median " << fmt(srt(srt.size() / 2), 5) << ", 10–90 % "
          << fmt(srt(srt.size() / 10), 5) << " … " << fmt(srt((9 * srt.size()) / 10), 5) << " cm⁻¹\n";
      app::CsvWriter w(app::join_path(dir, "site_shifts.csv"),
                       {"site energy shifts relative to the monomer energy in its solvent, cm⁻¹",
                        "dispersive = Q·Σ_k V_mk² (Q = " + fmt(D.Q_cm_nm2, 5) + " cm⁻¹·nm², calibrated to the continuum " +
                            fmt(D.continuum_cm, 5) + " cm⁻¹); total = dispersive + " + fmt(D.to_gas_cm, 5) + " + extra"},
                       {"x_nm", "y_nm", "z_nm", "r_xy_nm", "dispersive_cm", "total_shift_cm"});
      for (std::size_t i = 0; i < agg.size(); ++i)
        w.row({agg.positions[i].x(), agg.positions[i].y(), agg.positions[i].z(),
               std::hypot(agg.positions[i].x(), agg.positions[i].y()), D.dispersive_cm(static_cast<Eigen::Index>(i)),
               s(static_cast<Eigen::Index>(i))});
    }
    const int n_store = std::max(rc.nanoparticle.enabled ? rc.nanoparticle.realizations : 0,
                                 rc.holstein.enabled ? rc.holstein.realizations : 0);
    EnsembleSettings es = rc.ensemble;
    if (rc.nanoparticle.enabled) {
      es.extra_wmin_cm = rc.nanoparticle.wmin_cm;
      es.extra_wmax_cm = rc.nanoparticle.wmax_cm;
    }
    if (rc.dynamics.enabled) {
      es.dyn_initial = rc.dynamics.initial;
      for (int k = 0; k < rc.dynamics.n_times; ++k)
        es.dyn_times_ps.push_back(rc.dynamics.t_max_ps * k / std::max(1, rc.dynamics.n_times - 1));
    }
    const EnsembleResult er = run_ensemble(agg, rc.model, rc.grid, es, n_store);
    ShapeResult& sr = result.shapes.emplace_back();
    sr.name = cfg.name;
    sr.n_molecules = static_cast<int>(agg.size());
    sr.w = er.w;
    sr.abs = er.abs_mean;
    sr.em = er.em_mean;
    for (const auto& [k, st] : er.observables) sr.observables[k] = st.mean;
    out << "realizations: " << er.realizations << (er.converged ? " (converged)" : " (NOT converged within max)")
              << ", line-shape window residual: " << fmt(er.max_lineshape_tail, 3) << "\n";
    if (er.max_lineshape_tail > 1e-3)
      out << "  ! increase spectrum.t_max_ps: the line shape has not decayed within the window\n";
    {
      // Доля силы осциллятора вне спектрального окна (по представительной реализации)
      const auto& st = er.representative.states;
      double outside = 0.0;
      for (Eigen::Index n = 0; n < st.energies.size(); ++n)
        if (st.energies(n) < rc.grid.wmin_cm || st.energies(n) > rc.grid.wmax_cm) outside += st.mu2(n);
      outside /= st.mu2.sum();
      out << "exciton band: " << fmt(st.energies.minCoeff(), 6) << " … " << fmt(st.energies.maxCoeff(), 6)
                << " cm⁻¹; fraction of oscillator strength outside the spectral window: " << fmt(100.0 * outside, 3) << " %\n";
      if (outside > 0.05)
        out << "  ! more than 5 % of the oscillator strength lies outside [wmin, wmax]: peak positions and widths are unreliable — "
                     "widen the spectrum window\n";
    }

    {
      app::CsvWriter w(app::join_path(dir, "spectra.csv"),
                       {"configuration: " + cfg.name,
                        "abs: absorption cross section per aggregate, nm² (cumulant + Redfield); em: photons/(ps·cm⁻¹) per "
                        "unit pump (integral ≈ quantum yield)",
                        "_se — standard error of the ensemble mean, realizations: " +
                            std::to_string(er.realizations)},
                       {"w_cm", "abs_nm2", "abs_se", "em", "em_se"});
      for (std::size_t q = 0; q < er.w.size(); ++q)
        w.row({er.w[q], er.abs_mean[q], er.abs_se[q], er.em_mean[q], er.em_se[q]});
    }
    if (rc.holstein.enabled && er.engine && !er.stored.empty()) {
      // Вибронная модель Хольстейна (OPA) на тех же реализациях беспорядка
      const std::size_t nr = std::min<std::size_t>(er.stored.size(), static_cast<std::size_t>(rc.holstein.realizations));
      const std::vector<Realization> sub(er.stored.begin(), er.stored.begin() + static_cast<std::ptrdiff_t>(nr));
      std::vector<double> vm, vs;
      holstein_absorption(sub, rc.holstein.params, *er.engine, rc.model.n_medium, er.w, vm, vs);
      sr.vib_abs = vm;
      app::CsvWriter w(app::join_path(dir, "vibronic_absorption.csv"),
                       {"Holstein one-particle approximation: mode " + fmt(rc.holstein.params.omega_cm) + " cm⁻¹, S = " +
                            fmt(rc.holstein.params.huang_rhys, 3) + ", ν ≤ " + std::to_string(rc.holstein.params.nu_max),
                        "abs_nm2 — same units as spectra.csv; realizations " + std::to_string(nr)},
                       {"w_cm", "abs_nm2", "abs_se", "abs_electronic_nm2"});
      for (std::size_t q = 0; q < er.w.size(); ++q) w.row({er.w[q], vm[q], vs[q], er.abs_mean[q]});
      out << "vibronic (Holstein OPA) absorption: vibronic_absorption.csv, " << nr << " realizations, basis "
          << agg.size() * static_cast<std::size_t>(rc.holstein.params.nu_max + 1) << "\n";
    }
    {
      app::CsvWriter w(app::join_path(dir, "observables.csv"), {"scalar observables: mean, SE, SD"},
                       {"name", "mean", "se", "sd", "description"});
      for (const auto& [k, s] : er.observables)
        w.row_text({k, fmt(s.mean, 8), fmt(s.se(), 4), fmt(s.sd(), 4), observable_descriptions().at(k)});
    }
    const Realization& R = er.representative;
    {
      app::CsvWriter w(app::join_path(dir, "states_representative.csv"),
                       {"exciton states of one disorder realization (seed)"},
                       {"energy_cm", "mu2_over_mu2mono", "delocalization", "gamma_rad_ps", "hwhm_cm", "pump",
                        "population"});
      const double m2 = rc.model.coupling.mu_D * rc.model.coupling.mu_D;
      for (Eigen::Index n = 0; n < R.states.energies.size(); ++n)
        w.row({R.states.energies(n), R.states.mu2(n) / m2, 1.0 / R.states.ipr(n), R.gamma_rad_ps(n), R.hwhm_cm(n),
               R.pump(n), R.P(n)});
    }

    // Кинетика населённостей (секулярный Редфилд + распад)
    if (rc.dynamics.enabled && !er.dyn.empty()) {
      app::CsvWriter w(app::join_path(dir, "dynamics.csv"),
                       {"exciton population kinetics averaged over " + std::to_string(er.realizations) +
                            " disorder realizations (Redfield + decay), initial: " + rc.dynamics.initial,
                        "_se — standard error over realizations"},
                       {"t_ps", "total", "total_se", "lowest", "lowest_se", "brightest", "brightest_se",
                        "mean_energy_cm", "mean_energy_se"});
      for (std::size_t k = 0; k < er.dyn_t.size(); ++k) {
        std::vector<double> row = {er.dyn_t[k]};
        for (int q = 0; q < 4; ++q) {
          row.push_back(er.dyn[q][k].mean);
          row.push_back(er.dyn[q][k].se());
        }
        w.row(row);
      }
    }

    // HEOM-эталон
    if (rc.heom.enabled) {
      const int N = static_cast<int>(agg.size());
      if (N > rc.heom.max_sites) {
        out << "HEOM skipped: N = " << N << " > max_sites = " << rc.heom.max_sites << "\n";
      } else if (rc.model.bath.type != BathType::DrudeLorentz || !rc.model.bath.modes.empty()) {
        out << "HEOM skipped: only a Drude–Lorentz bath without modes is supported\n";
      } else {
        out << "HEOM: ADO = " << heom_ado_count(N, rc.heom.params.matsubara, rc.heom.params.depth) << "\n";
        if (rc.heom.absorption) {
          Eigen::MatrixX3d smu(N, 3);
          for (int i = 0; i < N; ++i) smu.row(i) = rc.model.coupling.mu_D * agg.dipoles[i].transpose();
          const auto ha = heom_absorption(R.H, smu, rc.heom.params, er.w);
          // Кумулянт для той же реализации (нормированная форма)
          const LineshapeEngine eng(rc.model.bath, 1.2 * (rc.grid.wmax_cm - rc.grid.wmin_cm) + 3000.0,
                                    rc.grid.t_max_ps / units::ps_per_tu);
          std::vector<double> cum(er.w.size(), 0.0);
          const Eigen::VectorXd Gout = -R.K_ps.diagonal() / units::rad_ps_per_cm;
          for (Eigen::Index k = 0; k < R.states.energies.size(); ++k)
            eng.add_line(cum, er.w, R.states.energies(k), R.states.mu2(k) / R.states.mu2.sum(), R.states.ipr(k),
                         Gout(k), false, 0);
          app::CsvWriter w(app::join_path(dir, "heom_absorption.csv"),
                           {"normalized absorption line shape of one realization: HEOM (reference) vs cumulant+Redfield",
                            "HEOM: K=" + std::to_string(rc.heom.params.matsubara) +
                                ", L=" + std::to_string(rc.heom.params.depth) + ", dt=" + fmt(ha.info.dt_fs) + " fs"},
                           {"w_cm", "heom", "cumulant"});
          for (std::size_t q = 0; q < er.w.size(); ++q) w.row({er.w[q], ha.lineshape[q], cum[q]});
        }
        if (rc.heom.populations) {
          const Eigen::Index n = R.states.energies.size();
          Eigen::MatrixXcd rho0 = Eigen::MatrixXcd::Zero(n, n);
          Eigen::VectorXd P0 = Eigen::VectorXd::Zero(n);
          if (rc.heom.initial.rfind("site:", 0) == 0) {
            const int s = std::stoi(rc.heom.initial.substr(5));
            rho0(s, s) = 1.0;
            for (Eigen::Index k = 0; k < n; ++k) P0(k) = R.states.c(s, k) * R.states.c(s, k);
          } else {
            const int k0 = pick_state(R.states, rc.heom.initial);
            rho0 = (R.states.c.col(k0) * R.states.c.col(k0).transpose()).cast<std::complex<double>>();
            P0(k0) = 1.0;
          }
          const auto hp = heom_populations(R.H, rho0, rc.heom.params, rc.heom.n_out);
          Eigen::VectorXd t(static_cast<Eigen::Index>(hp.t_ps.size()));
          for (std::size_t k = 0; k < hp.t_ps.size(); ++k) t(static_cast<Eigen::Index>(k)) = hp.t_ps[k];
          const Eigen::MatrixXd Pr =
              population_dynamics(R.K_ps, Eigen::VectorXd::Zero(n), R.db_energies, rc.model.bath.beta(), P0, t);
          std::vector<std::string> c = {"t_ps"};
          for (Eigen::Index k = 0; k < n; ++k) {
            c.push_back("heom_" + std::to_string(k));
            c.push_back("redfield_" + std::to_string(k));
          }
          app::CsvWriter w(app::join_path(dir, "heom_populations.csv"),
                           {"exciton state populations (ascending energy): HEOM vs secular Redfield",
                            "initial state: " + rc.heom.initial}, c);
          for (Eigen::Index k = 0; k < t.size(); ++k) {
            std::vector<double> row = {t(k)};
            for (Eigen::Index s = 0; s < n; ++s) {
              row.push_back(hp.exciton(s, k));
              row.push_back(Pr(s, k));
            }
            w.row(row);
          }
        }
      }
    }

    // Агрегат + наночастица
    std::vector<std::string> np_cols;
    if (rc.nanoparticle.enabled) {
      const auto& np = rc.nanoparticle;
      PlexcitonSystem sys(agg, np.spec, *material, rc.model.coupling.mu_D);
      std::vector<double> w(np.nw);
      {
        const double step = (np.wmax_cm - np.wmin_cm) / double(np.nw - 1);
        const double fw = er.observables.count("abs_fwhm_cm") ? er.observables.at("abs_fwhm_cm").mean : 0.0;
        if (fw > 0 && step > 0.2 * fw)
          out << "  ! NP frequency step " << fmt(step, 3) << " cm⁻¹ > 0.2·FWHM of the aggregate band (" << fmt(fw, 3)
              << " cm⁻¹): increase nanoparticle.nw\n";
      }
      for (std::size_t q = 0; q < np.nw; ++q) w[q] = np.wmin_cm + (np.wmax_cm - np.wmin_cm) * q / double(np.nw - 1);
      PlexcitonSpectra mean, se;
      if (sys.fast_l() > 0)
        out << "NP: fast solver — multipoles l ≤ " << sys.fast_l()
            << " exact (Woodbury), higher l diagonal in the exciton basis\n";
      // Отклик экситонов — та же кумулянтная теория, что и спектр агрегата без НЧ (единая форма линии)
      if (np.spec.size_sigma_rel <= 0.0) {
        sys.spectra(er.stored, er.engine.get(), np.homogeneous_hwhm_cm, w, mean, se);
      } else {
        // Распределение размеров ядра: узлы Гаусса–Эрмита; молекулы следуют за поверхностью (зазор сохраняется)
        const auto nodes = normal_quadrature(np.spec.size_nodes);
        out << "NP size distribution: σ_R/R = " << fmt(np.spec.size_sigma_rel, 3) << ", " << nodes.size()
            << " nodes, radius ";
        for (auto* v : {&mean, &se}) {
          v->w = w;
          for (auto* x : {&v->ext, &v->abs, &v->sca, &v->np_ext, &v->np_abs, &v->np_sca}) x->assign(w.size(), 0.0);
        }
        for (const auto& [x, wt] : nodes) {
          const double sc = std::max(0.3, 1.0 + np.spec.size_sigma_rel * x);
          NanoparticleSpec sp = np.spec;
          sp.a_xy_nm *= sc;
          sp.a_z_nm *= sc;
          plasmonics::Material mk = *material;
          mk.set_size_damping(np.kreibig_A, np.fermi_velocity_m_s, std::cbrt(sp.a_xy_nm * sp.a_xy_nm * sp.a_z_nm),
                              np.kreibig_omega_p_eV, np.kreibig_gamma_bulk_eV);
          const PlexcitonSystem sk(follow_core_scale(agg, np.spec, sc), sp, mk, rc.model.coupling.mu_D);
          PlexcitonSpectra mk_, sk_;
          sk.spectra(er.stored, er.engine.get(), np.homogeneous_hwhm_cm, w, mk_, sk_);
          out << fmt(sp.a_xy_nm, 3) << " (w " << fmt(wt, 2) << ") ";
          auto acc = [&](std::vector<double>& m, std::vector<double>& e, const std::vector<double>& mv,
                         const std::vector<double>& ev) {
            for (std::size_t q = 0; q < m.size(); ++q) {
              m[q] += wt * mv[q];
              e[q] += wt * wt * ev[q] * ev[q];
            }
          };
          acc(mean.ext, se.ext, mk_.ext, sk_.ext);
          acc(mean.abs, se.abs, mk_.abs, sk_.abs);
          acc(mean.sca, se.sca, mk_.sca, sk_.sca);
          acc(mean.np_ext, se.np_ext, mk_.np_ext, sk_.np_ext);
          acc(mean.np_abs, se.np_abs, mk_.np_abs, sk_.np_abs);
          acc(mean.np_sca, se.np_sca, mk_.np_sca, sk_.np_sca);
        }
        out << "nm\n";
        for (auto* x : {&se.ext, &se.abs, &se.sca, &se.np_ext, &se.np_abs, &se.np_sca})
          for (auto& v : *x) v = std::sqrt(v);
      }
      const auto mode = sys.dipole_mode_coupling(w);
      sr.mode_omega_cm = mode.omega_p_cm;
      sr.mode_gamma_cm = mode.gamma_p_cm;
      sr.g_micro_cm = mode.g_rms_cm;
      {
        app::CsvWriter wr(app::join_path(dir, "np_mode.csv"),
                          {"dipolar plasmon mode of the NP: α_a(ω) ≈ A/(ω_p² − ω² − iγω) + B fitted per axis;",
                           "g: collective coupling of that mode to all molecules, g² = (A/2ω_p)·Tr K_a (no peak fitting)"},
                          {"axis", "omega_p_cm", "gamma_p_cm", "g_cm"});
        const char* ax[3] = {"x", "y", "z"};
        for (int a = 0; a < 3; ++a)
          wr.row_text({ax[a], fmt(mode.omega_axis[a], 7), fmt(mode.gamma_axis[a], 5), fmt(mode.g_cm[a], 5)});
      }
      out << "NP dipolar mode: ω_p = " << fmt(mode.omega_p_cm, 6) << " cm⁻¹, γ_p = " << fmt(mode.gamma_p_cm, 4)
          << " cm⁻¹ (pole fit, rel. rms " << fmt(mode.fit_rel_rms, 2) << "); collective coupling g = "
          << fmt(mode.g_rms_cm, 4) << " cm⁻¹ (x/y/z: " << fmt(mode.g_cm[0], 3) << "/" << fmt(mode.g_cm[1], 3) << "/"
          << fmt(mode.g_cm[2], 3) << ")\n";
      {
        const auto& o = er.observables;
        const double gx = o.count("abs_fwhm_cm") ? o.at("abs_fwhm_cm").mean : 0.0;
        const double g = mode.g_rms_cm, gp = mode.gamma_p_cm;
        out << "  criteria (γ = FWHM): 2g = " << fmt(2 * g, 4) << " vs (γp+γx)/2 = " << fmt(0.5 * (gp + gx), 4)
            << (2 * g > 0.5 * (gp + gx) ? " → strong coupling expected" : " → not strong") << "; C = 4g²/(γpγx) = "
            << fmt(4 * g * g / std::max(gp * gx, 1e-300), 3) << "\n";
      }
      sr.has_np = true;
      sr.np_w = w;
      sr.ext = mean.ext;
      sr.abs_sys = mean.abs;
      sr.sca = mean.sca;
      sr.np_ext = mean.np_ext;
      sr.np_abs = mean.np_abs;
      sr.np_sca = mean.np_sca;
      const double trunc = sys.truncation_estimate(np.near_field_cm);
      out << "NP: realizations " << er.stored.size() << ", multipole truncation estimate " << fmt(trunc, 3) << "\n";
      if (trunc > 1e-2) out << "  ! increase nanoparticle.l_max or the aggregate–NP gap\n";
      {
        app::CsvWriter wr(app::join_path(dir, "plexciton.csv"),
                          {"cross sections of the aggregate+NP system and of the NP alone, nm² (orientation average)",
                           "multipole truncation estimate: " + fmt(trunc, 3)},
                          {"w_cm", "ext", "ext_se", "abs", "abs_se", "sca", "sca_se", "np_ext", "np_abs", "np_sca"});
        for (std::size_t q = 0; q < w.size(); ++q)
          wr.row({w[q], mean.ext[q], se.ext[q], mean.abs[q], se.abs[q], mean.sca[q], se.sca[q], mean.np_ext[q],
                  mean.np_abs[q], mean.np_sca[q]});
      }
      // Излучательная скорость мономера без НЧ (та же модель, что и для экситонов)
      double gamma0 = 0.0;
      {
        ExcitonStates m1;
        m1.energies = Eigen::VectorXd::Constant(1, rc.model.monomer_energy_cm);
        m1.c = Eigen::MatrixXd::Identity(1, 1);
        m1.mu = Eigen::MatrixX3d::Zero(1, 3);
        m1.mu(0, 2) = rc.model.coupling.mu_D;
        m1.mu2 = Eigen::VectorXd::Constant(1, rc.model.coupling.mu_D * rc.model.coupling.mu_D);
        m1.ipr = Eigen::VectorXd::Ones(1);
        gamma0 = radiative_rates_ps(m1, rc.model.decay)(0);
      }
      const auto nf = sys.near_field(np.near_field_cm, gamma0, rc.model.decay.nonradiative_ps);
      double maxE = 0.0, mingap = 1e300;
      {
        app::CsvWriter wr(app::join_path(dir, "molecules_near_field.csv"),
                          {"NP near field at the molecules at ν̃ = " + fmt(np.near_field_cm) + " cm⁻¹",
                           "enhancement = ⟨|E_loc/E₀|²⟩ over polarizations; shift = Re W_ii; decay = −2 Im W_ii (decay into the NP)",
                           "rad_enhancement = Γ_rad/Γ_rad,0 = |L_iᵀμ̂_i|² (Purcell, radiative part); nonrad_to_np = Γ₀ − 2Im W_ii − Γ_rad;",
                           "quantum_yield — isolated monomer at this position (Γ₀ = " + fmt(gamma0, 4) +
                               " ps⁻¹, intrinsic Γ_nr = " + fmt(rc.model.decay.nonradiative_ps, 4) + " ps⁻¹)"},
                          {"x_nm", "y_nm", "z_nm", "gap_nm", "field_enhancement", "shift_cm", "decay_to_np_ps",
                           "rad_enhancement", "nonrad_to_np_ps", "quantum_yield", "quantum_yield_free"});
        for (std::size_t i = 0; i < nf.size(); ++i) {
          wr.row({agg.positions[i].x(), agg.positions[i].y(), agg.positions[i].z(), nf[i].gap_nm,
                  nf[i].field_enhancement, nf[i].shift_cm, nf[i].decay_ps, nf[i].rad_enhancement,
                  nf[i].nonrad_to_np_ps, nf[i].quantum_yield, nf[i].quantum_yield_free});
          maxE = std::max(maxE, nf[i].field_enhancement);
          mingap = std::min(mingap, nf[i].gap_nm);
        }
      }
      np_cols = {fmt(maxE), fmt(mingap), fmt(trunc, 3)};
      if (np.emission && er.engine) {
        const auto em = sys.emission(er.stored, *er.engine, w, rc.model.decay.nonradiative_ps,
                                     rc.model.pump == PumpType::Laser, rc.model.excitation_cm);
        out << "hybrid emission: QY " << fmt(em.qy, 4) << " (without NP " << fmt(em.qy_free, 4)
            << "), mean excited-state lifetime " << fmt(em.lifetime_ps, 4) << " ps (without NP "
            << fmt(em.lifetime_free_ps, 4) << "), brightest state Γ_rad enhancement " << fmt(em.rad_enhancement_bright, 4)
            << "\n";
        app::CsvWriter wr(app::join_path(dir, "hybrid_emission.csv"),
                          {"PL of the aggregate next to the NP vs without it (same pump absorbed), photons/(ps·cm⁻¹)",
                           "QY " + fmt(em.qy, 4) + " vs " + fmt(em.qy_free, 4) + "; lifetime " + fmt(em.lifetime_ps, 4) +
                               " vs " + fmt(em.lifetime_free_ps, 4) + " ps",
                           "Γ_rad via effective dipole Σ c_in L_iᵀμ_i (radiation through the NP), Γ_nr,NP = absorption in the metal"},
                          {"w_cm", "pl_hybrid", "pl_free"});
        for (std::size_t q = 0; q < w.size(); ++q) wr.row({w[q], em.pl[q], em.pl_free[q]});
        sr.pl_hybrid = em.pl;
        sr.pl_free = em.pl_free;
        sr.qy_hybrid = em.qy;
        sr.qy_free = em.qy_free;
        np_cols.push_back(fmt(em.qy, 4));
        np_cols.push_back(fmt(em.qy_free, 4));
      } else {
        np_cols.push_back("");
        np_cols.push_back("");
      }
    }

    std::vector<std::string> row = {cfg.name, cfg.type, std::to_string(agg.size()), std::to_string(er.realizations),
                                    er.converged ? "1" : "0"};
    for (const auto& k : obs_names) {
      const auto it = er.observables.find(k);
      row.push_back(it != er.observables.end() ? fmt(it->second.mean, 8) : "");
      row.push_back(it != er.observables.end() ? fmt(it->second.se(), 4) : "");
    }
    row.insert(row.end(), np_cols.begin(), np_cols.end());
    summary.row_text(row);

    const auto& o = er.observables;
    out << "absorption: " << fmt(o.at("abs_peak_cm").mean, 7) << " ± " << fmt(o.at("abs_peak_cm").se(), 2)
              << " cm⁻¹, FWHM " << fmt(o.at("abs_fwhm_cm").mean, 4) << "; PL: " << fmt(o.at("em_peak_cm").mean, 7)
              << "; Stokes " << fmt(o.at("stokes_shift_cm").mean, 4) << "; QY " << fmt(o.at("quantum_yield").mean, 4)
              << "; superradiance " << fmt(o.at("superradiance_factor").mean, 4) << "\n";
  }
  const double secs =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  out << "\nDone in " << fmt(secs, 3) << " s. Results: " << rc.output_dir << "\n";
  return result;
}

}  // namespace ssp::app
