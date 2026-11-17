#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/modified_redfield.hpp"

#include "ssp/excitons/spectra.hpp"
#include "ssp/units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ssp::excitons {

std::vector<double> SpectrumGrid::points() const {
  if (nw < 2 || wmax_cm <= wmin_cm) throw std::invalid_argument("spectrum: invalid grid");
  std::vector<double> w(nw);
  for (std::size_t i = 0; i < nw; ++i) w[i] = wmin_cm + (wmax_cm - wmin_cm) * double(i) / double(nw - 1);
  return w;
}

void RunningStat::add(double x) {
  ++n;
  const double d = x - mean;
  mean += d / double(n);
  m2 += d * (x - mean);
}
double RunningStat::sd() const {
  if (se_fixed >= 0.0) return se_fixed * std::sqrt(double(std::max<long>(n, 1)));
  return n > 1 ? std::sqrt(m2 / double(n - 1)) : 0.0;
}
double RunningStat::se() const {
  if (se_fixed >= 0.0) return se_fixed;
  return n > 1 ? sd() / std::sqrt(double(n)) : 0.0;
}

const std::map<std::string, std::string>& observable_descriptions() {
  static const std::map<std::string, std::string> d = {
      {"abs_peak_cm", "absorption maximum of the ensemble-averaged spectrum, cm⁻¹"},
      {"abs_fwhm_cm", "FWHM of the main band of the ensemble-averaged absorption (homogeneous + inhomogeneous), cm⁻¹"},
      {"abs_shift_cm", "shift of the absorption maximum relative to the monomer (J<0, H>0), cm⁻¹"},
      {"abs_centroid_cm", "first moment of the averaged absorption within the spectral window, cm⁻¹"},
      {"abs_width_sd_cm", "standard deviation (2nd central moment) of the averaged absorption, cm⁻¹"},
      {"abs_n_peaks", "number of absorption maxima above 20 % of the main one (1 = single band)"},
      {"em_peak_cm", "fluorescence maximum of the ensemble-averaged spectrum, cm⁻¹"},
      {"em_fwhm_cm", "FWHM of the main band of the averaged fluorescence, cm⁻¹"},
      {"em_centroid_cm", "first moment of the averaged fluorescence, cm⁻¹"},
      {"stokes_shift_cm", "absorption maximum − fluorescence maximum (averaged spectra), cm⁻¹; unreliable if abs_n_peaks > 1"},
      {"stokes_centroid_cm", "absorption centroid − fluorescence centroid, cm⁻¹"},
      {"lowest_state_cm", "energy of the lowest exciton state (mean over realizations), cm⁻¹"},
      {"brightest_state_cm", "energy of the brightest state (mean over realizations), cm⁻¹"},
      {"superradiance_factor", "|μ|² of the brightest state / monomer μ² (mean over realizations)"},
      {"lowest_state_strength", "|μ|² of the lowest state / monomer μ² (mean over realizations)"},
      {"delocalization_brightest", "delocalization number 1/Σc⁴ of the brightest state (mean over realizations)"},
      {"delocalization_thermal", "thermal average of 1/Σc⁴ over steady-state populations"},
      {"quantum_yield", "fluorescence quantum yield"},
      {"radiative_rate_ps", "mean radiative rate ΣΓ_rad P/ΣP, ps⁻¹"},
      {"fluorescence_lifetime_ps", "fluorescence lifetime ΣP/Σpump, ps"},
      {"relaxation_time_ps", "intraband equilibration time (−1/λ₂ of the Redfield matrix), ps"},
      {"homogeneous_hwhm_brightest_cm", "Markovian (Redfield) homogeneous HWHM of the brightest transition, cm⁻¹ — estimate only"},
  };
  return d;
}

namespace {

// Положение максимума и FWHM главного пика (линейная интерполяция краёв).
void peak_and_fwhm(const std::vector<double>& w, const std::vector<double>& y, double& peak, double& fwhm) {
  const auto it = std::max_element(y.begin(), y.end());
  const std::size_t k = static_cast<std::size_t>(it - y.begin());
  const double ymax = *it;
  peak = w[k];
  if (k > 0 && k + 1 < w.size()) {  // параболическое уточнение максимума
    const double a = y[k - 1], b = y[k], c = y[k + 1];
    const double den = a - 2.0 * b + c;
    if (den < 0.0) peak = w[k] + 0.5 * (a - c) / den * (w[1] - w[0]);
  }
  const double half = 0.5 * ymax;
  std::size_t l = k, r = k;
  while (l > 0 && y[l] > half) --l;
  while (r + 1 < y.size() && y[r] > half) ++r;
  auto cross = [&](std::size_t i, std::size_t j) {
    const double t = (half - y[i]) / (y[j] - y[i]);
    return w[i] + t * (w[j] - w[i]);
  };
  const double wl = (y[l] <= half && l < k) ? cross(l, l + 1) : w[l];
  const double wr = (y[r] <= half && r > k) ? cross(r - 1, r) : w[r];
  fwhm = wr - wl;
}

// Число локальных максимумов выше rel·max.
int count_peaks(const std::vector<double>& y, double rel) {
  if (y.size() < 3) return 0;
  const double m = *std::max_element(y.begin(), y.end());
  int n = 0;
  for (std::size_t i = 1; i + 1 < y.size(); ++i)
    if (y[i] > y[i - 1] && y[i] >= y[i + 1] && y[i] > rel * m) ++n;
  return n;
}

void moments(const std::vector<double>& w, const std::vector<double>& y, double& c, double& sd) {
  double s0 = 0, s1 = 0, s2 = 0;
  for (std::size_t i = 0; i < w.size(); ++i) {
    s0 += y[i];
    s1 += y[i] * w[i];
  }
  c = s0 > 0 ? s1 / s0 : NAN;
  for (std::size_t i = 0; i < w.size(); ++i) s2 += y[i] * (w[i] - c) * (w[i] - c);
  sd = s0 > 0 ? std::sqrt(s2 / s0) : NAN;
}

}  // namespace

std::map<std::string, double> spectral_metrics(const std::vector<double>& w, const std::vector<double>& abs,
                                               const std::vector<double>& em, double monomer_cm) {
  std::map<std::string, double> o;
  double pk, fw, c, sd;
  peak_and_fwhm(w, abs, pk, fw);
  o["abs_peak_cm"] = pk;
  o["abs_fwhm_cm"] = fw;
  o["abs_shift_cm"] = pk - monomer_cm;
  moments(w, abs, c, sd);
  o["abs_centroid_cm"] = c;
  o["abs_width_sd_cm"] = sd;
  o["abs_n_peaks"] = count_peaks(abs, 0.2);
  double epk, efw, ec, esd;
  peak_and_fwhm(w, em, epk, efw);
  moments(w, em, ec, esd);
  o["em_peak_cm"] = epk;
  o["em_fwhm_cm"] = efw;
  o["em_centroid_cm"] = ec;
  o["stokes_shift_cm"] = pk - epk;
  o["stokes_centroid_cm"] = c - ec;
  return o;
}

Realization make_realization(const Aggregate& agg, const Eigen::MatrixXd& J, const ExcitonModel& m,
                             std::uint64_t seed, const Eigen::MatrixXd* chol, const LineshapeEngine* engine) {
  Realization R;
  std::mt19937_64 rng(seed);
  // Структурный беспорядок: своя геометрия и J для каждой реализации
  const bool structural = m.orientation_jitter_deg > 0.0 || m.position_jitter_nm > 0.0;
  Aggregate aj;
  if (structural) {
    aj = agg;
    std::normal_distribution<double> g(0.0, 1.0);
    const double sa = m.orientation_jitter_deg * units::pi / 180.0;
    for (std::size_t i = 0; i < aj.size(); ++i) {
      if (m.position_jitter_nm > 0.0) aj.positions[i] += m.position_jitter_nm * Vec3(g(rng), g(rng), g(rng));
      if (sa > 0.0) {
        Vec3 axis(g(rng), g(rng), g(rng));
        axis = (axis - axis.dot(aj.dipoles[i]) * aj.dipoles[i]).normalized();
        aj.dipoles[i] = Eigen::AngleAxisd(sa * g(rng), axis) * aj.dipoles[i];
      }
    }
  }
  const Aggregate& A = structural ? aj : agg;
  const Eigen::MatrixXd Jr = structural ? coupling_matrix(A, m.coupling) : J;
  const bool shifted = m.site_shift_cm.size() == static_cast<Eigen::Index>(A.size());
  if (chol) {
    std::normal_distribution<double> g(0.0, 1.0);
    Eigen::VectorXd z(Jr.rows());
    for (Eigen::Index i = 0; i < z.size(); ++i) z(i) = g(rng);
    R.H = Jr;
    R.H.diagonal() = Eigen::VectorXd::Constant(z.size(), m.monomer_energy_cm) + (*chol) * z;
  } else {
    R.H = frenkel_hamiltonian(Jr, m.monomer_energy_cm, m.disorder_sigma_cm, rng);
  }
  if (shifted) R.H.diagonal() += m.site_shift_cm;
  R.states = diagonalize(R.H, A, m.coupling.mu_D);
  const bool modified = m.relaxation == RelaxationTheory::ModifiedRedfield && engine;
  const Eigen::MatrixXd K =
      modified ? modified_redfield_rate_matrix(R.states, *engine) : redfield_rate_matrix(R.states, m.bath);
  R.db_energies = modified ? Eigen::VectorXd(R.states.energies - engine->lambda_total() * R.states.ipr)
                           : R.states.energies;
  R.K_ps = K * units::rad_ps_per_cm;
  R.hwhm_cm = redfield_homogeneous_hwhm(R.states, m.bath, K);
  DecayParams dp = m.decay;
  dp.mu_D = m.coupling.mu_D;
  dp.monomer_energy_cm = m.monomer_energy_cm;
  R.gamma_rad_ps = radiative_rates_ps(R.states, dp);
  return R;
}

EnsembleResult run_ensemble(const Aggregate& agg, const ExcitonModel& model_in, const SpectrumGrid& grid,
                            const EnsembleSettings& settings, int n_store) {
  ExcitonModel model = model_in;
  if (model.dispersion.enabled && model.site_shift_cm.size() != static_cast<Eigen::Index>(agg.size()))
    model.site_shift_cm = dispersion_site_shifts(agg, model.dispersion).shift_cm;
  EnsembleResult res;
  res.w = grid.points();
  const std::size_t nw = res.w.size();
  const Eigen::MatrixXd J = coupling_matrix(agg, model.coupling);
  const double N = static_cast<double>(agg.size());

  // Диапазон отстроек |ω − ε_n|: оценка спектра по Гершгорину + 6σ беспорядка.
  const double jmax = J.cwiseAbs().rowwise().sum().maxCoeff();
  const double smin = model.site_shift_cm.size() ? model.site_shift_cm.minCoeff() : 0.0;
  const double smax = model.site_shift_cm.size() ? model.site_shift_cm.maxCoeff() : 0.0;
  const double emin = model.monomer_energy_cm + smin - jmax - 6.0 * model.disorder_sigma_cm;
  const double emax = model.monomer_energy_cm + smax + jmax + 6.0 * model.disorder_sigma_cm;
  double wlo = grid.wmin_cm, whi = grid.wmax_cm;
  if (settings.extra_wmax_cm > settings.extra_wmin_cm) {
    wlo = std::min(wlo, settings.extra_wmin_cm);
    whi = std::max(whi, settings.extra_wmax_cm);
  }
  const double half_range = 1.2 * std::max(whi - emin, emax - wlo) + 4.0 * model.bath.reorganization() + 500.0;
  res.engine = std::make_shared<LineshapeEngine>(model.bath, half_range, grid.t_max_ps / units::ps_per_tu);
  const LineshapeEngine& engine = *res.engine;
  // Коррелированный энергетический беспорядок: множитель Холецкого ковариации σ² exp(−r/ξ)
  Eigen::MatrixXd chol;
  if (model.disorder_corr_length_nm > 0.0 && model.disorder_sigma_cm > 0.0) {
    const auto Nn = static_cast<Eigen::Index>(agg.size());
    Eigen::MatrixXd C(Nn, Nn);
    for (Eigen::Index i = 0; i < Nn; ++i)
      for (Eigen::Index j = 0; j < Nn; ++j)
        C(i, j) = model.disorder_sigma_cm * model.disorder_sigma_cm *
                  std::exp(-(agg.positions[i] - agg.positions[j]).norm() / model.disorder_corr_length_nm);
    C.diagonal().array() += 1e-12 * model.disorder_sigma_cm * model.disorder_sigma_cm;
    Eigen::LLT<Eigen::MatrixXd> llt(C);
    if (llt.info() != Eigen::Success) throw std::runtime_error("correlated disorder: covariance is not positive definite");
    chol = llt.matrixL();
  }
  const Eigen::MatrixXd* cholp = chol.size() ? &chol : nullptr;
  const double sigma_pref = absorption_cross_section_prefactor_nm2(model.n_medium);
  const double beta = model.bath.beta();

  std::vector<RunningStat> abs_stat(nw), em_stat(nw);
  const double nr = model.decay.nonradiative_ps;

  std::vector<std::vector<double>> all_abs, all_em;  // спектры всех реализаций (для бутстрэпа)
  res.dyn_t = settings.dyn_times_ps;
  int done = 0;
  while (done < settings.max_realizations) {
    const int batch = std::min(settings.batch, settings.max_realizations - done);
    std::vector<Realization> reals(batch);
    std::vector<std::vector<double>> abs_b(batch, std::vector<double>(nw, 0.0)), em_b(batch, std::vector<double>(nw, 0.0));
    std::vector<std::map<std::string, double>> obs_b(batch);
    std::vector<double> tail_b(batch, 0.0);
    std::vector<std::vector<std::vector<double>>> dyn_b(batch);

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < batch; ++b) {
      Realization R = make_realization(agg, J, model, settings.seed + static_cast<std::uint64_t>(done + b), cholp, &engine);
      const ExcitonStates& st = R.states;
      const Eigen::Index n = static_cast<Eigen::Index>(st.size());
      const Eigen::VectorXd Gout = -R.K_ps.diagonal() / units::rad_ps_per_cm;  // угл. см⁻¹

      // Поглощение
      double tail = 0.0;
      for (Eigen::Index k = 0; k < n; ++k) {
        engine.add_line(abs_b[b], res.w, st.energies(k), sigma_pref * st.mu2(k), st.ipr(k), Gout(k), false, 1);
        tail = std::max(tail, engine.tail(st.ipr(k), Gout(k)));
      }
      // Накачка
      R.pump = Eigen::VectorXd::Zero(n);
      if (model.pump == PumpType::OscillatorStrength) {
        R.pump = st.mu2 / st.mu2.sum();
      } else {
        const std::vector<double> wx = {model.excitation_cm};
        for (Eigen::Index k = 0; k < n; ++k) {
          std::vector<double> v = {0.0};
          engine.add_line(v, wx, st.energies(k), st.mu2(k), st.ipr(k), Gout(k), false, 1);
          R.pump(k) = std::max(v[0], 0.0);
        }
        const double s = R.pump.sum();
        if (s > 0.0) R.pump /= s;
        else R.pump = st.mu2 / st.mu2.sum();
      }
      // Стационарные населённости и флуоресценция
      const Eigen::VectorXd decay = R.gamma_rad_ps.array() + nr;
      R.P = steady_state(R.K_ps, decay, R.pump);
      for (Eigen::Index k = 0; k < n; ++k) {
        const double eps = st.energies(k);
        const double weight = R.P(k) * R.gamma_rad_ps(k) / (eps * eps * eps);
        engine.add_line(em_b[b], res.w, eps, weight, st.ipr(k), Gout(k), true, 3);
      }

      // Скалярные наблюдаемые
      // Неспектральные величины (свойства состояний) — среднее по реализациям;
      // спектральные — из усреднённого спектра после цикла.
      auto& o = obs_b[b];
      Eigen::Index ib;
      st.mu2.maxCoeff(&ib);
      const double mu2m = model.coupling.mu_D * model.coupling.mu_D;
      o["lowest_state_cm"] = st.energies(0);
      o["brightest_state_cm"] = st.energies(ib);
      o["superradiance_factor"] = st.mu2(ib) / mu2m;
      o["lowest_state_strength"] = st.mu2(0) / mu2m;
      o["delocalization_brightest"] = 1.0 / st.ipr(ib);
      const double Psum = R.P.sum();
      o["delocalization_thermal"] = (R.P.array() / st.ipr.array()).sum() / Psum;
      o["quantum_yield"] = R.gamma_rad_ps.dot(R.P) / R.pump.sum();
      o["radiative_rate_ps"] = R.gamma_rad_ps.dot(R.P) / Psum;
      o["fluorescence_lifetime_ps"] = Psum / R.pump.sum();
      o["homogeneous_hwhm_brightest_cm"] = R.hwhm_cm(ib);
      if (n > 1) {
        // −1/λ₂ симметризованной матрицы Редфилда (без распада)
        Eigen::VectorXd d(n);
        for (Eigen::Index k = 0; k < n; ++k)
          d(k) = std::isinf(beta) ? 1.0 : std::exp(-0.5 * beta * (st.energies(k) - st.energies(0)));
        Eigen::MatrixXd S = d.cwiseInverse().asDiagonal() * R.K_ps * d.asDiagonal();
        S = 0.5 * (S + S.transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S, Eigen::EigenvaluesOnly);
        const double l2 = es.eigenvalues()(n - 2);
        o["relaxation_time_ps"] = l2 < 0.0 ? -1.0 / l2 : std::numeric_limits<double>::infinity();
      }
      if (!settings.dyn_times_ps.empty()) {
        Eigen::VectorXd P0 = Eigen::VectorXd::Zero(n);
        if (settings.dyn_initial == "pump") P0 = R.pump;
        else if (settings.dyn_initial == "highest") P0(n - 1) = 1.0;
        else P0(ib) = 1.0;
        const Eigen::VectorXd t = Eigen::Map<const Eigen::VectorXd>(settings.dyn_times_ps.data(),
                                                                    static_cast<Eigen::Index>(settings.dyn_times_ps.size()));
        const Eigen::MatrixXd P = population_dynamics(R.K_ps, R.gamma_rad_ps.array() + nr, R.db_energies, beta, P0, t);
        auto& dv = dyn_b[b];
        dv.assign(4, std::vector<double>(static_cast<std::size_t>(t.size())));
        for (Eigen::Index k = 0; k < t.size(); ++k) {
          const double tot = P.col(k).sum();
          dv[0][k] = tot;
          dv[1][k] = P(0, k);
          dv[2][k] = P(ib, k);
          dv[3][k] = tot > 0 ? st.energies.dot(P.col(k)) / tot : NAN;
        }
      }
      tail_b[b] = tail;
      reals[b] = std::move(R);
    }

    for (int b = 0; b < batch; ++b) {
      for (std::size_t q = 0; q < nw; ++q) {
        abs_stat[q].add(abs_b[b][q]);
        em_stat[q].add(em_b[b][q]);
      }
      all_abs.push_back(abs_b[b]);
      all_em.push_back(em_b[b]);
      if (!dyn_b[b].empty()) {
        if (res.dyn.empty()) res.dyn.assign(4, std::vector<RunningStat>(dyn_b[b][0].size()));
        for (int k = 0; k < 4; ++k)
          for (std::size_t j = 0; j < dyn_b[b][k].size(); ++j)
            if (std::isfinite(dyn_b[b][k][j])) res.dyn[k][j].add(dyn_b[b][k][j]);
      }
      for (const auto& [k, v] : obs_b[b])
        if (std::isfinite(v)) res.observables[k].add(v);
      res.max_lineshape_tail = std::max(res.max_lineshape_tail, tail_b[b]);
      if (done + b == 0) res.representative = reals[b];
      if (static_cast<int>(res.stored.size()) < n_store) res.stored.push_back(reals[b]);
    }
    done += batch;

    if (done >= settings.min_realizations) {
      double amax = 0.0, emaxv = 0.0, ase = 0.0, ese = 0.0;
      for (std::size_t q = 0; q < nw; ++q) {
        amax = std::max(amax, abs_stat[q].mean);
        emaxv = std::max(emaxv, em_stat[q].mean);
        ase = std::max(ase, abs_stat[q].se());
        ese = std::max(ese, em_stat[q].se());
      }
      if (ase <= settings.rel_tolerance * amax && ese <= settings.rel_tolerance * emaxv) {
        res.converged = true;
        break;
      }
    }
  }
  res.realizations = done;
  res.abs_mean.resize(nw);
  res.abs_se.resize(nw);
  res.em_mean.resize(nw);
  res.em_se.resize(nw);
  for (std::size_t q = 0; q < nw; ++q) {
    res.abs_mean[q] = abs_stat[q].mean;
    res.abs_se[q] = abs_stat[q].se();
    res.em_mean[q] = em_stat[q].mean;
    res.em_se[q] = em_stat[q].se();
  }
  // Спектральные наблюдаемые — из усреднённого спектра; погрешность — бутстрэп по реализациям.
  {
    const auto m0 = spectral_metrics(res.w, res.abs_mean, res.em_mean, model.monomer_energy_cm);
    std::map<std::string, std::vector<double>> boot;
    std::mt19937_64 rng(settings.seed ^ 0x9e3779b97f4a7c15ULL);
    std::uniform_int_distribution<int> pick(0, done - 1);
    std::vector<double> a(nw), e(nw);
    for (int it = 0; it < settings.bootstrap && done > 1; ++it) {
      std::fill(a.begin(), a.end(), 0.0);
      std::fill(e.begin(), e.end(), 0.0);
      for (int r = 0; r < done; ++r) {
        const int k = pick(rng);
        for (std::size_t q = 0; q < nw; ++q) {
          a[q] += all_abs[k][q];
          e[q] += all_em[k][q];
        }
      }
      for (const auto& [k, v] : spectral_metrics(res.w, a, e, model.monomer_energy_cm)) boot[k].push_back(v);
    }
    for (const auto& [k, v] : m0) {
      RunningStat st;
      st.n = done;
      st.mean = v;
      double s1 = 0, s2 = 0;
      int nb = 0;
      for (double x : boot[k])
        if (std::isfinite(x)) {
          s1 += x;
          s2 += x * x;
          ++nb;
        }
      st.se_fixed = nb > 1 ? std::sqrt(std::max(0.0, s2 / nb - (s1 / nb) * (s1 / nb))) : 0.0;
      res.observables[k] = st;
    }
  }
  (void)N;
  return res;
}

}  // namespace ssp::excitons
