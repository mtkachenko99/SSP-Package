#include "config.hpp"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <set>
#include <stdexcept>

namespace ssp::app {

namespace {

// Обёртка над QJsonObject с учётом прочитанных ключей.
class Obj {
 public:
  Obj(const QJsonObject& o, std::string where) : o_(o), where_(std::move(where)) {}
  ~Obj() noexcept(false) {
    if (std::uncaught_exceptions() > 0) return;
    for (const auto& k : o_.keys())
      if (!used_.contains(k.toStdString()) && !k.startsWith("_"))
        throw std::invalid_argument(where_ + ": unknown key '" + k.toStdString() + "'");
  }
  bool has(const char* k) const { return o_.contains(k); }
  double num(const char* k, double def) {
    used_.insert(k);
    if (!o_.contains(k)) return def;
    if (!o_[k].isDouble()) throw std::invalid_argument(where_ + "." + k + ": number expected");
    return o_[k].toDouble();
  }
  int integer(const char* k, int def) { return static_cast<int>(num(k, def)); }
  bool flag(const char* k, bool def) {
    used_.insert(k);
    if (!o_.contains(k)) return def;
    if (!o_[k].isBool()) throw std::invalid_argument(where_ + "." + k + ": true/false expected");
    return o_[k].toBool();
  }
  std::string str(const char* k, const std::string& def) {
    used_.insert(k);
    if (!o_.contains(k)) return def;
    if (!o_[k].isString()) throw std::invalid_argument(where_ + "." + k + ": string expected");
    return o_[k].toString().toStdString();
  }
  QJsonObject obj(const char* k) {
    used_.insert(k);
    if (!o_.contains(k)) return {};
    if (!o_[k].isObject()) throw std::invalid_argument(where_ + "." + k + ": object expected");
    return o_[k].toObject();
  }
  QJsonArray arr(const char* k) {
    used_.insert(k);
    if (!o_.contains(k)) return {};
    if (!o_[k].isArray()) throw std::invalid_argument(where_ + "." + k + ": array expected");
    return o_[k].toArray();
  }
  QJsonValue raw(const char* k) {
    used_.insert(k);
    return o_.contains(k) ? o_[k] : QJsonValue(QJsonValue::Undefined);
  }
  std::string where(const char* k) const { return where_ + "." + k; }

 private:
  QJsonObject o_;
  std::string where_;
  std::set<std::string> used_;
};

excitons::Vec3 vec3(const QJsonArray& a, const std::string& where) {
  if (a.isEmpty()) return excitons::Vec3::Zero();
  if (a.size() != 3) throw std::invalid_argument(where + ": array of 3 numbers expected");
  return {a[0].toDouble(), a[1].toDouble(), a[2].toDouble()};
}

}  // namespace

RunConfig load_config(const std::string& path) {
  QFile f(QString::fromStdString(path));
  if (!f.open(QIODevice::ReadOnly)) throw std::runtime_error("cannot open config: " + path);
  const QByteArray data = f.readAll();
  return load_config_text(data.toStdString(),
                          QFileInfo(QString::fromStdString(path)).absolutePath().toStdString(), path);
}

RunConfig load_config_text(const std::string& json, const std::string& base_dir, const std::string& path) {
  QJsonParseError err;
  const QJsonDocument doc = QJsonDocument::fromJson(QByteArray::fromStdString(json), &err);
  if (err.error != QJsonParseError::NoError)
    throw std::runtime_error("JSON error: " + err.errorString().toStdString() + " (offset " +
                             std::to_string(err.offset) + ")");
  if (!doc.isObject()) throw std::runtime_error("the config must be a JSON object");

  RunConfig rc;
  rc.config_path = path;
  rc.config_dir = base_dir;
  Obj root(doc.object(), "config");
  rc.output_dir = root.str("output_dir", "ssp_results");
  // Относительный output_dir — относительно папки конфига (как и material_file)
  if (QFileInfo(QString::fromStdString(rc.output_dir)).isRelative() && !base_dir.empty())
    rc.output_dir = QDir(QString::fromStdString(base_dir)).absoluteFilePath(QString::fromStdString(rc.output_dir)).toStdString();

  // --- конфигурации (формы)
  for (const auto& v : root.arr("configurations")) {
    Obj c(v.toObject(), "configurations[]");
    Configuration cfg;
    cfg.type = c.str("geometry_type", "");
    cfg.name = c.str("name", cfg.type);
    if (cfg.type.empty()) throw std::invalid_argument("configurations[]: geometry_type is missing");
    const QJsonObject gp = c.obj("geometry_params");
    for (const auto& k : gp.keys()) {
      const QJsonValue x = gp[k];
      if (x.isBool()) cfg.params[k.toStdString()] = x.toBool() ? 1.0 : 0.0;
      else if (x.isDouble()) cfg.params[k.toStdString()] = x.toDouble();
      else throw std::invalid_argument("geometry_params." + k.toStdString() + ": number or true/false expected");
    }
    if (c.has("placement")) {
      Obj p(c.obj("placement"), "placement");
      cfg.placement.shift_nm = vec3(p.arr("shift_nm"), "placement.shift_nm");
      const QJsonArray e = p.arr("euler_deg");
      if (!e.isEmpty()) {
        if (e.size() != 3) throw std::invalid_argument("placement.euler_deg: 3 angles expected");
        for (int k = 0; k < 3; ++k) cfg.placement.euler_deg[k] = e[k].toDouble();
      }
      cfg.placement.center = p.flag("center", false);
    }
    rc.configurations.push_back(cfg);
  }
  if (rc.configurations.empty()) throw std::invalid_argument("config: the configurations list is empty");

  // --- экситонная модель
  {
    Obj e(root.obj("exciton"), "exciton");
    auto& m = rc.model;
    m.monomer_energy_cm = e.num("monomer_energy_cm", 17000.0);
    m.coupling.mu_D = e.num("transition_dipole_D", 10.0);
    m.disorder_sigma_cm = e.num("disorder_sigma_cm", 60.0);
    m.disorder_corr_length_nm = e.num("disorder_corr_length_nm", 0.0);
    m.orientation_jitter_deg = e.num("orientation_jitter_deg", 0.0);
    m.position_jitter_nm = e.num("position_jitter_nm", 0.0);
    bool use_coupling_charges = false;
    if (e.has("dispersion")) {
      Obj d(e.obj("dispersion"), "exciton.dispersion");
      auto& dp = m.dispersion;
      dp.enabled = d.flag("enabled", true);
      dp.F_eV = d.num("F_eV", 1.43);
      dp.n_aggregate = d.num("n_aggregate", 1.78);
      dp.n_solvent = d.num("n_solvent", m.n_medium);
      dp.dipole_length_nm = d.num("dipole_length_nm", 1.0);
      dp.cutoff_nm = d.num("cutoff_nm", 3.0);
      dp.softening_nm = d.num("softening_nm", 0.3);
      dp.extra_shift_cm = d.num("extra_shift_cm", 0.0);
      dp.orientation_samples = d.integer("orientation_samples", 16);
      // атомные переходные заряды для V_mk: свой файл или те же, что для связей
      const std::string cf = d.str("charges_file", "");
      if (!cf.empty()) {
        const QString q = QString::fromStdString(cf);
        dp.charges = excitons::load_transition_charges(
            (QFileInfo(q).isAbsolute() ? QFileInfo(q) : QFileInfo(QDir(QString::fromStdString(rc.config_dir)), q))
                .absoluteFilePath()
                .toStdString());
      } else {
        use_coupling_charges = d.flag("use_coupling_charges", false);  // заряды связи читаются ниже
      }
    }
    const std::string rel = e.str("relaxation", "redfield");
    if (rel == "redfield") m.relaxation = excitons::RelaxationTheory::Redfield;
    else if (rel == "modified_redfield") m.relaxation = excitons::RelaxationTheory::ModifiedRedfield;
    else throw std::invalid_argument("exciton.relaxation: redfield | modified_redfield");
    m.n_medium = e.num("n_medium", 1.0);
    Obj c(e.obj("coupling"), "exciton.coupling");
    const std::string model = c.str("model", "point_dipole");
    if (model == "point_dipole") m.coupling.model = excitons::CouplingModel::PointDipole;
    else if (model == "extended_dipole") m.coupling.model = excitons::CouplingModel::ExtendedDipole;
    else if (model == "legacy_J0") m.coupling.model = excitons::CouplingModel::LegacyJ0;
    else if (model == "transition_charges") m.coupling.model = excitons::CouplingModel::TransitionCharges;
    else throw std::invalid_argument("exciton.coupling.model: point_dipole | extended_dipole | legacy_J0 | transition_charges");
    const std::string qf = c.str("charges_file", "");
    if (!qf.empty()) {
      const QString q = QString::fromStdString(qf);
      m.coupling.charges = excitons::load_transition_charges(
          (QFileInfo(q).isAbsolute() ? QFileInfo(q) : QFileInfo(QDir(QString::fromStdString(rc.config_dir)), q))
              .absoluteFilePath()
              .toStdString());
    }
    if (m.coupling.model == excitons::CouplingModel::TransitionCharges && m.coupling.charges.empty())
      throw std::invalid_argument("exciton.coupling: transition_charges needs charges_file (x_nm y_nm z_nm q)");
    if (use_coupling_charges) {
      if (m.coupling.charges.empty())
        throw std::invalid_argument("exciton.dispersion.use_coupling_charges: no exciton.coupling.charges_file");
      m.dispersion.charges = m.coupling.charges;
    }
    m.coupling.screening = c.num("screening", 1.0 / (m.n_medium * m.n_medium));
    m.coupling.dipole_length_nm = c.num("dipole_length_nm", 0.8);
    m.coupling.J0_cm = c.num("J0_cm", 900.0);
    m.coupling.r0_nm = c.num("r0_nm", 1.0);
  }
  // --- баня
  {
    Obj b(root.obj("bath"), "bath");
    auto& bath = rc.model.bath;
    const std::string t = b.str("type", "drude_lorentz");
    if (t == "drude_lorentz") bath.type = excitons::BathType::DrudeLorentz;
    else if (t == "ohmic_exponential") bath.type = excitons::BathType::OhmicExponential;
    else throw std::invalid_argument("bath.type: drude_lorentz | ohmic_exponential");
    bath.lambda_cm = b.num("lambda_cm", 50.0);
    bath.cutoff_cm = b.num("cutoff_cm", 100.0);
    bath.temperature_K = b.num("temperature_K", 300.0);
    for (const auto& v : b.arr("modes")) {
      Obj mo(v.toObject(), "bath.modes[]");
      excitons::VibrationalMode vm;
      vm.omega_cm = mo.num("omega_cm", 0.0);
      vm.huang_rhys = mo.num("huang_rhys", 0.0);
      vm.gamma_cm = mo.num("gamma_cm", 10.0);
      if (vm.omega_cm <= 0.0) throw std::invalid_argument("bath.modes[].omega_cm must be > 0");
      bath.modes.push_back(vm);
    }
  }
  // --- распад и накачка
  {
    Obj d(root.obj("decay"), "decay");
    rc.model.decay.monomer_lifetime_ns = d.num("monomer_lifetime_ns", 0.0);
    rc.model.decay.n_medium = d.num("n_medium", rc.model.n_medium);
    rc.model.decay.nonradiative_ps = d.num("nonradiative_ps", 0.0);
    Obj p(root.obj("pump"), "pump");
    const std::string t = p.str("type", "oscillator_strength");
    if (t == "oscillator_strength") rc.model.pump = excitons::PumpType::OscillatorStrength;
    else if (t == "laser") rc.model.pump = excitons::PumpType::Laser;
    else throw std::invalid_argument("pump.type: oscillator_strength | laser");
    rc.model.excitation_cm = p.num("excitation_cm", rc.model.monomer_energy_cm);
  }
  // --- спектры и ансамбль
  {
    Obj s(root.obj("spectrum"), "spectrum");
    rc.grid.wmin_cm = s.num("wmin_cm", 15500.0);
    rc.grid.wmax_cm = s.num("wmax_cm", 18000.0);
    rc.grid.nw = static_cast<std::size_t>(s.integer("nw", 1800));
    rc.grid.t_max_ps = s.num("t_max_ps", 100.0);
    Obj e(root.obj("ensemble"), "ensemble");
    rc.ensemble.min_realizations = e.integer("min_realizations", 20);
    rc.ensemble.max_realizations = e.integer("max_realizations", 400);
    rc.ensemble.rel_tolerance = e.num("rel_tolerance", 0.01);
    rc.ensemble.seed = static_cast<std::uint64_t>(e.num("seed", 1234));
    rc.ensemble.batch = e.integer("batch", 16);
  }
  // --- кинетика
  {
    Obj d(root.obj("dynamics"), "dynamics");
    rc.dynamics.enabled = d.flag("enabled", true);
    rc.dynamics.t_max_ps = d.num("t_max_ps", 20.0);
    rc.dynamics.n_times = d.integer("n_times", 200);
    rc.dynamics.initial = d.str("initial", "brightest");
  }
  // --- наночастица
  {
    Obj n(root.obj("nanoparticle"), "nanoparticle");
    auto& np = rc.nanoparticle;
    np.enabled = n.flag("enabled", false);
    const std::string shape = n.str("shape", "sphere");
    const double radius = n.num("radius_nm", 15.0);
    np.spec.a_xy_nm = n.num("semi_axis_xy_nm", radius);
    np.spec.a_z_nm = n.num("semi_axis_z_nm", shape == "sphere" ? np.spec.a_xy_nm : radius);
    if (shape == "sphere") np.spec.a_z_nm = np.spec.a_xy_nm;
    else if (shape == "prolate" && np.spec.a_z_nm <= np.spec.a_xy_nm)
      throw std::invalid_argument("nanoparticle: a rod (prolate) needs semi_axis_z_nm > semi_axis_xy_nm; the given axes (z = " + std::to_string(np.spec.a_z_nm) + ", xy = " + std::to_string(np.spec.a_xy_nm) + " nm) describe a disk — choose 'oblate' or swap the axes");
    else if (shape == "oblate" && np.spec.a_z_nm >= np.spec.a_xy_nm)
      throw std::invalid_argument("nanoparticle: a disk (oblate) needs semi_axis_z_nm < semi_axis_xy_nm; the given axes (z = " + std::to_string(np.spec.a_z_nm) + ", xy = " + std::to_string(np.spec.a_xy_nm) + " nm) describe a rod — choose 'prolate' or swap the axes");
    else if (shape != "sphere" && shape != "prolate" && shape != "oblate")
      throw std::invalid_argument("nanoparticle.shape: sphere | prolate | oblate");
    np.spec.l_max = n.integer("l_max", 40);
    np.spec.n_medium = n.num("n_medium", 1.33);
    np.spec.dynamic_dipole = n.flag("dynamic_dipole", true);
    np.spec.shell_thickness_nm = n.num("shell_thickness_nm", 0.0);
    np.spec.shell_n = n.num("shell_n", 1.45);
    if (n.has("d_perp_nm")) {
      const QJsonValue v = n.raw("d_perp_nm");
      if (v.isArray() && v.toArray().size() == 2)
        np.spec.d_perp_nm = {v.toArray()[0].toDouble(), v.toArray()[1].toDouble()};
      else if (v.isDouble())
        np.spec.d_perp_nm = v.toDouble();
      else
        throw std::invalid_argument("nanoparticle.d_perp_nm: number or [re, im] (nm)");
    }
    np.spec.d_perp_nm += std::complex<double>(0.0, n.num("d_perp_im_nm", 0.0));  // мнимая часть отдельным ключом (GUI)
    if (shape != "sphere" && (np.spec.shell_thickness_nm > 0.0 || np.spec.d_perp_nm != std::complex<double>(0.0)))
      throw std::invalid_argument("nanoparticle: shell_thickness_nm and d_perp_nm are implemented for spheres only");
    {
      const QJsonValue v = n.raw("fast_l");
      if (v.isString() && v.toString() == "auto") np.spec.fast_l = -1;
      else if (v.isUndefined()) np.spec.fast_l = -1;
      else if (v.isDouble()) np.spec.fast_l = v.toInt();
      else throw std::invalid_argument("nanoparticle.fast_l: \"auto\" or integer (0 — exact solver)");
    }
    np.spec.size_sigma_rel = n.num("size_sigma_rel", 0.0);
    np.spec.size_nodes = n.integer("size_nodes", 5);
    np.emission = n.flag("emission", true);
    const std::string mat = n.str("material_file", "");
    if (!mat.empty()) {
      const QString m = QString::fromStdString(mat);
      np.material_file = (QFileInfo(m).isAbsolute() ? QFileInfo(m)
                                                    : QFileInfo(QDir(QString::fromStdString(rc.config_dir)), m))
                             .absoluteFilePath()
                             .toStdString();
    }
    if (n.has("drude")) {
      Obj d(n.obj("drude"), "nanoparticle.drude");
      np.use_drude = true;
      np.drude.eps_inf = d.num("eps_inf", 1.0);
      np.drude.omega_p_eV = d.num("omega_p_eV", 9.0);
      np.drude.gamma_eV = d.num("gamma_eV", 0.07);
    }
    np.kreibig_A = n.num("kreibig_A", 0.0);
    np.fermi_velocity_m_s = n.num("fermi_velocity_m_s", 1.40e6);
    np.kreibig_omega_p_eV = n.num("kreibig_omega_p_eV", 9.0);
    np.kreibig_gamma_bulk_eV = n.num("kreibig_gamma_bulk_eV", 0.07);
    np.realizations = n.integer("realizations", 30);
    np.wmin_cm = n.num("wmin_cm", rc.grid.wmin_cm);
    np.wmax_cm = n.num("wmax_cm", rc.grid.wmax_cm);
    np.nw = static_cast<std::size_t>(n.integer("nw", 600));
    np.homogeneous_hwhm_cm = n.num("homogeneous_hwhm_cm", -1.0);
    np.near_field_cm = n.num("near_field_cm", rc.model.monomer_energy_cm);
    if (np.enabled && np.material_file.empty() && !np.use_drude)
      throw std::invalid_argument("nanoparticle: set material_file or drude");
  }
  // --- скан (разбирается отдельно в scan.cpp; здесь только помечаем раздел как известный)
  (void)root.obj("scan");
  // --- слой на поверхности НЧ: полуоси берутся из раздела nanoparticle, если не заданы явно
  for (auto& c : rc.configurations) {
    if (c.type != "particle_shell") continue;
    // внешняя поверхность — с учётом лигандной оболочки
    if (!c.params.contains("a_xy")) c.params["a_xy"] = rc.nanoparticle.spec.outer_xy();
    if (!c.params.contains("a_z")) c.params["a_z"] = rc.nanoparticle.spec.outer_z();
  }
  // --- вибронная модель Хольстейна
  {
    Obj h(root.obj("holstein"), "holstein");
    auto& hc = rc.holstein;
    hc.enabled = h.flag("enabled", false);
    hc.params.omega_cm = h.num("omega_cm", 1300.0);
    hc.params.huang_rhys = h.num("huang_rhys", 0.6);
    hc.params.nu_max = h.integer("nu_max", 4);
    hc.realizations = h.integer("realizations", 20);
  }
  // --- HEOM
  {
    Obj h(root.obj("heom"), "heom");
    auto& hc = rc.heom;
    hc.enabled = h.flag("enabled", false);
    hc.params.matsubara = h.integer("matsubara", 1);
    hc.params.depth = h.integer("depth", 4);
    hc.params.terminator = h.flag("terminator", true);
    hc.params.dt_fs = h.num("dt_fs", 0.0);
    hc.params.t_max_ps = h.num("t_max_ps", 1.0);
    hc.max_sites = h.integer("max_sites", 10);
    hc.absorption = h.flag("absorption", true);
    hc.populations = h.flag("populations", true);
    hc.initial = h.str("initial", "brightest");
    hc.n_out = h.integer("n_out", 200);
    hc.params.lambda_cm = rc.model.bath.lambda_cm;
    hc.params.gamma_cm = rc.model.bath.cutoff_cm;
    hc.params.temperature_K = rc.model.bath.temperature_K;
  }
  return rc;
}

}  // namespace ssp::app
