#include "calibrate.hpp"

#include "output.hpp"

#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace ssp::app {

namespace {

std::string f6(double v, int p = 6) {
  std::ostringstream os;
  os.precision(p);
  os << v;
  return os.str();
}

std::string resolve(const std::string& base, const std::string& p) {
  const QString q = QString::fromStdString(p);
  return (QFileInfo(q).isAbsolute() ? QFileInfo(q) : QFileInfo(QDir(QString::fromStdString(base)), q))
      .absoluteFilePath()
      .toStdString();
}

// Таблица x, y[, σ]: разделители , ; табуляция пробел; строки без чисел (заголовки, #) пропускаются.
excitons::Spectrum read_table(const std::string& path, const std::string& unit, bool per_wavelength_density) {
  std::ifstream f(path);
  if (!f) throw std::invalid_argument("cannot open " + path);
  excitons::Spectrum s;
  std::string line;
  while (std::getline(f, line)) {
    for (char& c : line)
      if (c == ',' || c == ';' || c == '\t') c = ' ';
    std::istringstream is(line);
    double x, y, e;
    if (!(is >> x >> y)) continue;
    double nu = x;
    if (unit == "nm") nu = 1e7 / x;
    else if (unit == "eV") nu = x * 8065.544;
    else if (unit != "cm-1") throw std::invalid_argument("unit: cm-1 | nm | eV");
    // ФЛ на интервал длины волны → на интервал волнового числа: F_ν̃ = F_λ·λ²
    const double jac = (unit == "nm" && per_wavelength_density) ? x * x * 1e-7 : 1.0;
    s.w.push_back(nu);
    s.y.push_back(y * jac);
    if (is >> e) s.sigma.push_back(e * jac);
  }
  if (s.w.size() < 5) throw std::invalid_argument("too few data points in " + path);
  std::vector<std::size_t> idx(s.w.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return s.w[a] < s.w[b]; });
  excitons::Spectrum o;
  for (auto i : idx) {
    o.w.push_back(s.w[i]);
    o.y.push_back(s.y[i]);
    if (s.sigma.size() == s.w.size()) o.sigma.push_back(s.sigma[i]);
  }
  return o;
}

excitons::Spectrum crop(const excitons::Spectrum& s, double lo, double hi) {
  if (lo <= 0 && hi <= 0) return s;
  excitons::Spectrum o;
  for (std::size_t i = 0; i < s.w.size(); ++i)
    if ((lo <= 0 || s.w[i] >= lo) && (hi <= 0 || s.w[i] <= hi)) {
      o.w.push_back(s.w[i]);
      o.y.push_back(s.y[i]);
      if (s.sigma.size() == s.w.size()) o.sigma.push_back(s.sigma[i]);
    }
  return o;
}

}  // namespace

std::string config_dir_of(const std::string& config_path) {
  return QFileInfo(QString::fromStdString(config_path)).absolutePath().toStdString();
}

bool is_calibration_config(const std::string& json) {
  const auto doc = QJsonDocument::fromJson(QByteArray::fromStdString(json));
  return doc.isObject() && doc.object().contains("monomer_fit");
}

CalibrationOutput run_calibration(const std::string& json, const std::string& base_dir, const Log& log) {
  QJsonParseError err;
  const auto doc = QJsonDocument::fromJson(QByteArray::fromStdString(json), &err);
  if (!doc.isObject()) throw std::invalid_argument("JSON: " + err.errorString().toStdString());
  const QJsonObject root = doc.object();
  const QJsonObject c = root["monomer_fit"].toObject();
  for (const auto& k : c.keys())
    if (!QStringList{"absorption_file", "absorption_unit", "absolute_extinction", "emission_file", "emission_unit",
                     "emission_per_wavelength", "n_medium", "temperature_K", "wmin_cm", "wmax_cm", "start", "fit", "stokes_peak_cm",
                     "stokes_sigma_cm",
                     "output_dir"}
             .contains(k) &&
        !k.startsWith("_"))
      throw std::invalid_argument("monomer_fit: unknown key '" + k.toStdString() + "'");
  CalibrationOutput out;
  const double n_med = c["n_medium"].toDouble(1.33);
  out.abs = crop(read_table(resolve(base_dir, c["absorption_file"].toString().toStdString()),
                            c["absorption_unit"].toString("nm").toStdString(), false),
                 c["wmin_cm"].toDouble(0), c["wmax_cm"].toDouble(0));
  if (c.contains("emission_file") && !c["emission_file"].toString().isEmpty())
    out.em = crop(read_table(resolve(base_dir, c["emission_file"].toString().toStdString()),
                             c["emission_unit"].toString("nm").toStdString(), c["emission_per_wavelength"].toBool(true)),
                  c["wmin_cm"].toDouble(0), c["wmax_cm"].toDouble(0));
  excitons::MonomerFitParams p;
  p.temperature_K = c["temperature_K"].toDouble(300.0);
  const QJsonObject s = c["start"].toObject();
  p.E0_cm = s["E0_cm"].toDouble(18000.0);
  p.lambda_cm = s["lambda_cm"].toDouble(100.0);
  p.cutoff_cm = s["cutoff_cm"].toDouble(100.0);
  p.sigma_inh_cm = s["sigma_inh_cm"].toDouble(150.0);
  for (const auto& v : s["modes"].toArray()) {
    const QJsonObject m = v.toObject();
    excitons::VibrationalMode vm;
    vm.omega_cm = m["omega_cm"].toDouble(1200.0);
    vm.huang_rhys = m["huang_rhys"].toDouble(0.5);
    vm.gamma_cm = m["gamma_cm"].toDouble(20.0);
    p.modes.push_back(vm);
  }
  const QJsonObject fl = c["fit"].toObject();
  p.fit_E0 = fl["E0"].toBool(true);
  p.fit_lambda = fl["lambda"].toBool(true);
  p.fit_cutoff = fl["cutoff"].toBool(true);
  p.fit_sigma = fl["sigma_inh"].toBool(true);
  p.fit_mode_S = fl["mode_S"].toBool(true);
  p.fit_mode_omega = fl["mode_omega"].toBool(false);
  p.moment_start = fl["moment_start"].toBool(true);
  p.stokes_peak_cm = c["stokes_peak_cm"].toDouble(0.0);
  p.stokes_sigma_cm = c["stokes_sigma_cm"].toDouble(0.0);
  const bool absolute = c["absolute_extinction"].toBool(false);

  log("monomer calibration: " + std::to_string(out.abs.w.size()) + " absorption points, " +
      std::to_string(out.em.w.size()) + " emission points");
  out.fit = excitons::fit_monomer(out.abs, out.em, p, absolute, n_med);
  const auto& F = out.fit;
  log("fit " + std::string(F.converged ? "converged" : "NOT converged") + ", χ²_ν = " + f6(F.chi2_red, 4) +
      (F.chi2_red > 3 ? "  (model misses features of the data — see residuals)" : ""));
  for (std::size_t k = 0; k < F.names.size(); ++k)
    log("  " + F.names[k] + " = " + f6(F.values[k], 6) + " ± " + f6(F.errors[k], 3));
  double lt = F.p.lambda_cm;
  for (const auto& m : F.p.modes) lt += m.huang_rhys * m.omega_cm;
  log("  total reorganization λ_tot = " + f6(lt, 5) + " cm⁻¹ (Stokes shift of centroids 2λ_tot); peak Stokes shift " +
      f6(F.stokes_peak_cm, 5) + " cm⁻¹" +
      (p.stokes_peak_cm > 0 ? " (target " + f6(p.stokes_peak_cm, 5) + " ± " + f6(p.stokes_sigma_cm, 3) + ")" : ""));
  if (absolute)
    log("  μ = " + f6(F.mu_D, 4) + " D (∫ε/ν̃, n = " + f6(n_med, 3) + ", no local-field factor); τ_rad = " +
        f6(F.tau_rad_ns, 4) + " ns (model formula), Strickler–Berg " + f6(F.tau_sb_ns, 4) + " ns");

  out.output_dir = resolve(base_dir, c["output_dir"].toString("results/monomer_fit").toStdString());
  make_dirs(out.output_dir);
  {
    CsvWriter w(join_path(out.output_dir, "fit_parameters.csv"),
                {"monomer calibration: cumulant model (Drude + Brownian modes) + static Gaussian broadening",
                 "errors from (JᵀJ)⁻¹·χ²_ν; chi2_red = " + f6(F.chi2_red, 4)},
                {"name", "value", "error"});
    for (std::size_t k = 0; k < F.names.size(); ++k) w.row_text({F.names[k], f6(F.values[k], 8), f6(F.errors[k], 4)});
    w.row_text({"stokes_peak_cm (model)", f6(F.stokes_peak_cm, 6), ""});
    if (absolute) {
      w.row_text({"mu_D", f6(F.mu_D, 6), ""});
      w.row_text({"tau_rad_ns", f6(F.tau_rad_ns, 6), ""});
      w.row_text({"tau_strickler_berg_ns", f6(F.tau_sb_ns, 6), ""});
    }
  }
  auto write_spec = [&](const std::string& name, const excitons::Spectrum& s, const std::vector<double>& m) {
    CsvWriter w(join_path(out.output_dir, name), {"data and fitted model (same scale); resid = (data − model)/σ"},
                {"w_cm", "data", "model", "resid"});
    const double mx = s.y.empty() ? 1.0 : *std::max_element(s.y.begin(), s.y.end());
    for (std::size_t i = 0; i < s.w.size(); ++i) {
      const double sg = s.sigma.size() == s.w.size() ? s.sigma[i] : 0.02 * mx;
      w.row({s.w[i], s.y[i], m[i], (s.y[i] - m[i]) / sg});
    }
  };
  write_spec("fit_absorption.csv", out.abs, F.abs_model);
  if (!out.em.w.empty()) write_spec("fit_emission.csv", out.em, F.em_model);

  // Фрагмент конфига агрегата: параметры мономера → модель агрегата (предсказание без подгонки)
  QJsonObject ex, bath, cp;
  ex["monomer_energy_cm"] = F.p.E0_cm;
  if (absolute) ex["transition_dipole_D"] = F.mu_D;
  ex["n_medium"] = n_med;
  ex["_disorder_hint"] = QString("monomer static broadening σ_inh = %1 cm⁻¹ — upper bound for the site disorder")
                             .arg(F.p.sigma_inh_cm, 0, 'f', 1);
  bath["type"] = "drude_lorentz";
  bath["lambda_cm"] = F.p.lambda_cm;
  bath["cutoff_cm"] = F.p.cutoff_cm;
  bath["temperature_K"] = F.p.temperature_K;
  QJsonArray modes;
  for (const auto& m : F.p.modes)
    modes.append(QJsonObject{{"omega_cm", m.omega_cm}, {"huang_rhys", m.huang_rhys}, {"gamma_cm", m.gamma_cm}});
  bath["modes"] = modes;
  cp["exciton"] = ex;
  cp["bath"] = bath;
  out.config_snippet = QJsonDocument(cp).toJson(QJsonDocument::Indented).toStdString();
  {
    std::ofstream f(join_path(out.output_dir, "aggregate_config_snippet.json"));
    f << out.config_snippet;
  }
  log("results and config snippet: " + out.output_dir);
  return out;
}

}  // namespace ssp::app
