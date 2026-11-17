#include "scan.hpp"

#include "config.hpp"
#include "output.hpp"

#include <QDir>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QRegularExpression>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ssp::app {

namespace {

std::string fmtd(double v, int p = 6) {
  std::ostringstream os;
  os.precision(p);
  os << v;
  return os.str();
}

QJsonObject parse_root(const std::string& json) {
  QJsonParseError err;
  const QJsonDocument d = QJsonDocument::fromJson(QByteArray::fromStdString(json), &err);
  if (err.error != QJsonParseError::NoError || !d.isObject())
    throw std::runtime_error("JSON error: " + err.errorString().toStdString());
  return d.object();
}

void set_in(QJsonObject& o, const QStringList& keys, int i, double v) {
  if (i + 1 == keys.size()) {
    o[keys[i]] = v;
    return;
  }
  QJsonObject child = o[keys[i]].toObject();
  set_in(child, keys, i + 1, v);
  o[keys[i]] = child;
}

// Установка значения по пути. Поддерживаются "a.b.c", "geometry_params.X" (все формы),
// "configurations[name].a.b" и "configurations[3].a.b".
void set_path(QJsonObject& root, const std::string& path, double v) {
  const QString p = QString::fromStdString(path);
  static const QRegularExpression rx("^configurations\\[([^\\]]+)\\]\\.(.+)$");
  const auto m = rx.match(p);
  if (p.startsWith("geometry_params.") || p.startsWith("placement.") || m.hasMatch()) {
    QJsonArray arr = root["configurations"].toArray();
    bool hit = false;
    for (int i = 0; i < arr.size(); ++i) {
      QJsonObject c = arr[i].toObject();
      QString sub = p;
      if (m.hasMatch()) {
        bool isNum = false;
        const int idx = m.captured(1).toInt(&isNum);
        if (isNum ? idx != i : c["name"].toString() != m.captured(1)) continue;
        sub = m.captured(2);
      }
      set_in(c, sub.split('.'), 0, v);
      arr[i] = c;
      hit = true;
    }
    if (!hit) throw std::invalid_argument("scan: no configuration matches '" + path + "'");
    root["configurations"] = arr;
    return;
  }
  set_in(root, p.split('.'), 0, v);
}

void peak_fwhm(const std::vector<double>& w, const std::vector<double>& y, double& pk, double& fw) {
  if (y.empty()) {
    pk = fw = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  const std::size_t k = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
  pk = w[k];
  const double half = 0.5 * y[k];
  std::size_t l = k, r = k;
  while (l > 0 && y[l] > half) --l;
  while (r + 1 < y.size() && y[r] > half) ++r;
  auto cross = [&](std::size_t i, std::size_t j) { return w[i] + (half - y[i]) / (y[j] - y[i]) * (w[j] - w[i]); };
  fw = ((y[r] <= half && r > k) ? cross(r - 1, r) : w[r]) - ((y[l] <= half && l < k) ? cross(l, l + 1) : w[l]);
}

// Локальные максимумы выше rel·max, не более n штук (наибольшие), по возрастанию ω.
std::vector<double> local_peaks(const std::vector<double>& w, const std::vector<double>& y, double rel, int n) {
  std::vector<std::pair<double, double>> p;
  if (y.size() < 3) return {};
  const double ymax = *std::max_element(y.begin(), y.end());
  for (std::size_t i = 1; i + 1 < y.size(); ++i)
    if (y[i] > y[i - 1] && y[i] >= y[i + 1] && y[i] > rel * ymax) {
      // параболическое уточнение
      const double a = y[i - 1], b = y[i], c = y[i + 1], den = a - 2 * b + c;
      const double x = den < 0 ? w[i] + 0.5 * (a - c) / den * (w[i + 1] - w[i]) : w[i];
      p.emplace_back(y[i], x);
    }
  std::sort(p.begin(), p.end(), [](auto& a, auto& b) { return a.first > b.first; });
  std::vector<double> out;
  for (int k = 0; k < std::min<int>(n, static_cast<int>(p.size())); ++k) out.push_back(p[k].second);
  std::sort(out.begin(), out.end());
  return out;
}

void write_map(const std::string& path, const std::string& title, const std::vector<double>& w,
               const std::vector<double>& values, const std::vector<std::vector<double>>& z) {
  std::vector<std::string> cols = {"w_cm"};
  for (double v : values) cols.push_back("v=" + fmtd(v));
  CsvWriter wr(path, {title, "rows: wavenumber; columns: scan values"}, cols);
  for (std::size_t q = 0; q < w.size(); ++q) {
    std::vector<double> row = {w[q]};
    for (const auto& zz : z) row.push_back(zz[q]);
    wr.row(row);
  }
}

}  // namespace

void coupled_branches(double wp, double gp, double wx, double gx, double g, double& lower, double& upper) {
  using cd = std::complex<double>;
  const cd mean = 0.5 * (wp + wx) - cd(0.0, 0.25 * (gp + gx));
  const cd half = 0.5 * (wp - wx) - cd(0.0, 0.25 * (gp - gx));
  const cd root = std::sqrt(g * g + half * half);
  const double a = (mean + root).real(), b = (mean - root).real();
  lower = std::min(a, b);
  upper = std::max(a, b);
}

BranchFit fit_coupled_oscillators(const std::vector<double>& wp, const std::vector<double>& gp,
                                  const std::vector<double>& wx, const std::vector<double>& gx,
                                  const std::vector<std::vector<double>>& peaks, const std::string& name) {
  BranchFit f;
  f.spectrum = name;
  const std::size_t n = wp.size();
  auto cost = [&](double g, int* np) {
    double s = 0.0;
    int cnt = 0;
    for (std::size_t k = 0; k < n; ++k) {
      double lo, up;
      coupled_branches(wp[k], gp[k], wx[k], gx[k], g, lo, up);
      for (double pk : peaks[k]) {
        // Пики далеко от обеих ветвей (другие моды) не учитываются
        const double d = std::min(std::abs(pk - lo), std::abs(pk - up));
        if (d > 0.5 * (gp[k] + gx[k]) + 300.0) continue;
        s += d * d;
        ++cnt;
      }
    }
    if (np) *np = cnt;
    return cnt > 0 ? s / cnt : std::numeric_limits<double>::infinity();
  };
  double gmax = 50.0;
  for (std::size_t k = 0; k < n; ++k) gmax = std::max(gmax, std::abs(wp[k] - wx[k]) + gp[k] + gx[k]);
  double best = std::numeric_limits<double>::infinity(), gb = 0.0;
  const int G = 2000;
  for (int i = 0; i <= G; ++i) {
    const double g = gmax * i / G;
    const double c = cost(g, nullptr);
    if (c < best) {
      best = c;
      gb = g;
    }
  }
  double a = std::max(0.0, gb - gmax / G), b = gb + gmax / G;
  const double r = 0.5 * (std::sqrt(5.0) - 1.0);
  for (int it = 0; it < 80; ++it) {
    const double c1 = b - r * (b - a), c2 = a + r * (b - a);
    if (cost(c1, nullptr) < cost(c2, nullptr)) b = c2;
    else a = c1;
  }
  f.g_cm = 0.5 * (a + b);
  f.rms_cm = std::sqrt(cost(f.g_cm, &f.n_peaks));
  // Интервал: n·(cost(g) − cost_min)/cost_min ≤ 1
  {
    const double cmin = std::max(cost(f.g_cm, nullptr), 1e-12);
    const int n = std::max(f.n_peaks, 1);
    f.g_lo = f.g_cm;
    f.g_hi = f.g_cm;
    for (int i = 0; i <= G; ++i) {
      const double g = gmax * i / G;
      if (n * (cost(g, nullptr) - cmin) / cmin <= 1.0) {
        f.g_lo = std::min(f.g_lo, g);
        f.g_hi = std::max(f.g_hi, g);
      }
    }
  }
  // Разрешённость: в скольких точках есть пики у обеих ветвей
  for (std::size_t k = 0; k < n; ++k) {
    double lo, up;
    coupled_branches(wp[k], gp[k], wx[k], gx[k], f.g_cm, lo, up);
    const double tol = 0.25 * std::abs(up - lo) + 1.0;
    bool hl = false, hu = false;
    for (double pk : peaks[k]) {
      hl = hl || std::abs(pk - lo) < tol;
      hu = hu || std::abs(pk - up) < tol;
    }
    if (hl && hu && peaks[k].size() >= 2) ++f.two_peak_points;
  }
  f.resolved = f.two_peak_points >= 2;
  // Параметры в точке наименьшей расстройки
  std::size_t k0 = 0;
  for (std::size_t k = 1; k < n; ++k)
    if (std::abs(wp[k] - wx[k]) < std::abs(wp[k0] - wx[k0])) k0 = k;
  f.gamma_p = gp[k0];
  f.gamma_x = gx[k0];
  f.detuning_min = wp[k0] - wx[k0];
  const double c = 0.5 * (wp[k0] + wx[k0]);
  double below = -1e300, above = 1e300;
  for (double pk : peaks[k0]) {
    if (pk <= c) below = std::max(below, pk);
    else above = std::min(above, pk);
  }
  f.observed_splitting = (below > -1e299 && above < 1e299) ? above - below : 0.0;
  f.visible = 2.0 * f.g_cm > 0.5 * (f.gamma_p + f.gamma_x);
  f.anticrossing = f.g_cm > 0.25 * std::abs(f.gamma_p - f.gamma_x);
  f.cooperativity = 4.0 * f.g_cm * f.g_cm / std::max(f.gamma_p * f.gamma_x, 1e-300);
  return f;
}

ScanSpec parse_scan(const std::string& config_json) {
  ScanSpec s;
  const QJsonObject root = parse_root(config_json);
  if (!root.contains("scan")) return s;
  const QJsonObject o = root["scan"].toObject();
  for (const auto& k : o.keys())
    if (!QStringList{"enabled", "parameter", "values", "from", "to", "steps", "also_set"}.contains(k) && !k.startsWith("_"))
      throw std::invalid_argument("scan: unknown key '" + k.toStdString() + "'");
  s.enabled = o["enabled"].toBool(false);
  s.parameter = o["parameter"].toString().toStdString();
  for (const auto& v : o["values"].toArray()) s.values.push_back(v.toDouble());
  if (s.values.empty() && o.contains("from")) {
    const double a = o["from"].toDouble(), b = o["to"].toDouble();
    const int n = std::max(2, o["steps"].toInt(11));
    for (int i = 0; i < n; ++i) s.values.push_back(a + (b - a) * i / (n - 1));
  }
  if (o["also_set"].isString()) {
    for (const auto& t : o["also_set"].toString().split(',', Qt::SkipEmptyParts))
      if (!t.trimmed().isEmpty()) s.also_set.push_back(t.trimmed().toStdString());
  } else {
    for (const auto& v : o["also_set"].toArray()) s.also_set.push_back(v.toString().toStdString());
  }
  if (s.enabled && s.parameter.empty()) throw std::invalid_argument("scan: 'parameter' is empty");
  if (s.enabled && s.values.size() < 2) throw std::invalid_argument("scan: at least two values are needed");
  return s;
}

ScanResult run_scan(const std::string& config_json, const std::string& base_dir, const std::string& config_path,
                    const Log& log, const std::atomic<bool>* cancel, const Progress& progress) {
  const ScanSpec spec = parse_scan(config_json);
  if (!spec.enabled) throw std::invalid_argument("scan is not enabled");
  const RunConfig rc0 = load_config_text(config_json, base_dir, config_path);
  const std::string out = rc0.output_dir;
  make_dirs(out);
  {
    std::FILE* f = std::fopen(join_path(out, "config_used.json").c_str(), "w");
    if (f) {
      std::fputs(config_json.c_str(), f);
      std::fclose(f);
    }
  }
  const QJsonObject root0 = parse_root(config_json);
  ScanResult R;
  R.values = spec.values;
  const int nv = static_cast<int>(spec.values.size());
  {
    CsvWriter w(join_path(out, "scan_info.csv"), {"scan over " + spec.parameter}, {"index", "value", "folder"});
    for (int k = 0; k < nv; ++k) {
      char d[32];
      std::snprintf(d, sizeof d, "scan_%02d", k);
      w.row_text({std::to_string(k), fmtd(spec.values[k], 10), d});
    }
  }
  log("Scan over " + spec.parameter + ": " + std::to_string(nv) + " values");
  for (int k = 0; k < nv; ++k) {
    if (cancel && cancel->load()) {
      log("Scan stopped by user.");
      break;
    }
    QJsonObject root = root0;
    set_path(root, spec.parameter, spec.values[k]);
    for (const auto& p : spec.also_set) set_path(root, p, spec.values[k]);
    QJsonObject sc = root["scan"].toObject();
    sc["enabled"] = false;
    root["scan"] = sc;
    char d[32];
    std::snprintf(d, sizeof d, "scan_%02d", k);
    root["output_dir"] = QString::fromStdString(join_path(out, d));
    const std::string text = QJsonDocument(root).toJson(QJsonDocument::Indented).toStdString();
    const RunConfig rc = load_config_text(text, base_dir, config_path);
    log("\n##### scan " + std::to_string(k + 1) + "/" + std::to_string(nv) + ": " + spec.parameter + " = " +
        fmtd(spec.values[k]));
    const RunResult res = run(rc, text, log, cancel, [&](int i, int n, const std::string& name) {
      if (progress) progress(k * n + i, nv * n, d + std::string(": ") + name);
    });
    for (const auto& s : res.shapes) {
      auto it = std::find_if(R.shapes.begin(), R.shapes.end(), [&](const ScanShape& x) { return x.name == s.name; });
      if (it == R.shapes.end()) {
        R.shapes.push_back({});
        it = R.shapes.end() - 1;
        it->name = s.name;
        it->w_agg = s.w;
        it->has_np = s.has_np;
        it->w_np = s.np_w;
      }
      it->abs.push_back(s.abs);
      it->em.push_back(s.em);
      it->observables.push_back(s.observables);
      if (s.has_np) {
        it->ext.push_back(s.ext);
        it->abs_sys.push_back(s.abs_sys);
        it->sca.push_back(s.sca);
        it->np_abs.push_back(s.np_abs);
        double pk, fw;
        it->w_p.push_back(s.mode_omega_cm);      // полюс дипольной моды, а не FWHM поглощения НЧ
        it->gamma_p.push_back(s.mode_gamma_cm);
        it->g_micro.push_back(s.g_micro_cm);
        peak_fwhm(s.w, s.abs, pk, fw);
        it->w_x.push_back(pk);
        it->gamma_x.push_back(fw);
        it->peaks_abs.push_back(local_peaks(s.np_w, s.abs_sys, 0.05, 3));
        it->peaks_ext.push_back(local_peaks(s.np_w, s.ext, 0.05, 3));
        it->peaks_sca.push_back(local_peaks(s.np_w, s.sca, 0.05, 3));
      }
    }
  }

  // ---- сборка
  const std::size_t done = R.shapes.empty() ? 0 : R.shapes.front().abs.size();
  std::vector<double> vals(spec.values.begin(), spec.values.begin() + static_cast<long>(done));
  for (auto& sh : R.shapes) {
    const std::string dir = join_path(out, sh.name);
    make_dirs(dir);
    write_map(join_path(dir, "scan_abs.csv"), "aggregate absorption σ_abs (nm²) vs scan value", sh.w_agg, vals, sh.abs);
    write_map(join_path(dir, "scan_em.csv"), "aggregate fluorescence vs scan value", sh.w_agg, vals, sh.em);
    {
      std::vector<std::string> cols = {"value"};
      std::vector<std::string> keys;
      if (!sh.observables.empty())
        for (const auto& [k, v] : sh.observables[0]) {
          keys.push_back(k);
          cols.push_back(k);
        }
      CsvWriter w(join_path(dir, "scan_observables.csv"), {"observables vs " + spec.parameter}, cols);
      for (std::size_t i = 0; i < sh.observables.size(); ++i) {
        std::vector<double> row = {vals[i]};
        for (const auto& k : keys) row.push_back(sh.observables[i].count(k) ? sh.observables[i].at(k) : NAN);
        w.row(row);
      }
    }
    if (!sh.has_np || sh.ext.size() < 2) continue;
    write_map(join_path(dir, "scan_ext.csv"), "system extinction (nm²)", sh.w_np, vals, sh.ext);
    write_map(join_path(dir, "scan_abs_sys.csv"), "system absorption (nm²)", sh.w_np, vals, sh.abs_sys);
    write_map(join_path(dir, "scan_sca.csv"), "system scattering (nm²)", sh.w_np, vals, sh.sca);
    std::vector<std::vector<double>> diff(sh.abs_sys.size());
    for (std::size_t i = 0; i < diff.size(); ++i) {
      diff[i].resize(sh.w_np.size());
      for (std::size_t q = 0; q < sh.w_np.size(); ++q) diff[i][q] = sh.abs_sys[i][q] - sh.np_abs[i][q];
    }
    write_map(join_path(dir, "scan_abs_diff.csv"), "system − NP absorption (nm²)", sh.w_np, vals, diff);

    sh.fits = {fit_coupled_oscillators(sh.w_p, sh.gamma_p, sh.w_x, sh.gamma_x, sh.peaks_abs, "abs"),
               fit_coupled_oscillators(sh.w_p, sh.gamma_p, sh.w_x, sh.gamma_x, sh.peaks_ext, "ext"),
               fit_coupled_oscillators(sh.w_p, sh.gamma_p, sh.w_x, sh.gamma_x, sh.peaks_sca, "sca")};
    {
      // Микроскопическое g (без подгонки по пикам) в точке наименьшей расстройки
      BranchFit m = sh.fits[0];
      std::size_t k0 = 0;
      for (std::size_t k = 1; k < sh.w_p.size(); ++k)
        if (std::abs(sh.w_p[k] - sh.w_x[k]) < std::abs(sh.w_p[k0] - sh.w_x[k0])) k0 = k;
      m.spectrum = "micro";
      m.g_cm = m.g_lo = m.g_hi = sh.g_micro[k0];
      m.rms_cm = NAN;
      m.n_peaks = 0;
      m.resolved = true;
      m.two_peak_points = 0;
      m.visible = 2.0 * m.g_cm > 0.5 * (m.gamma_p + m.gamma_x);
      m.anticrossing = m.g_cm > 0.25 * std::abs(m.gamma_p - m.gamma_x);
      m.cooperativity = 4.0 * m.g_cm * m.g_cm / std::max(m.gamma_p * m.gamma_x, 1e-300);
      sh.fits.push_back(m);
    }
    {
      CsvWriter w(join_path(dir, "scan_branches.csv"),
                  {"w_p, gamma_p: dipolar NP mode (pole fit); w_x, gamma_x: aggregate alone (ensemble spectrum); hybrid peaks of the system spectra;",
                   "model_lower/upper: coupled-oscillator branches with g fitted to the absorption peaks"},
                  {"value", "w_p", "gamma_p", "w_x", "gamma_x", "abs_peak1", "abs_peak2", "abs_peak3", "ext_peak1",
                   "ext_peak2", "ext_peak3", "sca_peak1", "sca_peak2", "sca_peak3", "model_lower", "model_upper"});
      for (std::size_t i = 0; i < sh.w_p.size(); ++i) {
        std::vector<double> row = {vals[i], sh.w_p[i], sh.gamma_p[i], sh.w_x[i], sh.gamma_x[i]};
        for (const auto* P : {&sh.peaks_abs, &sh.peaks_ext, &sh.peaks_sca})
          for (int j = 0; j < 3; ++j) row.push_back(j < static_cast<int>((*P)[i].size()) ? (*P)[i][j] : NAN);
        double lo, up;
        coupled_branches(sh.w_p[i], sh.gamma_p[i], sh.w_x[i], sh.gamma_x[i], sh.fits[0].g_cm, lo, up);
        row.push_back(lo);
        row.push_back(up);
        w.row(row);
      }
    }
    {
      CsvWriter w(join_path(dir, "scan_coupling.csv"),
                  {"coupled-oscillator fit (γ = FWHM): visible splitting 2g > (γp+γx)/2; anticrossing g > |γp−γx|/4;",
                   "cooperativity C = 4g²/(γp γx); parameters taken at the smallest detuning"},
                  {"spectrum", "g_cm", "g_lo_cm", "g_hi_cm", "rabi_2g_cm", "fit_rms_cm", "n_peaks", "resolved",
                   "two_peak_points", "gamma_p_cm", "gamma_x_cm", "detuning_cm", "observed_splitting_cm",
                   "visible_splitting", "anticrossing", "cooperativity"});
      for (const auto& f : sh.fits)
        w.row_text({f.spectrum, fmtd(f.g_cm), fmtd(f.g_lo), fmtd(f.g_hi), fmtd(2 * f.g_cm), fmtd(f.rms_cm, 4),
                    std::to_string(f.n_peaks), f.resolved ? "yes" : "no", std::to_string(f.two_peak_points),
                    fmtd(f.gamma_p), fmtd(f.gamma_x), fmtd(f.detuning_min), fmtd(f.observed_splitting),
                    f.visible ? "yes" : "no", f.anticrossing ? "yes" : "no", fmtd(f.cooperativity, 4)});
    }
    const BranchFit& fa = sh.fits[0];
    const BranchFit& fm = sh.fits.back();
    log("\n" + sh.name + ": microscopic coupling g = " + fmtd(fm.g_cm, 4) + " cm⁻¹ (2g = " + fmtd(2 * fm.g_cm, 4) +
        "); γp = " + fmtd(fm.gamma_p, 4) + ", γx = " + fmtd(fm.gamma_x, 4) + " cm⁻¹ → " +
        (fm.visible ? "STRONG coupling expected (2g > (γp+γx)/2)" : "not strong") + ", C = " + fmtd(fm.cooperativity, 3));
    log(std::string("  peak-based fit (absorption): ") +
        (fa.resolved ? "g = " + fmtd(fa.g_cm, 4) + " [" + fmtd(fa.g_lo, 4) + ", " + fmtd(fa.g_hi, 4) + "] cm⁻¹"
                     : "second branch NOT resolved — peak-based g is undefined (only " +
                           std::to_string(fa.two_peak_points) + " scan points with two branch peaks)"));
  }
  log("\nScan results: " + out);
  return R;
}

}  // namespace ssp::app
