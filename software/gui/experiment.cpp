#include "experiment.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

double interp(const QVector<double>& x, const QVector<double>& y, double t) {
  if (x.isEmpty() || t < x.front() || t > x.back()) return std::numeric_limits<double>::quiet_NaN();
  const auto it = std::upper_bound(x.begin(), x.end(), t);
  const int i = std::min<int>(int(it - x.begin()), x.size() - 1);
  if (i == 0) return y[0];
  const double f = (t - x[i - 1]) / (x[i] - x[i - 1]);
  return (1.0 - f) * y[i - 1] + f * y[i];
}

struct Fit {
  double s = 0, b = 0, chi2 = std::numeric_limits<double>::infinity();
  int n = 0;
};

Fit linearFit(const QVector<double>& mx, const QVector<double>& my, const QVector<double>& ex,
              const QVector<double>& ey, const QVector<double>& es, double shift, bool offset) {
  // min Σ w(y − s f − b)²
  double Sw = 0, Sf = 0, Sy = 0, Sff = 0, Sfy = 0, Syy = 0;
  int n = 0;
  for (int i = 0; i < ex.size(); ++i) {
    const double f = interp(mx, my, ex[i] - shift);
    if (!std::isfinite(f)) continue;
    const double w = es.isEmpty() ? 1.0 : 1.0 / (es[i] * es[i]);
    Sw += w;
    Sf += w * f;
    Sy += w * ey[i];
    Sff += w * f * f;
    Sfy += w * f * ey[i];
    Syy += w * ey[i] * ey[i];
    ++n;
  }
  Fit r;
  r.n = n;
  if (n < 3) return r;
  if (offset) {
    const double det = Sff * Sw - Sf * Sf;
    if (std::abs(det) < 1e-300) return r;
    r.s = (Sfy * Sw - Sf * Sy) / det;
    r.b = (Sff * Sy - Sf * Sfy) / det;
  } else {
    if (Sff <= 0) return r;
    r.s = Sfy / Sff;
    r.b = 0;
  }
  r.chi2 = Syy - 2 * r.s * Sfy - 2 * r.b * Sy + r.s * r.s * Sff + 2 * r.s * r.b * Sf + r.b * r.b * Sw;
  return r;
}

}  // namespace

void convertToWavenumber(QVector<double>& x, QVector<double>& y, QVector<double>& s, const QString& unit) {
  if (unit == "nm") {
    for (double& v : x) v = 1e7 / v;
  } else if (unit == "eV") {
    for (double& v : x) v *= 8065.54429;
  }
  // упорядочить по возрастанию x
  QVector<int> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });
  QVector<double> x2, y2, s2;
  for (int i : idx) {
    x2 << x[i];
    y2 << y[i];
    if (!s.isEmpty()) s2 << s[i];
  }
  x.swap(x2);
  y.swap(y2);
  s.swap(s2);
}

CompareResult compareSpectra(const QVector<double>& mx, const QVector<double>& my, const QVector<double>& ex,
                             const QVector<double>& ey, const QVector<double>& es, const CompareOptions& opt) {
  CompareResult R;
  R.have_sigma = !es.isEmpty();
  double shift = 0.0;
  if (opt.fit_shift) {
    // Грубая сетка по δ, затем золотое сечение около минимума
    const int N = 121;
    double best = std::numeric_limits<double>::infinity(), bd = 0;
    for (int k = 0; k < N; ++k) {
      const double d = -opt.shift_range_cm + 2.0 * opt.shift_range_cm * k / (N - 1);
      const Fit f = linearFit(mx, my, ex, ey, es, d, opt.fit_offset);
      if (f.chi2 < best) {
        best = f.chi2;
        bd = d;
      }
    }
    const double h = 2.0 * opt.shift_range_cm / (N - 1);
    double a = bd - h, b = bd + h;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double c = b - g * (b - a), d = a + g * (b - a);
    for (int it = 0; it < 60; ++it) {
      if (linearFit(mx, my, ex, ey, es, c, opt.fit_offset).chi2 < linearFit(mx, my, ex, ey, es, d, opt.fit_offset).chi2) b = d;
      else a = c;
      c = b - g * (b - a);
      d = a + g * (b - a);
    }
    shift = 0.5 * (a + b);
  }
  const Fit f = linearFit(mx, my, ex, ey, es, shift, opt.fit_offset);
  if (f.n < 3 || !std::isfinite(f.chi2)) {
    R.message = "the experiment does not overlap the model range";
    return R;
  }
  R.scale = f.s;
  R.offset = f.b;
  R.shift = shift;
  R.n = f.n;
  R.dof = f.n - 1 - (opt.fit_offset ? 1 : 0) - (opt.fit_shift ? 1 : 0);
  double ss = 0.0;
  for (int i = 0; i < ex.size(); ++i) {
    const double m = interp(mx, my, ex[i] - shift);
    if (!std::isfinite(m)) continue;
    const double mod = f.s * m + f.b;
    const double r = ey[i] - mod;
    R.x << ex[i];
    R.y << ey[i];
    R.model << mod;
    R.residual << r;
    ss += r * r;
  }
  R.rmse = std::sqrt(ss / R.n);
  const double sig_est = std::sqrt(ss / std::max(1, R.dof));
  for (int i = 0, k = 0; i < ex.size(); ++i) {
    if (!std::isfinite(interp(mx, my, ex[i] - shift))) continue;
    R.normalized_residual << R.residual[k] / (R.have_sigma ? es[i] : sig_est);
    ++k;
  }
  R.chi2 = f.chi2;
  R.chi2_nu = R.have_sigma ? f.chi2 / std::max(1, R.dof) : 1.0;
  double num = 0, den = 0, c1 = 0;
  for (int i = 0; i < R.residual.size(); ++i) {
    den += R.residual[i] * R.residual[i];
    if (i > 0) {
      num += std::pow(R.residual[i] - R.residual[i - 1], 2);
      c1 += R.residual[i] * R.residual[i - 1];
    }
  }
  R.durbin_watson = den > 0 ? num / den : 0.0;
  R.lag1_autocorr = den > 0 ? c1 / den : 0.0;
  R.ok = true;
  return R;
}
