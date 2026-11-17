#include "ssp/plasmonics/particle.hpp"

#include "ssp/units.hpp"

#include <cmath>
#include <stdexcept>

namespace ssp::plasmonics {

namespace {

constexpr double kPi = ssp::units::pi;
const cd I(0.0, 1.0);

// Узлы и веса Гаусса–Лежандра (16 точек) на [−1, 1].
constexpr double kGLx[8] = {0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
                            0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499};
constexpr double kGLw[8] = {0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
                            0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};

// ln I_m(ξ), I_m = ∫₀^∞ cosh(mt)·(ξ + s·cosh t)^{−m−1} dt, s > 0 (вещественная форма), в лог-масштабе.
double log_heine_integral(int m, double xi, double s) {
  auto logf = [&](double t) {
    const double lc = m == 0 ? 0.0 : m * t + std::log1p(std::exp(-2.0 * m * t)) - std::log(2.0);
    return lc - (m + 1) * std::log(xi + s * std::cosh(t));
  };
  // Хвост ~ e^{−t}: интегрируем до точки, где лог-подынтегральное ниже максимума на 40.
  double fmax = logf(0.0);
  double t = 0.0;
  const double dt = 0.25;
  while (true) {
    t += dt;
    const double v = logf(t);
    fmax = std::max(fmax, v);
    if (v < fmax - 40.0 && t > 2.0) break;
    if (t > 2000.0) break;
  }
  const double T = t;
  const int panels = static_cast<int>(std::ceil(T / dt));
  double sum = 0.0;
  for (int p = 0; p < panels; ++p) {
    const double a = p * dt, c = 0.5 * dt, mid = a + c;
    for (int k = 0; k < 8; ++k) {
      sum += kGLw[k] * c * std::exp(logf(mid + c * kGLx[k]) - fmax);
      sum += kGLw[k] * c * std::exp(logf(mid - c * kGLx[k]) - fmax);
    }
  }
  return std::log(sum) + fmax;
}

double log_double_factorial_odd(int m) {  // ln (2m−1)!!
  double v = 0.0;
  for (int k = 1; k <= m; ++k) v += std::log(2.0 * k - 1.0);
  return v;
}

}  // namespace

void normalized_legendre(int L, int m, double eta, std::vector<double>& out) {
  out.assign(static_cast<std::size_t>(L + 1), 0.0);
  if (m > L) return;
  const double s = std::sqrt(std::max(0.0, 1.0 - eta * eta));
  double pmm = std::sqrt(0.5);
  for (int k = 1; k <= m; ++k) pmm *= std::sqrt((2.0 * k + 1.0) / (2.0 * k)) * s;
  out[m] = pmm;
  if (m + 1 <= L) out[m + 1] = std::sqrt(2.0 * m + 3.0) * eta * pmm;
  for (int l = m + 2; l <= L; ++l) {
    const double a = std::sqrt((4.0 * l * l - 1.0) / (double(l) * l - double(m) * m));
    const double b = std::sqrt((double(l - 1) * (l - 1) - double(m) * m) / (4.0 * (l - 1) * (l - 1) - 1.0));
    out[l] = a * (eta * out[l - 1] - b * out[l - 2]);
  }
}

QuasistaticParticle::QuasistaticParticle(double a_xy_nm, double a_z_nm, int l_max)
    : axy_(a_xy_nm), az_(a_z_nm), L_(l_max) {
  if (axy_ <= 0.0 || az_ <= 0.0) throw std::invalid_argument("particle: semi-axes must be > 0");
  if (L_ < 1 || L_ > 120) throw std::invalid_argument("particle: l_max must be in [1, 120]");
  const double amax = std::max(axy_, az_);
  if (std::abs(axy_ - az_) < 1e-9 * amax) {
    shape_ = ParticleShape::Sphere;
    xi0_ = axy_;
  } else if (az_ > axy_) {
    shape_ = ParticleShape::Prolate;
    f_ = std::sqrt(az_ * az_ - axy_ * axy_);
    xi0_ = az_ / f_;
  } else {
    shape_ = ParticleShape::Oblate;
    f_ = std::sqrt(axy_ * axy_ - az_ * az_);
    xi0_ = az_ / f_;
  }
  fd_h_ = 2e-3 * std::min(axy_, az_) / (L_ + 1);

  // Факторы деполяризации
  if (shape_ == ParticleShape::Sphere) {
    L_dep_ = Eigen::Vector3d::Constant(1.0 / 3.0);
  } else if (shape_ == ParticleShape::Prolate) {
    const double e = std::sqrt(1.0 - axy_ * axy_ / (az_ * az_));
    const double Lz = (1.0 - e * e) / (e * e) * (-1.0 + std::log((1.0 + e) / (1.0 - e)) / (2.0 * e));
    L_dep_ << 0.5 * (1.0 - Lz), 0.5 * (1.0 - Lz), Lz;
  } else {
    const double e = std::sqrt(1.0 - az_ * az_ / (axy_ * axy_));
    const double Lz = (1.0 - std::sqrt(1.0 - e * e) / e * std::asin(e)) / (e * e);
    L_dep_ << 0.5 * (1.0 - Lz), 0.5 * (1.0 - Lz), Lz;
  }

  // Члены разложения: нормировка и логарифмические производные на поверхности.
  lQ0_.assign(static_cast<std::size_t>(L_ + 1), {});
  std::vector<std::vector<cd>> lP0(static_cast<std::size_t>(L_ + 1));
  for (int m = 0; m <= L_; ++m) radial_log(xi0_, m, lP0[m], lQ0_[m]);

  const double h = (shape_ == ParticleShape::Sphere) ? 1e-4 * xi0_ : 1e-4 * std::max(xi0_, 0.05);
  std::vector<std::vector<cd>> Pp(5), Qp(5);
  for (int l = 1; l <= L_; ++l) {
    for (int m = 0; m <= l; ++m) {
      MultipoleTerm t;
      t.l = l;
      t.m = m;
      const double dm = (m == 0) ? 1.0 : 2.0;
      if (shape_ == ParticleShape::Sphere) {
        t.aP = double(l) / xi0_;
        t.aQ = -double(l + 1) / xi0_;
        t.w = dm * 2.0 / ((2.0 * l + 1.0) * xi0_);
      } else {
        // d ln P/dξ, d ln Q/dξ — разности 4-го порядка
        cd lp[5], lq[5];
        for (int k = -2; k <= 2; ++k) {
          if (k == 0) continue;
          std::vector<cd> P, Q;
          radial_log(xi0_ + k * h, m, P, Q);
          lp[k + 2] = P[l];
          lq[k + 2] = Q[l];
        }
        auto d4 = [&](const cd* v) { return (v[0] - 8.0 * v[1] + 8.0 * v[3] - v[4]) / (12.0 * h); };
        t.aP = d4(lp);
        t.aQ = d4(lq);
        // w = (2/f)(2−δ_m0)·κ·(l−m)!/(l+m)!·P(ξ₀)Q(ξ₀)
        const cd logw = std::log(2.0 / f_ * dm) + (std::lgamma(l - m + 1.0) - std::lgamma(l + m + 1.0)) +
                        lP0[m][l] + lQ0_[m][l];
        const cd kappa = (shape_ == ParticleShape::Prolate) ? cd(1.0) : I;
        t.w = kappa * std::exp(logw);
      }
      terms_.push_back(t);
    }
  }

  // Однородное поле: r_a = s_a P_1^m(ξ) P̄_1^m(η) trig(φ)  ⇒  u_a = s_a P_1^m(ξ₀).
  const double nz = std::sqrt(1.5), nx = std::sqrt(0.75);
  uterm_ = {term_index(1, 1), term_index(1, 1), term_index(1, 0)};
  usin_ = {false, true, false};
  if (shape_ == ParticleShape::Sphere) {
    ucoef_ = {xi0_ / nx, xi0_ / nx, xi0_ / nz};
  } else if (shape_ == ParticleShape::Prolate) {
    const cd ux = f_ / nx * std::sqrt(xi0_ * xi0_ - 1.0);
    ucoef_ = {ux, ux, f_ / nz * xi0_};
  } else {
    // s_a = −i f/N, P_1^0(iξ₀) = iξ₀, P_1^1(iξ₀) = i√(ξ₀²+1)
    const cd ux = f_ / nx * std::sqrt(xi0_ * xi0_ + 1.0);
    ucoef_ = {ux, ux, cd(f_ / nz * xi0_)};
  }
}

QuasistaticParticle::Coords QuasistaticParticle::coords(const Vec3& r) const {
  const double rho = std::hypot(r.x(), r.y());
  const double phi = std::atan2(r.y(), r.x());
  const double z = r.z();
  if (shape_ == ParticleShape::Sphere) {
    const double R = r.norm();
    return {R, R > 0.0 ? std::clamp(z / R, -1.0, 1.0) : 1.0, phi};
  }
  if (shape_ == ParticleShape::Prolate) {
    const double r1 = std::hypot(rho, z + f_), r2 = std::hypot(rho, z - f_);
    return {(r1 + r2) / (2.0 * f_), std::clamp((r1 - r2) / (2.0 * f_), -1.0, 1.0), phi};
  }
  const double f2 = f_ * f_;
  const double B = f2 - rho * rho - z * z;
  const double disc = std::sqrt(B * B + 4.0 * f2 * z * z);
  const double s = (B > 0.0) ? 2.0 * z * z / (B + disc) : (-B + disc) / (2.0 * f2);
  const double xi = std::sqrt(std::max(s, 0.0));
  double eta;
  if (xi > 1e-14) eta = z / (f_ * xi);
  else eta = std::sqrt(std::max(0.0, 1.0 - rho * rho / f2)) * (z >= 0.0 ? 1.0 : -1.0);
  return {xi, std::clamp(eta, -1.0, 1.0), phi};
}

bool QuasistaticParticle::inside(const Vec3& r) const {
  const double t = (r.x() * r.x() + r.y() * r.y()) / (axy_ * axy_) + r.z() * r.z() / (az_ * az_);
  return t < 1.0;
}

double QuasistaticParticle::surface_gap(const Vec3& r) const {
  const double R = r.norm();
  if (R == 0.0) return -std::min(axy_, az_);
  const Vec3 u = r / R;
  const double rs = 1.0 / std::sqrt((u.x() * u.x() + u.y() * u.y()) / (axy_ * axy_) + u.z() * u.z() / (az_ * az_));
  return R - rs;
}

void QuasistaticParticle::radial_log(double xi, int m, std::vector<cd>& lP, std::vector<cd>& lQ) const {
  // ln P_l^m(ζ), ln Q_l^m(ζ), l = m..L — в логарифмах, чтобы исключить пере/недополнение при больших l и ξ.
  lP.assign(static_cast<std::size_t>(L_ + 1), 0.0);
  lQ.assign(static_cast<std::size_t>(L_ + 1), 0.0);
  if (shape_ == ParticleShape::Sphere) {
    const double lx = std::log(xi);
    for (int l = m; l <= L_; ++l) {
      lP[l] = l * lx;
      lQ[l] = -(l + 1.0) * lx;
    }
    return;
  }
  const bool prolate = (shape_ == ParticleShape::Prolate);
  const cd zeta = prolate ? cd(xi) : I * xi;
  const double sr = prolate ? std::sqrt(std::max(xi * xi - 1.0, 0.0)) : std::sqrt(xi * xi + 1.0);
  const cd s = prolate ? cd(sr) : I * sr;

  // P (доминирующее решение): P_m^m = (2m−1)!!·s^m, отношения σ_l = P_l/P_{l−1} — восходящая рекуррентность.
  lP[m] = log_double_factorial_odd(m) + double(m) * std::log(s);
  cd sigma = 0.0;
  for (int l = m + 1; l <= L_; ++l) {
    sigma = (l == m + 1) ? (2.0 * m + 1.0) * zeta
                         : ((2.0 * l - 1.0) * zeta - double(l - 1 + m) / sigma) / double(l - m);
    lP[l] = lP[l - 1] + std::log(sigma);
  }

  // Q — минимальное решение. Отношения ρ_l = Q_l/Q_{l−1} — нисходящая непрерывная дробь
  // ρ_l = (l+m)/((2l+1)ζ − (l−m+1)ρ_{l+1}); |Q_m^m| = m!·I_m (Уиттекер–Ватсон; знак (−1)^m поглощён),
  // для сплюснутого фаза i^{−m−1}.
  const double rho = xi + sr;  // |ζ + s|
  const int extra = static_cast<int>(std::ceil(20.0 / std::log(rho))) + 20;
  const int top = L_ + std::min(extra, 20000);
  std::vector<cd> ratio(static_cast<std::size_t>(L_ + 1), 0.0);
  cd r_next = 0.0;
  for (int l = top; l > m; --l) {
    const cd r = double(l + m) / ((2.0 * l + 1.0) * zeta - double(l - m + 1) * r_next);
    if (l <= L_) ratio[l] = r;
    r_next = r;
  }
  lQ[m] = std::lgamma(m + 1.0) + log_heine_integral(m, xi, sr);
  if (!prolate) lQ[m] += cd(0.0, -(m + 1.0) * kPi / 2.0);
  for (int l = m + 1; l <= L_; ++l) lQ[l] = lQ[l - 1] + std::log(ratio[l]);
}

cd QuasistaticParticle::response(const MultipoleTerm& t, cd eps, double eps_m) const {
  return (eps - eps_m) * t.aP / (eps_m * t.aQ - eps * t.aP);
}

void QuasistaticParticle::values(const Vec3& r, std::vector<cd>& fc, std::vector<cd>& fs) const {
  const Coords c = coords(r);
  fc.assign(terms_.size(), 0.0);
  fs.assign(terms_.size(), 0.0);
  std::vector<cd> P, Q;
  std::vector<double> Pb;
  for (int m = 0; m <= L_; ++m) {
    radial_log(c.xi, m, P, Q);
    normalized_legendre(L_, m, c.eta, Pb);
    const double cm = std::cos(m * c.phi), sm = std::sin(m * c.phi);
    for (int l = std::max(m, 1); l <= L_; ++l) {
      const cd q = std::exp(Q[l] - lQ0_[m][l]);
      const int k = term_index(l, m);
      fc[k] = q * Pb[l] * cm;
      fs[k] = q * Pb[l] * sm;
    }
  }
}

void QuasistaticParticle::gradients(const Vec3& r, std::vector<Eigen::Vector3cd>& gc,
                                    std::vector<Eigen::Vector3cd>& gs) const {
  const std::size_t T = terms_.size();
  gc.assign(T, Eigen::Vector3cd::Zero());
  gs.assign(T, Eigen::Vector3cd::Zero());
  const double h = fd_h_;
  // На оси симметрии φ не определён, но f^{c,s} гладкие (∝ ρ^m) и точки шаблона лежат вне оси.
  const Vec3 r0 = r;
  std::vector<cd> fc[4], fs[4];
  const double cf[4] = {1.0, -8.0, 8.0, -1.0};
  const int off[4] = {-2, -1, 1, 2};
  for (int a = 0; a < 3; ++a) {
    for (int k = 0; k < 4; ++k) {
      Vec3 rp = r0;
      rp(a) += off[k] * h;
      values(rp, fc[k], fs[k]);
    }
    for (std::size_t t = 0; t < T; ++t) {
      cd dc = 0.0, ds = 0.0;
      for (int k = 0; k < 4; ++k) {
        dc += cf[k] * fc[k][t];
        ds += cf[k] * fs[k][t];
      }
      gc[t](a) = dc / (12.0 * h);
      gs[t](a) = ds / (12.0 * h);
    }
  }
}

cd QuasistaticParticle::induced_kernel(const Vec3& r, const Vec3& rp, cd eps, double eps_m) const {
  std::vector<cd> fc1, fs1, fc2, fs2;
  values(r, fc1, fs1);
  values(rp, fc2, fs2);
  cd K = 0.0;
  for (std::size_t t = 0; t < terms_.size(); ++t)
    K += terms_[t].w * response(terms_[t], eps, eps_m) * (fc1[t] * fc2[t] + fs1[t] * fs2[t]);
  return K;
}

cd QuasistaticParticle::free_kernel_expansion(const Vec3& r_in, const Vec3& r_out, int l_max) const {
  // 1/|r−r'| = Σ_{l,m} w_lm [P(ξ<)/P(ξ₀)][Q(ξ>)/Q(ξ₀)] P̄(η)P̄(η') cos mΔφ (включая l = 0).
  const Coords a = coords(r_in), b = coords(r_out);
  const int L = std::min(l_max, L_);
  cd S = 0.0;
  std::vector<cd> Pa, Qa, Pb, Qb, P0, Q0;
  std::vector<double> Na, Nb;
  for (int m = 0; m <= L; ++m) {
    radial_log(a.xi, m, Pa, Qa);
    radial_log(b.xi, m, Pb, Qb);
    radial_log(xi0_, m, P0, Q0);
    normalized_legendre(L_, m, a.eta, Na);
    normalized_legendre(L_, m, b.eta, Nb);
    const double dm = (m == 0) ? 1.0 : 2.0;
    for (int l = m; l <= L; ++l) {
      cd w;
      if (shape_ == ParticleShape::Sphere) {
        w = dm * 2.0 / ((2.0 * l + 1.0) * xi0_);
      } else {
        const cd logw = std::log(2.0 / f_ * dm) + (std::lgamma(l - m + 1.0) - std::lgamma(l + m + 1.0)) +
                        P0[l] + Q0[l];
        w = ((shape_ == ParticleShape::Prolate) ? cd(1.0) : I) * std::exp(logw);
      }
      S += w * std::exp(Pa[l] - P0[l] + Qb[l] - Q0[l]) * Na[l] * Nb[l] * std::cos(m * (a.phi - b.phi));
    }
  }
  return S;
}

Eigen::Vector3cd QuasistaticParticle::polarizability(cd eps, double eps_m) const {
  Eigen::Vector3cd a;
  const double v3 = axy_ * axy_ * az_ / 3.0;
  for (int k = 0; k < 3; ++k) a(k) = v3 * (eps - eps_m) / (eps_m + L_dep_(k) * (eps - eps_m));
  return a;
}

}  // namespace ssp::plasmonics
