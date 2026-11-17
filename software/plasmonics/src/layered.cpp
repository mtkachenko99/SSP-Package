#include "ssp/plasmonics/layered.hpp"

#include "ssp/units.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>

namespace ssp::plasmonics {

double SphereLayers::outer_radius() const {
  double r = core_radius_nm;
  for (double t : shell_thickness_nm) r += t;
  return r;
}

cd layered_response(int l, cd eps_core, const SphereLayers& s, double eps_m) {
  if (s.shell_thickness_nm.size() != s.shell_eps.size()) throw std::invalid_argument("layers: thickness/eps mismatch");
  // Неизвестные: P_0 (ядро), (P_k, Q_k) оболочек, Q_out. Уравнения: по 2 на каждую границу.
  const int K = static_cast<int>(s.shell_thickness_nm.size());
  const int n = 1 + 2 * K + 1;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(n, n);
  Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n);
  std::vector<cd> eps = {eps_core};
  for (const cd& e : s.shell_eps) eps.push_back(e);
  eps.push_back(eps_m);
  // Индексы: слой 0: P → 0; слой k (1..K): P → 2k−1, Q → 2k; внешний: Q → n−1, P = 1 (в правую часть).
  // Радиальные функции масштабированы на внешний радиус b: (r/b)^l и (b/r)^{l+1} — иначе при больших l
  // элементы матрицы различаются на десятки порядков. Тогда D̂_l = Q_out (с P_out = 1).
  auto idxP = [&](int k) { return k == 0 ? 0 : 2 * k - 1; };
  auto idxQ = [&](int k) { return k == K + 1 ? n - 1 : 2 * k; };
  const double bR = s.outer_radius();
  double R = s.core_radius_nm;
  const double L = l;
  for (int b = 0; b <= K; ++b) {  // граница между слоем b и b+1
    if (b > 0) R += s.shell_thickness_nm[b - 1];
    const cd d = (b == 0) ? s.d_perp_nm : cd(0.0);
    const double rp = std::pow(R / bR, l), rq = std::pow(bR / R, l + 1);
    const double dp = L / R * rp, dq = -(L + 1) / R * rq;
    const int e1 = 2 * b, e2 = 2 * b + 1;
    // (1) ε_in ∂φ_in − ε_out ∂φ_out = 0
    // (2) φ_out − φ_in + d(∂φ_out − ∂φ_in) = 0
    auto addInner = [&](int col, double val, double der) {
      A(e1, col) += eps[b] * der;
      A(e2, col) += -val - d * der;
    };
    auto addOuter = [&](int col, double val, double der, bool isP) {
      if (b + 1 == K + 1 && isP) {  // внешний P = 1 → в правую часть
        rhs(e1) -= -eps[b + 1] * der;
        rhs(e2) -= val + d * der;
        return;
      }
      A(e1, col) += -eps[b + 1] * der;
      A(e2, col) += val + d * der;
    };
    addInner(idxP(b), rp, dp);
    if (b > 0) addInner(idxQ(b), rq, dq);
    addOuter(b + 1 == K + 1 ? -1 : idxP(b + 1), rp, dp, true);
    addOuter(idxQ(b + 1), rq, dq, false);
  }
  const Eigen::VectorXcd x = A.partialPivLu().solve(rhs);
  return x(n - 1);
}

namespace {

// Риккати–Бессель для комплексного аргумента: ψ_n — через логарифмическую производную (нисходящая, устойчиво),
// χ_n — восходящая рекуррентность. Возвращает ψ_n, ψ'_n, χ_n, χ'_n для n = 0..N.
void riccati(cd z, int N, std::vector<cd>& psi, std::vector<cd>& dpsi, std::vector<cd>& chi, std::vector<cd>& dchi) {
  const int nmx = std::max(N, static_cast<int>(std::abs(z))) + 30;
  std::vector<cd> D(nmx + 1, 0.0);
  for (int n = nmx; n >= 1; --n) D[n - 1] = double(n) / z - 1.0 / (D[n] + double(n) / z);
  psi.assign(N + 1, 0.0);
  dpsi.assign(N + 1, 0.0);
  chi.assign(N + 1, 0.0);
  dchi.assign(N + 1, 0.0);
  psi[0] = std::sin(z);
  chi[0] = std::cos(z);
  dpsi[0] = std::cos(z);
  dchi[0] = -std::sin(z);
  for (int n = 1; n <= N; ++n) {
    psi[n] = psi[n - 1] / (D[n] + double(n) / z);
    dpsi[n] = D[n] * psi[n];
    const cd chim2 = (n >= 2) ? chi[n - 2] : -std::sin(z);  // χ_{−1} = −sin z
    chi[n] = (2.0 * n - 1.0) / z * chi[n - 1] - chim2;
    dchi[n] = chi[n - 1] - double(n) / z * chi[n];
  }
}

}  // namespace

MieResult mie_coated(double x, double y, cd m1, cd m2, int n_max) {
  MieResult r;
  r.x = y;
  const int N = n_max > 0 ? n_max : static_cast<int>(y + 4.0 * std::cbrt(y) + 2.0);
  std::vector<cd> p1x, dp1x, c1x, dc1x, p2x, dp2x, c2x, dc2x, p2y, dp2y, c2y, dc2y, py, dpy, cy, dcy;
  riccati(m1 * x, N, p1x, dp1x, c1x, dc1x);
  riccati(m2 * x, N, p2x, dp2x, c2x, dc2x);
  riccati(m2 * y, N, p2y, dp2y, c2y, dc2y);
  riccati(cd(y), N, py, dpy, cy, dcy);
  double qe = 0, qs = 0;
  for (int n = 1; n <= N; ++n) {
    const cd An = (m2 * p2x[n] * dp1x[n] - m1 * dp2x[n] * p1x[n]) / (m2 * c2x[n] * dp1x[n] - m1 * dc2x[n] * p1x[n]);
    const cd Bn = (m2 * p1x[n] * dp2x[n] - m1 * p2x[n] * dp1x[n]) / (m2 * dc2x[n] * p1x[n] - m1 * dp1x[n] * c2x[n]);
    const cd xi = py[n] - cd(0, 1) * cy[n], dxi = dpy[n] - cd(0, 1) * dcy[n];
    const cd fa = dp2y[n] - An * dc2y[n], ga = p2y[n] - An * c2y[n];
    const cd fb = dp2y[n] - Bn * dc2y[n], gb = p2y[n] - Bn * c2y[n];
    const cd an = (py[n] * fa - m2 * dpy[n] * ga) / (xi * fa - m2 * dxi * ga);
    const cd bn = (m2 * py[n] * fb - dpy[n] * gb) / (m2 * xi * fb - dxi * gb);
    r.a.push_back(an);
    r.b.push_back(bn);
    qe += (2.0 * n + 1.0) * (an.real() + bn.real());
    qs += (2.0 * n + 1.0) * (std::norm(an) + std::norm(bn));
  }
  r.q_ext = 2.0 / (y * y) * qe;
  r.q_sca = 2.0 / (y * y) * qs;
  r.q_abs = r.q_ext - r.q_sca;
  return r;
}

MieResult mie_coated_sphere(double a_nm, double b_nm, cd eps_core, cd eps_shell, double n_medium, double nu_cm,
                            int n_max) {
  const double k = 2.0 * units::pi * n_medium * nu_cm * 1e-7;
  MieResult r = mie_coated(k * a_nm, k * b_nm, std::sqrt(eps_core) / n_medium, std::sqrt(eps_shell) / n_medium, n_max);
  const double geo = units::pi * b_nm * b_nm;
  r.q_ext *= geo;
  r.q_sca *= geo;
  r.q_abs *= geo;
  return r;
}

std::vector<cd> mie_dynamic_factors(const SphereLayers& s, cd eps_core, double n_medium, double nu_cm, int l_max,
                                    int l_dyn) {
  std::vector<cd> f(static_cast<std::size_t>(l_max + 1), cd(1.0));
  const int L = std::min(l_max, l_dyn);
  const double k = 2.0 * units::pi * n_medium * nu_cm * 1e-7;
  const double b = s.outer_radius();
  // Квазистатический предел: a_l → κ_l x^{2l+1} D̂_l, κ_l = i(l+1)/[l(2l+1)!!(2l−1)!!], x = k b (B&H, малые x)
  const double x = k * b;
  // Однородная сфера — через тот же устойчивый расчёт (оболочка из материала ядра): восходящая рекуррентность ψ_n
  // в BHMIE теряет точность при l ≫ x.
  if (!s.homogeneous() && s.shell_thickness_nm.size() != 1)
    throw std::invalid_argument("Mie dressing: only one shell is supported");
  const cd m1 = std::sqrt(eps_core) / n_medium;
  const auto a = s.homogeneous() ? mie_coated(0.5 * x, x, m1, m1, L).a
                                 : mie_coated(k * s.core_radius_nm, x, m1, std::sqrt(s.shell_eps[0]) / n_medium, L).a;
  SphereLayers sq = s;
  sq.d_perp_nm = 0.0;
  const double eps_m = n_medium * n_medium;
  double logk = 0.0;  // log[x^{2l+1}/((2l+1)!!(2l−1)!!)]
  for (int l = 1; l <= L; ++l) {
    logk += 2.0 * std::log(x) - std::log(2.0 * l + 1.0) - std::log(2.0 * l - 1.0);
    const cd qs = cd(0.0, (l + 1.0) / l) * std::exp(logk + std::log(x)) * layered_response(l, eps_core, sq, eps_m);
    if (std::abs(qs) > 0 && static_cast<std::size_t>(l) <= a.size()) f[l] = a[l - 1] / qs;
  }
  return f;
}

}  // namespace ssp::plasmonics
