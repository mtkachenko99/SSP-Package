#include "ssp/plasmonics/dipole_sphere.hpp"

#include "ssp/units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ssp::plasmonics {

namespace {

using cl = std::complex<long double>;

// Риккати–Бессель ψ_n, ξ_n = ψ_n − iχ_n (= ρ h_n^(1)) и производные; ψ — через нисходящую лог. производную.
void riccati_l(cl z, int N, std::vector<cl>& psi, std::vector<cl>& dpsi, std::vector<cl>& xi, std::vector<cl>& dxi,
               std::vector<cl>& D) {
  const int nmx = std::max(N, static_cast<int>(std::abs(z))) + 60;
  std::vector<cl> Dl(nmx + 1, 0.0L);
  for (int n = nmx; n >= 1; --n) Dl[n - 1] = (long double)n / z - 1.0L / (Dl[n] + (long double)n / z);
  psi.assign(N + 1, 0.0L);
  dpsi.assign(N + 1, 0.0L);
  xi.assign(N + 1, 0.0L);
  dxi.assign(N + 1, 0.0L);
  D.assign(Dl.begin(), Dl.begin() + N + 1);
  std::vector<cl> chi(N + 1), dchi(N + 1);
  psi[0] = std::sin(z);
  dpsi[0] = std::cos(z);
  chi[0] = std::cos(z);
  dchi[0] = -std::sin(z);
  for (int n = 1; n <= N; ++n) {
    psi[n] = psi[n - 1] / (Dl[n] + (long double)n / z);
    dpsi[n] = Dl[n] * psi[n];
    const cl chim2 = n >= 2 ? chi[n - 2] : -std::sin(z);
    chi[n] = (2.0L * n - 1.0L) / z * chi[n - 1] - chim2;
    dchi[n] = chi[n - 1] - (long double)n / z * chi[n];
  }
  const cl I(0.0L, 1.0L);
  for (int n = 0; n <= N; ++n) {
    xi[n] = psi[n] - I * chi[n];
    dxi[n] = dpsi[n] - I * dchi[n];
  }
}

}  // namespace

DipoleSphereRates dipole_near_sphere(double a_nm, std::complex<double> eps, double n_medium, double nu_cm, double r_nm,
                                     int n_max) {
  if (r_nm <= a_nm) throw std::invalid_argument("dipole_near_sphere: the dipole must be outside the sphere");
  const long double k = 2.0L * units::pi * n_medium * nu_cm * 1e-7L;
  const long double x = k * a_nm, kr = k * r_nm;
  const cl m = std::sqrt(cl(eps.real(), eps.imag())) / (long double)n_medium;
  int N = n_max;
  if (N <= 0) {  // ближнее поле сходится как (a/r)^{2n}; дальнее — как у обычного ряда Ми
    const double nf = 16.0 * std::log(10.0) / (2.0 * std::log(r_nm / a_nm));
    // Γ_rad содержит и поле самого диполя (Σ по j_n²) — нужны n до ~kr
    N = static_cast<int>(std::max({double(x) + 4.0 * std::cbrt(double(x)) + 4.0, nf + 10.0,
                                   double(kr) + 4.0 * std::cbrt(double(kr)) + 10.0}));
    N = std::min(N, 4000);
  }
  std::vector<cl> px, dpx, xx, dxx, Dx, pm, dpm, xm, dxm, Dm, pr, dpr, xr, dxr, Dr;
  riccati_l(cl(x), N, px, dpx, xx, dxx, Dx);
  riccati_l(m * x, N, pm, dpm, xm, dxm, Dm);
  riccati_l(cl(kr), N, pr, dpr, xr, dxr, Dr);
  long double tp = 0, rp = 0, tl = 0, rl = 0;
  for (int n = 1; n <= N; ++n) {
    // BH: a_n = [D_n(mx)/m + n/x]ψ_n − ψ_{n−1} / ([D_n(mx)/m + n/x]ξ_n − ξ_{n−1}); b_n — с m·D_n
    const cl Da = Dm[n] / m + (long double)n / x, Db = m * Dm[n] + (long double)n / x;
    const cl an = (Da * px[n] - px[n - 1]) / (Da * xx[n] - xx[n - 1]);
    const cl bn = (Db * px[n] - px[n - 1]) / (Db * xx[n] - xx[n - 1]);
    // в точке диполя: j_n = ψ/kr, h_n = ξ/kr, (kr j_n)' = ψ', (kr h_n)' = ξ'
    const cl jn = pr[n] / kr, hn = xr[n] / kr, djn = dpr[n], dhn = dxr[n];
    const long double c = (long double)n * (n + 1) * (2 * n + 1);
    tp += c * (an * hn * hn / (kr * kr)).real();
    rp += c * std::norm(jn - an * hn) / (kr * kr);
    tl += (2.0L * n + 1.0L) * (bn * hn * hn + an * dhn * dhn / (kr * kr)).real();
    rl += (2.0L * n + 1.0L) * (std::norm(jn - bn * hn) + std::norm(djn - an * dhn) / (kr * kr));
  }
  DipoleSphereRates R;
  R.total_perp = double(1.0L - 1.5L * tp);
  R.rad_perp = double(1.5L * rp);
  R.total_par = double(1.0L - 0.75L * tl);
  R.rad_par = double(0.75L * rl);
  R.n_terms = N;
  return R;
}

}  // namespace ssp::plasmonics
