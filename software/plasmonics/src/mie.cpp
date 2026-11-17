#include "ssp/plasmonics/mie.hpp"

#include "ssp/units.hpp"

#include <cmath>
#include <stdexcept>

namespace ssp::plasmonics {

MieResult mie(double x, std::complex<double> m, int n_max) {
  using cd = std::complex<double>;
  if (x <= 0.0) throw std::invalid_argument("mie: x must be > 0");
  MieResult r;
  r.x = x;
  const int nstop = n_max > 0 ? n_max : static_cast<int>(x + 4.0 * std::cbrt(x) + 2.0);
  const cd mx = m * x;
  const int nmx = static_cast<int>(std::max<double>(nstop, std::abs(mx))) + 16;
  // Логарифмическая производная D_n(mx), нисходящая рекуррентность.
  std::vector<cd> D(nmx + 1, 0.0);
  for (int n = nmx; n >= 1; --n) D[n - 1] = double(n) / mx - 1.0 / (D[n] + double(n) / mx);
  // Функции Риккати–Бесселя ψ_n, χ_n (восходящие), ξ_n = ψ_n − iχ_n.
  double psi0 = std::cos(x), psi1 = std::sin(x);
  double chi0 = -std::sin(x), chi1 = std::cos(x);
  cd xi1(psi1, -chi1);
  double qext = 0.0, qsca = 0.0;
  for (int n = 1; n <= nstop; ++n) {
    const double fn = n;
    const double psi = (2.0 * fn - 1.0) * psi1 / x - psi0;
    const double chi = (2.0 * fn - 1.0) * chi1 / x - chi0;
    const cd xi(psi, -chi);
    const cd an = ((D[n] / m + fn / x) * psi - psi1) / ((D[n] / m + fn / x) * xi - xi1);
    const cd bn = ((m * D[n] + fn / x) * psi - psi1) / ((m * D[n] + fn / x) * xi - xi1);
    r.a.push_back(an);
    r.b.push_back(bn);
    qext += (2.0 * fn + 1.0) * (an.real() + bn.real());
    qsca += (2.0 * fn + 1.0) * (std::norm(an) + std::norm(bn));
    psi0 = psi1;
    psi1 = psi;
    chi0 = chi1;
    chi1 = chi;
    xi1 = cd(psi1, -chi1);
  }
  r.q_ext = 2.0 / (x * x) * qext;
  r.q_sca = 2.0 / (x * x) * qsca;
  r.q_abs = r.q_ext - r.q_sca;
  return r;
}

MieResult mie_sphere(double a_nm, std::complex<double> eps, double n_medium, double nu_cm) {
  const double k = 2.0 * units::pi * n_medium * nu_cm * 1e-7;  // нм⁻¹
  MieResult r = mie(k * a_nm, std::sqrt(eps) / n_medium);
  const double geo = units::pi * a_nm * a_nm;
  r.q_ext *= geo;  // далее q_* содержат сечения в нм²
  r.q_sca *= geo;
  r.q_abs *= geo;
  return r;
}

std::complex<double> mie_dipole_polarizability(const MieResult& r, double k_per_nm) {
  return std::complex<double>(0.0, 1.5) * r.a.at(0) / (k_per_nm * k_per_nm * k_per_nm);
}

}  // namespace ssp::plasmonics
