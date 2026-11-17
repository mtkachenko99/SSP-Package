#pragma once
// Точная электродинамика (теория Ми): скорости распада точечного диполя у однородной сферы
// (Ruppin 1982; Kim, Leung, Liu 1988; обзор — Novotny & Hecht, «Principles of Nano-Optics»).
// Нормировка — на излучательную скорость Γ₀ того же диполя в однородной среде n_medium.
//   радиальный (⊥):   Γ/Γ₀ = 1 − (3/2) Re Σ n(n+1)(2n+1) a_n [h_n(kr)/kr]²
//                    Γ_rad/Γ₀ = (3/2) Σ n(n+1)(2n+1) |j_n(kr) − a_n h_n(kr)|²/(kr)²
//   касательный (∥): Γ/Γ₀ = 1 − (3/4) Re Σ (2n+1)[b_n h_n² + a_n ((kr h_n)'/kr)²]
//                    Γ_rad/Γ₀ = (3/4) Σ (2n+1)[|j_n − b_n h_n|² + |(kr j_n)' − a_n (kr h_n)'|²/(kr)²]
// a_n, b_n — коэффициенты Ми (Bohren & Huffman, e^{−iωt}), h_n = h_n^(1). Знаки проверены законом сохранения
// энергии: для непоглощающей сферы Γ = Γ_rad. Расчёт в long double (h_n(kr) ~ (2n−1)!!/(kr)^{n+1} у поверхности).

#include <complex>

namespace ssp::plasmonics {

struct DipoleSphereRates {
  double total_perp = 1, rad_perp = 1;  // радиальная ориентация
  double total_par = 1, rad_par = 1;    // касательная
  int n_terms = 0;
};

// a_nm — радиус, r_nm — расстояние диполя от центра (> a), ε сферы, ν̃ (см⁻¹), n_max = 0 — автоматически.
DipoleSphereRates dipole_near_sphere(double a_nm, std::complex<double> eps, double n_medium, double nu_cm,
                                     double r_nm, int n_max = 0);

}  // namespace ssp::plasmonics
