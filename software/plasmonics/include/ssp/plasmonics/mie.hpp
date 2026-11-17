#pragma once
// Теория Ми для однородной сферы (Bohren & Huffman, BHMIE; время e^{−iωt}).

#include <complex>
#include <vector>

namespace ssp::plasmonics {

struct MieResult {
  std::vector<std::complex<double>> a, b;  // a_n, b_n, n = 1..n_max (индекс 0 ↔ n = 1)
  double q_ext = 0.0, q_sca = 0.0, q_abs = 0.0;
  double x = 0.0;
};

// x = k·a (k — волновое число в среде), m — относительный показатель преломления (n_частицы/n_среды).
MieResult mie(double x, std::complex<double> m, int n_max = 0);

// Сечения (нм²) для сферы радиуса a_nm, ε частицы, n_med среды, волновое число ν̃ (см⁻¹).
MieResult mie_sphere(double a_nm, std::complex<double> eps, double n_medium, double nu_cm);

// Динамическая дипольная поляризуемость (нм³, p = 4πε₀ε_m α E): α = (3i/2k³)·a₁.
std::complex<double> mie_dipole_polarizability(const MieResult& r, double k_per_nm);

}  // namespace ssp::plasmonics
