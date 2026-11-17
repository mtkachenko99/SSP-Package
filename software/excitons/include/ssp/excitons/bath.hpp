#pragma once
// Гармоническая баня на каждом сайте (некоррелированная): спектральная плотность J(ω) и производные величины.
// Соглашения (ħ = 1, ω в см⁻¹): λ = (1/π)∫J(ω)/ω dω,
// C(t) = (1/π)∫₀^∞ J(ω)[coth(βω/2)cos ωt − i sin ωt] dω,  C̃(ω) = ∫ e^{iωt}C(t)dt = 2J(ω)[n(ω)+1].

#include <complex>
#include <vector>

namespace ssp::excitons {

enum class BathType { DrudeLorentz, OhmicExponential };

// Недемпфированная колебательная мода (броуновский осциллятор): λ_k = S·ω_k.
struct VibrationalMode {
  double omega_cm = 0.0;
  double huang_rhys = 0.0;
  double gamma_cm = 10.0;
};

struct Bath {
  BathType type = BathType::DrudeLorentz;
  double lambda_cm = 50.0;   // энергия реорганизации обертона (Друде/омической части)
  double cutoff_cm = 100.0;  // ω_c (Друде: γ)
  std::vector<VibrationalMode> modes;
  double temperature_K = 300.0;

  double beta() const;                 // 1/(k_B T) в (см⁻¹)⁻¹; +inf при T = 0
  double J_over_w(double w) const;     // J(ω)/ω для ω ≥ 0 (конечна при ω → 0)
  double J(double w) const { return w * J_over_w(w); }
  double reorganization() const;       // полная λ (обертон + моды)
  // C̃(ω) = 2 J(|ω|)·ω/(1 − e^{−βω}); детальный баланс C̃(−ω) = e^{−βω}C̃(ω).
  double rate_spectrum(double w) const;
};

// Бозе-фактор n(ω) = 1/(e^{βω} − 1).
double bose(double w_cm, double beta);

// Функция уширения g(t) на равномерной сетке t_j = j·dt (единицы t_u), j = 0..n−1.
// Вычисляется как g(t) = ∫₀ᵗ(t−τ)C(τ)dτ, C(t) — через БПФ спектральной суммы на подсетке dt/oversample.
std::vector<std::complex<double>> lineshape_function(const Bath& bath, double dt, std::size_t n, int oversample = 8);

// Аналитическое g(t) для бани Друде–Лоренца (мацубаровский ряд) — для проверки.
std::complex<double> lineshape_function_drude(double lambda, double gamma, double beta, double t,
                                              int n_matsubara = 20000);

}  // namespace ssp::excitons
