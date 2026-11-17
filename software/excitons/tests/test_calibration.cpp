// Калибровка по мономеру: восстановление параметров по синтетическим спектрам и согласованность μ.
#include "ssp/excitons/monomer_fit.hpp"
#include "ssp/excitons/spectra.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>

using namespace ssp::excitons;

namespace {
MonomerFitParams truth() {
  MonomerFitParams p;
  p.E0_cm = 18500.0;
  p.lambda_cm = 120.0;
  p.cutoff_cm = 150.0;
  p.sigma_inh_cm = 200.0;
  VibrationalMode m;
  m.omega_cm = 1200.0;
  m.huang_rhys = 0.5;
  m.gamma_cm = 20.0;
  p.modes = {m};
  return p;
}
}  // namespace

TEST(Calibration, MonomerFitRecoversSyntheticParameters) {
  const MonomerFitParams T = truth();
  Spectrum a, e;
  for (double v = 15000.0; v <= 23000.0; v += 20.0) {
    a.w.push_back(v);
    e.w.push_back(v);
  }
  std::vector<double> A, E, A2, E2;
  monomer_model(T, a.w, A, E);
  const double ma = *std::max_element(A.begin(), A.end()), me = *std::max_element(E.begin(), E.end());
  std::mt19937_64 rng(4);
  std::normal_distribution<double> nz(0.0, 0.01);
  for (std::size_t i = 0; i < a.w.size(); ++i) {
    a.y.push_back(A[i] / ma * 7e4 + 7e4 * nz(rng));
    a.sigma.push_back(700.0);
    e.y.push_back(E[i] / me + nz(rng));
    e.sigma.push_back(0.01);
  }
  MonomerFitParams s = T;
  s.E0_cm += 300.0;
  s.lambda_cm = 60.0;
  s.cutoff_cm = 80.0;
  s.sigma_inh_cm = 300.0;
  s.modes[0].huang_rhys = 0.3;
  const auto R = fit_monomer(a, e, s, false, 1.33);
  const std::vector<double> want = {T.E0_cm, T.lambda_cm, T.cutoff_cm, T.sigma_inh_cm, T.modes[0].huang_rhys};
  ASSERT_EQ(R.values.size(), want.size());
  for (std::size_t k = 0; k < want.size(); ++k) {
    EXPECT_NEAR(R.values[k], want[k], std::max(4.0 * R.errors[k], 1e-3 * std::abs(want[k]))) << R.names[k];
    EXPECT_NEAR(R.values[k] / want[k], 1.0, 0.1) << R.names[k];
  }
  EXPECT_NEAR(R.chi2_red, 1.0, 0.2);
}

TEST(Calibration, TransitionDipoleFromExtinctionIsConsistentWithCrossSection) {
  // ε(ν̃) построено из σ_abs модели агрегата для μ = 8 Д → интеграл возвращает μ² = 64 Д²
  const MonomerFitParams T = truth();
  Spectrum eps;
  for (double v = 12000.0; v <= 30000.0; v += 10.0) eps.w.push_back(v);
  std::vector<double> A, E;
  monomer_model(T, eps.w, A, E);  // ν̃·L(ν̃), ∫L = 1
  const double P = absorption_cross_section_prefactor_nm2(1.33), NA = 6.02214076e23;
  for (std::size_t i = 0; i < eps.w.size(); ++i) {
    const double sigma_nm2 = P * 64.0 * A[i];
    eps.y.push_back(sigma_nm2 * 1e-14 * NA / (std::log(10.0) * 1000.0));
  }
  EXPECT_NEAR(mu2_from_extinction(eps, 1.33), 64.0, 0.3);
  // при n = 1 — стандартная формула μ² = 9.186·10⁻³∫ε/ν̃ dν̃ (Д²); в среде модель не вводит поправку локального
  // поля (σ ∝ 1/n, k_rad ∝ n), поэтому μ² ∝ n∫ε/ν̃ — согласованно со Стриклером–Бергом (k_r ∝ n²∫ε/ν̃)
  double integ = 0;
  for (std::size_t i = 1; i < eps.w.size(); ++i)
    integ += 0.5 * (eps.y[i] / eps.w[i] + eps.y[i - 1] / eps.w[i - 1]) * (eps.w[i] - eps.w[i - 1]);
  EXPECT_NEAR(9.186e-3 * integ / mu2_from_extinction(eps, 1.0), 1.0, 0.002);
  EXPECT_NEAR(mu2_from_extinction(eps, 1.33) / mu2_from_extinction(eps, 1.0), 1.33, 1e-9);
  // τ_rad по формуле модели = Стриклер–Берг (при той же ⟨ν̃⁻³⟩)
  Spectrum none;
  MonomerFitParams fixed = T;
  fixed.fit_E0 = fixed.fit_lambda = fixed.fit_cutoff = fixed.fit_sigma = fixed.fit_mode_S = false;
  const auto R = fit_monomer(eps, none, fixed, true, 1.33);
  EXPECT_NEAR(R.mu_D, 8.0, 0.02);
  EXPECT_NEAR(R.tau_rad_ns / R.tau_sb_ns, 1.0, 0.01);
}
