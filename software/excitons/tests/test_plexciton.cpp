// Связь агрегат–НЧ: предельные случаи.
#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/plexciton.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/plasmonics/dipole_sphere.hpp"
#include "ssp/units.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

using namespace ssp;
using namespace ssp::excitons;
using cd = std::complex<double>;

TEST(Plexciton, TransparentParticleReducesToBareAggregate) {
  // ε_НЧ = ε_m: W = 0, L = 1, α_NP = 0 ⇒ σ_ext — сумма лоренцианов экситонов.
  plasmonics::DrudeParams dp;
  dp.eps_inf = 1.33 * 1.33;
  dp.omega_p_eV = 0.0;
  const auto mat = plasmonics::Material::drude(dp);
  Aggregate a = build_geometry("ring", {{"N", 6}, {"radius", 2.0}});
  transform(a, Vec3(0, 0, 20.0), 0, 0, 0);
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = 10.0;
  spec.l_max = 10;
  spec.n_medium = 1.33;
  const PlexcitonSystem sys(a, spec, mat, 10.0);

  CouplingParams cp;
  const Eigen::MatrixXd J = coupling_matrix(a, cp);
  std::mt19937_64 rng(3);
  const Eigen::MatrixXd H = frenkel_hamiltonian(J, 17000.0, 0.0, rng);
  const ExcitonStates st = diagonalize(H, a, 10.0);
  std::vector<double> w = {16500.0, 17000.0, 17300.0};
  PlexcitonSpectra mean, se;
  Realization Rz;
  Rz.H = H;
  Rz.states = st;
  Rz.K_ps = Eigen::MatrixXd::Zero(st.size(), st.size());
  Rz.hwhm_cm = Eigen::VectorXd::Constant(st.size(), 50.0);
  sys.spectra({Rz}, nullptr, 50.0, w, mean, se);
  const double eps_m = 1.33 * 1.33;
  for (std::size_t q = 0; q < w.size(); ++q) {
    const double k = 2.0 * units::pi * 1.33 * w[q] * 1e-7;
    double im = 0.0;
    for (Eigen::Index n = 0; n < st.energies.size(); ++n)
      im += st.mu2(n) * 50.0 / ((st.energies(n) - w[q]) * (st.energies(n) - w[q]) + 2500.0);
    const double ext = 4.0 * units::pi * k * units::dd_cm_per_D2_nm3 / eps_m * im / 3.0;
    EXPECT_NEAR(mean.ext[q] / ext, 1.0, 1e-8);
    EXPECT_NEAR(mean.np_ext[q], 0.0, 1e-12);
  }
}

TEST(Plexciton, QuenchingRateGrowsTowardSurface) {
  plasmonics::DrudeParams dp;
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto mat = plasmonics::Material::drude(dp);
  Aggregate a;
  for (double d : {1.0, 2.0, 4.0, 8.0}) {
    a.positions.emplace_back(0, 0, 10.0 + d);
    a.dipoles.emplace_back(0, 0, 1);
  }
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = 10.0;
  spec.l_max = 60;
  const PlexcitonSystem sys(a, spec, mat, 10.0);
  const auto nf = sys.near_field(19000.0);
  for (std::size_t i = 1; i < nf.size(); ++i) EXPECT_GT(nf[i - 1].decay_ps, nf[i].decay_ps);
  EXPECT_GT(nf[0].decay_ps, 0.0);
  EXPECT_LT(sys.truncation_estimate(19000.0), 1e-2);
}

TEST(Plexciton, TransparentParticleMatchesCumulantSpectrum) {
  // Единая форма линии: с «прозрачной» НЧ спектр системы совпадает с кумулянтным спектром агрегата
  // (те же реализации беспорядка), т.е. σ_ext = 4πk·Im α  ↔  σ_abs(кумулянт).
  plasmonics::DrudeParams dp;
  dp.eps_inf = 1.33 * 1.33;
  dp.omega_p_eV = 0.0;
  const auto mat = plasmonics::Material::drude(dp);
  Aggregate a = build_geometry("ring", {{"N", 10}, {"radius", 2.0}});
  transform(a, Vec3(0, 0, 25.0), 0, 0, 0);
  ExcitonModel m;
  m.monomer_energy_cm = 17000.0;
  m.disorder_sigma_cm = 100.0;
  m.n_medium = 1.33;
  m.coupling.screening = 1.0 / (1.33 * 1.33);
  m.bath.lambda_cm = 60.0;
  m.bath.cutoff_cm = 100.0;
  m.decay.nonradiative_ps = 1e-3;
  SpectrumGrid g{15500.0, 18500.0, 600, 80.0};
  EnsembleSettings s;
  s.min_realizations = s.max_realizations = 16;
  s.extra_wmin_cm = 15500.0;
  s.extra_wmax_cm = 18500.0;
  const EnsembleResult er = run_ensemble(a, m, g, s, 16);
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = 10.0;
  spec.l_max = 6;
  spec.n_medium = 1.33;
  const PlexcitonSystem sys(a, spec, mat, m.coupling.mu_D);
  PlexcitonSpectra mean, se;
  sys.spectra(er.stored, er.engine.get(), -1.0, er.w, mean, se);
  double peak = 0.0, maxdiff = 0.0;
  for (std::size_t q = 0; q < er.w.size(); ++q) {
    peak = std::max(peak, er.abs_mean[q]);
    maxdiff = std::max(maxdiff, std::abs(mean.ext[q] - er.abs_mean[q]));
  }
  EXPECT_LT(maxdiff / peak, 0.01);
}

TEST(Plexciton, DipoleModeCouplingMatchesDrudeSphereAnalytics) {
  // Диполь на оси друдевской сферы (квазистатика): α = α_bg + A/(ω_r² − ω² − iγω),
  //   ω_r² = ω_p²/(ε∞ + 2ε_m), A = 3ε_m a³ ω_r²/(ε∞ + 2ε_m);  W^(l=1) = −α·4·5.034μ²/(ε_m r⁶) ⇒ g² = |K|·A/(2ω_r).
  const double eps_inf = 4.0, wp_eV = 9.0, gam_eV = 0.1, n = 1.33, em = n * n, a = 8.0, r = 12.0, mu = 10.0;
  plasmonics::DrudeParams dp;
  dp.eps_inf = eps_inf;
  dp.omega_p_eV = wp_eV;
  dp.gamma_eV = gam_eV;
  const auto mat = plasmonics::Material::drude(dp);
  Aggregate agg;
  agg.positions.emplace_back(0, 0, r);
  agg.dipoles.emplace_back(0, 0, 1);
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = a;
  spec.l_max = 4;
  spec.n_medium = n;
  spec.dynamic_dipole = false;
  const PlexcitonSystem sys(agg, spec, mat, mu);
  const double cm_per_eV = units::cm_per_eV;
  const double wr = wp_eV * cm_per_eV / std::sqrt(eps_inf + 2 * em);
  std::vector<double> w;
  for (int i = 0; i <= 400; ++i) w.push_back(wr - 3000.0 + 15.0 * i);
  const auto c = sys.dipole_mode_coupling(w);
  const double A = 3 * em * a * a * a * wr * wr / (eps_inf + 2 * em);
  const double K = 4.0 * units::dd_cm_per_D2_nm3 * mu * mu / (em * std::pow(r, 6));
  EXPECT_NEAR(c.omega_axis[2] / wr, 1.0, 2e-3);
  EXPECT_NEAR(c.gamma_axis[2] / (gam_eV * cm_per_eV), 1.0, 0.03);
  EXPECT_NEAR(c.g_cm[2] / std::sqrt(K * A / (2 * wr)), 1.0, 0.02);
  EXPECT_NEAR(c.g_cm[0], 0.0, 1e-6 * c.g_cm[2]);  // радиальный диполь не связан с поперечной (x) модой
}

TEST(Plexciton, PurcellRadialDipoleMatchesQuasistaticDipoleLimit) {
  // Далеко от сферы (вакуум): μ_eff = μ(1 + 2α/r³), α = a³(ε−1)/(ε+2); высшие мультиполи ~ (a/r)^5.
  plasmonics::DrudeParams dp;
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto mat = plasmonics::Material::drude(dp);
  Aggregate a;
  const double r = 40.0, R = 10.0;
  a.positions.emplace_back(0, 0, r);
  a.dipoles.emplace_back(0, 0, 1);
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = R;
  spec.l_max = 30;
  spec.n_medium = 1.0;
  spec.dynamic_dipole = false;
  const PlexcitonSystem sys(a, spec, mat, 10.0);
  const double nu = 15000.0;
  const auto nf = sys.near_field(nu, 0.1, 0.05);
  const std::complex<double> eps = mat.epsilon(nu);
  const std::complex<double> alpha = R * R * R * (eps - 1.0) / (eps + 2.0);
  const double expect = std::norm(1.0 + 2.0 * alpha / (r * r * r));
  EXPECT_NEAR(nf[0].rad_enhancement / expect, 1.0, 2e-3);
  EXPECT_GT(nf[0].quantum_yield, 0.0);
  EXPECT_NEAR(nf[0].quantum_yield_free, 0.1 / 0.15, 1e-12);
}

TEST(Plexciton, HybridEmissionTransparentAndQuenched) {
  plasmonics::DrudeParams dp;
  dp.eps_inf = 1.0;
  dp.omega_p_eV = 0.0;
  const auto clear = plasmonics::Material::drude(dp);
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto gold = plasmonics::Material::drude(dp);
  Aggregate a = build_geometry("ring", {{"N", 6}, {"radius", 2.0}});
  CouplingParams cp;
  const Eigen::MatrixXd J = coupling_matrix(a, cp);
  std::mt19937_64 rng(5);
  const Eigen::MatrixXd H = frenkel_hamiltonian(J, 17000.0, 0.0, rng);
  std::vector<double> w;
  for (int k = 0; k < 41; ++k) w.push_back(16000.0 + 50.0 * k);
  auto run = [&](const plasmonics::Material& m, double z) {
    Aggregate b = a;
    transform(b, Vec3(0, 0, z), 0, 0, 0);
    NanoparticleSpec spec;
    spec.a_xy_nm = spec.a_z_nm = 10.0;
    spec.l_max = 40;
    spec.n_medium = 1.0;
    const PlexcitonSystem sys(b, spec, m, 10.0);
    Realization Rz;
    Rz.H = H;
    Rz.states = diagonalize(H, b, 10.0);
    const auto n = Rz.states.size();
    Rz.K_ps = Eigen::MatrixXd::Zero(n, n);
    Rz.hwhm_cm = Eigen::VectorXd::Constant(n, 50.0);
    Rz.gamma_rad_ps = Rz.states.mu2 / 100.0 * 3e-4;  // ~0.3 нс⁻¹ на 10 Д
    Bath bth;
    bth.lambda_cm = 5.0;
    const LineshapeEngine eng(bth, 5000.0, 20.0);
    return sys.emission({Rz}, eng, w, 1e-4, false, 0.0);
  };
  const auto t = run(clear, 15.0);
  EXPECT_NEAR(t.qy / t.qy_free, 1.0, 1e-9);
  for (std::size_t q = 0; q < w.size(); ++q) EXPECT_NEAR(t.pl[q], t.pl_free[q], 1e-9 * (1e-300 + t.pl_free[q]) + 1e-30);
  const auto near = run(gold, 12.0), far = run(gold, 25.0);
  std::printf("QY free %.3f  far %.3f  near %.3f\n", near.qy_free, far.qy, near.qy);
  EXPECT_LT(near.qy, far.qy);
  EXPECT_LT(near.lifetime_ps, near.lifetime_free_ps);
}

TEST(Plexciton, FastSolverMatchesExact) {
  // Оболочка из ~300 молекул у сферы: низкие l точно + диагональ старших против полного решения
  plasmonics::DrudeParams dp;
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto gold = plasmonics::Material::drude(dp);
  Aggregate a = build_geometry("particle_shell", {{"a_xy", 10.0}, {"a_z", 10.0}, {"gap", 1.5}, {"N", 300},
                                                  {"dipole_mode", 1}, {"seed", 3}});
  CouplingParams cp;
  const Eigen::MatrixXd J = coupling_matrix(a, cp);
  std::mt19937_64 rng(7);
  const Eigen::MatrixXd H = frenkel_hamiltonian(J, 19000.0, 300.0, rng);
  Realization Rz;
  Rz.H = H;
  Rz.states = diagonalize(H, a, 8.0);
  const auto n = Rz.states.size();
  Rz.K_ps = Eigen::MatrixXd::Zero(n, n);
  Rz.hwhm_cm = Eigen::VectorXd::Constant(n, 150.0);
  std::vector<double> w;
  for (int k = 0; k < 31; ++k) w.push_back(16000.0 + 200.0 * k);
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = 10.0;
  spec.l_max = 30;
  spec.n_medium = 1.33;
  PlexcitonSpectra ex, fa, se;
  PlexcitonSystem(a, spec, gold, 8.0).spectra({Rz}, nullptr, -1.0, w, ex, se);
  spec.fast_l = 8;
  PlexcitonSystem(a, spec, gold, 8.0).spectra({Rz}, nullptr, -1.0, w, fa, se);
  double mx = 0, err = 0;
  for (std::size_t q = 0; q < w.size(); ++q) {
    mx = std::max(mx, std::abs(ex.ext[q] - ex.np_ext[q]));  // вклад молекул (и их влияния на НЧ)
    err = std::max(err, std::abs(ex.ext[q] - fa.ext[q]));
  }
  // Сходимость по порогу: fast_l = 16 точнее, чем 8
  spec.fast_l = 16;
  PlexcitonSpectra f16;
  PlexcitonSystem(a, spec, gold, 8.0).spectra({Rz}, nullptr, -1.0, w, f16, se);
  double err16 = 0;
  for (std::size_t q = 0; q < w.size(); ++q) err16 = std::max(err16, std::abs(ex.ext[q] - f16.ext[q]));
  std::printf("fast solver: max |Δσ_ext| = %.3g (l ≤ 8), %.3g (l ≤ 16) of molecular part %.3g nm²\n", err, err16, mx);
  EXPECT_LT(err / mx, 0.02);
  EXPECT_LT(err16, err);
}

TEST(Plexciton, RadialDipoleSelfFieldMatchesGerstenNitzanSeries) {
  // W_ii = pref·Σ_l (l+1)² D̂_l b^{2l+1}/r^{2l+4} — все мультиполи, вплоть до l ~ 80 (зазор 1.5 нм)
  plasmonics::DrudeParams dp;
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto gold = plasmonics::Material::drude(dp);
  for (const double shell : {0.0, 1.0}) {
    Aggregate a;
    const double R = 10.0, b = R + shell, r = b + 1.5;
    a.positions.emplace_back(0, 0, r);
    a.dipoles.emplace_back(0, 0, 1);
    NanoparticleSpec spec;
    spec.a_xy_nm = spec.a_z_nm = R;
    spec.shell_thickness_nm = shell;
    spec.shell_n = 1.5;
    spec.l_max = 90;
    spec.n_medium = 1.33;
    spec.dynamic_dipole = false;
    const PlexcitonSystem sys(a, spec, gold, 10.0);
    const double nu = 19000.0, em = 1.33 * 1.33;
    plasmonics::SphereLayers L;
    L.core_radius_nm = R;
    if (shell > 0) {
      L.shell_thickness_nm = {shell};
      L.shell_eps = {cd(2.25)};
    }
    cd sum = 0.0;
    for (int l = 1; l <= 90; ++l)
      sum += (l + 1.0) * (l + 1.0) * plasmonics::layered_response(l, gold.epsilon(nu), L, em) *
             std::pow(b / r, 2 * l + 1) / std::pow(r, 3);
    const cd W = units::dd_cm_per_D2_nm3 * 100.0 / em * sum;
    const auto nf = sys.near_field(nu);
    EXPECT_NEAR(nf[0].shift_cm / W.real(), 1.0, 1e-6) << "shell " << shell;
    EXPECT_NEAR(nf[0].decay_ps / (-2.0 * W.imag() * units::rad_ps_per_cm), 1.0, 1e-6) << "shell " << shell;
  }
}

TEST(Plexciton, DecayRatesMatchExactMieForSmallSpheres) {
  // Модель (квазистатические мультиполи + Ми для l ≤ 8 + запаздывание дипольной связи) против точной теории Ми
  // для диполя у сферы R = 10 нм в воде: зазоры 0.5–20 нм, обе ориентации, λ = 520 (у резонанса) и 650 нм.
  // Допуски — измеренная точность модели: Γ_tot 3 %, Γ_rad 10 % (у резонанса друдевской сферы).
  plasmonics::DrudeParams dp;
  dp.eps_inf = 9.5;
  dp.omega_p_eV = 9.0;
  dp.gamma_eV = 0.07;
  const auto gold = plasmonics::Material::drude(dp);
  const double R = 10.0, nm = 1.33;
  const std::vector<double> gaps = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};
  Aggregate a;
  for (double g : gaps) {
    a.positions.emplace_back(0, 0, R + g);
    a.dipoles.emplace_back(0, 0, 1);
    a.positions.emplace_back(0, 0, -(R + g));
    a.dipoles.emplace_back(1, 0, 0);
  }
  NanoparticleSpec spec;
  spec.a_xy_nm = spec.a_z_nm = R;
  spec.l_max = 120;
  spec.n_medium = nm;
  const PlexcitonSystem sys(a, spec, gold, 10.0);
  for (const double lam : {520.0, 650.0}) {
    const double nu = 1e7 / lam;
    ExcitonStates m1;
    m1.energies = Eigen::VectorXd::Constant(1, nu);
    m1.mu2 = Eigen::VectorXd::Constant(1, 100.0);
    DecayParams d;
    d.n_medium = nm;
    d.mu_D = 10.0;
    d.monomer_energy_cm = nu;
    const double g0 = radiative_rates_ps(m1, d)(0);
    const auto nf = sys.near_field(nu, g0, 0.0);
    for (std::size_t k = 0; k < gaps.size(); ++k) {
      const auto ex = plasmonics::dipole_near_sphere(R, gold.epsilon(nu), nm, nu, R + gaps[k]);
      for (int o = 0; o < 2; ++o) {
        const auto& q = nf[2 * k + static_cast<std::size_t>(o)];
        const double et = o == 0 ? ex.total_perp : ex.total_par, er = o == 0 ? ex.rad_perp : ex.rad_par;
        EXPECT_NEAR((1.0 + q.decay_ps / g0) / et, 1.0, 0.03) << lam << " nm, gap " << gaps[k] << ", o " << o;
        if (er > 0.2) EXPECT_NEAR(q.rad_enhancement / er, 1.0, 0.1) << lam << " nm, gap " << gaps[k] << ", o " << o;
        else EXPECT_NEAR(q.rad_enhancement, er, 0.1);
      }
    }
  }
}
