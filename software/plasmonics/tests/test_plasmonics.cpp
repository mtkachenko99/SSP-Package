// Верификация модуля плазмоники: Ми, квазистатика сферы и сфероидов, материалы.
#include "ssp/plasmonics/material.hpp"
#include "ssp/plasmonics/mie.hpp"
#include "ssp/plasmonics/particle.hpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace ssp::plasmonics;

TEST(Mie, RayleighLimitMatchesQuasistaticPolarizability) {
  // x ≪ 1: α = (3i/2k³)a₁ → a³(m²−1)/(m²+2)
  const cd m(1.5, 0.1);
  const double a = 1.0, k = 1e-3;
  const MieResult r = mie(k * a, m);
  const cd alpha = mie_dipole_polarizability(r, k);
  const cd expect = a * a * a * (m * m - 1.0) / (m * m + 2.0);
  EXPECT_NEAR(std::abs(alpha - expect) / std::abs(expect), 0.0, 1e-5);
}

TEST(Mie, LosslessSphereHasNoAbsorption) {
  const MieResult r = mie(5.0, cd(1.5, 0.0));
  EXPECT_NEAR(r.q_abs / r.q_ext, 0.0, 1e-12);
  EXPECT_GT(r.q_ext, 2.0);  // крупная частица: Q_ext порядка 2–4
}

TEST(Mie, LargeSphereExtinctionParadox) {
  // Q_ext → 2 при x → ∞ (сильно поглощающая сфера); сходимость медленная, ~x^{−2/3}
  const double q1 = mie(400.0, cd(1.5, 1.0)).q_ext, q2 = mie(3200.0, cd(1.5, 1.0)).q_ext;
  EXPECT_GT(q1, 2.0);
  EXPECT_LT(q2 - 2.0, 0.5 * (q1 - 2.0));  // остаток уменьшается ~4 раза при x×8
  EXPECT_NEAR(q2, 2.0, 0.015);
}

TEST(Legendre, NormalizedOrthonormal) {
  // ∫₋₁¹ P̄_l^m P̄_l'^m dη = δ_ll' (численно, Гаусс–Лежандр на равномерной сетке трапеций)
  const int L = 12, m = 3, n = 20001;
  std::vector<double> buf;
  std::vector<std::vector<double>> P(n);
  for (int k = 0; k < n; ++k) {
    normalized_legendre(L, m, -1.0 + 2.0 * k / (n - 1), buf);
    P[k] = buf;
  }
  for (int l1 = m; l1 <= L; ++l1)
    for (int l2 = m; l2 <= L; ++l2) {
      double s = 0.0;
      for (int k = 0; k < n; ++k) s += (k == 0 || k == n - 1 ? 0.5 : 1.0) * P[k][l1] * P[k][l2];
      s *= 2.0 / (n - 1);
      EXPECT_NEAR(s, l1 == l2 ? 1.0 : 0.0, 2e-4);
    }
}

namespace {
void check_free_expansion(const QuasistaticParticle& p, const Vec3& a, const Vec3& b, double tol) {
  const cd S = p.free_kernel_expansion(a, b, p.l_max());
  const double exact = 1.0 / (a - b).norm();
  EXPECT_NEAR(S.real(), exact, tol * exact) << "форма " << int(p.shape());
  EXPECT_NEAR(S.imag(), 0.0, tol * exact);
}
}  // namespace

TEST(Quasistatic, SphereLaplaceExpansion) {
  QuasistaticParticle p(10.0, 10.0, 60);
  check_free_expansion(p, Vec3(3.0, 11.0, 2.0), Vec3(-15.0, 4.0, 9.0), 1e-9);
}

TEST(Quasistatic, ProlateLaplaceExpansion) {
  QuasistaticParticle p(8.0, 20.0, 60);
  check_free_expansion(p, Vec3(9.0, 1.0, 5.0), Vec3(-4.0, 12.0, 25.0), 1e-8);
  check_free_expansion(p, Vec3(0.5, 0.3, 21.0), Vec3(14.0, -6.0, -3.0), 1e-8);
}

TEST(Quasistatic, OblateLaplaceExpansion) {
  QuasistaticParticle p(20.0, 8.0, 60);
  check_free_expansion(p, Vec3(9.0, 1.0, 9.0), Vec3(-24.0, 12.0, 5.0), 1e-8);
  check_free_expansion(p, Vec3(21.0, 0.5, 0.2), Vec3(3.0, -6.0, -19.0), 1e-8);
}

TEST(Quasistatic, SphereKernelMatchesImageSeries) {
  // K = −Σ α_l a^{2l+1} P_l(cos γ)/(r r')^{l+1}
  const double a = 10.0, em = 1.77;
  const cd eps(-5.0, 1.2);
  QuasistaticParticle p(a, a, 80);
  const Vec3 r1(12.0, 1.0, 3.0), r2(-2.0, 13.0, 1.0);
  const double c = r1.normalized().dot(r2.normalized());
  cd K = 0.0;
  double Pm = 1.0, P = c;  // P_0, P_1
  for (int l = 1; l <= 80; ++l) {
    if (l > 1) {
      const double Pn = ((2.0 * l - 1.0) * c * P - (l - 1.0) * Pm) / l;
      Pm = P;
      P = Pn;
    }
    const cd al = double(l) * (eps - em) / (double(l) * eps + (l + 1.0) * em);
    K -= al * std::pow(a, 2 * l + 1) * P / std::pow(r1.norm() * r2.norm(), l + 1);
  }
  const cd Kp = p.induced_kernel(r1, r2, eps, em);
  EXPECT_NEAR(std::abs(Kp - K) / std::abs(K), 0.0, 1e-10);
}

TEST(Quasistatic, SpheroidTendsToSphere) {
  const double em = 1.77;
  const cd eps(-6.0, 0.8);
  QuasistaticParticle s(10.0, 10.0, 30), pr(10.0, 10.02, 30), ob(10.0, 9.98, 30);
  const Vec3 r1(3.0, 11.0, 4.0), r2(-6.0, 2.0, 12.0);
  const cd Ks = s.induced_kernel(r1, r2, eps, em);
  EXPECT_NEAR(std::abs(pr.induced_kernel(r1, r2, eps, em) - Ks) / std::abs(Ks), 0.0, 1e-2);
  EXPECT_NEAR(std::abs(ob.induced_kernel(r1, r2, eps, em) - Ks) / std::abs(Ks), 0.0, 1e-2);
}

namespace {
// Индуцированное поле однородного E₀ вдоль оси a на большом расстоянии → поле диполя α.
void check_dipole_far_field(const QuasistaticParticle& p, int axis) {
  const double em = 1.77;
  const cd eps(-4.0, 0.5);
  const Eigen::Vector3cd alpha = p.polarizability(eps, em);
  const auto& terms = p.terms();
  const int t = p.uniform_term()[axis];
  const cd D = p.response(terms[t], eps, em);
  const cd u = p.uniform_coeff()[axis];
  // точка на оси a, далеко: E_ind = 2α/r³ ê_a
  Vec3 r = Vec3::Zero();
  r(axis) = 600.0;
  std::vector<Eigen::Vector3cd> gc, gs;
  p.gradients(r, gc, gs);
  const Eigen::Vector3cd g = p.uniform_is_sine()[axis] ? gs[t] : gc[t];
  const cd E = u * D * g(axis);
  const cd expect = 2.0 * alpha(axis) / std::pow(600.0, 3);
  EXPECT_NEAR(std::abs(E - expect) / std::abs(expect), 0.0, 5e-3) << "ось " << axis << " форма " << int(p.shape());
}
}  // namespace

TEST(Quasistatic, DipoleResponseMatchesDepolarizationFactors) {
  for (auto p : {QuasistaticParticle(10.0, 10.0, 4), QuasistaticParticle(8.0, 20.0, 4),
                 QuasistaticParticle(20.0, 8.0, 4)})
    for (int a = 0; a < 3; ++a) check_dipole_far_field(p, a);
}

TEST(Quasistatic, SphereSurfaceFieldOnAxis) {
  // E на оси z вне сферы: E₀(1 + 2α/r³), α = a³(ε−ε_m)/(ε+2ε_m)
  const double a = 10.0, em = 1.77, r = 12.0;
  const cd eps(-5.0, 0.3);
  QuasistaticParticle p(a, a, 4);
  std::vector<Eigen::Vector3cd> gc, gs;
  p.gradients(Vec3(0, 0, r), gc, gs);
  const int t = p.uniform_term()[2];
  const cd E = 1.0 + p.uniform_coeff()[2] * p.response(p.terms()[t], eps, em) * gc[t](2);
  const cd alpha = a * a * a * (eps - em) / (eps + 2.0 * em);
  EXPECT_NEAR(std::abs(E - (1.0 + 2.0 * alpha / (r * r * r))), 0.0, 1e-7);
}

TEST(Quasistatic, DipoleSelfInteractionGerstenNitzan) {
  // Самодействие диполя у сферы (квазистатика): W = (μ·∇)(μ·∇')K_ind,
  //   W_⊥ = −Σ α_l a^{2l+1}(l+1)²/r^{2l+4},  W_∥ = −½Σ α_l a^{2l+1} l(l+1)/r^{2l+4}.
  const double a = 10.0, r = 12.0, em = 1.77;
  const cd eps(-6.0, 0.7);
  const int L = 120;
  QuasistaticParticle p(a, a, L);
  std::vector<Eigen::Vector3cd> gc, gs;
  p.gradients(Vec3(0, 0, r), gc, gs);
  auto self = [&](const Vec3& mu) {
    cd W = 0.0;
    for (std::size_t t = 0; t < p.terms().size(); ++t) {
      const cd D = p.response(p.terms()[t], eps, em);
      const cd vc = mu.cast<cd>().dot(gc[t]), vs = mu.cast<cd>().dot(gs[t]);
      W += p.terms()[t].w * D * (vc * vc + vs * vs);
    }
    return W;
  };
  cd Wp = 0.0, Wt = 0.0;
  for (int l = 1; l <= L; ++l) {
    const cd al = double(l) * (eps - em) / (double(l) * eps + (l + 1.0) * em);
    const double g = std::pow(a, 2 * l + 1) / std::pow(r, 2 * l + 4);
    Wp -= al * g * (l + 1.0) * (l + 1.0);
    Wt -= 0.5 * al * g * double(l) * (l + 1.0);
  }
  // Точность ограничена конечными разностями 4-го порядка в градиентах (~1e-6).
  EXPECT_NEAR(std::abs(self(Vec3(0, 0, 1)) - Wp) / std::abs(Wp), 0.0, 1e-5);
  EXPECT_NEAR(std::abs(self(Vec3(1, 0, 0)) - Wt) / std::abs(Wt), 0.0, 1e-5);
  // Проводящая сфера притягивает диполь (сдвиг вниз), поглощающая — даёт распад: Im W < 0.
  EXPECT_LT(Wp.imag(), 0.0);
}

TEST(Material, JohnsonChristyGold) {
  const Material au = Material::from_refractiveindex_yaml(std::string(SSP_DATA_DIR) + "/main/Au/nk/Johnson.yml");
  // λ = 0.5209 мкм: n = 0.62, k = 2.081 (табличная точка J&C)
  const cd e = au.epsilon(1e4 / 0.5209);
  const cd expect = cd(0.62, 2.081) * cd(0.62, 2.081);
  EXPECT_NEAR(std::abs(e - expect), 0.0, 1e-9);
  EXPECT_THROW(au.epsilon(1e4 / 5.0), std::out_of_range);
}

// ---------------------------------------------------------------- слоистая сфера, d-параметры, покрытая Ми
#include "ssp/plasmonics/layered.hpp"
#include "ssp/plasmonics/dipole_sphere.hpp"

TEST(Layered, ReducesToHomogeneousSphere) {
  const cd e(-7.0, 1.1);
  const double em = 1.77;
  SphereLayers s;
  s.core_radius_nm = 10.0;
  for (int l = 1; l <= 60; ++l) {  // большие l — проверка обусловленности системы
    const cd hom = -double(l) * (e - em) / (double(l) * e + (l + 1.0) * em);
    EXPECT_NEAR(std::abs(layered_response(l, e, s, em) - hom), 0.0, 1e-12) << "l = " << l;
    // оболочка с ε = ε_m: отклик ядра, нормированный на внешний радиус
    SphereLayers s1 = s;
    s1.shell_thickness_nm = {2.0};
    s1.shell_eps = {cd(em)};
    EXPECT_NEAR(std::abs(layered_response(l, e, s1, em) / (hom * std::pow(10.0 / 12.0, 2 * l + 1)) - 1.0), 0.0, 1e-6);
    // оболочка с ε = ε_ядра: однородная сфера радиуса 12
    SphereLayers s2 = s;
    s2.shell_thickness_nm = {2.0};
    s2.shell_eps = {e};
    EXPECT_NEAR(std::abs(layered_response(l, e, s2, em) - hom), 0.0, 1e-12);
  }
}

TEST(Layered, FeibelmanSpillOutRedshiftsLinearly) {
  // Друде (ε∞ = 1), сфера 2 нм в вакууме: d⊥ > 0 (spill-out) — красный сдвиг дипольного резонанса, линейный по d⊥.
  const double wp = 9.0, gam = 0.05, em = 1.0;
  auto peak = [&](cd d) {
    SphereLayers s;
    s.core_radius_nm = 2.0;
    s.d_perp_nm = d;
    double best = 0, wbest = 0;
    for (double w = 4.0; w < 6.5; w += 1e-4) {
      const cd e = 1.0 - wp * wp / cd(w * w, gam * w);
      const double v = (-layered_response(1, e, s, em)).imag();
      if (v > best) {
        best = v;
        wbest = w;
      }
    }
    return wbest;
  };
  const double w0 = peak(0.0), w1 = peak(0.05), w2 = peak(0.10);
  EXPECT_NEAR(w0, wp / std::sqrt(3.0), 2e-3);
  EXPECT_LT(w1, w0);
  EXPECT_NEAR((w0 - w2) / (w0 - w1), 2.0, 0.1);
}

TEST(CoatedMie, LimitsMatchHomogeneousMie) {
  const cd mc(0.4, 3.0), ms(1.45, 0.0);
  const double x = 0.6, y = 0.75;
  const MieResult c1 = mie_coated(x, y, mc, cd(1.0));  // оболочка = среда
  const MieResult h1 = mie(x, mc);
  EXPECT_NEAR(c1.q_ext * y * y, h1.q_ext * x * x, 1e-9 * h1.q_ext * x * x);
  const MieResult c2 = mie_coated(x, y, mc, mc);  // оболочка = ядро
  const MieResult h2 = mie(y, mc);
  EXPECT_NEAR(c2.q_ext, h2.q_ext, 1e-9 * h2.q_ext);
  EXPECT_NEAR(c2.q_sca, h2.q_sca, 1e-9 * h2.q_sca);
  (void)ms;
}

TEST(MieDressing, DipoleFactorReproducesMieA1) {
  SphereLayers s;
  s.core_radius_nm = 25.0;
  const cd e(-6.0, 2.0);
  const double n = 1.33, nu = 18500.0, k = 2.0 * M_PI * n * nu * 1e-7;
  const auto f = mie_dynamic_factors(s, e, n, nu, 4);
  const cd a_qs = -layered_response(1, e, s, n * n) * std::pow(25.0, 3);
  const cd a_mie = mie_dipole_polarizability(mie_sphere(25.0, e, n, nu), k);
  EXPECT_NEAR(std::abs(f[1] * a_qs - a_mie) / std::abs(a_mie), 0.0, 1e-6);
  EXPECT_GT(std::abs(f[1] - 1.0), 0.05);  // для 50 нм сферы динамика существенна
}

TEST(MieDressing, HigherMultipoleFactorsTendToOneForSmallSpheres) {
  // f_l = a_l / a_l^{qs} → 1 при x → 0 для всех l (в т.ч. с оболочкой), отклонение ~ x²
  const cd e(-7.0, 1.1);
  for (const bool shell : {false, true}) {
    SphereLayers s;
    s.core_radius_nm = 2.0;
    if (shell) {
      s.shell_thickness_nm = {0.5};
      s.shell_eps = {cd(2.1)};
    }
    const auto f = mie_dynamic_factors(s, e, 1.33, 12000.0, 8);
    for (int l = 1; l <= 8; ++l) EXPECT_NEAR(std::abs(f[l] - 1.0), 0.0, 2e-2) << "l = " << l << " shell " << shell;
  }
}

TEST(DipoleSphere, LosslessSphereConservesEnergyAndFarLimit) {
  // Непоглощающая сфера: вся энергия излучается (Γ = Γ_rad); далеко от сферы Γ → Γ₀
  for (const double r : {21.0, 25.0, 40.0, 120.0}) {
    const auto R = dipole_near_sphere(20.0, {2.25, 0.0}, 1.33, 16000.0, r);
    EXPECT_NEAR(R.total_perp / R.rad_perp, 1.0, 1e-8) << r;
    EXPECT_NEAR(R.total_par / R.rad_par, 1.0, 1e-8) << r;
  }
  const auto F = dipole_near_sphere(20.0, {-10.0, 1.0}, 1.33, 16000.0, 3000.0);
  EXPECT_NEAR(F.total_perp, 1.0, 2e-3);
  EXPECT_NEAR(F.rad_par, 1.0, 2e-3);
}

TEST(DipoleSphere, SmallSphereRadiativeRateMatchesDipoleLimit) {
  // a ≪ λ, r = 4a: Γ_rad⊥/Γ₀ ≈ |1 + 2α/r³|², Γ_rad∥/Γ₀ ≈ |1 − α/r³|², α = a³(m²−1)/(m²+2)
  const cd eps(-12.0, 1.2);
  const double a = 2.0, r = 8.0, nm = 1.0;
  const auto R = dipole_near_sphere(a, eps, nm, 15000.0, r);
  const cd al = a * a * a * (eps - 1.0) / (eps + 2.0);
  EXPECT_NEAR(R.rad_perp / std::norm(1.0 + 2.0 * al / (r * r * r)), 1.0, 1e-2);
  EXPECT_NEAR(R.rad_par / std::norm(1.0 - al / (r * r * r)), 1.0, 1e-2);
  EXPECT_GT(R.total_perp, R.rad_perp);  // поглощение в металле
}
