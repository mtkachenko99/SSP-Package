// Слой на поверхности НЧ и анализ ветвей скана.
#include "scan.hpp"

#include "ssp/excitons/geometry.hpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace ssp;

TEST(ParticleShell, ConstantGapAlongNormalAndUniformCoverage) {
  const double a = 10.0, c = 20.0, gap = 1.5;
  const auto agg = excitons::build_geometry("particle_shell", {{"a_xy", a}, {"a_z", c}, {"gap", gap}, {"N", 2000}});
  ASSERT_EQ(agg.size(), 2000u);
  int upper = 0;
  for (std::size_t i = 0; i < agg.size(); ++i) {
    const auto& r = agg.positions[i];
    // Расстояние до поверхности по нормали: найти ближайшую точку сфероида (минимизация по параметрическим углам)
    // ближайшая точка лежит в той же меридиональной плоскости: грубый поиск по θ ∈ [0, π] + уточнение
    const double ph = std::atan2(r.y(), r.x());
    auto dist = [&](double th) {
      const excitons::Vec3 s(a * std::sin(th) * std::cos(ph), a * std::sin(th) * std::sin(ph), c * std::cos(th));
      return (r - s).norm();
    };
    double tb = 0.0, best = 1e300;
    for (int k = 0; k <= 4000; ++k) {
      const double th = M_PI * k / 4000;
      if (dist(th) < best) {
        best = dist(th);
        tb = th;
      }
    }
    double lo = std::max(0.0, tb - M_PI / 4000), hi = std::min(M_PI, tb + M_PI / 4000);
    for (int it = 0; it < 60; ++it) {
      const double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
      if (dist(m1) < dist(m2)) hi = m2;
      else lo = m1;
    }
    best = dist(0.5 * (lo + hi));
    EXPECT_NEAR(best, gap, 1e-6);
    EXPECT_NEAR(agg.dipoles[i].norm(), 1.0, 1e-12);
    if (r.z() > 0) ++upper;
  }
  EXPECT_NEAR(upper, 1000, 5);  // симметрия по площади
  EXPECT_GT(excitons::min_pair_distance(agg), 0.5);
}

TEST(ParticleShell, CapCoversOnlyPolarRegion) {
  const auto agg = excitons::build_geometry("particle_shell", {{"N", 300}, {"theta_max_deg", 40.0}});
  for (const auto& r : agg.positions) EXPECT_GT(r.z() / r.norm(), std::cos(40.5 * M_PI / 180.0));
}

TEST(CoupledOscillators, FitRecoversCouplingFromSyntheticBranches) {
  const double g = 180.0, wx = 17000.0, gp = 900.0, gx = 80.0;
  std::vector<double> wp, gpv, wxv, gxv;
  std::vector<std::vector<double>> peaks;
  for (int k = 0; k < 15; ++k) {
    const double w = 16000.0 + 150.0 * k;
    double lo, up;
    app::coupled_branches(w, gp, wx, gx, g, lo, up);
    wp.push_back(w);
    gpv.push_back(gp);
    wxv.push_back(wx);
    gxv.push_back(gx);
    peaks.push_back({lo, up});
  }
  const auto f = app::fit_coupled_oscillators(wp, gpv, wxv, gxv, peaks, "abs");
  EXPECT_NEAR(f.g_cm, g, 0.5);
  EXPECT_LT(f.rms_cm, 0.1);
  EXPECT_FALSE(f.visible);  // 2g = 360 < (900 + 80)/2 = 490
  EXPECT_FALSE(f.anticrossing);  // g = 180 < |γp − γx|/4 = 205: исключительная точка не пройдена
  EXPECT_NEAR(f.cooperativity, 4.0 * g * g / (gp * gx), 1e-6 * f.cooperativity);
}

#include "ssp/excitons/ensemble.hpp"

TEST(ParticleShell, DensityAndDomainOrder) {
  const double a = 20.0, gap = 1.0, rho = 0.5;
  const auto agg = excitons::build_geometry("particle_shell", {{"a_xy", a}, {"a_z", a}, {"gap", gap}, {"density_per_nm2", rho}});
  const double area = 4.0 * M_PI * (a + gap) * (a + gap);
  EXPECT_NEAR(double(agg.size()), rho * area, 0.01 * rho * area);
  // ориентационный порядок ближайших соседей: ⟨|d_i·d_j|⟩ для r < 2 нм
  auto order = [](const excitons::Aggregate& g) {
    double s = 0;
    int n = 0;
    for (std::size_t i = 0; i < g.size(); ++i)
      for (std::size_t j = i + 1; j < g.size(); ++j)
        if ((g.positions[i] - g.positions[j]).norm() < 2.0) {
          s += std::abs(g.dipoles[i].dot(g.dipoles[j]));
          ++n;
        }
    return s / n;
  };
  const auto dom = excitons::build_geometry(
      "particle_shell", {{"a_xy", a}, {"a_z", a}, {"gap", gap}, {"density_per_nm2", rho}, {"domains", 30}, {"domain_jitter_deg", 5}});
  EXPECT_LT(order(agg), 0.75);  // случайная касательная ориентация: ⟨|cos|⟩ ≈ 2/π
  EXPECT_GT(order(dom), 0.9);  // внутри доменов почти параллельны (пары через границы доменов снижают среднее)
}

TEST(Disorder, CorrelatedSiteEnergiesFollowExponentialCovariance) {
  excitons::Aggregate agg = excitons::build_geometry("linear_chain", {{"N", 6}, {"spacing", 1.0}});
  excitons::ExcitonModel m;
  m.disorder_sigma_cm = 100.0;
  m.disorder_corr_length_nm = 2.0;
  const auto N = static_cast<Eigen::Index>(agg.size());
  Eigen::MatrixXd C(N, N);
  for (Eigen::Index i = 0; i < N; ++i)
    for (Eigen::Index j = 0; j < N; ++j) C(i, j) = 1e4 * std::exp(-std::abs(double(i - j)) / 2.0);
  const Eigen::MatrixXd L = Eigen::LLT<Eigen::MatrixXd>(C).matrixL();
  const Eigen::MatrixXd J = excitons::coupling_matrix(agg, m.coupling);
  double s01 = 0, s03 = 0, s00 = 0;
  const int R = 6000;
  for (int r = 0; r < R; ++r) {
    const auto Rz = excitons::make_realization(agg, J, m, 1000 + r, &L);
    const Eigen::VectorXd d = Rz.H.diagonal().array() - m.monomer_energy_cm;
    s00 += d(0) * d(0);
    s01 += d(0) * d(1);
    s03 += d(0) * d(3);
  }
  EXPECT_NEAR(s00 / R, 1e4, 600);
  EXPECT_NEAR(s01 / s00, std::exp(-0.5), 0.04);
  EXPECT_NEAR(s03 / s00, std::exp(-1.5), 0.04);
}
