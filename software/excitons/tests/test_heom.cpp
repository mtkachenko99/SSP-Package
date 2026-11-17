// HEOM: точное совпадение с кумулянтом для мономера и согласие с Редфилдом в режиме слабой связи.
#include "ssp/excitons/bath.hpp"
#include "ssp/excitons/heom.hpp"
#include "ssp/excitons/modified_redfield.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/excitons/spectra.hpp"
#include "ssp/units.hpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace ssp;
using namespace ssp::excitons;

TEST(Heom, MonomerAbsorptionEqualsCumulant) {
  HeomParams p;
  p.lambda_cm = 40.0;
  p.gamma_cm = 150.0;
  p.temperature_K = 300.0;
  p.matsubara = 3;
  p.depth = 10;
  p.t_max_ps = 1.5;
  Eigen::MatrixXd H(1, 1);
  H << 16000.0;
  Eigen::MatrixX3d mu(1, 3);
  mu << 1.0, 0.0, 0.0;
  std::vector<double> w;
  for (double x = 15000.0; x <= 17000.0; x += 5.0) w.push_back(x);
  const HeomAbsorption h = heom_absorption(H, mu, p, w);

  Bath b;
  b.lambda_cm = p.lambda_cm;
  b.cutoff_cm = p.gamma_cm;
  const LineshapeEngine eng(b, 5000.0, 20.0);  // длинное окно: разрешение ≪ ширины линии
  std::vector<double> c(w.size(), 0.0);
  eng.add_line(c, w, 16000.0, 1.0, 1.0, 0.0, false, 0);
  double peak = 0.0, maxdiff = 0.0;
  for (std::size_t i = 0; i < w.size(); ++i) {
    peak = std::max(peak, c[i]);
    maxdiff = std::max(maxdiff, std::abs(c[i] - h.lineshape[i]));
  }
  EXPECT_LT(maxdiff / peak, 1e-2);
}

TEST(Heom, DimerPopulationsWeakCouplingMatchRedfield) {
  // Слабая связь с баней (λ ≪ J, Δ): HEOM и секулярный Редфилд дают одинаковое равновесие и близкую кинетику.
  Eigen::MatrixXd H(2, 2);
  H << 17100.0, 120.0, 120.0, 16900.0;
  HeomParams p;
  p.lambda_cm = 5.0;
  p.gamma_cm = 100.0;
  p.temperature_K = 300.0;
  p.matsubara = 2;
  p.depth = 6;
  p.t_max_ps = 3.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  const Eigen::MatrixXd U = es.eigenvectors();
  Eigen::MatrixXcd rho0 = (U.col(1) * U.col(1).transpose()).cast<std::complex<double>>();
  const HeomPopulations hp = heom_populations(H, rho0, p, 60);

  ExcitonStates st;
  st.energies = es.eigenvalues();
  st.c = U;
  Bath b;
  b.lambda_cm = p.lambda_cm;
  b.cutoff_cm = p.gamma_cm;
  const Eigen::MatrixXd K = redfield_rate_matrix(st, b) * units::rad_ps_per_cm;
  Eigen::VectorXd t(static_cast<Eigen::Index>(hp.t_ps.size()));
  for (std::size_t k = 0; k < hp.t_ps.size(); ++k) t(static_cast<Eigen::Index>(k)) = hp.t_ps[k];
  Eigen::VectorXd P0(2);
  P0 << 0.0, 1.0;
  const Eigen::MatrixXd Pr = population_dynamics(K, Eigen::VectorXd::Zero(2), st.energies, b.beta(), P0, t);
  double maxdiff = 0.0;
  for (Eigen::Index k = 0; k < t.size(); ++k) maxdiff = std::max(maxdiff, std::abs(hp.exciton(0, k) - Pr(0, k)));
  EXPECT_LT(maxdiff, 0.03);
  EXPECT_NEAR(hp.exciton.col(t.size() - 1).sum(), 1.0, 1e-6);
}

namespace {
ExcitonStates dimer_states(const Eigen::MatrixXd& H) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  ExcitonStates st;
  st.energies = es.eigenvalues();
  st.c = es.eigenvectors();
  st.ipr = st.c.array().pow(4).colwise().sum().transpose();
  return st;
}
}  // namespace

TEST(ModifiedRedfield, WeakCouplingReducesToRedfield) {
  Eigen::MatrixXd H(2, 2);
  H << 17100.0, 120.0, 120.0, 16900.0;
  const ExcitonStates st = dimer_states(H);
  Bath b;
  b.lambda_cm = 1.0;
  b.cutoff_cm = 100.0;
  const LineshapeEngine eng(b, 4000.0, 40.0);
  const Eigen::MatrixXd Km = modified_redfield_rate_matrix(st, eng), Kr = redfield_rate_matrix(st, b);
  EXPECT_NEAR(Km(0, 1) / Kr(0, 1), 1.0, 0.03);  // вниз
  EXPECT_NEAR(Km(1, 0) / Kr(1, 0), 1.0, 0.03);  // вверх
}

TEST(ModifiedRedfield, DetailedBalanceWithReorganizedEnergies) {
  Eigen::MatrixXd H(2, 2);
  H << 17300.0, 80.0, 80.0, 16900.0;
  const ExcitonStates st = dimer_states(H);
  Bath b;
  b.lambda_cm = 80.0;
  b.cutoff_cm = 100.0;
  const LineshapeEngine eng(b, 6000.0, 40.0);
  const Eigen::MatrixXd K = modified_redfield_rate_matrix(st, eng);
  const double lam = b.reorganization();
  const double dF = (st.energies(1) - lam * st.ipr(1)) - (st.energies(0) - lam * st.ipr(0));
  EXPECT_NEAR(K(0, 1) / K(1, 0) / std::exp(b.beta() * dF), 1.0, 0.03);
}

TEST(ModifiedRedfield, DimerIntermediateCouplingCloserToHeomThanRedfield) {
  // λ ~ J: секулярный Редфилд сильно недооценивает скорость (перенормировка λ не учтена); модифицированный ближе
  // к HEOM (немного быстрее на t < 0.5 пс — марковское приближение) и даёт верное равновесие
  Eigen::MatrixXd H(2, 2);
  H << 17250.0, 60.0, 60.0, 16950.0;
  HeomParams p;
  p.lambda_cm = 60.0;
  p.gamma_cm = 100.0;
  p.temperature_K = 300.0;
  p.matsubara = 1;
  p.depth = 6;
  p.t_max_ps = 1.5;
  const ExcitonStates st = dimer_states(H);
  Eigen::MatrixXcd rho0 = (st.c.col(1) * st.c.col(1).transpose()).cast<std::complex<double>>();
  const HeomPopulations hp = heom_populations(H, rho0, p, 40);
  Bath b;
  b.lambda_cm = p.lambda_cm;
  b.cutoff_cm = p.gamma_cm;
  const LineshapeEngine eng(b, 6000.0, 40.0);
  Eigen::VectorXd t(static_cast<Eigen::Index>(hp.t_ps.size()));
  for (std::size_t k = 0; k < hp.t_ps.size(); ++k) t(static_cast<Eigen::Index>(k)) = hp.t_ps[k];
  Eigen::VectorXd P0(2);
  P0 << 0.0, 1.0;
  const Eigen::VectorXd Fe = st.energies - b.reorganization() * st.ipr;
  const Eigen::MatrixXd Pm = population_dynamics(modified_redfield_rate_matrix(st, eng) * units::rad_ps_per_cm,
                                                 Eigen::VectorXd::Zero(2), Fe, b.beta(), P0, t);
  const Eigen::MatrixXd Pr = population_dynamics(redfield_rate_matrix(st, b) * units::rad_ps_per_cm,
                                                 Eigen::VectorXd::Zero(2), st.energies, b.beta(), P0, t);
  double em = 0, er = 0;
  for (Eigen::Index k = 0; k < t.size(); ++k) {
    em = std::max(em, std::abs(hp.exciton(0, k) - Pm(0, k)));
    er = std::max(er, std::abs(hp.exciton(0, k) - Pr(0, k)));
  }
  std::printf("dimer λ = J: max |ΔP| vs HEOM — modified Redfield %.3f, Redfield %.3f\n", em, er);
  EXPECT_LT(em, 0.7 * er);
  EXPECT_LT(em, 0.2);
  EXPECT_NEAR(Pm(0, t.size() - 1), hp.exciton(0, t.size() - 1), 0.02);
}
