#include "ssp/excitons/relaxation.hpp"

#include "ssp/units.hpp"

#include <cmath>
#include <stdexcept>

namespace ssp::excitons {

Eigen::MatrixXd redfield_rate_matrix(const ExcitonStates& st, const Bath& bath) {
  const auto N = static_cast<Eigen::Index>(st.size());
  const Eigen::MatrixXd p = st.c.array().square().matrix();  // p(i,n) = c_in²
  const Eigen::MatrixXd overlap = p.transpose() * p;         // Σ_i c_in² c_im²
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(N, N);
  for (Eigen::Index n = 0; n < N; ++n)
    for (Eigen::Index m = 0; m < N; ++m)
      if (m != n) K(m, n) = overlap(m, n) * bath.rate_spectrum(st.energies(n) - st.energies(m));
  for (Eigen::Index n = 0; n < N; ++n) K(n, n) = -K.col(n).sum();
  return K;
}

Eigen::VectorXd redfield_homogeneous_hwhm(const ExcitonStates& st, const Bath& bath, const Eigen::MatrixXd& K) {
  const double C0 = bath.rate_spectrum(0.0);
  Eigen::VectorXd h(K.cols());
  for (Eigen::Index n = 0; n < K.cols(); ++n) h(n) = -0.5 * K(n, n) + 0.5 * st.ipr(n) * C0;
  return h;
}

Eigen::VectorXd radiative_rates_ps(const ExcitonStates& st, const DecayParams& p) {
  using namespace units;
  Eigen::VectorXd g(st.size());
  for (Eigen::Index n = 0; n < g.size(); ++n) {
    const double eps = st.energies(n);
    if (p.monomer_lifetime_ns > 0.0) {
      const double ratio = st.mu2(n) / (p.mu_D * p.mu_D);
      g(n) = ratio * std::pow(eps / p.monomer_energy_cm, 3) / (p.monomer_lifetime_ns * 1e3);
    } else {
      const double w = rad_s_per_cm * eps;
      const double mu2 = st.mu2(n) * debye * debye;
      g(n) = p.n_medium * w * w * w * mu2 / (3.0 * pi * eps0 * hbar * c * c * c) * 1e-12;
    }
  }
  return g;
}

Eigen::VectorXd steady_state(const Eigen::MatrixXd& K_ps, const Eigen::VectorXd& decay_ps, const Eigen::VectorXd& pump) {
  Eigen::MatrixXd A = K_ps;
  A.diagonal() -= decay_ps;
  if (decay_ps.minCoeff() <= 0.0) throw std::invalid_argument("steady_state: all decay rates must be > 0");
  Eigen::VectorXd P = A.partialPivLu().solve(-pump);
  return P;
}

Eigen::MatrixXd population_dynamics(const Eigen::MatrixXd& K_ps, const Eigen::VectorXd& decay_ps,
                                    const Eigen::VectorXd& energies_cm, double beta, const Eigen::VectorXd& P0,
                                    const Eigen::VectorXd& times_ps) {
  // Детальный баланс K(m,n)e^{−βε_n} = K(n,m)e^{−βε_m} ⇒ S = D⁻¹(K − Γ)D симметрична, D = diag(e^{−βε/2}).
  const Eigen::Index N = K_ps.rows();
  Eigen::MatrixXd A = K_ps;
  A.diagonal() -= decay_ps;
  Eigen::VectorXd d(N);
  const double e0 = energies_cm.minCoeff();
  for (Eigen::Index n = 0; n < N; ++n) d(n) = std::isinf(beta) ? 1.0 : std::exp(-0.5 * beta * (energies_cm(n) - e0));
  Eigen::MatrixXd S = d.cwiseInverse().asDiagonal() * A * d.asDiagonal();
  S = 0.5 * (S + S.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  const Eigen::VectorXd y0 = es.eigenvectors().transpose() * d.cwiseInverse().asDiagonal() * P0;
  Eigen::MatrixXd out(N, times_ps.size());
  for (Eigen::Index k = 0; k < times_ps.size(); ++k) {
    const Eigen::VectorXd y = (es.eigenvalues() * times_ps(k)).array().exp() * y0.array();
    out.col(k) = d.asDiagonal() * (es.eigenvectors() * y);
  }
  return out;
}

}  // namespace ssp::excitons
