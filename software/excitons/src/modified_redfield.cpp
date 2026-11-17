#include "ssp/excitons/modified_redfield.hpp"

#include <cmath>
#include <complex>
#include <vector>

namespace ssp::excitons {

Eigen::MatrixXd modified_redfield_rate_matrix(const ExcitonStates& st, const LineshapeEngine& engine) {
  using cd = std::complex<double>;
  const auto& g = engine.g();
  const double dt = engine.dt(), lam = engine.lambda_total();
  const std::size_t nt = g.size();
  // ġ, g̈ — конечные разности (g(0) = ġ(0) = 0)
  std::vector<cd> g1(nt), g2(nt);
  for (std::size_t j = 0; j < nt; ++j) {
    if (j == 0) {
      g1[j] = 0.0;
      g2[j] = 2.0 * g[1] / (dt * dt);
    } else if (j + 1 < nt) {
      g1[j] = (g[j + 1] - g[j - 1]) / (2.0 * dt);
      g2[j] = (g[j + 1] - 2.0 * g[j] + g[j - 1]) / (dt * dt);
    } else {
      g1[j] = (g[j] - g[j - 1]) / dt;
      g2[j] = g2[j - 1];
    }
  }
  const Eigen::MatrixXd& C = st.c;
  const Eigen::Index n = C.cols();
  const Eigen::MatrixXd C2 = C.array().square().matrix(), C3 = C.array().cube().matrix();
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
#pragma omp parallel for schedule(dynamic)
  for (Eigen::Index M = 0; M < n; ++M)
    for (Eigen::Index N = 0; N < n; ++N) {
      if (M == N) continue;
      const double a = st.ipr(M), b = st.ipr(N);
      const double p = C2.col(M).dot(C2.col(N));
      const double q1 = C.col(N).dot(C3.col(M)), q2 = C3.col(N).dot(C.col(M));
      const double Om = st.energies(M) - st.energies(N) - 2.0 * lam * (a - p);
      const double Lam = a + b - 2.0 * p;
      double acc = 0.0;
      for (std::size_t j = 0; j < nt; ++j) {
        const double t = dt * static_cast<double>(j);
        const cd e = std::exp(cd(0.0, Om * t) - Lam * g[j]);
        const cd br = (q1 - q2) * g1[j] + cd(0.0, 2.0 * lam * q1);
        const double f = (e * (p * g2[j] - br * br)).real();
        acc += (j == 0 ? 0.5 : 1.0) * f;
        if (j > 16 && std::abs(e) < 1e-7) break;
      }
      K(N, M) = std::max(0.0, 2.0 * acc * dt);
    }
  for (Eigen::Index M = 0; M < n; ++M) K(M, M) = -K.col(M).sum();
  return K;
}

}  // namespace ssp::excitons
