#include "ssp/excitons/vibronic.hpp"

#include "ssp/units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ssp::excitons {

VibronicStates holstein_opa(const Eigen::MatrixXd& H, const Eigen::MatrixX3d& mu_sites, const HolsteinParams& p) {
  if (p.nu_max < 0 || p.omega_cm <= 0.0 || p.huang_rhys < 0.0) throw std::invalid_argument("holstein: invalid parameters");
  const Eigen::Index N = H.rows(), V = p.nu_max + 1, D = N * V;
  Eigen::VectorXd F(V);
  for (Eigen::Index v = 0; v < V; ++v)
    F(v) = std::exp(-0.5 * p.huang_rhys + 0.5 * v * std::log(std::max(p.huang_rhys, 1e-300)) - 0.5 * std::lgamma(v + 1.0));
  if (p.huang_rhys == 0.0) {
    F.setZero();
    F(0) = 1.0;
  }
  const Eigen::MatrixXd FF = F * F.transpose();
  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(D, D);  // индекс n·V + ν
  for (Eigen::Index n = 0; n < N; ++n)
    for (Eigen::Index m = 0; m < N; ++m) {
      if (n == m) {
        for (Eigen::Index v = 0; v < V; ++v)
          Hv(n * V + v, n * V + v) = H(n, n) - p.huang_rhys * p.omega_cm + double(v) * p.omega_cm;
      } else {
        Hv.block(n * V, m * V, V, V) = H(n, m) * FF;
      }
    }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hv);
  const Eigen::MatrixXd& C = es.eigenvectors();
  Eigen::MatrixXd Mv(D, 3);
  for (Eigen::Index n = 0; n < N; ++n)
    for (Eigen::Index v = 0; v < V; ++v) Mv.row(n * V + v) = F(v) * mu_sites.row(n);
  const Eigen::MatrixXd mu = C.transpose() * Mv;
  VibronicStates out;
  out.energies = es.eigenvalues();
  out.mu2 = mu.rowwise().squaredNorm();
  out.ipr.resize(D);
  for (Eigen::Index k = 0; k < D; ++k) {
    double s = 0.0;
    for (Eigen::Index n = 0; n < N; ++n) {
      const double w = C.col(k).segment(n * V, V).squaredNorm();
      s += w * w;
    }
    out.ipr(k) = s;
  }
  return out;
}

Eigen::MatrixX3d site_dipoles(const ExcitonStates& st) { return st.c * st.mu; }

void holstein_absorption(const std::vector<Realization>& reals, const HolsteinParams& p, const LineshapeEngine& engine,
                         double n_medium, const std::vector<double>& w, std::vector<double>& mean,
                         std::vector<double>& se) {
  const std::size_t nw = w.size();
  mean.assign(nw, 0.0);
  se.assign(nw, 0.0);
  if (reals.empty()) return;
  const double pref = absorption_cross_section_prefactor_nm2(n_medium);
  std::vector<std::vector<double>> per(reals.size(), std::vector<double>(nw, 0.0));
#pragma omp parallel for schedule(dynamic)
  for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(reals.size()); ++r) {
    const Realization& R = reals[static_cast<std::size_t>(r)];
    const VibronicStates vs = holstein_opa(R.H, site_dipoles(R.states), p);
    // Время жизни вибронного состояния — как у электронного экситона, ближайшего по энергии (скорость ухода по
    // Редфилду); внутримолекулярная колебательная релаксация отдельно не учитывается.
    const Eigen::VectorXd Gout = -R.K_ps.diagonal() / units::rad_ps_per_cm;
    const Eigen::VectorXd& Ee = R.states.energies;
    for (Eigen::Index k = 0; k < vs.energies.size(); ++k) {
      if (vs.mu2(k) <= 1e-10 * vs.mu2.maxCoeff()) continue;
      const double e = vs.energies(k);
      const auto it = std::lower_bound(Ee.data(), Ee.data() + Ee.size(), e);
      Eigen::Index j = std::clamp<Eigen::Index>(it - Ee.data(), 0, Ee.size() - 1);
      if (j > 0 && std::abs(Ee(j - 1) - e) < std::abs(Ee(j) - e)) --j;
      engine.add_line(per[static_cast<std::size_t>(r)], w, e, pref * vs.mu2(k), vs.ipr(k), Gout(j), false, 1);
    }
  }
  const double Rn = static_cast<double>(reals.size());
  for (std::size_t q = 0; q < nw; ++q) {
    double s1 = 0, s2 = 0;
    for (const auto& v : per) {
      s1 += v[q];
      s2 += v[q] * v[q];
    }
    mean[q] = s1 / Rn;
    se[q] = reals.size() > 1 ? std::sqrt(std::max(0.0, s2 / Rn - mean[q] * mean[q]) / (Rn - 1.0)) : 0.0;
  }
}

}  // namespace ssp::excitons
