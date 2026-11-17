#include "ssp/excitons/heom.hpp"

#include "ssp/units.hpp"

#include <cmath>
#include <map>
#include <stdexcept>

namespace ssp::excitons {

namespace {

using cd = std::complex<double>;
const cd I(0.0, 1.0);

// Индексы ADO: векторы n[j*(K+1)+k] ≥ 0 с суммой ≤ L; соседи n ± e_s.
struct Hierarchy {
  int M = 0;  // число мод N·(K+1)
  std::vector<std::vector<int>> nvec;
  std::vector<std::vector<int>> plus, minus;  // индекс соседа или −1

  Hierarchy(int n_modes, int depth) : M(n_modes) {
    std::map<std::vector<int>, int> index;
    std::vector<int> cur(M, 0);
    // Порядок: по возрастанию уровня, затем лексикографически.
    std::vector<std::vector<int>> level = {cur};
    nvec.push_back(cur);
    index[cur] = 0;
    for (int L = 1; L <= depth; ++L) {
      std::vector<std::vector<int>> next;
      for (const auto& v : level)
        for (int s = 0; s < M; ++s) {
          auto u = v;
          ++u[s];
          if (!index.contains(u)) {
            index[u] = static_cast<int>(nvec.size());
            nvec.push_back(u);
            next.push_back(u);
          }
        }
      level.swap(next);
    }
    plus.assign(nvec.size(), std::vector<int>(M, -1));
    minus.assign(nvec.size(), std::vector<int>(M, -1));
    for (std::size_t a = 0; a < nvec.size(); ++a)
      for (int s = 0; s < M; ++s) {
        auto u = nvec[a];
        ++u[s];
        if (auto it = index.find(u); it != index.end()) plus[a][s] = it->second;
        if (nvec[a][s] > 0) {
          auto d = nvec[a];
          --d[s];
          minus[a][s] = index.at(d);
        }
      }
  }
};

struct BathCoeffs {
  std::vector<cd> c;      // c_k
  std::vector<double> nu; // ν_k
  double delta = 0.0;     // терминатор Σ_{k>K} c_k/ν_k
};

BathCoeffs drude_coeffs(const HeomParams& p) {
  if (p.temperature_K <= 0.0) throw std::invalid_argument("HEOM: T > 0 required");
  const double beta = 1.0 / (units::kB_cm * p.temperature_K);
  const double lam = p.lambda_cm, g = p.gamma_cm;
  BathCoeffs b;
  b.nu.push_back(g);
  b.c.push_back(lam * g * cd(1.0 / std::tan(0.5 * beta * g), -1.0));
  double sum_kept = 0.0;
  for (int k = 1; k <= p.matsubara; ++k) {
    const double nu = 2.0 * units::pi * k / beta;
    const double ck = 4.0 * lam * g / beta * nu / (nu * nu - g * g);
    b.nu.push_back(nu);
    b.c.push_back(ck);
    sum_kept += ck / nu;
  }
  if (p.terminator) b.delta = 2.0 * lam / (beta * g) - lam / std::tan(0.5 * beta * g) - sum_kept;
  return b;
}

double auto_dt_tu(const Eigen::MatrixXd& Hs, const BathCoeffs& b, int depth) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hs, Eigen::EigenvaluesOnly);
  const double spread = es.eigenvalues().cwiseAbs().maxCoeff();
  double numax = 0.0;
  for (double v : b.nu) numax = std::max(numax, v);
  double cmax = 0.0;
  for (const cd& c : b.c) cmax = std::max(cmax, std::abs(c));
  const double rate = spread + depth * numax + std::sqrt(cmax) * std::sqrt(double(depth) + 1.0) * 2.0;
  return 0.5 / rate;  // RK4 устойчив до ~2.8/rate; берём с запасом
}

}  // namespace

std::size_t heom_ado_count(int n_sites, int matsubara, int depth) {
  const int M = n_sites * (matsubara + 1);
  double c = 1.0;  // C(M + L, L)
  for (int i = 1; i <= depth; ++i) c = c * (M + i) / i;
  return static_cast<std::size_t>(std::llround(c));
}

HeomAbsorption heom_absorption(const Eigen::MatrixXd& H_cm, const Eigen::MatrixX3d& site_mu, const HeomParams& p,
                               const std::vector<double>& w_cm) {
  const int N = static_cast<int>(H_cm.rows());
  const int K1 = p.matsubara + 1;
  const BathCoeffs b = drude_coeffs(p);
  const Hierarchy h(N * K1, p.depth);
  const std::size_t A = h.nvec.size();

  HeomAbsorption out;
  out.e0_cm = H_cm.diagonal().mean();
  const Eigen::MatrixXd Hs = H_cm - out.e0_cm * Eigen::MatrixXd::Identity(N, N);
  const double dt = p.dt_fs > 0.0 ? p.dt_fs * 1e-3 / units::ps_per_tu : auto_dt_tu(Hs, b, p.depth);
  const std::size_t steps = static_cast<std::size_t>(std::ceil(p.t_max_ps / units::ps_per_tu / dt));
  out.info = {A, dt * units::ps_per_tu * 1e3, steps};

  std::vector<double> damp(A, 0.0);
  for (std::size_t a = 0; a < A; ++a)
    for (int s = 0; s < h.M; ++s) damp[a] += h.nvec[a][s] * b.nu[s % K1];

  // Состояние: A векторов длины N, хранится как матрица N × A.
  auto rhs = [&](const Eigen::MatrixXcd& X, Eigen::MatrixXcd& dX) {
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t a = 0; a < static_cast<std::ptrdiff_t>(A); ++a) {
      Eigen::VectorXcd d = -I * (Hs * X.col(a)) - (damp[a] + b.delta) * X.col(a);
      for (int s = 0; s < h.M; ++s) {
        const int j = s / K1, k = s % K1;
        if (const int ap = h.plus[a][s]; ap >= 0) d(j) += -I * X(j, ap);
        if (const int am = h.minus[a][s]; am >= 0) d(j) += -I * double(h.nvec[a][s]) * b.c[k] * X(j, am);
      }
      dX.col(a) = d;
    }
  };

  out.response.assign(steps + 1, 0.0);
  out.t_ps.resize(steps + 1);
  for (std::size_t q = 0; q <= steps; ++q) out.t_ps[q] = q * dt * units::ps_per_tu;
  for (int pol = 0; pol < 3; ++pol) {
    const Eigen::VectorXd mu = site_mu.col(pol);
    if (mu.squaredNorm() == 0.0) continue;
    Eigen::MatrixXcd X = Eigen::MatrixXcd::Zero(N, static_cast<Eigen::Index>(A));
    X.col(0) = mu.cast<cd>();
    Eigen::MatrixXcd k1(N, A), k2(N, A), k3(N, A), k4(N, A);
    out.response[0] += mu.dot(mu);
    for (std::size_t q = 1; q <= steps; ++q) {
      rhs(X, k1);
      rhs(X + 0.5 * dt * k1, k2);
      rhs(X + 0.5 * dt * k2, k3);
      rhs(X + dt * k3, k4);
      X += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
      out.response[q] += mu.cast<cd>().dot(X.col(0));
    }
  }
  // L(ω) = (1/π) Re∫ e^{i(ω−e0)t} R(t)/R(0) dt (трапеции), t в единицах t_u.
  const cd R0 = out.response[0];
  out.w_cm = w_cm;
  out.lineshape.assign(w_cm.size(), 0.0);
  for (std::size_t iw = 0; iw < w_cm.size(); ++iw) {
    const double dw = w_cm[iw] - out.e0_cm;
    cd s = 0.0;
    for (std::size_t q = 0; q <= steps; ++q) {
      const double wt = (q == 0 || q == steps) ? 0.5 : 1.0;
      s += wt * std::exp(I * (dw * q * dt)) * out.response[q] / R0;
    }
    out.lineshape[iw] = dt * s.real() / units::pi;
  }
  return out;
}

HeomPopulations heom_populations(const Eigen::MatrixXd& H_cm, const Eigen::MatrixXcd& rho0, const HeomParams& p,
                                 int n_out) {
  const int N = static_cast<int>(H_cm.rows());
  const int K1 = p.matsubara + 1;
  const BathCoeffs b = drude_coeffs(p);
  const Hierarchy h(N * K1, p.depth);
  const std::size_t A = h.nvec.size();

  const double e0 = H_cm.diagonal().mean();
  const Eigen::MatrixXd Hs = H_cm - e0 * Eigen::MatrixXd::Identity(N, N);
  const double dt = p.dt_fs > 0.0 ? p.dt_fs * 1e-3 / units::ps_per_tu : auto_dt_tu(Hs, b, p.depth);
  const std::size_t steps = static_cast<std::size_t>(std::ceil(p.t_max_ps / units::ps_per_tu / dt));
  HeomPopulations out;
  out.info = {A, dt * units::ps_per_tu * 1e3, steps};

  std::vector<double> damp(A, 0.0);
  for (std::size_t a = 0; a < A; ++a)
    for (int s = 0; s < h.M; ++s) damp[a] += h.nvec[a][s] * b.nu[s % K1];

  using Mat = Eigen::MatrixXcd;
  auto rhs = [&](const std::vector<Mat>& X, std::vector<Mat>& dX) {
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t a = 0; a < static_cast<std::ptrdiff_t>(A); ++a) {
      const Mat& r = X[a];
      Mat d = -I * (Hs * r - r * Hs) - damp[a] * r;
      if (b.delta != 0.0) {
        // −δ Σ_j [V_j,[V_j,ρ]]: элемент (u,v) умножается на число сайтов j, для которых δ_uj ≠ δ_vj, т.е. 2(1−δ_uv)
        for (int u = 0; u < N; ++u)
          for (int v = 0; v < N; ++v)
            if (u != v) d(u, v) -= 2.0 * b.delta * r(u, v);
      }
      for (int s = 0; s < h.M; ++s) {
        const int j = s / K1, k = s % K1;
        if (const int ap = h.plus[a][s]; ap >= 0) {
          // −i[V_j, ρ+]
          const Mat& rp = X[ap];
          d.row(j) += -I * rp.row(j);
          d.col(j) -= -I * rp.col(j);
        }
        if (const int am = h.minus[a][s]; am >= 0) {
          // −i n (c V_j ρ− − c* ρ− V_j)
          const Mat& rm = X[am];
          const double nn = h.nvec[a][s];
          d.row(j) += -I * nn * b.c[k] * rm.row(j);
          d.col(j) -= -I * nn * std::conj(b.c[k]) * rm.col(j);
        }
      }
      dX[a] = d;
    }
  };

  std::vector<Mat> X(A, Mat::Zero(N, N)), k1 = X, k2 = X, k3 = X, k4 = X, tmp = X;
  X[0] = rho0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_cm);
  const Eigen::MatrixXd U = es.eigenvectors();

  const std::size_t every = std::max<std::size_t>(1, steps / std::max(1, n_out));
  std::vector<double> ts;
  std::vector<Eigen::VectorXd> site, exc;
  auto record = [&](std::size_t q) {
    ts.push_back(q * dt * units::ps_per_tu);
    site.push_back(X[0].diagonal().real());
    exc.push_back((U.transpose() * X[0] * U).diagonal().real());
  };
  record(0);
  auto axpy = [&](const std::vector<Mat>& base, const std::vector<Mat>& k, double f, std::vector<Mat>& o) {
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t a = 0; a < static_cast<std::ptrdiff_t>(A); ++a) o[a] = base[a] + f * k[a];
  };
  for (std::size_t q = 1; q <= steps; ++q) {
    rhs(X, k1);
    axpy(X, k1, 0.5 * dt, tmp);
    rhs(tmp, k2);
    axpy(X, k2, 0.5 * dt, tmp);
    rhs(tmp, k3);
    axpy(X, k3, dt, tmp);
    rhs(tmp, k4);
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t a = 0; a < static_cast<std::ptrdiff_t>(A); ++a)
      X[a] += dt / 6.0 * (k1[a] + 2.0 * k2[a] + 2.0 * k3[a] + k4[a]);
    if (q % every == 0 || q == steps) record(q);
  }
  out.t_ps = ts;
  out.site.resize(N, static_cast<Eigen::Index>(ts.size()));
  out.exciton.resize(N, static_cast<Eigen::Index>(ts.size()));
  for (std::size_t k = 0; k < ts.size(); ++k) {
    out.site.col(static_cast<Eigen::Index>(k)) = site[k];
    out.exciton.col(static_cast<Eigen::Index>(k)) = exc[k];
  }
  return out;
}

}  // namespace ssp::excitons
