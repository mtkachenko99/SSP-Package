#include "ssp/excitons/hamiltonian.hpp"

#include "ssp/units.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ssp::excitons {

double pair_coupling(const Vec3& ri, const Vec3& mui, const Vec3& rj, const Vec3& muj, const CouplingParams& p) {
  const Vec3 rij = rj - ri;
  const double r = rij.norm();
  if (r < 1e-12) return 0.0;
  switch (p.model) {
    case CouplingModel::LegacyJ0: {
      const Vec3 n = rij / r;
      const double kappa = mui.dot(muj) - 3.0 * mui.dot(n) * muj.dot(n);
      return p.J0_cm * std::pow(p.r0_nm / r, 3) * kappa;
    }
    case CouplingModel::PointDipole: {
      const Vec3 n = rij / r;
      const double kappa = mui.dot(muj) - 3.0 * mui.dot(n) * muj.dot(n);
      return p.screening * units::dd_cm_per_D2_nm3 * p.mu_D * p.mu_D * kappa / (r * r * r);
    }
    case CouplingModel::TransitionCharges:
      throw std::invalid_argument("transition charges: use coupling_matrix (needs molecular frames)");
    case CouplingModel::ExtendedDipole: {
      // Заряды q = μ/l в точках r ± (l/2)μ̂. Энергия Σ q_a q_b / r_ab в единицах (Д/нм)²·нм.
      const double l = p.dipole_length_nm;
      if (l <= 0.0) throw std::invalid_argument("extended dipole: dipole_length_nm must be > 0");
      const double q = p.mu_D / l;  // Д/нм
      const Vec3 a[2] = {ri + 0.5 * l * mui, ri - 0.5 * l * mui};
      const Vec3 b[2] = {rj + 0.5 * l * muj, rj - 0.5 * l * muj};
      const double sa[2] = {1.0, -1.0};
      double s = 0.0;
      for (int u = 0; u < 2; ++u)
        for (int v = 0; v < 2; ++v) s += sa[u] * sa[v] / (a[u] - b[v]).norm();
      return p.screening * units::dd_cm_per_D2_nm3 * q * q * s;
    }
  }
  return 0.0;
}

std::vector<Eigen::Vector4d> load_transition_charges(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open transition charges file: " + path);
  std::vector<Eigen::Vector4d> q;
  std::string line;
  while (std::getline(in, line)) {
    for (char& ch : line)
      if (ch == ',' || ch == ';' || ch == '\t') ch = ' ';
    std::istringstream ss(line);
    double x, y, z, c;
    if (ss >> x >> y >> z >> c) q.emplace_back(x, y, z, c);
  }
  if (q.size() < 2) throw std::runtime_error("transition charges: fewer than two charges in " + path);
  return q;
}

Eigen::MatrixXd coupling_matrix(const Aggregate& agg, const CouplingParams& p) {
  const auto N = static_cast<Eigen::Index>(agg.size());
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
  if (p.model == CouplingModel::TransitionCharges) {
    if (p.charges.size() < 2) throw std::invalid_argument("transition_charges: no charges loaded");
    // Нормировка: |Σ q r| = μ (в Д), заряды — в Д/нм
    Vec3 d = Vec3::Zero();
    for (const auto& c : p.charges) d += c(3) * c.head<3>();
    const double scale = p.mu_D / d.norm();
    // Ориентация: дипольный момент зарядов поворачивается на ось x молекулы
    const Eigen::Matrix3d align = Eigen::Quaterniond::FromTwoVectors(d.normalized(), Vec3::UnitX()).toRotationMatrix();
    std::vector<std::vector<Vec3>> pos(static_cast<std::size_t>(N));
    for (Eigen::Index i = 0; i < N; ++i) {
      const Eigen::Matrix3d F = agg.frame(static_cast<std::size_t>(i));
      for (const auto& c : p.charges) pos[i].push_back(agg.positions[i] + F * (align * Vec3(c.head<3>())));
    }
    for (Eigen::Index i = 0; i < N; ++i)
      for (Eigen::Index j = i + 1; j < N; ++j) {
        double s = 0.0;
        for (std::size_t a = 0; a < p.charges.size(); ++a)
          for (std::size_t b = 0; b < p.charges.size(); ++b)
            s += p.charges[a](3) * p.charges[b](3) / (pos[i][a] - pos[j][b]).norm();
        J(i, j) = J(j, i) = p.screening * units::dd_cm_per_D2_nm3 * scale * scale * s;
      }
    return J;
  }
  for (Eigen::Index i = 0; i < N; ++i)
    for (Eigen::Index j = i + 1; j < N; ++j) {
      const double v = pair_coupling(agg.positions[i], agg.dipoles[i], agg.positions[j], agg.dipoles[j], p);
      J(i, j) = v;
      J(j, i) = v;
    }
  return J;
}

Eigen::MatrixXd frenkel_hamiltonian(const Eigen::MatrixXd& J, double monomer_energy_cm, double disorder_sigma_cm,
                                    std::mt19937_64& rng) {
  Eigen::MatrixXd H = J;
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (Eigen::Index i = 0; i < H.rows(); ++i)
    H(i, i) = monomer_energy_cm + (disorder_sigma_cm > 0.0 ? disorder_sigma_cm * gauss(rng) : 0.0);
  return H;
}

ExcitonStates diagonalize(const Eigen::MatrixXd& H, const Aggregate& agg, double mu_D) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  if (es.info() != Eigen::Success) throw std::runtime_error("Hamiltonian diagonalization failed");
  ExcitonStates s;
  s.energies = es.eigenvalues();
  s.c = es.eigenvectors();
  const auto N = H.rows();
  Eigen::MatrixX3d site_mu(N, 3);
  for (Eigen::Index i = 0; i < N; ++i) site_mu.row(i) = mu_D * agg.dipoles[static_cast<std::size_t>(i)].transpose();
  s.mu = s.c.transpose() * site_mu;
  s.mu2 = s.mu.rowwise().squaredNorm();
  s.ipr = s.c.array().pow(4).colwise().sum().transpose();
  return s;
}

namespace {
// Σ_k V_mk² по соседям (ячеечный поиск по центрам молекул). Молекула m — набор точечных зарядов pts[m] (общие
// значения q); V_mk = Σ_ab q_a q_b /√(r_ab² + a²).
Eigen::VectorXd sum_v2(const Aggregate& a, const std::vector<std::vector<Vec3>>& pts, const std::vector<double>& q,
                       double cutoff, double soft) {
  const std::size_t N = a.size();
  Vec3 lo = a.positions[0], hi = a.positions[0];
  for (const auto& r : a.positions) {
    lo = lo.cwiseMin(r);
    hi = hi.cwiseMax(r);
  }
  const Eigen::Vector3i nc = ((hi - lo) / cutoff).array().floor().cast<int>() + 1;
  auto cell = [&](const Vec3& r) -> Eigen::Vector3i {  // явный тип: иначе возвращается выражение Eigen на локальную
    const Eigen::Vector3i c = ((r - lo) / cutoff).array().floor().cast<int>();
    return c.cwiseMin(nc - Eigen::Vector3i::Ones()).cwiseMax(Eigen::Vector3i::Zero());
  };
  std::vector<std::vector<std::size_t>> cells(static_cast<std::size_t>(nc.prod()));
  auto idx = [&](const Eigen::Vector3i& c) { return static_cast<std::size_t>((c.z() * nc.y() + c.y()) * nc.x() + c.x()); };
  for (std::size_t i = 0; i < N; ++i) cells[idx(cell(a.positions[i]))].push_back(i);
  const double s2 = soft * soft;
  Eigen::VectorXd S = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(N));
#pragma omp parallel for schedule(dynamic, 64)
  for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(N); ++ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    const Eigen::Vector3i c = cell(a.positions[i]);
    double s = 0.0;
    for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          const Eigen::Vector3i cc = c + Eigen::Vector3i(dx, dy, dz);
          if ((cc.array() < 0).any() || (cc.array() >= nc.array()).any()) continue;
          for (const std::size_t k : cells[idx(cc)]) {
            if (k == i || (a.positions[k] - a.positions[i]).norm() > cutoff) continue;
            double v = 0.0;
            for (std::size_t u = 0; u < q.size(); ++u)
              for (std::size_t w = 0; w < q.size(); ++w)
                v += q[u] * q[w] / std::sqrt((pts[i][u] - pts[k][w]).squaredNorm() + s2);
            s += v * v;
          }
        }
    S(static_cast<Eigen::Index>(i)) = s;
  }
  return S;
}
}  // namespace

DispersionResult dispersion_site_shifts(const Aggregate& agg, const DispersionParams& p) {
  DispersionResult R;
  const auto N = static_cast<Eigen::Index>(agg.size());
  R.shift_cm = R.dispersive_cm = Eigen::VectorXd::Zero(N);
  if (!p.enabled || N == 0) return R;
  if (p.n_aggregate < 1.0 || p.n_solvent < 1.0 || p.F_eV < 0.0 || p.dipole_length_nm <= 0.0 || p.cutoff_nm <= 0.0)
    throw std::invalid_argument("dispersion: invalid parameters");
  const double cm_per_eV = 8065.544;
  auto bayliss = [](double n) { return (n * n - 1.0) / (2.0 * n * n + 1.0); };
  R.continuum_cm = -p.F_eV * bayliss(p.n_aggregate) * cm_per_eV;
  R.to_gas_cm = p.F_eV * bayliss(p.n_solvent) * cm_per_eV;
  // Распределение заряда молекулы в её системе (x — диполь, z — нормаль): протяжённый диполь ±1 на ±l/2 или
  // атомные переходные заряды, нормированные на |Σqr| = l.
  std::vector<Vec3> local;
  std::vector<double> q;
  if (p.charges.empty()) {
    local = {Vec3(0.5 * p.dipole_length_nm, 0, 0), Vec3(-0.5 * p.dipole_length_nm, 0, 0)};
    q = {1.0, -1.0};
  } else {
    Vec3 d = Vec3::Zero();
    for (const auto& c : p.charges) d += c(3) * c.head<3>();
    if (d.norm() < 1e-12) throw std::invalid_argument("dispersion: transition charges have zero dipole");
    const Eigen::Matrix3d align = Eigen::Quaterniond::FromTwoVectors(d.normalized(), Vec3::UnitX()).toRotationMatrix();
    for (const auto& c : p.charges) {
      local.push_back(align * Vec3(c.head<3>()));
      q.push_back(c(3) * p.dipole_length_nm / d.norm());
    }
  }
  auto place = [&](const std::vector<Eigen::Matrix3d>& frames) {
    std::vector<std::vector<Vec3>> pts(agg.size());
    for (std::size_t i = 0; i < agg.size(); ++i)
      for (const auto& r : local) pts[i].push_back(agg.positions[i] + frames[i] * r);
    return pts;
  };
  // калибровка: случайные ориентации (случайные повороты всей молекулы) на тех же положениях
  std::mt19937_64 rng(p.seed);
  std::normal_distribution<double> g(0.0, 1.0);
  double avg = 0.0;
  const int K = std::max(1, p.orientation_samples);
  for (int k = 0; k < K; ++k) {
    std::vector<Eigen::Matrix3d> fr(agg.size());
    for (auto& F : fr) F = Eigen::Quaterniond(g(rng), g(rng), g(rng), g(rng)).normalized().toRotationMatrix();
    avg += sum_v2(agg, place(fr), q, p.cutoff_nm, p.softening_nm).mean();
  }
  avg /= K;
  if (!(avg > 0.0)) throw std::invalid_argument("dispersion: no neighbours within the cutoff");
  R.Q_cm_nm2 = R.continuum_cm / avg;
  std::vector<Eigen::Matrix3d> fr(agg.size());
  for (std::size_t i = 0; i < agg.size(); ++i) fr[i] = agg.frame(i);
  const Eigen::VectorXd S = sum_v2(agg, place(fr), q, p.cutoff_nm, p.softening_nm);
  R.dispersive_cm = R.Q_cm_nm2 * S;
  R.shift_cm = R.dispersive_cm.array() + R.to_gas_cm + p.extra_shift_cm;
  return R;
}

}  // namespace ssp::excitons
