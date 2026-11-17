#include "ssp/excitons/huckel.hpp"

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace ssp::excitons {

std::vector<XyzAtom> read_xyz(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open xyz file: " + path);
  std::vector<XyzAtom> a;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream is(line);
    std::string el;
    double x, y, z;
    if (!(is >> el >> x >> y >> z)) continue;  // строка числа атомов и комментарий пропускаются
    if (el.empty() || !std::isalpha(static_cast<unsigned char>(el[0]))) continue;
    a.push_back({el, Vec3(x, y, z) * 0.1});
  }
  if (a.empty()) throw std::runtime_error("no atoms in " + path);
  return a;
}

HuckelResult huckel_transition_charges(const std::vector<XyzAtom>& atoms, const HuckelOptions& opt) {
  const std::size_t n = atoms.size();
  auto heavy = [&](std::size_t i) { return atoms[i].element != "H"; };
  // соседи (включая H — для числа связей)
  std::vector<std::vector<std::size_t>> nb(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      const double d = (atoms[i].r - atoms[j].r).norm();
      const double lim = (heavy(i) && heavy(j)) ? opt.bond_max_nm : 0.125;
      if (d < lim) {
        nb[i].push_back(j);
        nb[j].push_back(i);
      }
    }
  // кандидаты: sp2-C (3 соседа) и N, связанные с ними
  std::vector<char> cand(n, 0);
  for (std::size_t i = 0; i < n; ++i)
    if (atoms[i].element == "C" && nb[i].size() == 3) cand[i] = 1;
  for (std::size_t i = 0; i < n; ++i)
    if (atoms[i].element == "N")
      for (auto j : nb[i])
        if (cand[j] && atoms[j].element == "C") cand[i] = 1;
  // наибольшая связная компонента кандидатов
  std::vector<int> comp(n, -1);
  int best = -1;
  std::size_t bestsz = 0;
  for (std::size_t s = 0, c = 0; s < n; ++s) {
    if (!cand[s] || comp[s] >= 0) continue;
    std::vector<std::size_t> st = {s};
    comp[s] = static_cast<int>(c);
    std::size_t sz = 0;
    while (!st.empty()) {
      const auto i = st.back();
      st.pop_back();
      ++sz;
      for (auto j : nb[i])
        if (cand[j] && comp[j] < 0) {
          comp[j] = static_cast<int>(c);
          st.push_back(j);
        }
    }
    if (sz > bestsz) {
      bestsz = sz;
      best = static_cast<int>(c);
    }
    ++c;
  }
  HuckelResult R;
  for (std::size_t i = 0; i < n; ++i)
    if (comp[i] == best) R.pi_atoms.push_back(i);
  const auto m = static_cast<Eigen::Index>(R.pi_atoms.size());
  if (m < 3) throw std::invalid_argument("huckel: no conjugated π-system found");
  std::map<std::size_t, Eigen::Index> loc;
  for (Eigen::Index k = 0; k < m; ++k) loc[R.pi_atoms[static_cast<std::size_t>(k)]] = k;
  // матрица Хюккеля в единицах |β| (β < 0): H = −(h на диагонали + k·β-связи); энергия по возрастанию
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);
  int electrons = -opt.pi_charge;
  for (Eigen::Index k = 0; k < m; ++k) {
    const auto i = R.pi_atoms[static_cast<std::size_t>(k)];
    const bool N = atoms[i].element == "N";
    const bool pyrrole = N && nb[i].size() == 3;
    electrons += N ? (pyrrole ? 2 : 1) : 1;
    H(k, k) = N ? -(pyrrole ? opt.h_N_pyrrole : opt.h_N_pyridine) : 0.0;
    for (auto j : nb[i]) {
      const auto it = loc.find(j);
      if (it == loc.end()) continue;
      const bool Nj = atoms[j].element == "N";
      double kk;
      if (N || Nj) {
        const bool pyr = (N && pyrrole) || (Nj && nb[j].size() == 3);
        kk = pyr ? opt.k_CN_pyrrole : opt.k_CN_pyridine;
      } else {
        // β_CC по длине связи: 1.0 при 0.140 нм, экспоненциально (0.134 → 1.1, 0.146 → 0.9)
        kk = std::exp(-((atoms[i].r - atoms[j].r).norm() - 0.140) / 0.0569);
      }
      H(k, it->second) = -kk;
    }
  }
  if (electrons % 2 != 0 || electrons <= 0 || electrons / 2 >= m)
    throw std::invalid_argument("huckel: odd or invalid number of π electrons (" + std::to_string(electrons) +
                                "), check pi_charge");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  R.energies = es.eigenvalues();
  R.homo = electrons / 2 - 1;
  const Eigen::VectorXd cH = es.eigenvectors().col(R.homo), cL = es.eigenvectors().col(R.homo + 1);
  // заряды и система молекулы
  Vec3 cen = Vec3::Zero();
  for (auto i : R.pi_atoms) cen += atoms[i].r;
  cen /= static_cast<double>(m);
  Eigen::MatrixXd P(m, 3);
  Eigen::VectorXd q(m);
  Vec3 d = Vec3::Zero();
  for (Eigen::Index k = 0; k < m; ++k) {
    P.row(k) = (atoms[R.pi_atoms[static_cast<std::size_t>(k)]].r - cen).transpose();
    q(k) = std::sqrt(2.0) * cH(k) * cL(k);
    d += q(k) * P.row(k).transpose();
  }
  R.dipole_e_nm = d.norm();
  if (R.dipole_e_nm < 1e-9) throw std::invalid_argument("huckel: HOMO→LUMO transition is dipole-forbidden");
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinV);
  Vec3 z = svd.matrixV().col(2);  // нормаль к плоскости π-системы
  const Vec3 x = d.normalized();
  z = (z - z.dot(x) * x).normalized();
  const Vec3 y = z.cross(x);
  for (Eigen::Index k = 0; k < m; ++k) {
    const Vec3 r = P.row(k).transpose();
    R.charges.emplace_back(r.dot(x), r.dot(y), r.dot(z), q(k));
  }
  return R;
}

void write_transition_charges(const std::string& path, const HuckelResult& h, const std::string& comment) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  f << "# " << comment << "\n"
    << "# Hückel HOMO→LUMO transition charges (q in e, unnormalized; |Σqr| = " << h.dipole_e_nm << " e·nm)\n"
    << "# molecular frame: x — transition dipole, z — normal of the π-plane, origin — π-centroid\n"
    << "# x_nm y_nm z_nm q\n";
  f.precision(8);
  for (const auto& c : h.charges) f << c(0) << " " << c(1) << " " << c(2) << " " << c(3) << "\n";
}

}  // namespace ssp::excitons
