#pragma once
// Переходные заряды (TrEsp) по атомной геометрии в π-электронном приближении Хюккеля: для цианинов S₁ — почти
// чистый переход HOMO→LUMO, q_i ∝ c_i^HOMO c_i^LUMO на атомах π-системы. Это форма переходной плотности
// без TD-DFT; абсолютная величина нормируется на измеренный μ при использовании (coupling_matrix).
// π-система: атомы C с тремя соседями и N, связанные с ними (наибольшая связная компонента). Электроны: C — 1,
// N с тремя соседями (пиррольный) — 2, N с двумя (пиридиновый) — 1, минус заряд π-системы.
// Параметры (Streitwieser): α_N = α + h_N β, β_CN = k_CN β; β_CC зависит от длины связи.

#include "ssp/excitons/geometry.hpp"

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace ssp::excitons {

struct XyzAtom {
  std::string element;
  Vec3 r;  // нм
};

std::vector<XyzAtom> read_xyz(const std::string& path);  // XYZ в Å → нм

struct HuckelOptions {
  double h_N_pyrrole = 1.5, k_CN_pyrrole = 0.8;
  double h_N_pyridine = 0.5, k_CN_pyridine = 1.0;
  int pi_charge = 1;           // заряд π-системы (катионный цианин: +1)
  double bond_max_nm = 0.19;   // порог связи тяжёлых атомов (C–Cl ≈ 0.175 нм; несвязанные ≥ 0.24 нм)
};

struct HuckelResult {
  std::vector<std::size_t> pi_atoms;      // индексы атомов π-системы в исходном списке
  std::vector<Eigen::Vector4d> charges;   // (x, y, z, q) в системе молекулы: x — диполь, z — нормаль плоскости
  Eigen::VectorXd energies;               // в единицах |β|, по возрастанию энергии
  int homo = -1;
  double dipole_e_nm = 0.0;               // |Σ q r| до нормировки (q в долях e)
};

HuckelResult huckel_transition_charges(const std::vector<XyzAtom>& atoms, const HuckelOptions& opt = {});
void write_transition_charges(const std::string& path, const HuckelResult& h, const std::string& comment);

}  // namespace ssp::excitons
