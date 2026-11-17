#pragma once
// Гамильтониан Френкеля одноэкситонного блока и его диагонализация.

#include "ssp/excitons/geometry.hpp"

#include <Eigen/Dense>

#include <random>
#include <vector>

namespace ssp::excitons {

enum class CouplingModel {
  PointDipole,     // J = f_s·5.034·μ²·κ/r³
  ExtendedDipole,  // заряды ±μ/l на расстоянии l (модель Куна)
  LegacyJ0,        // J = J0·(r0/r)³·κ — как в testcalc.ipynb
  TransitionCharges  // TrEsp: Σ q_a q_b / r_ab по переходным зарядам молекул (заряды — в системе молекулы)
};

struct CouplingParams {
  CouplingModel model = CouplingModel::PointDipole;
  double mu_D = 10.0;               // переходный диполь мономера, Д
  double screening = 1.0;           // множитель экранирования f_s (обычно 1/n²)
  double dipole_length_nm = 0.8;    // для ExtendedDipole
  double J0_cm = 900.0;             // для LegacyJ0
  double r0_nm = 1.0;               // для LegacyJ0
  // TrEsp: (x, y, z, q) в системе молекулы (нм; x — вдоль переходного диполя); нормируются на μ_D
  std::vector<Eigen::Vector4d> charges;
};

// Загрузка переходных зарядов из CSV/текста: x_nm y_nm z_nm q (лишние столбцы и строки-комментарии игнорируются).
std::vector<Eigen::Vector4d> load_transition_charges(const std::string& path);

// Матрица связей J (см⁻¹), нули на диагонали.
Eigen::MatrixXd coupling_matrix(const Aggregate& agg, const CouplingParams& p);

// Энергия связи пары (см⁻¹) по выбранной модели.
double pair_coupling(const Vec3& ri, const Vec3& mui, const Vec3& rj, const Vec3& muj, const CouplingParams& p);

// Зависящий от окружения дисперсионный (ван-дер-ваальсов) сдвиг энергий сайтов (Megow et al., PCCP 17, 6741
// (2015); arXiv:1412.3554): ΔE_m = Q·Σ_k V_mk², V_mk — кулоновская энергия двух протяжённых диполей с единичными
// зарядами (длина dipole_length_nm, вдоль диполей молекул). Q калибруется по континууму: при случайных взаимных
// ориентациях (те же положения) ⟨Σ_k V_mk²⟩ должно давать сдвиг Байлиса −F(n²−1)/(2n²+1) с n = n_aggregate.
// Энергия сайта: E_m = E_mono(растворитель) + F·f(n_solvent) [возврат к газовой фазе] + ΔE_m + extra_shift_cm.
struct DispersionParams {
  bool enabled = false;
  double F_eV = 1.43;              // коэффициент Байлиса красителя
  double n_aggregate = 1.78;       // показатель преломления плотного слоя красителя (длинноволновый предел)
  double n_solvent = 1.33;         // растворитель, в котором измерена энергия мономера
  double dipole_length_nm = 1.0;   // длина протяжённого диполя для V_mk
  double cutoff_nm = 3.0;          // радиус суммирования (Megow: 2 нм дают > 99 %)
  double softening_nm = 0.3;       // размытие зарядов: 1/√(r² + a²) — переходная плотность не точечная; без него
                                   // почти совпадающие заряды плотно упакованных соседей дают расходимость
  double extra_shift_cm = 0.0;     // остаточный контакт с растворителем у молекул агрегата
  int orientation_samples = 16;
  unsigned seed = 1;
  // Атомные переходные заряды (x, y, z, q в системе молекулы; x — диполь, z — нормаль): если заданы, V_mk
  // считается по ним вместо протяжённого диполя. Нормируются на |Σqr| = dipole_length_nm (единичный «заряд»
  // в тех же единицах, что и у протяжённого диполя), так что Q сохраняет смысл.
  std::vector<Eigen::Vector4d> charges;
};
struct DispersionResult {
  Eigen::VectorXd shift_cm;        // полный сдвиг сайта относительно энергии мономера в растворителе
  Eigen::VectorXd dispersive_cm;   // ΔE_m (только дисперсионная часть)
  double Q_cm_nm2 = 0.0;           // калибровочный множитель (см⁻¹·нм²)
  double continuum_cm = 0.0;       // −F(n_agg²−1)/(2n_agg²+1)
  double to_gas_cm = 0.0;          // +F(n_s²−1)/(2n_s²+1)
};
DispersionResult dispersion_site_shifts(const Aggregate& agg, const DispersionParams& p);

// H = diag(E_mono + δ) + J, δ ~ N(0, sigma).
Eigen::MatrixXd frenkel_hamiltonian(const Eigen::MatrixXd& J, double monomer_energy_cm, double disorder_sigma_cm,
                                    std::mt19937_64& rng);

struct ExcitonStates {
  Eigen::VectorXd energies;     // ε_n, см⁻¹ (по возрастанию)
  Eigen::MatrixXd c;            // c(i, n): амплитуда сайта i в состоянии n
  Eigen::MatrixX3d mu;          // μ_n (Д), строка n
  Eigen::VectorXd mu2;          // |μ_n|² (Д²)
  Eigen::VectorXd ipr;          // s_n = Σ_i c_in⁴
  std::size_t size() const { return static_cast<std::size_t>(energies.size()); }
};

// Диагонализация: собственные состояния и переходные диполи μ_n = μ·Σ_i c_in μ̂_i.
ExcitonStates diagonalize(const Eigen::MatrixXd& H, const Aggregate& agg, double mu_D);

}  // namespace ssp::excitons
