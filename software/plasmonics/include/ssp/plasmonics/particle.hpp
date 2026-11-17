#pragma once
// Квазистатический отклик однородной частицы (сфера, вытянутый и сплюснутый сфероид с осью z) на произвольные
// внешние потенциалы: точное разложение по сферическим/сфероидальным гармоникам до l ≤ L_max.
//
// Наведённый потенциал точечного заряда q' в r' (в среде ε_m), вычисленный в r:
//   φ_ind(r) = q'/(4πε₀ε_m) · K(r, r'),   K(r, r') = Σ_{l≥1,m} w_lm D̂_lm(ε) [f^c_lm(r)f^c_lm(r') + f^s_lm(r)f^s_lm(r')]
//   f^{c,s}_lm(r) = q_lm(ξ) P̄_l^m(η) {cos, sin}(mφ),  q_lm(ξ) = Q_l^m(ξ)/Q_l^m(ξ₀)
//   D̂_lm(ε) = (ε − ε_m) a_P / (ε_m a_Q − ε a_P),  a_P = d ln P_l^m/dξ |ξ₀,  a_Q = d ln Q_l^m/dξ |ξ₀.
// Для сферы ξ = r, P = r^l, Q = r^{−l−1} и K = −Σ α_l a^{2l+1} P_l(cos γ)/(r r')^{l+1}.
// Энергия взаимодействия диполей через частицу: +(μ_i·∇)(μ_j·∇')K/(4πε₀ε_m). Однородное поле E₀ вдоль оси a наводит
//   E_ind(r) = E₀ u_a D̂_{1m} ∇f_{1m}(r).

#include <Eigen/Dense>

#include <array>
#include <complex>
#include <vector>

namespace ssp::plasmonics {

using cd = std::complex<double>;
using Vec3 = Eigen::Vector3d;

enum class ParticleShape { Sphere, Prolate, Oblate };

struct MultipoleTerm {
  int l = 0, m = 0;
  cd w;       // нормировка разложения
  cd aP, aQ;  // логарифмические производные радиальных функций на поверхности
};

class QuasistaticParticle {
 public:
  // Полуоси: a_xy (в плоскости xy) и a_z (вдоль оси симметрии z). Форма выбирается по их соотношению;
  // при |a_xy − a_z| < 1e-9·a используется сфера.
  QuasistaticParticle(double a_xy_nm, double a_z_nm, int l_max);

  ParticleShape shape() const { return shape_; }
  int l_max() const { return L_; }
  double a_xy() const { return axy_; }
  double a_z() const { return az_; }
  double focal() const { return f_; }

  bool inside(const Vec3& r) const;
  // Расстояние до поверхности (приближённо, по направлению от центра; точно для сферы).
  double surface_gap(const Vec3& r) const;

  const std::vector<MultipoleTerm>& terms() const { return terms_; }
  cd response(const MultipoleTerm& t, cd eps, double eps_m) const;

  // f^{c,s} для всех членов (порядок terms()).
  void values(const Vec3& r, std::vector<cd>& fc, std::vector<cd>& fs) const;
  // ∇f^{c,s} (центральные разности 4-го порядка).
  void gradients(const Vec3& r, std::vector<Eigen::Vector3cd>& gc, std::vector<Eigen::Vector3cd>& gs) const;

  // Индукционный потенциал K(r, r') при заданных ε (для проверок).
  cd induced_kernel(const Vec3& r, const Vec3& rp, cd eps, double eps_m) const;
  // Разложение свободного 1/|r − r'| по тем же функциям (для проверки нормировки), r ближе к частице.
  cd free_kernel_expansion(const Vec3& r_inner, const Vec3& r_outer, int l_max) const;

  // Однородное поле: E_ind(r) = Σ_a E₀_a u_a D̂ ∇f — индексы членов (l=1) и коэффициенты u_a для x, y, z.
  const std::array<int, 3>& uniform_term() const { return uterm_; }
  const std::array<bool, 3>& uniform_is_sine() const { return usin_; }
  const std::array<cd, 3>& uniform_coeff() const { return ucoef_; }

  // Квазистатическая поляризуемость вдоль x, y, z (нм³, p = 4πε₀ε_m α E): (abc/3)(ε−ε_m)/(ε_m + L(ε−ε_m)).
  Eigen::Vector3cd polarizability(cd eps, double eps_m) const;
  Eigen::Vector3d depolarization() const { return L_dep_; }
  Eigen::Vector3d semi_axes() const { return Eigen::Vector3d(axy_, axy_, az_); }

 private:
  ParticleShape shape_;
  double axy_, az_, f_ = 0.0, xi0_ = 0.0;
  int L_;
  std::vector<MultipoleTerm> terms_;
  std::array<int, 3> uterm_{};
  std::array<bool, 3> usin_{};
  std::array<cd, 3> ucoef_{};
  Eigen::Vector3d L_dep_;
  double fd_h_ = 1e-4;

  struct Coords {
    double xi, eta, phi;
  };
  Coords coords(const Vec3& r) const;
  // ln P_l^m, ln Q_l^m (l = m..L) при радиальной координате ξ (комплексный аргумент для сплюснутого).
  void radial_log(double xi, int m, std::vector<cd>& lP, std::vector<cd>& lQ) const;
  std::vector<std::vector<cd>> lQ0_;  // ln Q_l^m(ξ₀) для нормировки q
  int term_index(int l, int m) const { return l * (l + 1) / 2 + m - 1; }
};

// Нормированные присоединённые функции Лежандра P̄_l^m(η) (без фазы Кондона–Шортли), l = m..L.
void normalized_legendre(int L, int m, double eta, std::vector<double>& out);

}  // namespace ssp::plasmonics
