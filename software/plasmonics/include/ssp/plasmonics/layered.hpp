#pragma once
// Слоистая сфера: металлическое ядро + диэлектрические оболочки (лиганды/спейсер) + граничные условия Фейбельмана.
//
// Квазистатика, мультиполь порядка l. В слое k потенциал φ_k = P_k r^l + Q_k r^{−l−1}; ядро Q_0 = 0; снаружи P = 1
// (внешний источник), Q_out — наведённый коэффициент. На границе радиуса R между слоями (вн. → нар.):
//   ε_in ∂φ_in = ε_out ∂φ_out                     (d_∥ = 0: нормальная компонента D непрерывна)
//   φ_out − φ_in = −d⊥ (∂φ_out − ∂φ_in)          (Фейбельман, d⊥ > 0 — центроид заряда снаружи; только на поверхности ядра)
// Нормированный отклик D̂_l = Q_out / b^{2l+1} (b — внешний радиус); для однородной сферы D̂_l = −l(ε−ε_m)/(lε+(l+1)ε_m).
// Дипольная поляризуемость (нм³): α = −D̂_1·b³.

#include "ssp/plasmonics/mie.hpp"

#include <complex>
#include <vector>

namespace ssp::plasmonics {

using cd = std::complex<double>;

struct SphereLayers {
  double core_radius_nm = 10.0;
  std::vector<double> shell_thickness_nm;  // изнутри наружу
  std::vector<cd> shell_eps;               // ε оболочек (не зависят от частоты)
  cd d_perp_nm = 0.0;                      // параметр Фейбельмана d⊥ на поверхности ядра
  double outer_radius() const;
  bool homogeneous() const { return shell_thickness_nm.empty(); }
};

// D̂_l для слоистой сферы (см. выше).
cd layered_response(int l, cd eps_core, const SphereLayers& s, double eps_m);

// Ми для сферы с одной оболочкой (Bohren & Huffman, BHCOAT): x = k·a (ядро), y = k·b (оболочка), m1, m2 — относительные
// показатели ядра и оболочки. Сечения в нм² — через mie_coated_sphere.
MieResult mie_coated(double x, double y, cd m1, cd m2, int n_max = 0);
MieResult mie_coated_sphere(double a_nm, double b_nm, cd eps_core, cd eps_shell, double n_medium, double nu_cm,
                            int n_max = 0);

// Динамический множитель мультиполя l: f_l = a_l(x)/[a_l(x_s)(x/x_s)^{2l+1}], x_s ≪ x — отношение точного
// коэффициента Ми к его квазистатическому пределу (радиационное затухание, динамическая деполяризация).
// Для l > l_dyn (f_l − 1 ~ x² и мельче) возвращается 1.
std::vector<cd> mie_dynamic_factors(const SphereLayers& s, cd eps_core, double n_medium, double nu_cm, int l_max,
                                    int l_dyn = 8);

}  // namespace ssp::plasmonics
