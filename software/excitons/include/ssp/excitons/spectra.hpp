#pragma once
// Секулярные кумулянтные формы линий экситонных переходов (поглощение и флуоресценция).
//   L_abs,n(ω) = (1/π) Re∫₀^T exp[i(ω − ε_n)t − s_n g(t) − Γ_n t/2] dt
//   L_em,n(ω)  = (1/π) Re∫₀^T exp[i(ω − ε_n + 2 s_n λ)t − s_n g*(t) − Γ_n t/2] dt
// s_n = Σ_i c_in⁴, Γ_n — полная скорость ухода из состояния (угл. см⁻¹). Площадь каждой линии = 1.

#include "ssp/excitons/bath.hpp"

#include <complex>
#include <vector>

namespace ssp::excitons {

class LineshapeEngine {
 public:
  // half_range_cm — максимальное |ω − ε_n|, которое нужно разрешить; t_max_tu — длина временного окна.
  LineshapeEngine(const Bath& bath, double half_range_cm, double t_max_tu);
  ~LineshapeEngine();
  LineshapeEngine(const LineshapeEngine&) = delete;
  LineshapeEngine& operator=(const LineshapeEngine&) = delete;

  // Добавить к out[k] вклад weight·(w_k)^power·L_n(w_k) на сетке w (см⁻¹).
  void add_line(std::vector<double>& out, const std::vector<double>& w, double center_cm, double weight,
                double s, double Gamma_cm, bool emission, int power) const;

  // Комплексный отклик перехода X(ω) = ∫₀^T exp[i(ω − ε)t − s g(t) − Γt/2] dt на сетке w (поглощение).
  // Re X = π·L(ω); Im X — дисперсионная часть (Крамерс–Крониг). Вне разрешённого диапазона — асимптотика i/(ω − ε).
  // Резольвента экситона в линейном отклике: (ε − ω − iΓ/2)⁻¹ ↔ i·X(ω).
  std::vector<std::complex<double>> complex_response(const std::vector<double>& w, double center_cm, double s,
                                                     double Gamma_cm) const;

  // |F(T)| — остаток подынтегральной функции в конце окна (контроль обрезания).
  double tail(double s, double Gamma_cm) const;

  double dt() const { return dt_; }
  std::size_t n() const { return n_; }
  double lambda_total() const { return lambda_; }
  const std::vector<std::complex<double>>& g() const { return g_; }

 private:
  double dt_ = 0.0;
  std::size_t n_ = 0;
  double lambda_ = 0.0;
  std::vector<std::complex<double>> g_;
  void* plan_ = nullptr;  // fftw_plan
};

// Абсолютное сечение поглощения (нм²) одного перехода с |μ|² (Д²) и нормированной линией L (см):
// σ = π ν̃ |μ|² L / (3 ε₀ ħ c n_med).
double absorption_cross_section_prefactor_nm2(double n_medium);

}  // namespace ssp::excitons
