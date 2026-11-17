#include "ssp/excitons/bath.hpp"

#include "ssp/units.hpp"

#include <fftw3.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ssp::excitons {

namespace {
constexpr double kPi = ssp::units::pi;
}

double Bath::beta() const {
  if (temperature_K <= 0.0) return std::numeric_limits<double>::infinity();
  return 1.0 / (units::kB_cm * temperature_K);
}

double Bath::J_over_w(double w) const {
  double v = 0.0;
  const double wc = cutoff_cm;
  if (lambda_cm > 0.0) {
    if (type == BathType::DrudeLorentz) v += 2.0 * lambda_cm * wc / (w * w + wc * wc);
    else v += (kPi * lambda_cm / wc) * std::exp(-w / wc);
  }
  for (const auto& m : modes) {
    const double lk = m.huang_rhys * m.omega_cm;
    const double d = m.omega_cm * m.omega_cm - w * w;
    v += 2.0 * lk * m.omega_cm * m.omega_cm * m.gamma_cm / (d * d + m.gamma_cm * m.gamma_cm * w * w);
  }
  return v;
}

double Bath::reorganization() const {
  double l = lambda_cm;
  for (const auto& m : modes) l += m.huang_rhys * m.omega_cm;
  return l;
}

double Bath::rate_spectrum(double w) const {
  const double b = beta();
  const double jw = J_over_w(std::abs(w));
  if (std::isinf(b)) return w > 0.0 ? 2.0 * jw * w : 0.0;
  if (std::abs(b * w) < 1e-12) return 2.0 * jw / b;
  return 2.0 * jw * w / (-std::expm1(-b * w));
}

double bose(double w, double beta) {
  if (w <= 0.0) return 0.0;
  if (std::isinf(beta)) return 0.0;
  const double x = beta * w;
  if (x > 700.0) return 0.0;
  return 1.0 / std::expm1(x);
}

std::vector<std::complex<double>> lineshape_function(const Bath& bath, double dt, std::size_t n, int oversample) {
  if (dt <= 0.0 || n < 2) throw std::invalid_argument("lineshape_function: invalid grid");
  // Друдевская часть: точная мацубаровская формула (у J ~ 1/ω хвост, C(t) сингулярна при t → 0 и плохо
  // интегрируется численно). Остальное (омическая часть, колебательные моды) — через БПФ.
  if (bath.type == BathType::DrudeLorentz && bath.lambda_cm > 0.0) {
    Bath rest = bath;
    rest.lambda_cm = 0.0;
    std::vector<std::complex<double>> g(n, 0.0);
    if (!rest.modes.empty()) g = lineshape_function(rest, dt, n, oversample);
    const double b = bath.beta();
    if (std::isinf(b)) throw std::invalid_argument("lineshape_function: a Drude bath requires T > 0");
    for (std::size_t j = 1; j < n; ++j)
      g[j] += lineshape_function_drude(bath.lambda_cm, bath.cutoff_cm, b, dt * static_cast<double>(j), 2000);
    return g;
  }
  bool any = bath.lambda_cm > 0.0;
  for (const auto& m : bath.modes) any = any || m.huang_rhys > 0.0;
  if (!any) return std::vector<std::complex<double>>(n, 0.0);
  // Подсетка τ_j = j·h, h = dt/oversample; период БПФ ≥ 2·T, чтобы избежать наложения.
  const double h = dt / oversample;
  const double T = dt * static_cast<double>(n - 1);
  std::size_t N = 1;
  while (static_cast<double>(N) * h < 2.0 * T + 4.0 * h) N <<= 1;
  const double dw = 2.0 * kPi / (static_cast<double>(N) * h);
  const double b = bath.beta();

  // a_k = J(ω)coth(βω/2), b_k = J(ω) на ω_k = k·dω; веса трапеций.
  fftw_complex* A = fftw_alloc_complex(N);
  fftw_complex* B = fftw_alloc_complex(N);
  for (std::size_t k = 0; k < N; ++k) {
    const double w = static_cast<double>(k) * dw;
    const double wt = (k == 0 || k == N - 1) ? 0.5 : 1.0;
    double a;
    if (k == 0) a = std::isinf(b) ? 0.0 : 2.0 * bath.J_over_w(0.0) / b;  // lim J·coth = 2J'(0)/β
    else a = bath.J(w) * (std::isinf(b) ? 1.0 : 1.0 / std::tanh(0.5 * b * w));
    A[k][0] = wt * a;
    A[k][1] = 0.0;
    B[k][0] = wt * bath.J(w);
    B[k][1] = 0.0;
  }
  fftw_plan pa = fftw_plan_dft_1d(static_cast<int>(N), A, A, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan pb = fftw_plan_dft_1d(static_cast<int>(N), B, B, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(pa);
  fftw_execute(pb);

  // C(τ_j) = (dω/π)[Re A_j − i Im B_j]; g = ∫∫C: двойное интегрирование трапециями.
  const std::size_t M = static_cast<std::size_t>(n - 1) * oversample + 1;
  std::vector<std::complex<double>> g_fine(M);
  std::complex<double> G1 = 0.0, g = 0.0;
  std::complex<double> C_prev(dw / kPi * A[0][0], -dw / kPi * B[0][1]);
  std::complex<double> G1_prev = 0.0;
  g_fine[0] = 0.0;
  for (std::size_t j = 1; j < M; ++j) {
    const std::complex<double> Cj(dw / kPi * A[j][0], -dw / kPi * B[j][1]);
    G1 = G1_prev + 0.5 * h * (C_prev + Cj);
    g = g + 0.5 * h * (G1_prev + G1);
    g_fine[j] = g;
    C_prev = Cj;
    G1_prev = G1;
  }
  fftw_destroy_plan(pa);
  fftw_destroy_plan(pb);
  fftw_free(A);
  fftw_free(B);

  std::vector<std::complex<double>> out(n);
  for (std::size_t j = 0; j < n; ++j) out[j] = g_fine[j * oversample];
  return out;
}

std::complex<double> lineshape_function_drude(double lambda, double gamma, double beta, double t, int n_matsubara) {
  // g(t) = (λ/γ)(cot(βγ/2) − i)(e^{−γt} + γt − 1) + (4λγ/β) Σ_k (e^{−ν_k t} + ν_k t − 1)/(ν_k(ν_k² − γ²)).
  using cd = std::complex<double>;
  const double x = std::exp(-gamma * t) + gamma * t - 1.0;
  cd g = (lambda / gamma) * cd(1.0 / std::tan(0.5 * beta * gamma), -1.0) * x;
  double s = 0.0;
  for (int k = 1; k <= n_matsubara; ++k) {
    const double nu = 2.0 * kPi * k / beta;
    s += (std::exp(-nu * t) + nu * t - 1.0) / (nu * (nu * nu - gamma * gamma));
  }
  // Хвост ряда: (ν t − 1)/ν³ при ν t ≫ 1, Σ_{k>K} 1/k² ≈ 1/K, Σ 1/k³ ≈ 1/(2K²).
  const double c = 2.0 * kPi / beta;
  const double K = n_matsubara;
  s += t / (c * c) * (1.0 / K - 0.5 / (K * K)) - 1.0 / (c * c * c) * (0.5 / (K * K));
  g += 4.0 * lambda * gamma / beta * s;
  return g;
}

}  // namespace ssp::excitons
