#include "ssp/excitons/spectra.hpp"

#include "ssp/units.hpp"

#include <fftw3.h>

#include <cmath>
#include <mutex>
#include <stdexcept>

namespace ssp::excitons {

namespace {
std::mutex& fftw_planner_mutex() {
  static std::mutex m;
  return m;
}
}  // namespace

LineshapeEngine::LineshapeEngine(const Bath& bath, double half_range_cm, double t_max_tu) {
  if (half_range_cm <= 0.0 || t_max_tu <= 0.0) throw std::invalid_argument("LineshapeEngine: invalid parameters");
  dt_ = units::pi / half_range_cm;
  n_ = 1;
  while (static_cast<double>(n_) * dt_ < t_max_tu) n_ <<= 1;
  lambda_ = bath.reorganization();
  g_ = lineshape_function(bath, dt_, n_);
  std::lock_guard<std::mutex> lock(fftw_planner_mutex());
  fftw_complex* buf = fftw_alloc_complex(n_);
  plan_ = fftw_plan_dft_1d(static_cast<int>(n_), buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_free(buf);
}

LineshapeEngine::~LineshapeEngine() {
  std::lock_guard<std::mutex> lock(fftw_planner_mutex());
  if (plan_) fftw_destroy_plan(static_cast<fftw_plan>(plan_));
}

double LineshapeEngine::tail(double s, double Gamma_cm) const {
  const double t = dt_ * static_cast<double>(n_ - 1);
  return std::exp(-s * g_.back().real() - 0.5 * Gamma_cm * t);
}

void LineshapeEngine::add_line(std::vector<double>& out, const std::vector<double>& w, double center_cm,
                               double weight, double s, double Gamma_cm, bool emission, int power) const {
  fftw_complex* X = fftw_alloc_complex(n_);
  const double shift = emission ? 2.0 * s * lambda_ : 0.0;
  for (std::size_t j = 0; j < n_; ++j) {
    const double t = dt_ * static_cast<double>(j);
    const std::complex<double> gj = emission ? std::conj(g_[j]) : g_[j];
    std::complex<double> F = std::exp(-s * gj - 0.5 * Gamma_cm * t + std::complex<double>(0.0, shift * t));
    if (j == 0) F *= 0.5;  // трапеция
    X[j][0] = F.real();
    X[j][1] = F.imag();
  }
  fftw_execute_dft(static_cast<fftw_plan>(plan_), X, X);
  // Δ_k = k·dΔ, k ∈ [−n/2, n/2); L(Δ) = dt·Re X/π.
  const double dD = 2.0 * units::pi / (static_cast<double>(n_) * dt_);
  const double half = 0.5 * static_cast<double>(n_);
  auto L_at = [&](long k) {
    const std::size_t idx = static_cast<std::size_t>((k % static_cast<long>(n_) + static_cast<long>(n_)) %
                                                     static_cast<long>(n_));
    return dt_ * X[idx][0] / units::pi;
  };
  for (std::size_t q = 0; q < w.size(); ++q) {
    const double x = (w[q] - center_cm) / dD;
    if (x <= -half || x >= half - 1.0) continue;  // вне разрешённого диапазона
    const long k0 = static_cast<long>(std::floor(x));
    const double f = x - static_cast<double>(k0);
    const double L = (1.0 - f) * L_at(k0) + f * L_at(k0 + 1);
    out[q] += weight * std::pow(w[q], power) * L;
  }
  fftw_free(X);
}

std::vector<std::complex<double>> LineshapeEngine::complex_response(const std::vector<double>& w, double center_cm,
                                                                    double s, double Gamma_cm) const {
  fftw_complex* X = fftw_alloc_complex(n_);
  for (std::size_t j = 0; j < n_; ++j) {
    const double t = dt_ * static_cast<double>(j);
    std::complex<double> F = std::exp(-s * g_[j] - 0.5 * Gamma_cm * t);
    if (j == 0) F *= 0.5;
    X[j][0] = F.real();
    X[j][1] = F.imag();
  }
  fftw_execute_dft(static_cast<fftw_plan>(plan_), X, X);
  const double dD = 2.0 * units::pi / (static_cast<double>(n_) * dt_);
  const double half = 0.5 * static_cast<double>(n_);
  auto at = [&](long k) {
    const std::size_t idx = static_cast<std::size_t>((k % static_cast<long>(n_) + static_cast<long>(n_)) %
                                                     static_cast<long>(n_));
    return std::complex<double>(dt_ * X[idx][0], dt_ * X[idx][1]);
  };
  std::vector<std::complex<double>> out(w.size());
  for (std::size_t q = 0; q < w.size(); ++q) {
    const double d = w[q] - center_cm;
    const double x = d / dD;
    if (x <= -half || x >= half - 1.0) {
      out[q] = std::complex<double>(0.0, 1.0 / d);
      continue;
    }
    const long k0 = static_cast<long>(std::floor(x));
    const double f = x - static_cast<double>(k0);
    out[q] = (1.0 - f) * at(k0) + f * at(k0 + 1);
  }
  fftw_free(X);
  return out;
}

double absorption_cross_section_prefactor_nm2(double n_medium) {
  using namespace units;
  // σ[м²] = π·ν̃[см⁻¹]·μ²[Д²]·debye²·L[см] / (3 ε₀ ħ c n). Здесь возвращается всё, кроме ν̃·μ²·L.
  return pi * debye * debye / (3.0 * eps0 * hbar * c * n_medium) * 1e18;
}

}  // namespace ssp::excitons
