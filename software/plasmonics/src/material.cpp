#include "ssp/plasmonics/material.hpp"

#include "ssp/units.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ssp::plasmonics {

namespace {
std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r");
  if (a == std::string::npos) return {};
  const auto b = s.find_last_not_of(" \t\r");
  return s.substr(a, b - a + 1);
}
std::size_t indent(const std::string& s) { return s.find_first_not_of(' '); }
}  // namespace

Material Material::from_refractiveindex_yaml(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open material file: " + path);
  Material m;
  m.tabulated_ = true;
  m.name_ = path;
  std::string line;
  enum class St { None, Refs, SeekData, Data } st = St::None;
  bool have_nk = false;
  std::size_t data_indent = 0;
  while (std::getline(in, line)) {
    const std::string t = trim(line);
    if (st == St::Data) {
      if (t.empty()) continue;
      if (indent(line) <= data_indent) break;  // конец блока данных
      std::istringstream ss(t);
      double l, n, k;
      if (!(ss >> l >> n >> k)) throw std::runtime_error("cannot parse data line: " + line);
      m.lam_.push_back(l);
      m.n_.push_back(n);
      m.k_.push_back(k);
      continue;
    }
    if (t.rfind("REFERENCES:", 0) == 0) { st = St::Refs; continue; }
    if (st == St::Refs) {
      if (!line.empty() && line[0] != ' ') st = St::None;
      else { m.reference_ += t + " "; continue; }
    }
    if (t.rfind("- type:", 0) == 0) {
      have_nk = (trim(t.substr(7)) == "tabulated nk");
      st = have_nk ? St::SeekData : St::None;
      continue;
    }
    if (st == St::SeekData && t.rfind("data:", 0) == 0) {
      st = St::Data;
      data_indent = indent(line);
    }
  }
  if (m.lam_.size() < 2) throw std::runtime_error("no 'tabulated nk' block in the file: " + path);
  // Сортировка по длине волны (на случай неупорядоченных таблиц).
  std::vector<std::size_t> idx(m.lam_.size());
  for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return m.lam_[a] < m.lam_[b]; });
  auto reorder = [&](std::vector<double>& v) {
    std::vector<double> o(v.size());
    for (std::size_t i = 0; i < idx.size(); ++i) o[i] = v[idx[i]];
    v.swap(o);
  };
  reorder(m.lam_);
  reorder(m.n_);
  reorder(m.k_);
  return m;
}

Material Material::drude(const DrudeParams& p, const std::string& name) {
  Material m;
  m.tabulated_ = false;
  m.drude_ = p;
  m.name_ = name;
  return m;
}

void Material::set_size_damping(double A, double vF, double R_eff_nm, double omega_p_eV, double gamma_bulk_eV) {
  if (A <= 0.0) {
    kreibig_extra_eV_ = 0.0;
    return;
  }
  // ħ·A·v_F/R в эВ
  kreibig_extra_eV_ = A * units::hbar * vF / (R_eff_nm * 1e-9) / units::e_charge;
  kreibig_wp_eV_ = omega_p_eV;
  kreibig_gb_eV_ = gamma_bulk_eV;
}

std::complex<double> Material::epsilon(double nu_cm) const {
  using cd = std::complex<double>;
  const double E_eV = nu_cm / units::cm_per_eV;
  cd eps;
  if (tabulated_) {
    const double lam = 1e4 / nu_cm;  // мкм
    if (lam < lam_.front() || lam > lam_.back())
      throw std::out_of_range("wavelength " + std::to_string(lam * 1e3) + " nm is outside the material table " + name_);
    // Монотонная кубическая интерполяция (PCHIP, Фритч–Карлсон) n(λ) и k(λ): гладкая, без ложных экстремумов.
    const auto it = std::upper_bound(lam_.begin(), lam_.end(), lam);
    const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(it - lam_.begin()), lam_.size() - 1);
    const std::size_t i0 = i - 1;
    auto pchip = [&](const std::vector<double>& y) {
      const std::size_t n = lam_.size();
      auto slope = [&](std::size_t j) { return (y[j + 1] - y[j]) / (lam_[j + 1] - lam_[j]); };
      auto deriv = [&](std::size_t j) {
        if (j == 0) return slope(0);
        if (j == n - 1) return slope(n - 2);
        const double d0 = slope(j - 1), d1 = slope(j);
        if (d0 * d1 <= 0.0) return 0.0;
        const double h0 = lam_[j] - lam_[j - 1], h1 = lam_[j + 1] - lam_[j];
        const double w1 = 2 * h1 + h0, w2 = h1 + 2 * h0;
        return (w1 + w2) / (w1 / d0 + w2 / d1);
      };
      const double h = lam_[i] - lam_[i0], t = (lam - lam_[i0]) / h;
      const double h00 = (1 + 2 * t) * (1 - t) * (1 - t), h10 = t * (1 - t) * (1 - t), h01 = t * t * (3 - 2 * t),
                   h11 = t * t * (t - 1);
      return h00 * y[i0] + h10 * h * deriv(i0) + h01 * y[i] + h11 * h * deriv(i);
    };
    const cd N(pchip(n_), pchip(k_));
    eps = N * N;
  } else {
    const double wp = drude_.omega_p_eV, g = drude_.gamma_eV;
    eps = drude_.eps_inf - wp * wp / cd(E_eV * E_eV, g * E_eV);
  }
  if (kreibig_extra_eV_ > 0.0) {
    const double wp2 = kreibig_wp_eV_ * kreibig_wp_eV_;
    eps += wp2 / cd(E_eV * E_eV, kreibig_gb_eV_ * E_eV) -
           wp2 / cd(E_eV * E_eV, (kreibig_gb_eV_ + kreibig_extra_eV_) * E_eV);
  }
  return eps;
}

}  // namespace ssp::plasmonics
