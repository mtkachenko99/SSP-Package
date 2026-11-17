#include "ssp/excitons/monomer_fit.hpp"

#include "ssp/excitons/spectra.hpp"
#include "ssp/units.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <tuple>
#include <stdexcept>

namespace ssp::excitons {

namespace {

constexpr double kNA = 6.02214076e23;

// g(t) моды с S = 1 (g линейна по S) — БПФ на длинном окне (разрешение спектральной плотности ≪ γ),
// мелкий шаг; кэшируется по (ω, γ, T).
struct UnitG {
  double dt = 0.0;
  std::vector<std::complex<double>> g;
};
const UnitG& unit_mode_g(double omega, double gamma, double temperature) {
  static std::mutex mtx;
  static std::map<std::tuple<double, double, double>, UnitG> cache;
  std::lock_guard<std::mutex> lock(mtx);
  const auto key = std::make_tuple(omega, gamma, temperature);
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  Bath b;
  b.type = BathType::DrudeLorentz;
  b.lambda_cm = 0.0;
  b.temperature_K = temperature;
  VibrationalMode m;
  m.omega_cm = omega;
  m.huang_rhys = 1.0;
  m.gamma_cm = gamma;
  b.modes = {m};
  UnitG u;
  u.dt = units::pi / (40.0 * std::max(omega, 100.0));  // ~80 точек на период
  const double T = std::clamp(60.0 / std::max(gamma, 1.0), 2.0, 20.0);
  // окно БПФ задаёт разрешение по ω, поэтому оно длинное (T ≫ 1/γ)
  u.g = lineshape_function(b, u.dt, static_cast<std::size_t>(T / u.dt) + 2, 2);
  return cache.emplace(key, std::move(u)).first->second;
}

double trapz(const std::vector<double>& x, const std::vector<double>& y) {
  double s = 0.0;
  for (std::size_t i = 1; i < x.size(); ++i) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return s;
}


// Вектор свободных параметров (положительные — в логарифме) ↔ MonomerFitParams
struct Packer {
  const MonomerFitParams& base;
  std::vector<std::string> names;
  explicit Packer(const MonomerFitParams& b) : base(b) {
    if (b.fit_E0) names.push_back("E0_cm");
    if (b.fit_lambda) names.push_back("lambda_cm");
    if (b.fit_cutoff) names.push_back("cutoff_cm");
    if (b.fit_sigma) names.push_back("sigma_inh_cm");
    for (std::size_t k = 0; k < b.modes.size(); ++k) {
      if (b.fit_mode_S) names.push_back("mode" + std::to_string(k) + "_S");
      if (b.fit_mode_omega) names.push_back("mode" + std::to_string(k) + "_omega_cm");
    }
  }
  Eigen::VectorXd pack(const MonomerFitParams& p) const {
    Eigen::VectorXd x(static_cast<Eigen::Index>(names.size()));
    Eigen::Index i = 0;
    if (base.fit_E0) x(i++) = p.E0_cm;
    if (base.fit_lambda) x(i++) = std::log(std::max(p.lambda_cm, 1e-3));
    if (base.fit_cutoff) x(i++) = std::log(std::max(p.cutoff_cm, 1e-3));
    if (base.fit_sigma) x(i++) = std::log(std::max(p.sigma_inh_cm, 1e-3));
    for (const auto& m : p.modes) {
      if (base.fit_mode_S) x(i++) = std::log(std::max(m.huang_rhys, 1e-4));
      if (base.fit_mode_omega) x(i++) = m.omega_cm;
    }
    return x;
  }
  MonomerFitParams unpack(const Eigen::VectorXd& x) const {
    MonomerFitParams p = base;
    Eigen::Index i = 0;
    if (base.fit_E0) p.E0_cm = x(i++);
    // физические границы (за ними модель «замораживается», якобиан → 0)
    if (base.fit_lambda) p.lambda_cm = std::clamp(std::exp(x(i++)), 1.0, 5000.0);
    if (base.fit_cutoff) p.cutoff_cm = std::clamp(std::exp(x(i++)), 5.0, 5000.0);
    if (base.fit_sigma) p.sigma_inh_cm = std::clamp(std::exp(x(i++)), 1.0, 3000.0);
    for (auto& m : p.modes) {
      if (base.fit_mode_S) m.huang_rhys = std::clamp(std::exp(x(i++)), 1e-3, 5.0);
      if (base.fit_mode_omega) m.omega_cm = x(i++);
    }
    return p;
  }
  // значение и погрешность в «физических» единицах
  void physical(const Eigen::VectorXd& x, const Eigen::VectorXd& sx, std::vector<double>& v,
                std::vector<double>& e) const {
    v.clear();
    e.clear();
    const MonomerFitParams p = unpack(x);  // значения — с учётом границ
    std::vector<double> all;
    if (base.fit_E0) all.push_back(p.E0_cm);
    if (base.fit_lambda) all.push_back(p.lambda_cm);
    if (base.fit_cutoff) all.push_back(p.cutoff_cm);
    if (base.fit_sigma) all.push_back(p.sigma_inh_cm);
    for (const auto& m : p.modes) {
      if (base.fit_mode_S) all.push_back(m.huang_rhys);
      if (base.fit_mode_omega) all.push_back(m.omega_cm);
    }
    for (Eigen::Index i = 0; i < x.size(); ++i) {
      const auto& n = names[static_cast<std::size_t>(i)];
      const bool lin = n == "E0_cm" || n.find("_omega_cm") != std::string::npos;
      v.push_back(all[static_cast<std::size_t>(i)]);
      e.push_back(lin ? sx(i) : all[static_cast<std::size_t>(i)] * sx(i));
    }
  }
};

}  // namespace

void monomer_model(const MonomerFitParams& p, const std::vector<double>& w, std::vector<double>& abs,
                   std::vector<double>& em) {
  Bath b;
  b.type = BathType::DrudeLorentz;
  b.lambda_cm = p.lambda_cm;
  b.cutoff_cm = p.cutoff_cm;
  b.modes = p.modes;
  b.temperature_K = p.temperature_K;
  const double wmin = *std::min_element(w.begin(), w.end()), wmax = *std::max_element(w.begin(), w.end());
  // шаг по времени не зависит от параметров (иначе дискретизация вносит шум в численный якобиан) и в 4 раза
  // мельче предела Найквиста — погрешность трапеций на краях сетки ≪ 1 %
  const double half = 4.0 * ((wmax - wmin) + 4000.0);
  // Неоднородное уширение G_σ — свёртка по частоте = множитель exp(−σ²t²/2) во времени, поэтому форма линии
  // считается прямым преобразованием Фурье сразу в точках w (без сетки, свёртки и интерполяции):
  //   A(ω) = ω·(1/π) Re Σ_j exp[i(ω − E0)t_j − g(t_j) − σ²t_j²/2] dt,  F(ω) = ω³·… с e^{2iλt} и g*.
  const double sg = p.sigma_inh_cm;
  const double dt = units::pi / half;
  const double beta = b.beta();
  const double kT = std::isinf(beta) ? 0.0 : 1.0 / beta;
  // окно: пока exp(−Re g − σ²t²/2) > 10⁻⁹ (оценка по высокотемпературной дефазировке Друде и σ)
  const double rate = std::max(2.0 * b.lambda_cm * std::max(kT, 0.5 * b.cutoff_cm) / b.cutoff_cm, 1.0);
  double tmax = std::min(30.0 / rate, 7.0 / std::max(sg, 1.0));
  tmax = std::clamp(tmax, 50.0 * dt, 6.0);
  std::size_t nt = static_cast<std::size_t>(tmax / dt) + 1;
  std::vector<std::complex<double>> g(nt, 0.0);
  for (std::size_t j = 1; j < nt; ++j)
    g[j] = lineshape_function_drude(b.lambda_cm, b.cutoff_cm, beta, dt * static_cast<double>(j), 2000);
  double lam = b.lambda_cm;
  for (const auto& m : p.modes) {
    const auto& gu = unit_mode_g(m.omega_cm, m.gamma_cm, p.temperature_K);
    const double h = gu.dt;
    for (std::size_t j = 1; j < nt; ++j) {
      const double x = dt * static_cast<double>(j) / h;
      const std::size_t k = static_cast<std::size_t>(x);
      if (k + 1 >= gu.g.size()) {
        g[j] += m.huang_rhys * gu.g.back();
        continue;
      }
      const double f = x - static_cast<double>(k);
      g[j] += m.huang_rhys * ((1.0 - f) * gu.g[k] + f * gu.g[k + 1]);
    }
    lam += m.huang_rhys * m.omega_cm;
  }
  std::vector<std::complex<double>> fa(nt), fe(nt);
  for (std::size_t j = 0; j < nt; ++j) {
    const double t = dt * static_cast<double>(j), wj = (j == 0 ? 0.5 : 1.0) * dt / units::pi;
    const double env = std::exp(-0.5 * sg * sg * t * t);
    fa[j] = wj * env * std::exp(-g[j]);
    fe[j] = wj * env * std::exp(std::complex<double>(0.0, 2.0 * lam * t) - std::conj(g[j]));
  }
  while (nt > 2 && std::abs(fa[nt - 1]) < 1e-12 * dt) --nt;
  abs.assign(w.size(), 0.0);
  em.assign(w.size(), 0.0);
#pragma omp parallel for schedule(static) if (w.size() > 256)
  for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(w.size()); ++i) {
    const double x = w[static_cast<std::size_t>(i)];
    const std::complex<double> step_ph = std::exp(std::complex<double>(0.0, (x - p.E0_cm) * dt));
    std::complex<double> ph = 1.0, sa = 0.0, se = 0.0;
    for (std::size_t j = 0; j < nt; ++j) {
      sa += fa[j] * ph;
      se += fe[j] * ph;
      ph *= step_ph;
    }
    abs[static_cast<std::size_t>(i)] = x * sa.real();
    em[static_cast<std::size_t>(i)] = x * x * x * se.real();
  }
}
double monomer_stokes_peak(const MonomerFitParams& p) {
  std::vector<double> u, A, E;
  for (double x = p.E0_cm - 6000.0; x <= p.E0_cm + 3000.0; x += 2.0) u.push_back(x);
  monomer_model(p, u, A, E);
  // максимум с параболической интерполяцией по трём точкам — гладкая функция параметров (для якобиана)
  auto peak = [&](const std::vector<double>& y) {
    const std::size_t i = static_cast<std::size_t>(std::max_element(y.begin(), y.end()) - y.begin());
    if (i == 0 || i + 1 >= y.size()) return u[i];
    const double d = y[i - 1] - 2.0 * y[i] + y[i + 1];
    return u[i] + (std::abs(d) > 0 ? 0.5 * (y[i - 1] - y[i + 1]) / d * 2.0 : 0.0);
  };
  return peak(A) - peak(E);
}

double mu2_from_extinction(const Spectrum& eps, double n_medium) {
  std::vector<double> f(eps.w.size());
  for (std::size_t i = 0; i < f.size(); ++i) {
    const double sigma_nm2 = eps.y[i] * std::log(10.0) * 1000.0 / kNA * 1e14;  // см² → нм²
    f[i] = sigma_nm2 / eps.w[i];
  }
  return trapz(eps.w, f) / absorption_cross_section_prefactor_nm2(n_medium);
}

MonomerFitResult fit_monomer(const Spectrum& abs, const Spectrum& em, const MonomerFitParams& start,
                             bool absolute_extinction, double n_medium, int max_iter) {
  if (abs.w.size() < 5) throw std::invalid_argument("monomer fit: absorption spectrum needs ≥ 5 points");
  MonomerFitParams st = start;
  if (start.moment_start) {
    auto moments = [](const Spectrum& sp, int power, double& m1, double& var) {
      double s0 = 0, s1 = 0, s2 = 0;
      const double thr = 0.05 * *std::max_element(sp.y.begin(), sp.y.end());  // крылья с шумом — не учитывать
      for (std::size_t i = 1; i < sp.w.size(); ++i) {
        const double dw = sp.w[i] - sp.w[i - 1], x = 0.5 * (sp.w[i] + sp.w[i - 1]);
        const double yy = 0.5 * (sp.y[i] + sp.y[i - 1]);
        const double y = (yy > thr ? yy : 0.0) / std::pow(x, power);
        s0 += y * dw;
        s1 += y * x * dw;
        s2 += y * x * x * dw;
      }
      m1 = s1 / s0;
      var = s2 / s0 - m1 * m1;
    };
    double ma = 0, va = 0;
    moments(abs, 1, ma, va);
    st.E0_cm = ma;
    double lmodes = 0, vmodes = 0;
    Bath tb;
    tb.temperature_K = start.temperature_K;
    const double beta = tb.beta();
    for (const auto& m : st.modes) {
      lmodes += m.huang_rhys * m.omega_cm;
      vmodes += m.huang_rhys * m.omega_cm * m.omega_cm / std::tanh(0.5 * beta * m.omega_cm);
    }
    if (em.w.size() > 4) {
      double me = 0, ve = 0;
      moments(em, 3, me, ve);
      st.lambda_cm = std::max(10.0, 0.5 * (ma - me) - lmodes);
    }
    const double kT = std::isinf(beta) ? 0.0 : 1.0 / beta;
    const double vhom = 2.0 * st.lambda_cm * kT + vmodes;  // высокотемпературная оценка для Друде
    st.sigma_inh_cm = std::sqrt(std::max(100.0, va - vhom));
  }
  const Packer pk(st);
  const std::size_t na = abs.w.size(), ne = em.w.size();
  auto weights = [](const Spectrum& s) {
    std::vector<double> sg = s.sigma;
    if (sg.size() != s.y.size()) {
      const double mx = s.y.empty() ? 1.0 : *std::max_element(s.y.begin(), s.y.end());
      sg.assign(s.y.size(), 0.02 * mx);
    }
    for (auto& v : sg) v = std::max(v, 1e-300);
    return sg;
  };
  const auto sa = weights(abs), se = weights(em);
  // остатки при оптимальных масштабах
  auto residuals = [&](const Eigen::VectorXd& x, double* sca, double* sce, std::vector<double>* ma,
                       std::vector<double>* me) {
    const auto p = pk.unpack(x);
    std::vector<double> wa = abs.w, A, E, Ad, Ed;
    monomer_model(p, wa, A, Ed);
    std::vector<double> EA;
    if (ne) monomer_model(p, em.w, Ad, E);
    auto scale = [](const std::vector<double>& m, const std::vector<double>& y, const std::vector<double>& s) {
      double num = 0, den = 0;
      for (std::size_t i = 0; i < m.size(); ++i) {
        num += m[i] * y[i] / (s[i] * s[i]);
        den += m[i] * m[i] / (s[i] * s[i]);
      }
      return den > 0 ? num / den : 0.0;
    };
    const double ka = scale(A, abs.y, sa), ke = ne ? scale(E, em.y, se) : 0.0;
    const bool stokes = start.stokes_peak_cm > 0.0 && start.stokes_sigma_cm > 0.0;
    Eigen::VectorXd r(static_cast<Eigen::Index>(na + ne + (stokes ? 1 : 0)));
    for (std::size_t i = 0; i < na; ++i) r(static_cast<Eigen::Index>(i)) = (ka * A[i] - abs.y[i]) / sa[i];
    for (std::size_t i = 0; i < ne; ++i) r(static_cast<Eigen::Index>(na + i)) = (ke * E[i] - em.y[i]) / se[i];
    // сдвиг пиков — с параболической интерполяцией максимумов (иначе ступенчатый якобиан)
    if (stokes) r(static_cast<Eigen::Index>(na + ne)) = (monomer_stokes_peak(p) - start.stokes_peak_cm) / start.stokes_sigma_cm;
    if (sca) *sca = ka;
    if (sce) *sce = ke;
    if (ma) {
      ma->resize(na);
      for (std::size_t i = 0; i < na; ++i) (*ma)[i] = ka * A[i];
    }
    if (me) {
      me->resize(ne);
      for (std::size_t i = 0; i < ne; ++i) (*me)[i] = ke * E[i];
    }
    return r;
  };
  Eigen::VectorXd x = pk.pack(st);
  const Eigen::Index np = x.size();
  Eigen::VectorXd r = residuals(x, nullptr, nullptr, nullptr, nullptr);
  double chi = r.squaredNorm();
  MonomerFitResult out;
  Eigen::MatrixXd Jm;
  auto jacobian = [&](const Eigen::VectorXd& x0, const Eigen::VectorXd& r0) {
    Eigen::MatrixXd J(r0.size(), np);
    for (Eigen::Index j = 0; j < np; ++j) {
      const bool lin = pk.names[static_cast<std::size_t>(j)] == "E0_cm" ||
                       pk.names[static_cast<std::size_t>(j)].find("_omega_cm") != std::string::npos;
      const double h = lin ? 1.0 : 1e-3;
      Eigen::VectorXd xp = x0;
      xp(j) += h;
      J.col(j) = (residuals(xp, nullptr, nullptr, nullptr, nullptr) - r0) / h;
    }
    return J;
  };
  // Левенберг–Марквардт из точки x (на месте); возвращает число итераций, converged — признак остановки
  auto lm = [&](Eigen::VectorXd& xx, Eigen::VectorXd& rr, double& cc, int iters, bool& conv) {
    double lmu = 1e-2;
    conv = false;
    int k = 0;
    for (; k < iters && np > 0; ++k) {
      const Eigen::MatrixXd J = jacobian(xx, rr);
      const Eigen::MatrixXd JtJ = J.transpose() * J;
      const Eigen::VectorXd g = J.transpose() * rr;
      bool improved = false;
      double rel = 1.0;
      for (int tries = 0; tries < 12; ++tries) {
        Eigen::MatrixXd Aug = JtJ;
        Aug.diagonal() += lmu * JtJ.diagonal().cwiseMax(1e-12);
        const Eigen::VectorXd xn = xx - Aug.ldlt().solve(g);
        const Eigen::VectorXd rn = residuals(xn, nullptr, nullptr, nullptr, nullptr);
        const double cn = rn.squaredNorm();
        if (std::isfinite(cn) && cn < cc) {
          rel = (cc - cn) / cc;
          xx = xn;
          rr = rn;
          cc = cn;
          lmu = std::max(lmu / 3.0, 1e-9);
          improved = true;
          break;
        }
        lmu *= 4.0;
      }
      if (!improved || rel < 1e-7) {
        conv = true;
        break;
      }
    }
    return k;
  };
  int it = 0;
  bool conv = false;
  if (start.moment_start && np > 0) {
    // Мультистарт: моменты фиксируют E0, λ_tot и σ, но не разделяют вклад Друде (ω_c) и мод (S) —
    // 9 коротких запусков по сетке ω_c × S, лучший доводится до сходимости.
    Eigen::VectorXd best = x;
    double bestc = chi;
    for (const double wc : {50.0, 150.0, 400.0})
      for (const double sf : {0.5, 1.0, 2.0}) {
        MonomerFitParams c = st;
        c.cutoff_cm = wc;
        double lm_modes = 0;
        for (auto& m : c.modes) {
          m.huang_rhys = std::clamp(m.huang_rhys * sf, 1e-3, 5.0);
          lm_modes += m.huang_rhys * m.omega_cm;
        }
        // λ Друде — так, чтобы λ_tot (стоксов сдвиг) остался прежним
        double lm0 = 0;
        for (const auto& m : st.modes) lm0 += m.huang_rhys * m.omega_cm;
        c.lambda_cm = std::max(10.0, st.lambda_cm + lm0 - lm_modes);
        Eigen::VectorXd xc = pk.pack(c), rc = residuals(xc, nullptr, nullptr, nullptr, nullptr);
        double cc = rc.squaredNorm();
        bool cv = false;
        lm(xc, rc, cc, 25, cv);
        if (std::isfinite(cc) && cc < bestc) {
          bestc = cc;
          best = xc;
        }
      }
    x = best;
    r = residuals(x, nullptr, nullptr, nullptr, nullptr);
    chi = r.squaredNorm();
  }
  it = lm(x, r, chi, max_iter, conv);
  out.converged = conv;
  if (np > 0) Jm = jacobian(x, r);
  const double dof = std::max<double>(1.0, double(r.size()) - double(np) - (ne ? 2.0 : 1.0));
  out.stokes_peak_cm = monomer_stokes_peak(pk.unpack(x));
  out.chi2_red = chi / dof;
  Eigen::VectorXd sx = Eigen::VectorXd::Zero(np);
  if (np > 0) {
    // псевдообратная: параметр на границе (нулевой столбец якобиана) получает погрешность 0, а не NaN
    const Eigen::MatrixXd cov =
        (Jm.transpose() * Jm).completeOrthogonalDecomposition().pseudoInverse() * out.chi2_red;
    sx = cov.diagonal().cwiseMax(0.0).cwiseSqrt();
  }
  out.p = pk.unpack(x);
  out.names = pk.names;
  pk.physical(x, sx, out.values, out.errors);
  out.iterations = it;
  residuals(x, &out.abs_scale, &out.em_scale, &out.abs_model, &out.em_model);

  // ⟨ν̃⁻³⟩^{−1/3} по модельному спектру ФЛ на широкой сетке
  {
    std::vector<double> u, A, E;
    for (double v = out.p.E0_cm - 8000.0; v <= out.p.E0_cm + 4000.0; v += 5.0) u.push_back(v);
    monomer_model(out.p, u, A, E);
    double s0 = 0, s3 = 0;
    for (std::size_t i = 0; i < u.size(); ++i) {
      s0 += E[i];
      s3 += E[i] / (u[i] * u[i] * u[i]);
    }
    out.nu_em_cm = s3 > 0 ? std::cbrt(s0 / s3) : out.p.E0_cm;
  }
  if (absolute_extinction) {
    const double mu2 = mu2_from_extinction(abs, n_medium);
    out.mu_D = std::sqrt(std::max(mu2, 0.0));
    using namespace units;
    const double w = rad_s_per_cm * out.nu_em_cm;
    const double g = n_medium * w * w * w * mu2 * debye * debye / (3.0 * pi * eps0 * hbar * c * c * c);  // с⁻¹
    out.tau_rad_ns = g > 0 ? 1e9 / g : 0.0;
    std::vector<double> f(abs.w.size());
    for (std::size_t i = 0; i < f.size(); ++i) f[i] = abs.y[i] / abs.w[i];
    const double kr = 2.880e-9 * n_medium * n_medium * std::pow(out.nu_em_cm, 3) * trapz(abs.w, f);
    out.tau_sb_ns = kr > 0 ? 1e9 / kr : 0.0;
  }
  return out;
}

}  // namespace ssp::excitons
