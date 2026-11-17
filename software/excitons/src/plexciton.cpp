#include "ssp/excitons/plexciton.hpp"

#include "ssp/excitons/relaxation.hpp"
#include "ssp/plasmonics/mie.hpp"
#include "ssp/units.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ssp::excitons {

using cd = std::complex<double>;

namespace {
// Поле единичного диполя μ̂ в точке r (нм) в начале координат: квазистатическое ĝ = (3n(n·μ̂) − μ̂)/r³ и
// запаздывающее g(k) = e^{ikr}[(1 − ikr)(3n(n·μ̂) − μ̂) + (kr)²(μ̂ − n(n·μ̂))]/r³ (n = −r/|r| — направление
// от диполя к центру; формула чётна по n).
void dipole_field_at_origin(const Eigen::Vector3d& r, const Eigen::Vector3d& mu, double k, Eigen::Vector3d& g0,
                            Eigen::Vector3cd& gk) {
  const double R = r.norm();
  const Eigen::Vector3d n = r / R;
  const Eigen::Vector3d q = 3.0 * n * n.dot(mu) - mu, t = mu - n * n.dot(mu);
  const double r3 = R * R * R;
  g0 = q / r3;
  const std::complex<double> ph = std::exp(std::complex<double>(0.0, k * R));
  gk = ph * ((1.0 - std::complex<double>(0.0, k * R)) * q.cast<std::complex<double>>() +
             (k * R) * (k * R) * t.cast<std::complex<double>>()) / r3;
}
std::complex<double> absorbing_part(std::complex<double> al, double k, bool dynamic) {
  return dynamic ? std::complex<double>(al.real(), al.imag() - 2.0 / 3.0 * k * k * k * std::norm(al)) : al;
}
}  // namespace

PlexcitonSystem::PlexcitonSystem(const Aggregate& agg, const NanoparticleSpec& spec, const plasmonics::Material& mat,
                                 double mu_D)
    : agg_(agg), spec_(spec), mat_(mat), mu_D_(mu_D), eps_m_(spec.n_medium * spec.n_medium) {
  const bool sphere = std::abs(spec.a_xy_nm - spec.a_z_nm) < 1e-9 * std::max(spec.a_xy_nm, spec.a_z_nm);
  if (!sphere && (spec.shell_thickness_nm > 0.0 || spec.d_perp_nm != cd(0.0)))
    throw std::invalid_argument("nanoparticle: ligand shell and Feibelman d⊥ are implemented for spheres only");
  particle_ = std::make_unique<plasmonics::QuasistaticParticle>(spec.outer_xy(), spec.outer_z(), spec.l_max);
  if (sphere) {
    layered_ = true;
    layers_.core_radius_nm = spec.a_xy_nm;
    if (spec.shell_thickness_nm > 0.0) {
      layers_.shell_thickness_nm = {spec.shell_thickness_nm};
      layers_.shell_eps = {cd(spec.shell_n * spec.shell_n)};
    }
    layers_.d_perp_nm = spec.d_perp_nm;
  }
  const auto N = static_cast<Eigen::Index>(agg_.size());
  for (Eigen::Index i = 0; i < N; ++i) {
    const double gap = particle_->surface_gap(agg_.positions[i]);
    if (particle_->inside(agg_.positions[i]) || gap < 0.05) {
      std::ostringstream os;
      os << "molecule " << i << " is inside the nanoparticle or on its surface (gap " << gap
         << " nm); move the aggregate (placement)";
      throw std::invalid_argument(os.str());
    }
  }
  // fast_l < 0 — автоматически: для N ≥ 600 точно до l, при котором (b/(b+gap_min))^{2l} < 10⁻², дальше — диагональ
  if (spec_.fast_l < 0) {
    spec_.fast_l = 0;
    if (N >= 600) {
      double gmin = 1e300;
      for (Eigen::Index i = 0; i < N; ++i) gmin = std::min(gmin, particle_->surface_gap(agg_.positions[i]));
      const double b = std::max(spec.outer_xy(), spec.outer_z());
      const int l = static_cast<int>(std::ceil(std::log(1e-2) / (2.0 * std::log(b / (b + gmin)))));
      // ранг низкой части 2T_c = l(l+3) не больше N/4 — иначе Вудбери не быстрее точного решения
      int lr = 1;
      while ((lr + 1) * (lr + 4) <= N / 4) ++lr;
      const int lf = std::min(l, lr);
      if (lf < particle_->l_max() && lf >= 4) spec_.fast_l = lf;
    }
  }
  const auto& terms = particle_->terms();
  const auto T = static_cast<Eigen::Index>(terms.size());
  Vc_.resize(N, T);
  Vs_.resize(N, T);
  g10c_.resize(N);
  g11c_.resize(N);
  g11s_.resize(N);
  const auto& ut = particle_->uniform_term();
#pragma omp parallel for schedule(dynamic)
  for (Eigen::Index i = 0; i < N; ++i) {
    std::vector<Eigen::Vector3cd> gc, gs;
    particle_->gradients(agg_.positions[i], gc, gs);
    const Eigen::Vector3cd d = agg_.dipoles[i].cast<cd>();
    for (Eigen::Index t = 0; t < T; ++t) {
      Vc_(i, t) = d.dot(gc[t]);  // dot() без сопряжения для вещественного d
      Vs_(i, t) = d.dot(gs[t]);
    }
    g11c_[i] = gc[ut[0]];
    g11s_[i] = gs[ut[1]];
    g10c_[i] = gc[ut[2]];
  }
  // Для сферы — свёртка по m заранее (если хватает памяти, ≤ 1.5 ГБ)
  const int L = particle_->l_max();
  if (particle_->shape() == plasmonics::ParticleShape::Sphere && spec_.fast_l == 0 &&
      double(N) * double(N) * L * sizeof(double) < 1.5e9) {
    Gl_.assign(static_cast<std::size_t>(L), Eigen::MatrixXd::Zero(N, N));
#pragma omp parallel for schedule(dynamic)
    for (int l = 1; l <= L; ++l) {
      const Eigen::Index t0 = Eigen::Index(l) * (l + 1) / 2 - 1, nt = l + 1;  // члены (l, m = 0..l)
      Eigen::VectorXd w(nt);
      for (Eigen::Index m = 0; m < nt; ++m) w(m) = terms[t0 + m].w.real();
      const Eigen::MatrixXd C = Vc_.middleCols(t0, nt).real(), S = Vs_.middleCols(t0, nt).real();
      Gl_[l - 1] = C * w.asDiagonal() * C.transpose() + S * w.asDiagonal() * S.transpose();
    }
  }
}

PlexcitonSystem::Freq PlexcitonSystem::at(double nu_cm, int l_cut, bool with_W) const {
  Freq F;
  const cd eps = mat_.epsilon(nu_cm);
  F.k = 2.0 * units::pi * spec_.n_medium * nu_cm * 1e-7;
  const auto& terms = particle_->terms();
  const auto T = static_cast<Eigen::Index>(terms.size());
  Eigen::VectorXcd d(T);
  cd D10 = 0.0, D11 = 0.0;
  if (layered_) {
    // Сфера: D̂_l по слоям (оболочка, d⊥), динамика каждого мультиполя — множитель Ми (без d⊥)
    const int L = particle_->l_max();
    std::vector<cd> Dl(static_cast<std::size_t>(L + 1));
    const auto f = spec_.dynamic_dipole ? plasmonics::mie_dynamic_factors(layers_, eps, spec_.n_medium, nu_cm, L)
                                        : std::vector<cd>(static_cast<std::size_t>(L + 1), cd(1.0));
    for (int l = 1; l <= L; ++l) Dl[l] = plasmonics::layered_response(l, eps, layers_, eps_m_) * f[l];
    const double b = layers_.outer_radius();
    F.alpha_np.setConstant(-Dl[1] * b * b * b);
    for (Eigen::Index t = 0; t < T; ++t) {
      const auto& tm = terms[t];
      d(t) = (l_cut > 0 && tm.l > l_cut) ? cd(0.0) : tm.w * Dl[tm.l];
    }
    D10 = D11 = Dl[1];
  } else {
    const Eigen::Vector3cd a_qs = particle_->polarizability(eps, eps_m_);
    Eigen::Vector3cd a_dyn = a_qs;
    if (spec_.dynamic_dipole) {
      // MLWA (Meier–Wokaun): α' = α/(1 − (k²/l_a)α − i(2/3)k³α)
      const Eigen::Vector3d ax = particle_->semi_axes();
      for (int a = 0; a < 3; ++a)
        a_dyn(a) = a_qs(a) / (1.0 - F.k * F.k / ax(a) * a_qs(a) - cd(0.0, 2.0 / 3.0) * std::pow(F.k, 3) * a_qs(a));
    }
    F.alpha_np = a_dyn;
    auto ratio = [](cd num, cd den) { return std::abs(den) > 1e-300 ? num / den : cd(1.0); };
    const cd s_xy = ratio(a_dyn(0), a_qs(0)), s_z = ratio(a_dyn(2), a_qs(2));
    for (Eigen::Index t = 0; t < T; ++t) {
      const auto& tm = terms[t];
      cd D = particle_->response(tm, eps, eps_m_);
      if (tm.l == 1) {
        D *= (tm.m == 0) ? s_z : s_xy;
        (tm.m == 0 ? D10 : D11) = D;
      }
      d(t) = (l_cut > 0 && tm.l > l_cut) ? cd(0.0) : tm.w * D;
    }
  }
  // Энергия взаимодействия через наведённое поле: W_ij = +(μ_i·∇)(μ_j·∇')K_ind/(4πε₀ε_m) (как и свободная J).
  const double pref = units::dd_cm_per_D2_nm3 * mu_D_ * mu_D_ / eps_m_;
  F.d = d;
  // Только поглощение в металле: из дипольных членов убирается радиационное затухание НЧ (оптическая теорема
  // Im α_abs = Im α − (2/3)k³|α|²); для l ≥ 2 оно ~k^{2l+1} и пренебрежимо.
  F.d_abs = d;
  if (spec_.dynamic_dipole) {  // в квазистатике радиационного затухания в α нет
    auto ratio = [&](cd al) {
      if (std::abs(al) < 1e-300) return cd(1.0);
      const cd aa(al.real(), al.imag() - 2.0 / 3.0 * std::pow(F.k, 3) * std::norm(al));
      return aa / al;
    };
    const cd rz = ratio(F.alpha_np(2)), rxy = ratio(F.alpha_np(0));
    for (Eigen::Index t = 0; t < T; ++t)
      if (terms[t].l == 1) F.d_abs(t) *= terms[t].m == 0 ? rz : rxy;
  }
  if (!with_W) {
  } else if (!Gl_.empty()) {
    const auto N = static_cast<Eigen::Index>(agg_.size());
    Eigen::MatrixXd Wr = Eigen::MatrixXd::Zero(N, N), Wi = Eigen::MatrixXd::Zero(N, N);
    for (int l = 1; l <= particle_->l_max(); ++l) {
      if (l_cut > 0 && l > l_cut) break;
      const Eigen::Index t0 = Eigen::Index(l) * (l + 1) / 2 - 1;
      const cd Dl = d(t0) / terms[t0].w;  // D̂_l (одинаков для всех m)
      Wr.noalias() += Dl.real() * Gl_[l - 1];
      Wi.noalias() += Dl.imag() * Gl_[l - 1];
    }
    F.W = pref * (Wr.cast<cd>() + cd(0.0, 1.0) * Wi.cast<cd>());
  } else {
    F.W = pref * (Vc_ * d.asDiagonal() * Vc_.transpose() + Vs_ * d.asDiagonal() * Vs_.transpose());
  }

  const auto& u = particle_->uniform_coeff();
  const auto N = agg_.size();
  F.L.resize(N);
  for (std::size_t i = 0; i < N; ++i) {
    Eigen::Matrix3cd L = Eigen::Matrix3cd::Identity();
    L.col(0) += u[0] * D11 * g11c_[i];
    L.col(1) += u[1] * D11 * g11s_[i];
    L.col(2) += u[2] * D10 * g10c_[i];
    F.L[i] = L;
  }
  return F;
}

void PlexcitonSystem::project(const Eigen::MatrixXd& U, Eigen::MatrixXcd& Pc, Eigen::MatrixXcd& Ps) const {
  // V вещественны для сферы и сфероида (действительные функции f), но храним комплексными
  const Eigen::MatrixXd Ut = U.transpose();
  Pc = (Ut * Vc_.real()).cast<cd>() + cd(0.0, 1.0) * (Ut * Vc_.imag()).cast<cd>();
  Ps = (Ut * Vs_.real()).cast<cd>() + cd(0.0, 1.0) * (Ut * Vs_.imag()).cast<cd>();
}

void PlexcitonSystem::spectra(const std::vector<Realization>& reals, const LineshapeEngine* engine,
                              double fixed_hwhm_cm, const std::vector<double>& w, PlexcitonSpectra& mean,
                              PlexcitonSpectra& se) const {
  const std::size_t nw = w.size(), R = reals.size();
  if (R == 0) throw std::invalid_argument("plexciton: no realizations");
  for (auto* s : {&mean, &se}) {
    s->w = w;
    for (auto* v : {&s->ext, &s->abs, &s->sca, &s->np_ext, &s->np_abs, &s->np_sca}) v->assign(nw, 0.0);
  }
  const auto N = static_cast<Eigen::Index>(agg_.size());
  const double c_mol = units::dd_cm_per_D2_nm3 * mu_D_ * mu_D_ / eps_m_;
  const bool sphere = particle_->shape() == plasmonics::ParticleShape::Sphere;
  const bool cumulant = engine != nullptr && fixed_hwhm_cm < 0.0;
  constexpr double pi = units::pi;

  // Суммы по реализациям (для среднего и SE) по каждой частоте
  std::vector<double> s1e(nw, 0), s2e(nw, 0), s1a(nw, 0), s2a(nw, 0), s1s(nw, 0), s2s(nw, 0);
  std::vector<double> np_ext(nw), np_sca(nw), np_ext_d(nw), np_sca_d(nw);
  std::vector<char> np_done(nw, 0);

  // Диагональ экситонного блока a_n(ω) = 1/резольвента: −i/X_n(ω) (кумулянт) или ε_n − ω − iγ_n.
  // Хранится для группы реализаций, чтобы ограничить память (≈ 2·10⁷ комплексных чисел).
  const std::size_t per = static_cast<std::size_t>(N) * nw;
  const std::size_t chunk = std::max<std::size_t>(1, std::min<std::size_t>(R, 20000000 / std::max<std::size_t>(per, 1)));
  for (std::size_t r0 = 0; r0 < R; r0 += chunk) {
    const std::size_t r1 = std::min(R, r0 + chunk);
    std::vector<std::vector<std::complex<double>>> diag(r1 - r0, std::vector<std::complex<double>>(per));
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (std::ptrdiff_t rr = 0; rr < static_cast<std::ptrdiff_t>(r1 - r0); ++rr)
      for (Eigen::Index n = 0; n < N; ++n) {
        const Realization& Rz = reals[r0 + static_cast<std::size_t>(rr)];
        const double eps = Rz.states.energies(n);
        auto& dg = diag[static_cast<std::size_t>(rr)];
        if (cumulant) {
          const double G = -Rz.K_ps(n, n) / units::rad_ps_per_cm;  // скорость ухода, угл. см⁻¹
          const auto X = engine->complex_response(w, eps, Rz.states.ipr(n), G);
          for (std::size_t q = 0; q < nw; ++q) dg[static_cast<std::size_t>(n) * nw + q] = cd(0.0, -1.0) / X[q];
        } else {
          const double hw = fixed_hwhm_cm >= 0.0 ? fixed_hwhm_cm : Rz.hwhm_cm(n);
          for (std::size_t q = 0; q < nw; ++q) dg[static_cast<std::size_t>(n) * nw + q] = cd(eps - w[q], -hw);
        }
      }

    // Быстрый решатель: проекции членов l ≤ fast_l (низкий ранг) и квадраты проекций старших членов (диагональ)
    const bool fast = spec_.fast_l > 0 && spec_.fast_l < particle_->l_max();
    const Eigen::Index T = Vc_.cols();
    const Eigen::Index Tc = fast ? Eigen::Index(spec_.fast_l) * (spec_.fast_l + 3) / 2 : 0;
    std::vector<Eigen::MatrixXcd> Blow(fast ? r1 - r0 : 0), Shigh(fast ? r1 - r0 : 0);
    if (fast) {
#pragma omp parallel for schedule(dynamic)
      for (std::ptrdiff_t rr = 0; rr < static_cast<std::ptrdiff_t>(r1 - r0); ++rr) {
        Eigen::MatrixXcd Pc, Ps;
        project(reals[r0 + static_cast<std::size_t>(rr)].states.c, Pc, Ps);
        auto& B = Blow[static_cast<std::size_t>(rr)];
        B.resize(N, 2 * Tc);
        B << Pc.leftCols(Tc), Ps.leftCols(Tc);
        Shigh[static_cast<std::size_t>(rr)] =
            Pc.rightCols(T - Tc).array().square().matrix() + Ps.rightCols(T - Tc).array().square().matrix();
      }
    }

#pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t q = 0; q < static_cast<std::ptrdiff_t>(nw); ++q) {
      const double nu = w[q];
      const Freq F = at(nu, -1, !fast);
      Eigen::MatrixXcd M(N, 3);
      for (Eigen::Index i = 0; i < N; ++i) M.row(i) = (F.L[i].transpose() * agg_.dipoles[i].cast<cd>()).transpose();
      auto cross = [&](const Eigen::Matrix3cd& a, double& ext, double& sca) {
        ext = 4.0 * pi * F.k * a.trace().imag() / 3.0;
        sca = 8.0 * pi / 3.0 * std::pow(F.k, 4) * a.cwiseAbs2().sum() / 3.0;
      };
      const Eigen::Matrix3cd a_np = F.alpha_np.asDiagonal();
      if (!np_done[q]) {
        cross(a_np, np_ext_d[q], np_sca_d[q]);
        np_ext[q] = np_ext_d[q];
        np_sca[q] = np_sca_d[q];
        if (sphere && spec_.d_perp_nm == cd(0.0)) {  // НЧ отдельно — полная теория Ми (все мультиполи)
          const auto mr = spec_.shell_thickness_nm > 0.0
                              ? plasmonics::mie_coated_sphere(spec_.a_xy_nm, spec_.outer_xy(), mat_.epsilon(nu),
                                                              cd(spec_.shell_n * spec_.shell_n), spec_.n_medium, nu)
                              : plasmonics::mie_sphere(spec_.a_xy_nm, mat_.epsilon(nu), spec_.n_medium, nu);
          np_ext[q] = mr.q_ext;
          np_sca[q] = mr.q_sca;
        }
        np_done[q] = 1;
      }
      for (std::size_t rr = 0; rr < r1 - r0; ++rr) {
        const Realization& Rz = reals[r0 + rr];
        const Eigen::MatrixXd& U = Rz.states.c;
        const Eigen::MatrixXcd Mx = U.transpose().cast<cd>() * M;
        Eigen::MatrixXcd X;
        if (!fast) {
          // В базисе экситонов: A' = diag(a_n) + Uᵀ W U, M' = Uᵀ M
          Eigen::MatrixXcd A = U.transpose().cast<cd>() * F.W * U.cast<cd>();
          for (Eigen::Index n = 0; n < N; ++n) A(n, n) += diag[rr][static_cast<std::size_t>(n) * nw + static_cast<std::size_t>(q)];
          X = A.partialPivLu().solve(Mx);
        } else {
          // A' = A₀ + pref·B·diag(δ)·Bᵀ, A₀ = diag(a_n + pref·Σ_{t>Tc} d_t S_nt); Вудбери:
          // X = Y − Z (I + pref·δ∘BᵀZ)⁻¹ (pref·δ∘BᵀY), Y = A₀⁻¹M', Z = A₀⁻¹B
          const Eigen::MatrixXcd& B = Blow[rr];
          Eigen::VectorXcd dl(2 * Tc);
          dl << F.d.head(Tc), F.d.head(Tc);
          dl *= c_mol;
          const Eigen::VectorXcd a0 = Shigh[rr] * (c_mol * F.d.tail(T - Tc));
          Eigen::VectorXcd inv(N);
          for (Eigen::Index n = 0; n < N; ++n)
            inv(n) = 1.0 / (a0(n) + diag[rr][static_cast<std::size_t>(n) * nw + static_cast<std::size_t>(q)]);
          const Eigen::MatrixXcd Y = inv.asDiagonal() * Mx, Z = inv.asDiagonal() * B;
          Eigen::MatrixXcd S = dl.asDiagonal() * (B.transpose() * Z);
          S.diagonal().array() += 1.0;
          X = Y - Z * S.partialPivLu().solve(dl.asDiagonal() * (B.transpose() * Y));
        }
        const Eigen::Matrix3cd a_tot = a_np + c_mol * (Mx.transpose() * X);
        double e, s;
        cross(a_tot, e, s);
        e += np_ext[q] - np_ext_d[q];  // поправка на высшие мультиполи НЧ (сфера)
        s += np_sca[q] - np_sca_d[q];
        const double a = e - s;
        s1e[q] += e;
        s2e[q] += e * e;
        s1s[q] += s;
        s2s[q] += s * s;
        s1a[q] += a;
        s2a[q] += a * a;
      }
    }
  }
  auto fin = [&](double s1, double s2, double& m, double& err) {
    m = s1 / R;
    err = R > 1 ? std::sqrt(std::max(0.0, (s2 / R - m * m)) * R / (R - 1.0) / R) : 0.0;
  };
  for (std::size_t q = 0; q < nw; ++q) {
    fin(s1e[q], s2e[q], mean.ext[q], se.ext[q]);
    fin(s1a[q], s2a[q], mean.abs[q], se.abs[q]);
    fin(s1s[q], s2s[q], mean.sca[q], se.sca[q]);
    mean.np_ext[q] = np_ext[q];
    mean.np_sca[q] = np_sca[q];
    mean.np_abs[q] = np_ext[q] - np_sca[q];
  }
}

std::vector<PlexcitonSystem::MoleculeNearField> PlexcitonSystem::near_field(double nu_cm, double gamma0_ps,
                                                                            double gamma_nr_ps) const {
  const Freq F = at(nu_cm, -1, false);
  const double pref = units::dd_cm_per_D2_nm3 * mu_D_ * mu_D_ / eps_m_;
  const Eigen::MatrixXcd S2 = Vc_.array().square().matrix() + Vs_.array().square().matrix();
  const Eigen::VectorXcd Wd = pref * (S2 * F.d);          // W_ii
  const Eigen::VectorXcd Wa = pref * (S2 * F.d_abs);      // W_ii без радиационного затухания НЧ
  std::vector<MoleculeNearField> out(agg_.size());
  for (std::size_t i = 0; i < agg_.size(); ++i) {
    auto& o = out[i];
    const auto ii = static_cast<Eigen::Index>(i);
    o.field_enhancement = F.L[i].cwiseAbs2().sum() / 3.0;
    o.shift_cm = Wd(ii).real();
    const Eigen::Vector3cd me = F.L[i].transpose() * agg_.dipoles[i].cast<cd>();
    o.rad_enhancement = me.squaredNorm();
    // Γ_tot = Γ₀|μ_eff/μ|² + Γ_abs: излучение молекулы вместе с наведённым диполем НЧ (с интерференцией) и
    // поглощение в металле (все мультиполи). −2Im W вместо этого теряло бы интерференционный член.
    // Дипольный (l = 1) вклад в поглощение — с запаздывающим полем молекулы в центре НЧ (|g(k)|² вместо ĝ²)
    double ret = 0.0;
    if (spec_.dynamic_dipole) {  // в чисто квазистатическом режиме — без запаздывания
      Eigen::Vector3d g0;
      Eigen::Vector3cd gk;
      dipole_field_at_origin(agg_.positions[i], agg_.dipoles[i], F.k, g0, gk);
      for (int a = 0; a < 3; ++a)
        ret += absorbing_part(F.alpha_np(a), F.k, spec_.dynamic_dipole).imag() * (std::norm(gk(a)) - g0(a) * g0(a));
    }
    o.nonrad_to_np_ps = std::max(0.0, (-2.0 * Wa(ii).imag() + 2.0 * pref * ret) * units::rad_ps_per_cm);
    o.decay_ps = o.nonrad_to_np_ps + (gamma0_ps > 0.0 ? gamma0_ps * (o.rad_enhancement - 1.0) : 0.0);
    o.gap_nm = particle_->surface_gap(agg_.positions[i]);
    if (gamma0_ps > 0.0) {
      const double grad = gamma0_ps * o.rad_enhancement;
      o.quantum_yield = grad / (grad + o.nonrad_to_np_ps + gamma_nr_ps);
      o.quantum_yield_free = gamma0_ps / (gamma0_ps + gamma_nr_ps);
    }
  }
  return out;
}

PlexcitonSystem::HybridEmission PlexcitonSystem::emission(const std::vector<Realization>& reals,
                                                          const LineshapeEngine& engine, const std::vector<double>& w,
                                                          double nonradiative_ps, bool laser_pump,
                                                          double excitation_cm, int nodes) const {
  HybridEmission H;
  H.w = w;
  H.pl.assign(w.size(), 0.0);
  H.pl_free.assign(w.size(), 0.0);
  if (reals.empty()) return H;
  const auto N = static_cast<Eigen::Index>(agg_.size());
  // Узлы по энергии, охватывающие все экситонные состояния всех реализаций
  double emin = 1e300, emax = -1e300;
  for (const auto& R : reals) {
    emin = std::min(emin, R.states.energies.minCoeff());
    emax = std::max(emax, R.states.energies.maxCoeff());
  }
  const int K = std::max(2, nodes);
  std::vector<double> wn(K);
  for (int k = 0; k < K; ++k) wn[k] = emin + (emax - emin) * k / (K - 1);
  std::vector<Freq> Fn(K);
#pragma omp parallel for schedule(dynamic)
  for (int k = 0; k < K; ++k) Fn[k] = at(wn[k], -1, false);
  Freq Fx;
  if (laser_pump) Fx = at(excitation_cm, -1, false);

  double qy = 0, qyf = 0, rb = 0, lt = 0, ltf = 0;
  for (const auto& R : reals) {
    const Eigen::MatrixXd& C = R.states.c;
    const Eigen::Index n = R.states.energies.size();
    // На узлах: W̃_nn и |μ_eff,n|²/|μ_n|²
    Eigen::MatrixXd Wim(n, K), Wre(n, K), Renh(n, K);
    // Диагональ W̃_nn = pref·Σ_t d_t (Ṽc² + Ṽs²) — точно, без построения N×N матрицы W
    Eigen::MatrixXcd Pc, Ps;
    project(C, Pc, Ps);
    const Eigen::MatrixXcd S2 = Pc.array().square().matrix() + Ps.array().square().matrix();
    const double pref = units::dd_cm_per_D2_nm3 * mu_D_ * mu_D_ / eps_m_;
    // Статические и запаздывающие поля сайтов в центре НЧ (для дипольной поправки на запаздывание)
    Eigen::MatrixXd G0(N, 3);
    for (Eigen::Index i = 0; i < N; ++i) {
      Eigen::Vector3d g0;
      Eigen::Vector3cd gk;
      dipole_field_at_origin(agg_.positions[i], agg_.dipoles[i], 0.0, g0, gk);
      G0.row(i) = g0.transpose();
    }
    const Eigen::MatrixXd P0 = C.transpose() * G0;  // n × 3
    for (int k = 0; k < K; ++k) {
      Eigen::MatrixXcd Gk(N, 3);
      for (Eigen::Index i = 0; i < N; ++i) {
        Eigen::Vector3d g0;
        Eigen::Vector3cd gk;
        dipole_field_at_origin(agg_.positions[i], agg_.dipoles[i], Fn[k].k, g0, gk);
        Gk.row(i) = gk.transpose();
      }
      const Eigen::MatrixXcd Pk = C.transpose().cast<cd>() * Gk;
      Eigen::VectorXcd Wd = pref * (S2 * Fn[k].d_abs);  // только поглощение в металле
      for (int a = 0; a < 3 && spec_.dynamic_dipole; ++a) {
        const double ia = absorbing_part(Fn[k].alpha_np(a), Fn[k].k, spec_.dynamic_dipole).imag();
        // −2Im W_nn += 2·pref·Im α_abs·(|P_k|² − P_0²) ⇔ Im W_nn −= pref·Im α_abs·(...)
        Wd -= cd(0.0, 1.0) * (pref * ia * (Pk.col(a).cwiseAbs2() - P0.col(a).cwiseAbs2())).cast<cd>();
      }
      Eigen::MatrixXcd M(N, 3);
      for (Eigen::Index i = 0; i < N; ++i) M.row(i) = (Fn[k].L[i].transpose() * agg_.dipoles[i].cast<cd>()).transpose();
      const Eigen::MatrixXcd Meff = C.transpose().cast<cd>() * M;
      for (Eigen::Index q = 0; q < n; ++q) {
        const cd wnn = Wd(q);
        Wim(q, k) = wnn.imag();
        Wre(q, k) = wnn.real();
        const double mu2 = R.states.mu2(q) / (mu_D_ * mu_D_);
        Renh(q, k) = mu2 > 1e-12 ? Meff.row(q).squaredNorm() / mu2 : 1.0;
      }
    }
    auto interp = [&](const Eigen::MatrixXd& T, Eigen::Index q, double e) {
      const double x = std::clamp((e - wn[0]) / (wn[K - 1] - wn[0] + 1e-300) * (K - 1), 0.0, double(K - 1));
      const int k0 = std::min(K - 2, int(x));
      const double f = x - k0;
      return (1 - f) * T(q, k0) + f * T(q, k0 + 1);
    };
    Eigen::VectorXd grad(n), gnr(n), pump(n), pumpf(n);
    const Eigen::VectorXd Gout = -R.K_ps.diagonal() / units::rad_ps_per_cm;
    Eigen::MatrixXcd Mx;
    if (laser_pump) {
      Eigen::MatrixXcd M(N, 3);
      for (Eigen::Index i = 0; i < N; ++i) M.row(i) = (Fx.L[i].transpose() * agg_.dipoles[i].cast<cd>()).transpose();
      Mx = C.transpose().cast<cd>() * M;
    }
    for (Eigen::Index q = 0; q < n; ++q) {
      const double e = R.states.energies(q), g0 = R.gamma_rad_ps(q);
      grad(q) = g0 * interp(Renh, q, e);
      gnr(q) = std::max(0.0, -2.0 * interp(Wim, q, e) * units::rad_ps_per_cm);
      double L = 1.0;
      if (laser_pump) {
        std::vector<double> v = {0.0};
        engine.add_line(v, {excitation_cm}, e, 1.0, R.states.ipr(q), Gout(q), false, 0);
        L = std::max(v[0], 0.0);
        pump(q) = Mx.row(q).squaredNorm() * mu_D_ * mu_D_ * L;
        pumpf(q) = R.states.mu2(q) * L;
      } else {
        pump(q) = interp(Renh, q, e) * R.states.mu2(q);
        pumpf(q) = R.states.mu2(q);
      }
    }
    // Нормировка накачки: одинаковое число поглощённых фотонов для гибрида и свободного агрегата
    pump /= std::max(pump.sum(), 1e-300);
    pumpf /= std::max(pumpf.sum(), 1e-300);
    const Eigen::VectorXd dec = grad + gnr + Eigen::VectorXd::Constant(n, nonradiative_ps);
    const Eigen::VectorXd decf = R.gamma_rad_ps + Eigen::VectorXd::Constant(n, nonradiative_ps);
    const Eigen::VectorXd P = steady_state(R.K_ps, dec, pump), Pf = steady_state(R.K_ps, decf, pumpf);
    qy += grad.dot(P);
    qyf += R.gamma_rad_ps.dot(Pf);
    lt += P.sum();
    ltf += Pf.sum();
    Eigen::Index ib;
    R.states.mu2.maxCoeff(&ib);
    rb += grad(ib) / std::max(R.gamma_rad_ps(ib), 1e-300);
    for (Eigen::Index q = 0; q < n; ++q) {
      const double e = R.states.energies(q);
      engine.add_line(H.pl, w, e, P(q) * grad(q) / (e * e * e), R.states.ipr(q), Gout(q), true, 3);
      engine.add_line(H.pl_free, w, e, Pf(q) * R.gamma_rad_ps(q) / (e * e * e), R.states.ipr(q), Gout(q), true, 3);
    }
  }
  const double Rn = static_cast<double>(reals.size());
  for (auto& v : H.pl) v /= Rn;
  for (auto& v : H.pl_free) v /= Rn;
  H.qy = qy / Rn;
  H.qy_free = qyf / Rn;
  H.rad_enhancement_bright = rb / Rn;
  H.lifetime_ps = lt / Rn;
  H.lifetime_free_ps = ltf / Rn;
  return H;
}

PlexcitonSystem::DipoleModeCoupling PlexcitonSystem::dipole_mode_coupling(const std::vector<double>& w) const {
  DipoleModeCoupling out;
  if (w.size() < 8) return out;
  // α по осям на сетке (динамическая дипольная поляризуемость, нм³)
  std::vector<Eigen::Vector3cd> alpha(w.size());
  for (std::size_t q = 0; q < w.size(); ++q) alpha[q] = at(w[q], -1, false).alpha_np;
  const auto& terms = particle_->terms();
  const auto& ut = particle_->uniform_term();
  const auto& us = particle_->uniform_is_sine();
  const double pref = units::dd_cm_per_D2_nm3 * mu_D_ * mu_D_ / eps_m_;
  double g2sum = 0.0;
  int best_axis = 0;
  for (int a = 0; a < 3; ++a) {
    // Подгонка лоренцевым полюсом вокруг максимума Im α
    std::size_t k = 0;
    for (std::size_t q = 1; q < w.size(); ++q)
      if (alpha[q](a).imag() > alpha[k](a).imag()) k = q;
    const double imax = alpha[k](a).imag();
    std::size_t lo = k, hi = k;
    while (lo > 0 && alpha[lo](a).imag() > 0.3 * imax) --lo;
    while (hi + 1 < w.size() && alpha[hi](a).imag() > 0.3 * imax) ++hi;
    double best = std::numeric_limits<double>::infinity(), bw = w[k], bg = 0, bA = 0;
    const double span = std::max(w[hi] - w[lo], 4.0 * (w[1] - w[0]));
    for (int i = 0; i <= 120; ++i) {
      const double w0 = w[lo] + (w[hi] - w[lo]) * i / 120.0;
      for (int j = 1; j <= 120; ++j) {
        const double gm = 3.0 * span * j / 120.0;
        // линейный МНК по A (вещ.) и B (компл.): α_q ≈ A·f_q + B
        Eigen::MatrixXd Mx(2 * (hi - lo + 1), 3);
        Eigen::VectorXd y(2 * (hi - lo + 1));
        for (std::size_t q = lo; q <= hi; ++q) {
          const cd f = 1.0 / cd(w0 * w0 - w[q] * w[q], -gm * w[q]);
          const Eigen::Index r = 2 * static_cast<Eigen::Index>(q - lo);
          Mx.row(r) << f.real(), 1.0, 0.0;
          Mx.row(r + 1) << f.imag(), 0.0, 1.0;
          y(r) = alpha[q](a).real();
          y(r + 1) = alpha[q](a).imag();
        }
        const Eigen::Vector3d c = Mx.colPivHouseholderQr().solve(y);
        const double res = (Mx * c - y).squaredNorm();
        if (res < best) {
          best = res;
          bw = w0;
          bg = gm;
          bA = c(0);
        }
      }
    }
    double norm = 0.0;
    for (std::size_t q = lo; q <= hi; ++q) norm += std::norm(alpha[q](a));
    // K_a: W^(l=1,a)_ii = α_a·K_a(i,i); отношение D̂/α одинаково для квазистатики и динамики
    const cd eps = mat_.epsilon(bw);
    const int t = ut[a];
    cd ratio;
    if (layered_) {
      const double b = layers_.outer_radius();
      ratio = terms[t].w * (-1.0 / (b * b * b));  // α = −D̂_1 b³ ⇒ w·D̂/α = −w/b³
    } else {
      const Eigen::Vector3cd aqs = particle_->polarizability(eps, eps_m_);
      ratio = terms[t].w * particle_->response(terms[t], eps, eps_m_) / aqs(a);
    }
    double tr = 0.0;
    for (Eigen::Index i = 0; i < Vc_.rows(); ++i) {
      const cd v = us[a] ? Vs_(i, t) : Vc_(i, t);
      tr += (pref * ratio * v * v).real();
    }
    const double g2 = std::abs(tr) * bA / (2.0 * bw);
    out.g_cm[a] = std::sqrt(std::max(0.0, g2));
    out.omega_axis[a] = bw;
    out.gamma_axis[a] = bg;
    g2sum += std::max(0.0, g2);
    if (a == 0 || out.g_cm[a] > out.g_cm[best_axis]) {
      best_axis = a;
      out.omega_p_cm = bw;
      out.gamma_p_cm = bg;
      out.residue = bA;
      out.fit_rel_rms = std::sqrt(best / std::max(norm, 1e-300));
    }
  }
  out.g_rms_cm = std::sqrt(g2sum / 3.0);
  return out;
}

double PlexcitonSystem::truncation_estimate(double nu_cm) const {
  const int L = particle_->l_max();
  const Freq full = at(nu_cm);
  const Freq part = at(nu_cm, std::max(1, (3 * L) / 4));
  const double num = (full.W - part.W).cwiseAbs().maxCoeff();
  const double den = full.W.cwiseAbs().maxCoeff();
  return den > 0.0 ? num / den : 0.0;
}

std::vector<std::pair<double, double>> normal_quadrature(int n) {
  if (n <= 1) return {{0.0, 1.0}};
  Eigen::MatrixXd Jm = Eigen::MatrixXd::Zero(n, n);
  for (int k = 1; k < n; ++k) Jm(k - 1, k) = Jm(k, k - 1) = std::sqrt(double(k));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Jm);
  std::vector<std::pair<double, double>> out;
  for (int k = 0; k < n; ++k) out.emplace_back(es.eigenvalues()(k), es.eigenvectors()(0, k) * es.eigenvectors()(0, k));
  return out;
}

Aggregate follow_core_scale(const Aggregate& agg, const NanoparticleSpec& spec, double s) {
  Aggregate a = agg;
  for (auto& r : a.positions) {
    const double n = r.norm();
    if (n < 1e-12) continue;
    const Vec3 u = r / n;
    const double rho = 1.0 / std::sqrt((u.x() * u.x() + u.y() * u.y()) / (spec.a_xy_nm * spec.a_xy_nm) +
                                       u.z() * u.z() / (spec.a_z_nm * spec.a_z_nm));
    r += (s - 1.0) * rho * u;
  }
  return a;
}

}  // namespace ssp::excitons
