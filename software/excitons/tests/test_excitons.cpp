// Верификация экситонного ядра: единицы, связь, геометрия, формы линий, Редфилд, кинетика.
#include "ssp/excitons/bath.hpp"
#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/geometry.hpp"
#include "ssp/excitons/hamiltonian.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/excitons/vibronic.hpp"
#include "ssp/excitons/huckel.hpp"

#include <fstream>
#include "ssp/excitons/spectra.hpp"
#include "ssp/units.hpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace ssp;
using namespace ssp::excitons;

TEST(Units, DipoleDipoleConstant) {
  // (1 Д)²/(4πε₀ нм³) ≈ 5.0341 см⁻¹
  EXPECT_NEAR(units::dd_cm_per_D2_nm3, 5.0341, 1e-3);
  EXPECT_NEAR(units::kB_cm, 0.6950348, 1e-6);
  EXPECT_NEAR(units::ps_per_tu, 5.30884, 1e-4);
}

TEST(Coupling, PointDipoleGeometryFactors) {
  CouplingParams p;
  p.mu_D = 10.0;
  const Vec3 x(1, 0, 0), z(0, 0, 1);
  // Параллельные «бок о бок» (H): +μ²/r³; «голова к хвосту» (J): −2μ²/r³
  EXPECT_NEAR(pair_coupling(Vec3::Zero(), z, Vec3(1, 0, 0), z, p), 503.41, 0.1);
  EXPECT_NEAR(pair_coupling(Vec3::Zero(), x, Vec3(1, 0, 0), x, p), -1006.8, 0.2);
}

TEST(Coupling, ExtendedDipoleTendsToPointDipole) {
  CouplingParams pd, ed;
  pd.mu_D = ed.mu_D = 8.0;
  ed.model = CouplingModel::ExtendedDipole;
  ed.dipole_length_nm = 0.5;
  const Vec3 a = Vec3(1, 2, 0.5).normalized(), b = Vec3(-0.3, 1, 2).normalized();
  const Vec3 r(30.0, 5.0, -10.0);
  const double jp = pair_coupling(Vec3::Zero(), a, r, b, pd), je = pair_coupling(Vec3::Zero(), a, r, b, ed);
  EXPECT_NEAR(je / jp, 1.0, 1e-3);
  // На малых расстояниях модели различаются заметно
  const Vec3 r1(0.9, 0.0, 0.0);
  EXPECT_GT(std::abs(pair_coupling(Vec3::Zero(), a, r1, b, ed) / pair_coupling(Vec3::Zero(), a, r1, b, pd) - 1.0), 0.02);
}

TEST(Geometry, AllTypesBuildWithDefaults) {
  for (const auto& t : geometry_types()) {
    const Aggregate a = build_geometry(t, {});
    EXPECT_GT(a.size(), 0u) << t;
    for (const auto& d : a.dipoles) EXPECT_NEAR(d.norm(), 1.0, 1e-12) << t;
  }
  EXPECT_THROW(build_geometry("ring", {{"radius_typo", 3.0}}), std::invalid_argument);
}

TEST(Geometry, LegacyNotebookCylinderMatches) {
  // Совпадение с build_cylinder ноутбука: молекула (m=1, n=2)
  const Aggregate a = build_geometry("cylinder", {{"n_rings", 8}, {"molecules_per_ring", 10}, {"radius", 3.8},
                                                  {"pitch", 1.0}, {"helical_twist_deg", 18.0},
                                                  {"dipole_tilt_deg", 22.0}});
  const double phi = 2 * units::pi * 2 / 10 + 18.0 * units::pi / 180.0;
  const Vec3 r = a.positions[10 + 2];
  EXPECT_NEAR(r.x(), 3.8 * std::cos(phi), 1e-12);
  EXPECT_NEAR(r.y(), 3.8 * std::sin(phi), 1e-12);
  EXPECT_NEAR(r.z(), 1.0, 1e-12);
}

TEST(Hamiltonian, RingIsCirculantWithAnalyticSpectrum) {
  // Кольцо с касательными диполями: циркулянт, E_k = ε + Σ_j J_0j cos(2πkj/N)
  const int N = 12;
  const Aggregate a = build_geometry("ring", {{"N", N}, {"radius", 3.0}, {"tilt_deg", 0.0}});
  CouplingParams p;
  const Eigen::MatrixXd J = coupling_matrix(a, p);
  std::mt19937_64 rng(1);
  const Eigen::MatrixXd H = frenkel_hamiltonian(J, 17000.0, 0.0, rng);
  const ExcitonStates st = diagonalize(H, a, p.mu_D);
  std::vector<double> ek(N);
  for (int k = 0; k < N; ++k) {
    double e = 17000.0;
    for (int j = 1; j < N; ++j) e += J(0, j) * std::cos(2.0 * units::pi * k * j / N);
    ek[k] = e;
  }
  std::sort(ek.begin(), ek.end());
  for (int k = 0; k < N; ++k) EXPECT_NEAR(st.energies(k), ek[k], 1e-9);
  // Сумма сил осцилляторов сохраняется: Σ|μ_n|² = N μ²
  EXPECT_NEAR(st.mu2.sum(), N * p.mu_D * p.mu_D, 1e-9);
}

TEST(Bath, LineshapeFunctionMatchesAnalyticDrude) {
  Bath b;
  b.lambda_cm = 60.0;
  b.cutoff_cm = 120.0;
  b.temperature_K = 300.0;
  const double dt = 2e-4;
  const std::size_t n = 20000;  // T = 4 t_u ≈ 21 пс
  const auto g = lineshape_function(b, dt, n);
  for (std::size_t j : {100u, 1000u, 5000u, 19999u}) {
    const auto ga = lineshape_function_drude(b.lambda_cm, b.cutoff_cm, b.beta(), dt * j);
    EXPECT_NEAR(std::abs(g[j] - ga) / std::abs(ga), 0.0, 2e-3) << "t = " << dt * j;
  }
}

TEST(Bath, DetailedBalance) {
  Bath b;
  b.temperature_K = 250.0;
  for (double w : {10.0, 150.0, 800.0})
    EXPECT_NEAR(b.rate_spectrum(-w) / b.rate_spectrum(w), std::exp(-b.beta() * w), 1e-12);
  // Реорганизация: (1/π)∫J/ω = λ
  double s = 0.0;
  const double h = 0.01;
  for (double w = h / 2; w < 2e5; w += h) s += b.J_over_w(w) * h;
  EXPECT_NEAR(s / units::pi, b.lambda_cm, 0.01 * b.lambda_cm);
}

TEST(Lineshape, MonomerAreaAndKmsDetailedBalance) {
  // Площадь линии = 1; флуоресценция и поглощение мономера связаны соотношением КМШ:
  //   F(ω)/A(ω) = exp[−β(ω − ω₀₀)],  ω₀₀ = ε − λ  ⇒  d ln(F/A)/dω = −β.
  Bath b;
  b.lambda_cm = 80.0;
  b.cutoff_cm = 100.0;
  const LineshapeEngine eng(b, 8000.0, 60.0);
  std::vector<double> w(6001);
  for (std::size_t i = 0; i < w.size(); ++i) w[i] = 13000.0 + i;  // шаг 1 см⁻¹
  std::vector<double> A(w.size(), 0.0), F(w.size(), 0.0);
  eng.add_line(A, w, 16000.0, 1.0, 1.0, 0.0, false, 0);
  eng.add_line(F, w, 16000.0, 1.0, 1.0, 0.0, true, 0);
  double sa = 0, sf = 0;
  for (std::size_t i = 0; i < w.size(); ++i) {
    sa += A[i];
    sf += F[i];
  }
  EXPECT_NEAR(sa, 1.0, 5e-3);
  EXPECT_NEAR(sf, 1.0, 5e-3);
  for (double x : {15800.0, 15920.0, 16050.0}) {
    const std::size_t i = static_cast<std::size_t>(x - 13000.0), j = i + 40;
    const double slope = (std::log(F[j] / A[j]) - std::log(F[i] / A[i])) / 40.0;
    EXPECT_NEAR(slope / (-b.beta()), 1.0, 1e-3) << "ω = " << x;
    const double w00 = w[i] + std::log(F[i] / A[i]) / b.beta();
    EXPECT_NEAR(w00, 16000.0 - b.lambda_cm, 0.5);
  }
}

TEST(Lineshape, KmsForOhmicBathWithModeViaFft) {
  // Численный путь (БПФ C(t) и двойное интегрирование): омическая баня + колебательная мода.
  Bath b;
  b.type = BathType::OhmicExponential;
  b.lambda_cm = 60.0;
  b.cutoff_cm = 150.0;
  b.modes.push_back({600.0, 0.15, 30.0});
  const LineshapeEngine eng(b, 6000.0, 20.0);
  std::vector<double> w(4001);
  for (std::size_t i = 0; i < w.size(); ++i) w[i] = 14000.0 + i;
  std::vector<double> A(w.size(), 0.0), F(w.size(), 0.0);
  eng.add_line(A, w, 16000.0, 1.0, 1.0, 0.0, false, 0);
  eng.add_line(F, w, 16000.0, 1.0, 1.0, 0.0, true, 0);
  double sa = 0.0;
  for (double v : A) sa += v;
  EXPECT_NEAR(sa, 1.0, 5e-3);
  const double l = b.reorganization();
  for (double x : {15700.0, 15900.0, 16050.0}) {
    const std::size_t i = static_cast<std::size_t>(x - 14000.0), j = i + 40;
    const double slope = (std::log(F[j] / A[j]) - std::log(F[i] / A[i])) / 40.0;
    EXPECT_NEAR(slope / (-b.beta()), 1.0, 2e-3) << "ω = " << x;
    EXPECT_NEAR(w[i] + std::log(F[i] / A[i]) / b.beta(), 16000.0 - l, 1.0);
  }
}

TEST(Redfield, StationaryPopulationsAreBoltzmann) {
  const Aggregate a = build_geometry("linear_chain", {{"N", 8}, {"spacing", 1.2}});
  CouplingParams p;
  const Eigen::MatrixXd J = coupling_matrix(a, p);
  std::mt19937_64 rng(5);
  const Eigen::MatrixXd H = frenkel_hamiltonian(J, 17000.0, 50.0, rng);
  const ExcitonStates st = diagonalize(H, a, p.mu_D);
  Bath b;
  b.temperature_K = 200.0;
  const Eigen::MatrixXd K = redfield_rate_matrix(st, b) * units::rad_ps_per_cm;
  // Очень медленный распад: населённости → Больцман
  const Eigen::VectorXd decay = Eigen::VectorXd::Constant(st.size(), 1e-9);
  const Eigen::VectorXd P = steady_state(K, decay, Eigen::VectorXd::Constant(st.size(), 1.0 / st.size()));
  Eigen::VectorXd B(st.size());
  for (Eigen::Index n = 0; n < B.size(); ++n) B(n) = std::exp(-b.beta() * (st.energies(n) - st.energies(0)));
  B /= B.sum();
  const Eigen::VectorXd Pn = P / P.sum();
  for (Eigen::Index n = 0; n < B.size(); ++n) EXPECT_NEAR(Pn(n), B(n), 1e-6);
  // Динамика сохраняет населённость при нулевом распаде и сходится к Больцману
  Eigen::VectorXd P0 = Eigen::VectorXd::Zero(st.size());
  P0(st.size() - 1) = 1.0;
  Eigen::VectorXd t(2);
  t << 0.0, 1e4;
  const Eigen::MatrixXd Pt = population_dynamics(K, Eigen::VectorXd::Zero(st.size()), st.energies, b.beta(), P0, t);
  EXPECT_NEAR(Pt.col(1).sum(), 1.0, 1e-9);
  for (Eigen::Index n = 0; n < B.size(); ++n) EXPECT_NEAR(Pt(n, 1), B(n), 1e-6);
}

TEST(Radiative, MonomerRateFromFormula) {
  // μ = 10 Д, 17000 см⁻¹, n = 1: Γ = ω³μ²/(3πε₀ħc³) ≈ 1.54·10⁸ с⁻¹
  ExcitonStates st;
  st.energies = Eigen::VectorXd::Constant(1, 17000.0);
  st.mu2 = Eigen::VectorXd::Constant(1, 100.0);
  DecayParams d;
  const double g = radiative_rates_ps(st, d)(0) * 1e12;
  EXPECT_NEAR(g / 1.54e8, 1.0, 0.01);
}

TEST(Ensemble, SuperradiantJAggregateNarrowsAndShifts) {
  // Линейная J-цепочка: полоса смещена в красную область, сужена (обменное сужение) относительно мономера;
  // без беспорядка стоксов сдвиг J-полосы меньше, чем у мономера (λ_n = λ/N_del).
  ExcitonModel m;
  m.monomer_energy_cm = 17000.0;
  m.disorder_sigma_cm = 150.0;
  m.coupling.mu_D = 10.0;
  m.bath.lambda_cm = 30.0;
  m.bath.cutoff_cm = 100.0;
  m.decay.nonradiative_ps = 1e-3;
  SpectrumGrid g{15000.0, 18500.0, 1400, 60.0};
  EnsembleSettings s;
  s.min_realizations = 32;
  s.max_realizations = 64;
  const Aggregate chain = build_geometry("linear_chain", {{"N", 30}, {"spacing", 1.0}, {"dipole_angle_deg", 0.0}});
  const Aggregate mono = build_geometry("linear_chain", {{"N", 1}});
  const EnsembleResult rc = run_ensemble(chain, m, g, s);
  const EnsembleResult rm = run_ensemble(mono, m, g, s);
  EXPECT_LT(rc.observables.at("abs_shift_cm").mean, -1000.0);
  EXPECT_LT(rc.observables.at("abs_fwhm_cm").mean, 0.6 * rm.observables.at("abs_fwhm_cm").mean);
  EXPECT_GT(rc.observables.at("superradiance_factor").mean, 5.0);
  ExcitonModel m0 = m;
  m0.disorder_sigma_cm = 0.0;
  EnsembleSettings s0 = s;
  s0.min_realizations = s0.max_realizations = s0.batch = 1;
  const EnsembleResult c0 = run_ensemble(chain, m0, g, s0), m00 = run_ensemble(mono, m0, g, s0);
  EXPECT_LT(c0.observables.at("stokes_shift_cm").mean, 0.5 * m00.observables.at("stokes_shift_cm").mean);
}

TEST(Coupling, TransitionChargesReproduceExtendedAndPointDipole) {
  Aggregate a = build_geometry("linear_chain", {{"N", 3}, {"spacing", 1.1}, {"dipole_angle_deg", 30.0}});
  CouplingParams ed, tq, pd;
  ed.model = CouplingModel::ExtendedDipole;
  ed.dipole_length_nm = 0.6;
  ed.mu_D = tq.mu_D = pd.mu_D = 9.0;
  tq.model = CouplingModel::TransitionCharges;
  // Два заряда ±1 на ±0.3 нм вдоль оси y файла: после выравнивания — тот же протяжённый диполь
  tq.charges = {Eigen::Vector4d(0.0, 0.3, 0.0, 1.0), Eigen::Vector4d(0.0, -0.3, 0.0, -1.0)};
  const Eigen::MatrixXd Je = coupling_matrix(a, ed), Jt = coupling_matrix(a, tq);
  EXPECT_NEAR((Je - Jt).cwiseAbs().maxCoeff(), 0.0, 1e-9 * Je.cwiseAbs().maxCoeff());
  // Несимметричное распределение (с квадрупольной частью) на больших расстояниях → точечный диполь
  tq.charges = {Eigen::Vector4d(0.4, 0.0, 0.0, 0.7), Eigen::Vector4d(-0.2, 0.1, 0.0, -0.4),
                Eigen::Vector4d(-0.3, -0.1, 0.05, -0.3)};
  Aggregate far = build_geometry("linear_chain", {{"N", 2}, {"spacing", 40.0}, {"dipole_angle_deg", 50.0}});
  EXPECT_NEAR(coupling_matrix(far, tq)(0, 1) / coupling_matrix(far, pd)(0, 1), 1.0, 2e-3);
}

TEST(Holstein, MonomerGivesPoissonProgression) {
  Eigen::MatrixXd H = Eigen::MatrixXd::Constant(1, 1, 18000.0);
  Eigen::MatrixX3d mu = Eigen::MatrixX3d::Zero(1, 3);
  mu(0, 0) = 10.0;
  HolsteinParams p;
  p.omega_cm = 1300.0;
  p.huang_rhys = 0.8;
  p.nu_max = 8;
  const VibronicStates vs = holstein_opa(H, mu, p);
  double first = 0.0;
  for (int v = 0; v <= p.nu_max; ++v) {
    EXPECT_NEAR(vs.energies(v), 18000.0 - 0.8 * 1300.0 + v * 1300.0, 1e-8);
    EXPECT_NEAR(vs.mu2(v) / 100.0, std::exp(-0.8) * std::pow(0.8, v) / std::tgamma(v + 1.0), 1e-12);
    first += vs.mu2(v) * vs.energies(v) / 100.0;
  }
  EXPECT_NEAR(first, 18000.0, 1.0);  // первый момент = вертикальная энергия
}

TEST(Holstein, JAggregateEnhancesZeroZeroBand) {
  // Цепочка с отрицательной связью (J): отношение 0-0/0-1 больше, чем у мономера (e^{−S}/S·e^{−S} → 1/S); H — меньше.
  HolsteinParams p;
  p.omega_cm = 1300.0;
  p.huang_rhys = 0.6;
  p.nu_max = 4;
  auto ratio = [&](double Jnn) {
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(12, 12);
    H.diagonal().setConstant(18000.0);
    for (int i = 0; i + 1 < 12; ++i) H(i, i + 1) = H(i + 1, i) = Jnn;
    Eigen::MatrixX3d mu = Eigen::MatrixX3d::Zero(12, 3);
    mu.col(0).setConstant(10.0);
    const VibronicStates vs = holstein_opa(H, mu, p);
    double a00 = 0, a01 = 0, tot = 0;
    for (Eigen::Index k = 0; k < vs.energies.size(); ++k) {
      tot += vs.mu2(k);
      if (vs.energies(k) < 18000.0 - 0.6 * 1300.0 + 650.0) a00 += vs.mu2(k);
      else if (vs.energies(k) < 18000.0 - 0.6 * 1300.0 + 1950.0) a01 += vs.mu2(k);
    }
    EXPECT_NEAR(tot / (12 * 100.0), 1.0, 1e-3);  // правило сумм (хвост ν > 4 ~ 4·10⁻⁴)
    return a00 / a01;
  };
  const double mono = 1.0 / 0.6;
  EXPECT_GT(ratio(-300.0), 1.3 * mono);
  EXPECT_LT(ratio(+300.0), 0.7 * mono);
}

TEST(Dispersion, RandomOrientationsReproduceContinuumAndTubeWallsDiffer) {
  // Калибровка Q: при случайных ориентациях средний сдвиг = континуальный −F(n²−1)/(2n²+1) (другая выборка
  // ориентаций — совпадение в пределах статистики); в трубке внутренняя стенка смещена в красную сторону сильнее.
  const Aggregate a = build_geometry("brickwork_tube", {{"length", 15.0}});
  DispersionParams p;
  p.enabled = true;
  p.n_solvent = 1.33;
  const auto D = dispersion_site_shifts(a, p);
  EXPECT_NEAR(D.continuum_cm, -1.43 * (1.78 * 1.78 - 1) / (2 * 1.78 * 1.78 + 1) * 8065.544, 1e-6);
  Aggregate r = a;
  std::mt19937_64 rng(99);
  std::normal_distribution<double> g(0.0, 1.0);
  for (auto& d : r.dipoles) d = Vec3(g(rng), g(rng), g(rng)).normalized();
  DispersionParams q = p;
  q.orientation_samples = 1;  // Q — по одной случайной выборке, применяем к другой
  const double Q = dispersion_site_shifts(a, p).Q_cm_nm2;
  const auto Dr = dispersion_site_shifts(r, q);
  EXPECT_NEAR(Q * (Dr.dispersive_cm.array() / Dr.Q_cm_nm2).mean() / D.continuum_cm, 1.0, 0.05);
  double in = 0, out = 0;
  int ni = 0, no = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double z = a.positions[i].z();
    if (z < 4.0 || z > 11.0) continue;  // без концов трубки
    const double rr = std::hypot(a.positions[i].x(), a.positions[i].y());
    (rr < 4.2 ? in : out) += D.dispersive_cm(static_cast<Eigen::Index>(i));
    ++(rr < 4.2 ? ni : no);
  }
  EXPECT_LT(in / ni, out / no);
  EXPECT_LT(in / ni, 0.0);
}

TEST(Geometry, SeamlessBrickworkTubeHasNoClosePairs) {
  const Aggregate a = build_geometry("brickwork_tube", {{"length", 15.0}});
  EXPECT_GT(min_pair_distance(a), 0.45);  // не ближе шага стопки (0.5 нм) с поправкой на кривизну
  int nin = 0;
  for (const auto& r : a.positions) nin += std::hypot(r.x(), r.y()) < 4.2;
  const double ratio = double(a.size() - nin) / nin;
  EXPECT_NEAR(ratio, 4.65 / 3.79, 0.1);  // одинаковая плотность на стенках: N_out/N_in ≈ R_out/R_in
  ASSERT_EQ(a.normals.size(), a.size());
  for (std::size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a.dipoles[i].dot(a.normals[i]), 0.0, 1e-12);
}

TEST(Dispersion, TwoPointChargesReproduceExtendedDipole) {
  const Aggregate a = build_geometry("brickwork_tube", {{"length", 10.0}});
  DispersionParams p;
  p.enabled = true;
  p.orientation_samples = 4;
  const auto D1 = dispersion_site_shifts(a, p);
  p.charges = {Eigen::Vector4d(0.5, 0, 0, 1.0), Eigen::Vector4d(-0.5, 0, 0, -1.0)};  // l = 1 нм
  const auto D2 = dispersion_site_shifts(a, p);
  // одинаковые S для упорядоченной структуры; Q — по другим (случайным поворотам) выборкам, поэтому сравниваем ΣV²
  EXPECT_LT(((D1.dispersive_cm / D1.Q_cm_nm2) - (D2.dispersive_cm / D2.Q_cm_nm2)).cwiseAbs().maxCoeff(), 1e-9);
  EXPECT_NEAR(D1.Q_cm_nm2 / D2.Q_cm_nm2, 1.0, 0.03);
}

TEST(Huckel, StreptocyanineTransitionChargesAlongChain) {
  // Пентаметинстрептоцианин H2N–(CH)5–NH2⁺ в плоскости xy (зигзаг): π-система из 7 атомов, 8 электронов;
  // переходный диполь HOMO→LUMO — вдоль цепи N…N, сумма зарядов 0, нормаль ⟂ плоскости.
  const double b = 1.39, h = 1.08;
  std::vector<Vec3> heavy;
  for (int k = 0; k < 7; ++k) heavy.emplace_back(k * b * std::cos(M_PI / 6), (k % 2) * b * std::sin(M_PI / 6), 0.0);
  const std::string path = ::testing::TempDir() + "strepto.xyz";
  {
    std::ofstream f(path);
    f << "16\nstreptocyanine\n";
    for (int k = 0; k < 7; ++k) {
      const Vec3& r = heavy[k];
      f << (k == 0 || k == 6 ? "N " : "C ") << r.x() << " " << r.y() << " " << r.z() << "\n";
      const double sy = (k % 2) ? 1.0 : -1.0;  // H — наружу от зигзага
      if (k == 0 || k == 6) {
        const double sx = k == 0 ? -1.0 : 1.0;
        f << "H " << r.x() + sx * h * std::cos(M_PI / 6) << " " << r.y() + sy * h * std::sin(M_PI / 6) << " 0\n";
        f << "H " << r.x() << " " << r.y() + sy * h << " 0\n";
      } else {
        f << "H " << r.x() << " " << r.y() + sy * h << " 0\n";
      }
    }
  }
  const auto atoms = read_xyz(path);
  ASSERT_EQ(atoms.size(), 16u);  // 7 тяжёлых + 5 H у C + 4 H у N
  const auto H = huckel_transition_charges(atoms);
  ASSERT_EQ(H.pi_atoms.size(), 7u);
  EXPECT_EQ(H.homo, 3);
  double qs = 0, ymax = 0, zmax = 0;
  for (const auto& c : H.charges) {
    qs += c(3);
    ymax = std::max(ymax, std::abs(c(1)));
    zmax = std::max(zmax, std::abs(c(2)));
  }
  EXPECT_NEAR(qs, 0.0, 1e-10);
  EXPECT_LT(zmax, 1e-9);   // плоская молекула: z — нормаль
  EXPECT_LT(ymax, 0.08);   // x — вдоль цепи (зигзаг ±0.035 нм)
  EXPECT_GT(H.dipole_e_nm, 0.05);
}
