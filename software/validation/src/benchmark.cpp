// ssp_benchmark — валидация моделей против независимых эталонов.
//   ssp_benchmark purcell <material.yml|drude> <out.csv>
//     Скорости распада диполя у Au-сферы: модель программы (квазистатические мультиполи + множители Ми,
//     PlexcitonSystem::near_field) против точной теории Ми (dipole_near_sphere). Радиусы, зазоры, ориентации, λ.
#include "ssp/excitons/plexciton.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/plasmonics/dipole_sphere.hpp"
#include "ssp/plasmonics/material.hpp"
#include "ssp/plasmonics/mie.hpp"

#include <array>
#include <cmath>
#include <map>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace ssp;
using namespace ssp::excitons;

namespace {

plasmonics::Material load_material(const std::string& m) {
  if (m == "drude") {
    plasmonics::DrudeParams dp;
    dp.eps_inf = 9.5;
    dp.omega_p_eV = 9.0;
    dp.gamma_eV = 0.07;
    return plasmonics::Material::drude(dp);
  }
  return plasmonics::Material::from_refractiveindex_yaml(m);
}

int purcell(const std::string& mat_path, const std::string& out_path) {
  const auto mat = load_material(mat_path);
  const double n_med = 1.33, mu = 10.0;
  // Γ₀ того же диполя в однородной среде (та же формула, что для экситонов)
  ExcitonStates m1;
  m1.energies = Eigen::VectorXd::Constant(1, 1.0);
  m1.mu2 = Eigen::VectorXd::Constant(1, mu * mu);
  std::ofstream f(out_path);
  f << "# dipole near a gold sphere in water (n = 1.33): model (quasistatic multipoles + Mie factors) vs exact Mie\n"
    << "# gamma ratios are normalized to the free radiative rate in the medium; qy assumes intrinsic QY = 1\n"
    << "radius_nm,gap_nm,lambda_nm,orientation,exact_total,model_total,exact_rad,model_rad,exact_qy,model_qy,"
       "err_total,err_rad\n";
  // сводка по радиусам: max |ошибка Γ_tot| (отн.), max |ошибка Γ_rad| — отн. при Γ_rad > 0.2Γ₀, иначе абс. в Γ₀
  std::map<double, std::array<double, 2>> worst;
  for (const double R : {5.0, 10.0, 20.0, 40.0}) {
    for (const double lam : {520.0, 600.0, 650.0}) {
      const double nu = 1e7 / lam;
      Aggregate a;
      std::vector<double> gaps = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};
      for (double g : gaps) {
        a.positions.emplace_back(0, 0, R + g);
        a.dipoles.emplace_back(0, 0, 1);  // радиальный
        a.positions.emplace_back(0, 0, -(R + g));
        a.dipoles.emplace_back(1, 0, 0);  // касательный
      }
      NanoparticleSpec spec;
      spec.a_xy_nm = spec.a_z_nm = R;
      spec.l_max = 120;
      spec.n_medium = n_med;
      const PlexcitonSystem sys(a, spec, mat, mu);
      DecayParams dp;
      dp.n_medium = n_med;
      dp.mu_D = mu;
      m1.energies(0) = nu;
      dp.monomer_energy_cm = nu;
      const double g0 = radiative_rates_ps(m1, dp)(0);
      const auto nf = sys.near_field(nu, g0, 0.0);
      for (std::size_t k = 0; k < gaps.size(); ++k) {
        const auto ex = plasmonics::dipole_near_sphere(R, mat.epsilon(nu), n_med, nu, R + gaps[k]);
        for (int o = 0; o < 2; ++o) {
          const auto& q = nf[2 * k + static_cast<std::size_t>(o)];
          const double mt = 1.0 + q.decay_ps / g0, mr = q.rad_enhancement;
          const double et = o == 0 ? ex.total_perp : ex.total_par, er = o == 0 ? ex.rad_perp : ex.rad_par;
          const double e1 = mt / et - 1.0, e2 = mr / er - 1.0;
          f << R << "," << gaps[k] << "," << lam << "," << (o == 0 ? "radial" : "tangential") << "," << et << "," << mt
            << "," << er << "," << mr << "," << er / et << "," << mr / mt << "," << e1 << "," << e2 << "\n";
          auto& wr = worst[R];
          wr[0] = std::max(wr[0], std::abs(e1));
          wr[1] = std::max(wr[1], er > 0.2 ? std::abs(e2) : std::abs(mr - er));
        }
      }
    }
  }
  std::printf("radius_nm  max|err Γ_tot|  max|err Γ_rad| (rel. if Γ_rad > 0.2Γ₀, else abs. in Γ₀)\n");
  for (const auto& [R, w] : worst) std::printf("%9.0f  %13.3f  %13.3f\n", R, w[0], w[1]);
  std::printf("table: %s\n", out_path.c_str());
  return 0;
}

// Положение пика и молярный коэффициент экстинкции Au-сфер в воде против паспортных данных (Sigma-Aldrich:
// λmax 5 нм 510–521, 10 нм 515–523, 50 нм 531–535 нм; ε(5 нм) = 1.10e7, ε(10 нм) = 1.01e8 M⁻¹см⁻¹).
int extinction(const std::vector<std::string>& files, const std::string& out_path) {
  struct Ref {
    double d, lo, hi, eps;
  };
  const std::vector<Ref> refs = {{5, 510, 521, 1.10e7}, {10, 515, 523, 1.01e8}, {20, 0, 0, 0}, {50, 531, 535, 0},
                                 {80, 0, 0, 0}, {100, 0, 0, 0}};
  std::ofstream f(out_path);
  f << "# gold spheres in water (n = 1.33), Mie; reference: Sigma-Aldrich datasheets (lambda_max ranges, molar "
       "extinction at lambda_max)\n"
    << "material,kreibig_A,diameter_nm,lambda_max_nm,eps_max_M-1cm-1,ref_lo_nm,ref_hi_nm,ref_eps,in_range,eps_ratio\n";
  constexpr double NA = 6.02214076e23;
  for (const auto& file : files) {
    for (const double A : {0.0, 1.0}) {
      std::printf("%s, Kreibig A = %.0f\n", file.substr(file.find_last_of('/') + 1).c_str(), A);
      for (const auto& r : refs) {
        auto mat = plasmonics::Material::from_refractiveindex_yaml(file);
        mat.set_size_damping(A, 1.40e6, 0.5 * r.d, 9.0, 0.07);
        double best = 0, lbest = 0;
        for (double lam = 450.0; lam <= 700.0; lam += 0.25) {
          const double nu = 1e7 / lam;
          const double s = plasmonics::mie_sphere(0.5 * r.d, mat.epsilon(nu), 1.33, nu).q_ext;  // нм²
          if (s > best) {
            best = s;
            lbest = lam;
          }
        }
        const double eps = best * 1e-14 * NA / (std::log(10.0) * 1000.0);  // σ (см²) → ε (M⁻¹см⁻¹)
        const bool inr = r.lo > 0 && lbest >= r.lo && lbest <= r.hi;
        f << file.substr(file.find_last_of('/') + 1) << "," << A << "," << r.d << "," << lbest << "," << eps << ","
          << r.lo << "," << r.hi << "," << r.eps << "," << (r.lo > 0 ? (inr ? "yes" : "NO") : "") << ","
          << (r.eps > 0 ? eps / r.eps : 0.0) << "\n";
        std::printf("  d = %5.0f nm: λmax %6.2f nm %s, ε %.3g%s\n", r.d, lbest,
                    r.lo > 0 ? (inr ? "(in ref. range)" : "(OUT of ref. range)") : "", eps,
                    r.eps > 0 ? (" (ref " + std::to_string(r.eps) + ", ratio " + std::to_string(eps / r.eps) + ")").c_str()
                              : "");
      }
    }
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::string cmd = argc > 1 ? argv[1] : "";
    if (cmd == "purcell" && argc == 4) return purcell(argv[2], argv[3]);
    if (cmd == "extinction" && argc >= 4) {
      std::vector<std::string> files(argv + 2, argv + argc - 1);
      return extinction(files, argv[argc - 1]);
    }
    std::cerr << "usage: ssp_benchmark purcell <material.yml|drude> <out.csv>\n"
                 "       ssp_benchmark extinction <material.yml>... <out.csv>\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
