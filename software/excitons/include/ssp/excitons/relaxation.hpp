#pragma once
// Секулярный Редфилд для населённостей экситонных состояний и кинетика с распадом.

#include "ssp/excitons/bath.hpp"
#include "ssp/excitons/hamiltonian.hpp"

#include <Eigen/Dense>

namespace ssp::excitons {

// Матрица скоростей K (угл. см⁻¹): K(m,n) = k_{n→m} = Σ_i c_in² c_im² C̃(ε_n − ε_m), K(n,n) = −Σ_m k_{n→m}.
Eigen::MatrixXd redfield_rate_matrix(const ExcitonStates& st, const Bath& bath);

// Полуширина (HWHM, см⁻¹) оптического перехода g→n по Редфилду: ½Σ_m k_{n→m} + ½ s_n C̃(0).
Eigen::VectorXd redfield_homogeneous_hwhm(const ExcitonStates& st, const Bath& bath, const Eigen::MatrixXd& K);

struct DecayParams {
  // Излучательный распад: если monomer_lifetime_ns > 0 — по нему (Γ_n = |μ_n|²/μ²·(ε_n/E_mono)³/τ),
  // иначе по формуле Γ = n_med ω³|μ|²/(3πε₀ħc³).
  double monomer_lifetime_ns = 0.0;
  double n_medium = 1.0;
  double monomer_energy_cm = 17000.0;
  double mu_D = 10.0;
  double nonradiative_ps = 0.0;  // Γ_nr, пс⁻¹ (одинаковый для всех состояний)
};

// Излучательные скорости (пс⁻¹) экситонных состояний.
Eigen::VectorXd radiative_rates_ps(const ExcitonStates& st, const DecayParams& p);

// Стационарные населённости: (K_ps − diag(Γ_rad + Γ_nr)) P = −pump. K_ps — в пс⁻¹.
Eigen::VectorXd steady_state(const Eigen::MatrixXd& K_ps, const Eigen::VectorXd& decay_ps, const Eigen::VectorXd& pump);

// Решение P(t) = exp[(K − Γ)t] P0 на сетке времён (пс) через симметризацию детальным балансом.
// Возвращает матрицу n_states × n_times.
Eigen::MatrixXd population_dynamics(const Eigen::MatrixXd& K_ps, const Eigen::VectorXd& decay_ps,
                                    const Eigen::VectorXd& energies_cm, double beta, const Eigen::VectorXd& P0,
                                    const Eigen::VectorXd& times_ps);

}  // namespace ssp::excitons
