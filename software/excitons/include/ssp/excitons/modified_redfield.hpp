#pragma once
// Модифицированная теория Редфилда (Zhang, Meier, Chernyak, Mukamel 1998; Yang & Fleming 2002):
// диагональная часть экситон-фононного взаимодействия учитывается непертурбативно (через g(t)), недиагональная —
// во втором порядке. Скорость перехода M → N (некоррелированные одинаковые бани на сайтах, Mukamel-соглашение g):
//   k_{M→N} = 2 Re ∫₀^∞ dt exp[iΩt − Λ g(t)] · { p g̈(t) − [(q₁ − q₂) ġ(t) + 2iλ q₁]² },
//   Ω = ε_M − ε_N − 2λ(a − p),  Λ = a + b − 2p,
//   a = Σc_iM⁴, b = Σc_iN⁴, p = Σc_iM²c_iN², q₁ = Σc_iN c_iM³, q₂ = Σc_iN³ c_iM, λ — полная реорганизация.
// Это перекрытие флуоресценции донора M с поглощением акцептора N; при λ → 0 переходит в обычный Редфилд,
// детальный баланс — по свободным энергиям ε_n − λ s_n (проверено численно; знак λ-члена — для g с Im g → −λt).

#include "ssp/excitons/hamiltonian.hpp"
#include "ssp/excitons/spectra.hpp"

#include <Eigen/Dense>

namespace ssp::excitons {

// Матрица скоростей (угл. см⁻¹) в формате redfield_rate_matrix: K(N,M) = k_{M→N}, K(M,M) = −Σ_N k_{M→N}.
// Используются g(t) и шаг по времени LineshapeEngine; интегрирование до затухания |exp(−Λ g)| < 1e-7 или конца окна.
Eigen::MatrixXd modified_redfield_rate_matrix(const ExcitonStates& st, const LineshapeEngine& engine);

}  // namespace ssp::excitons
