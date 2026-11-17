#pragma once
// Физические константы (CODATA 2018) и пересчёт единиц проекта.
// Соглашение: энергии — см⁻¹, длины — нм, дипольные моменты — Дебай, скорости — пс⁻¹.
// Внутри динамических решателей ħ = 1: энергия в см⁻¹ трактуется как угловая частота
// в единицах 2πc·(1 см⁻¹); соответствующая единица времени t_u = 1/(2πc·1 см⁻¹).

#include <numbers>

namespace ssp::units {

inline constexpr double pi = std::numbers::pi;

// СИ
inline constexpr double h = 6.62607015e-34;          // Дж·с
inline constexpr double hbar = h / (2.0 * pi);       // Дж·с
inline constexpr double c = 2.99792458e8;            // м/с
inline constexpr double eps0 = 8.8541878128e-12;     // Ф/м
inline constexpr double kB = 1.380649e-23;           // Дж/К
inline constexpr double e_charge = 1.602176634e-19;  // Кл
inline constexpr double debye = 3.33564095198e-30;   // Кл·м

// см⁻¹ ↔ Дж, эВ
inline constexpr double J_per_cm = h * c * 100.0;
inline constexpr double cm_per_eV = e_charge / J_per_cm;  // ≈ 8065.54
inline constexpr double kB_cm = kB / J_per_cm;             // ≈ 0.695035 см⁻¹/К

// Угловая частота, соответствующая 1 см⁻¹: 2πc·100 рад/с
inline constexpr double rad_s_per_cm = 2.0 * pi * c * 100.0;
inline constexpr double rad_ps_per_cm = rad_s_per_cm * 1e-12;  // ≈ 0.188365
// Единица времени при ħ = 1 и энергиях в см⁻¹
inline constexpr double ps_per_tu = 1.0 / rad_ps_per_cm;       // ≈ 5.3088 пс

// Диполь-дипольная энергия: (1 Д)²/(4πε₀·(1 нм)³) в см⁻¹ ≈ 5.0341
inline constexpr double dd_cm_per_D2_nm3 = debye * debye / (4.0 * pi * eps0 * 1e-27) / J_per_cm;

// Длина волны (нм) ↔ волновое число (см⁻¹)
inline constexpr double nm_to_cm(double lambda_nm) { return 1e7 / lambda_nm; }
inline constexpr double cm_to_nm(double nu_cm) { return 1e7 / nu_cm; }

}  // namespace ssp::units
