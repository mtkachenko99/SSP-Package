#pragma once
// Сравнение модельного спектра с экспериментом (методология M6–M7):
//   y_i ≈ s·f(x_i − δ) + b,  s — масштаб, b — базовая линия (опц.), δ — общий сдвиг по энергии (опц.).
// s, b — взвешенный линейный МНК при каждом δ; δ — одномерная минимизация χ² (сетка + золотое сечение).

#include <QString>
#include <QVector>

struct CompareOptions {
  bool fit_offset = false;
  bool fit_shift = false;
  double shift_range_cm = 300.0;
};

struct CompareResult {
  bool ok = false;
  QString message;
  double scale = 1.0, offset = 0.0, shift = 0.0;
  double chi2 = 0.0, chi2_nu = 0.0, rmse = 0.0;
  double durbin_watson = 0.0;  // ≈2 — остатки некоррелированы; ≪2 — структура (неучтённый эффект)
  double lag1_autocorr = 0.0;
  int n = 0, dof = 0;
  bool have_sigma = false;
  QVector<double> x, y, model, residual, normalized_residual;
};

// mx, my — модель на равномерной или неравномерной возрастающей сетке (см⁻¹); ex, ey, es — эксперимент
// (см⁻¹; es пустой — веса 1, χ²_ν оценивается по разбросу остатков).
CompareResult compareSpectra(const QVector<double>& mx, const QVector<double>& my, const QVector<double>& ex,
                             const QVector<double>& ey, const QVector<double>& es, const CompareOptions& opt);

// Перевод оси x эксперимента в см⁻¹: unit = "cm" | "nm" | "eV". Для нм порядок точек обращается.
void convertToWavenumber(QVector<double>& x, QVector<double>& y, QVector<double>& s, const QString& unit);
