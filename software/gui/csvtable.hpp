#pragma once
// Чтение CSV-результатов ssp_aggregate (строки «#» — комментарии) и экспериментальных таблиц.

#include <QString>
#include <QStringList>
#include <QVector>

struct CsvTable {
  QStringList comments;
  QStringList header;
  QVector<QStringList> text;     // сырые ячейки
  QVector<QVector<double>> cols; // числовые столбцы (NaN для нечисловых)

  bool load(const QString& path, QString* error = nullptr);
  int column(const QString& name) const { return header.indexOf(name); }
  QVector<double> col(const QString& name) const;
  int rows() const { return text.size(); }
};

// Экспериментальная таблица: 2–3 числовых столбца (x, y[, σ]), разделители — запятая, ;, табуляция, пробелы.
bool loadXYTable(const QString& path, QVector<double>& x, QVector<double>& y, QVector<double>& s, QString* error);
