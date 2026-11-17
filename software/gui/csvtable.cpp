#include "csvtable.hpp"

#include <QFile>
#include <QRegularExpression>
#include <QTextStream>

#include <cmath>
#include <limits>

bool CsvTable::load(const QString& path, QString* error) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
    if (error) *error = "cannot open " + path;
    return false;
  }
  QTextStream ts(&f);
  *this = CsvTable();
  while (!ts.atEnd()) {
    const QString line = ts.readLine();
    if (line.startsWith('#')) {
      comments << line.mid(1).trimmed();
      continue;
    }
    if (line.trimmed().isEmpty()) continue;
    // Поле в кавычках может содержать запятые (столбец description)
    QStringList cells;
    QString cur;
    bool quoted = false;
    for (QChar c : line) {
      if (c == '"') quoted = !quoted;
      else if (c == ',' && !quoted) {
        cells << cur;
        cur.clear();
      } else cur += c;
    }
    cells << cur;
    if (header.isEmpty()) header = cells;
    else text << cells;
  }
  cols.resize(header.size());
  for (int c = 0; c < header.size(); ++c) {
    cols[c].resize(text.size());
    for (int r = 0; r < text.size(); ++r) {
      bool ok = false;
      const double v = c < text[r].size() ? text[r][c].toDouble(&ok) : 0.0;
      cols[c][r] = ok ? v : std::numeric_limits<double>::quiet_NaN();
    }
  }
  return true;
}

QVector<double> CsvTable::col(const QString& name) const {
  const int c = column(name);
  return c >= 0 ? cols[c] : QVector<double>();
}

bool loadXYTable(const QString& path, QVector<double>& x, QVector<double>& y, QVector<double>& s, QString* error) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
    if (error) *error = "cannot open " + path;
    return false;
  }
  x.clear();
  y.clear();
  s.clear();
  QTextStream ts(&f);
  const QRegularExpression sep("[,;\\t ]+");
  bool has_sigma = true;
  while (!ts.atEnd()) {
    const QString line = ts.readLine().trimmed();
    if (line.isEmpty() || line.startsWith('#')) continue;
    const QStringList p = line.split(sep, Qt::SkipEmptyParts);
    if (p.size() < 2) continue;
    bool ok1, ok2, ok3 = false;
    const double a = p[0].toDouble(&ok1), b = p[1].toDouble(&ok2);
    if (!ok1 || !ok2) continue;  // заголовок
    const double c = p.size() >= 3 ? p[2].toDouble(&ok3) : 0.0;
    x << a;
    y << b;
    s << (ok3 ? c : 0.0);
    has_sigma = has_sigma && ok3 && c > 0.0;
  }
  if (!has_sigma) s.clear();
  if (x.size() < 3) {
    if (error) *error = "the file has fewer than three numeric rows";
    return false;
  }
  return true;
}
