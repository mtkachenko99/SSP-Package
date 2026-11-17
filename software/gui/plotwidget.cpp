#include "plotwidget.hpp"

#include <QApplication>
#include <QContextMenuEvent>
#include <QFile>
#include <QFileDialog>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QSvgGenerator>
#include <QTextStream>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// «Красивые» деления оси: шаг 1, 2, 5 × 10^k.
QVector<double> niceTicks(double a, double b, int target, double& step) {
  QVector<double> t;
  if (!(b > a)) return t;
  const double raw = (b - a) / std::max(target, 1);
  const double mag = std::pow(10.0, std::floor(std::log10(raw)));
  const double f = raw / mag;
  step = (f < 1.5 ? 1.0 : f < 3.5 ? 2.0 : f < 7.5 ? 5.0 : 10.0) * mag;
  for (double v = std::ceil(a / step) * step; v <= b + 1e-9 * step; v += step) t.push_back(std::abs(v) < 1e-12 * step ? 0.0 : v);
  return t;
}

QString tickLabel(double v, double step) {
  const double av = std::abs(v);
  if (av != 0.0 && (av >= 1e5 || av < 1e-3)) return QString::number(v, 'g', 3);
  const int dec = std::max(0, -static_cast<int>(std::floor(std::log10(step) + 1e-9)));
  return QString::number(v, 'f', std::min(dec, 6));
}

}  // namespace

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent) {
  setMouseTracking(true);
  setMinimumSize(320, 180);
  setFocusPolicy(Qt::ClickFocus);
}

QColor PlotWidget::paletteColor(int i) {
  // Палитра Okabe–Ito (различима при дальтонизме) + доп. цвета
  static const QColor c[] = {QColor(0, 114, 178), QColor(213, 94, 0), QColor(0, 158, 115), QColor(204, 121, 167),
                             QColor(230, 159, 0), QColor(86, 180, 233), QColor(120, 120, 120), QColor(160, 50, 50)};
  return c[i % 8];
}

void PlotWidget::clear() {
  series_.clear();
  refLines_.clear();
  vLines_.clear();
  vLabels_.clear();
  update();
}

QColor PlotWidget::rdbu(double t) {
  // Расходящаяся палитра (синий — 0, белый — 0.5, красный — 1) для знакопеременных величин
  static const double c[5][3] = {{33, 102, 172}, {146, 197, 222}, {247, 247, 247}, {244, 165, 130}, {178, 24, 43}};
  t = std::clamp(t, 0.0, 1.0) * 4.0;
  const int i = std::min(3, int(t));
  const double f = t - i;
  return QColor(int(c[i][0] + f * (c[i + 1][0] - c[i][0])), int(c[i][1] + f * (c[i + 1][1] - c[i][1])),
                int(c[i][2] + f * (c[i + 1][2] - c[i][2])));
}

void PlotWidget::setHeatmap(const QVector<double>& x, const QVector<double>& y, const QVector<QVector<double>>& z,
                            const QString& label) {
  hx_ = x;
  hy_ = y;
  hz_ = z;
  hlabel_ = label;
  hzmin_ = std::numeric_limits<double>::infinity();
  hzmax_ = -hzmin_;
  for (const auto& row : z)
    for (double v : row)
      if (std::isfinite(v)) {
        hzmin_ = std::min(hzmin_, v);
        hzmax_ = std::max(hzmax_, v);
      }
  if (!(hzmax_ > hzmin_)) hzmax_ = hzmin_ + 1;
  // Знакопеременные данные: симметричная шкала вокруг нуля и расходящаяся палитра
  hdiverging_ = hzmin_ < 0.0 && hzmax_ > 0.0;
  if (hdiverging_) {
    const double m = std::max(-hzmin_, hzmax_);
    hzmin_ = -m;
    hzmax_ = m;
  }
  update();
}

void PlotWidget::clearHeatmap() {
  hx_.clear();
  hy_.clear();
  hz_.clear();
  update();
}

QColor PlotWidget::viridis(double t) {
  static const double c[9][3] = {{68, 1, 84},   {71, 44, 122},  {59, 81, 139},  {44, 113, 142}, {33, 144, 141},
                                 {39, 173, 129}, {92, 200, 99}, {170, 220, 50}, {253, 231, 37}};
  t = std::clamp(t, 0.0, 1.0) * 8.0;
  const int i = std::min(7, int(t));
  const double f = t - i;
  return QColor(int(c[i][0] + f * (c[i + 1][0] - c[i][0])), int(c[i][1] + f * (c[i + 1][1] - c[i][1])),
                int(c[i][2] + f * (c[i + 1][2] - c[i][2])));
}

int PlotWidget::addSeries(const Series& s) {
  series_.push_back(s);
  if (!series_.back().color.isValid()) series_.back().color = paletteColor(series_.size() - 1);
  return series_.size() - 1;
}

void PlotWidget::setTitle(const QString& t) {
  title_ = t;
  update();
}

void PlotWidget::setLabels(const QString& x, const QString& y) {
  xlabel_ = x;
  ylabel_ = y;
  update();
}

void PlotWidget::autoscale() {
  double a = std::numeric_limits<double>::infinity(), b = -a;
  for (const auto& s : series_) {
    if (!s.visible) continue;
    for (double v : s.x) {
      a = std::min(a, v);
      b = std::max(b, v);
    }
  }
  if (!hx_.isEmpty()) {
    a = std::min(a, hx_.front());
    b = std::max(b, hx_.back());
  }
  if (!std::isfinite(a)) {
    a = 0;
    b = 1;
  }
  if (b <= a) b = a + 1;
  // Для точечных рядов — поля по x, чтобы крайние точки не сливались с рамкой
  bool points = false;
  for (const auto& s : series_) points = points || (s.visible && s.style == Style::Points);
  if (points && hx_.isEmpty()) {
    const double pad = 0.03 * (b - a);
    a -= pad;
    b += pad;
  }
  x0_ = a;
  x1_ = b;
  autoscaleY();
  update();
  emit xRangeChanged(x0_, x1_);
}

void PlotWidget::autoscaleY() {
  if (hy_.size() >= 2) {
    const double h0 = 0.5 * (hy_[1] - hy_[0]), h1 = 0.5 * (hy_.back() - hy_[hy_.size() - 2]);
    y0_ = hy_.front() - h0;
    y1_ = hy_.back() + h1;
    return;
  }
  double a = std::numeric_limits<double>::infinity(), b = -a;
  for (const auto& s : series_) {
    if (!s.visible) continue;
    for (int i = 0; i < s.x.size(); ++i) {
      if (s.x[i] < x0_ || s.x[i] > x1_ || !std::isfinite(s.y[i])) continue;
      double lo = s.y[i], hi = s.y[i];
      if (!s.lo.isEmpty()) {
        lo = std::min(lo, s.lo[i]);
        hi = std::max(hi, s.hi[i]);
      }
      a = std::min(a, lo);
      b = std::max(b, hi);
    }
  }
  for (double r : refLines_) {
    a = std::min(a, r);
    b = std::max(b, r);
  }
  if (zeroLine_) {
    a = std::min(a, 0.0);
    b = std::max(b, 0.0);
  }
  if (!std::isfinite(a)) {
    a = 0;
    b = 1;
  }
  const double pad = (b > a) ? 0.06 * (b - a) : std::max(1e-12, std::abs(a) * 0.1 + 1e-12);
  y0_ = a - pad;
  y1_ = b + pad;
}

void PlotWidget::setXRange(double a, double b, bool emitSignal) {
  if (!(b > a)) return;
  x0_ = a;
  x1_ = b;
  autoscaleY();
  update();
  if (emitSignal) emit xRangeChanged(a, b);
}

QRectF PlotWidget::plotRect(const QRect& area) const {
  const QFontMetrics fm(font());
  const double left = fm.horizontalAdvance("−00000.0") + fm.height() + 10;
  const double top = (title_.isEmpty() ? 8 : fm.height() + 10) + (nmAxis_ ? 2 * fm.height() + 6 : 0);
  const double bottom = 2 * fm.height() + 14;
  const double right = hz_.isEmpty() ? 14 : 14 + 28 + fm.horizontalAdvance("0.000e+00") + fm.height();
  return QRectF(area.left() + left, area.top() + top, area.width() - left - right, area.height() - top - bottom);
}

QPointF PlotWidget::toPixel(const QRectF& r, double x, double y) const {
  return {r.left() + (x - x0_) / (x1_ - x0_) * r.width(), r.bottom() - (y - y0_) / (y1_ - y0_) * r.height()};
}
double PlotWidget::xFromPixel(const QRectF& r, double px) const { return x0_ + (px - r.left()) / r.width() * (x1_ - x0_); }
double PlotWidget::yFromPixel(const QRectF& r, double py) const { return y0_ + (r.bottom() - py) / r.height() * (y1_ - y0_); }

void PlotWidget::paintEvent(QPaintEvent*) {
  QPainter p(this);
  render(p, rect());
  if (dragging_) {
    p.setPen(QPen(palette().color(QPalette::Highlight), 1, Qt::DashLine));
    QColor f = palette().color(QPalette::Highlight);
    f.setAlpha(40);
    p.setBrush(f);
    const QRectF r = plotRect(rect());
    p.drawRect(QRectF(QPointF(dragStart_.x(), r.top()), QPointF(dragNow_.x(), r.bottom())));
  }
}

void PlotWidget::render(QPainter& p, const QRect& area) {
  p.setRenderHint(QPainter::Antialiasing, true);
  const QColor bg = palette().color(QPalette::Base), fg = palette().color(QPalette::Text);
  QColor grid = fg;
  grid.setAlpha(35);
  p.fillRect(area, bg);
  const QRectF r = plotRect(area);
  const QFontMetrics fm(p.font());

  double xs, ys;
  const auto xt = niceTicks(x0_, x1_, std::max(3, int(r.width() / 90)), xs);
  const auto yt = niceTicks(y0_, y1_, std::max(3, int(r.height() / 45)), ys);
  p.setPen(QPen(grid, 1));
  for (double v : xt) p.drawLine(toPixel(r, v, y0_), toPixel(r, v, y1_));
  for (double v : yt) p.drawLine(toPixel(r, x0_, v), toPixel(r, x1_, v));

  p.save();
  p.setClipRect(r);
  if (!hz_.isEmpty() && hx_.size() >= 2 && hy_.size() >= 1) {
    p.setPen(Qt::NoPen);
    auto edge = [](const QVector<double>& v, int i) {
      if (v.size() == 1) return i == 0 ? v[0] - 0.5 : v[0] + 0.5;
      if (i == 0) return v[0] - 0.5 * (v[1] - v[0]);
      if (i == v.size()) return v.back() + 0.5 * (v.back() - v[v.size() - 2]);
      return 0.5 * (v[i - 1] + v[i]);
    };
    for (int iy = 0; iy < hy_.size() && iy < hz_.size(); ++iy) {
      const double ya = edge(hy_, iy), yb = edge(hy_, iy + 1);
      for (int ix = 0; ix < hx_.size() && ix < hz_[iy].size(); ++ix) {
        const double v = hz_[iy][ix];
        if (!std::isfinite(v)) continue;
        const QPointF p0 = toPixel(r, edge(hx_, ix), ya), p1 = toPixel(r, edge(hx_, ix + 1), yb);
        const double tt = (v - hzmin_) / (hzmax_ - hzmin_);
        p.setBrush(hdiverging_ ? rdbu(tt) : viridis(tt));
        p.drawRect(QRectF(p0, p1).normalized().adjusted(-0.5, -0.5, 0.5, 0.5));
      }
    }
  }
  if (zeroLine_) {
    p.setPen(QPen(fg, 1));
    p.drawLine(toPixel(r, x0_, 0), toPixel(r, x1_, 0));
  }
  for (double v : refLines_) {
    QColor c = fg;
    c.setAlpha(120);
    p.setPen(QPen(c, 1, Qt::DashLine));
    p.drawLine(toPixel(r, x0_, v), toPixel(r, x1_, v));
  }
  for (int i = 0; i < vLines_.size(); ++i) {
    QColor c = fg;
    c.setAlpha(150);
    p.setPen(QPen(c, 1, Qt::DashDotLine));
    const QPointF a = toPixel(r, vLines_[i], y0_), b = toPixel(r, vLines_[i], y1_);
    p.drawLine(a, b);
    if (i < vLabels_.size()) p.drawText(QPointF(a.x() + 3, a.y() - 4), vLabels_[i]);
  }
  for (const auto& s : series_) {
    if (!s.visible || s.x.isEmpty()) continue;
    if (!s.lo.isEmpty()) {
      QPainterPath band;
      band.moveTo(toPixel(r, s.x[0], s.hi[0]));
      for (int i = 1; i < s.x.size(); ++i) band.lineTo(toPixel(r, s.x[i], s.hi[i]));
      for (int i = s.x.size() - 1; i >= 0; --i) band.lineTo(toPixel(r, s.x[i], s.lo[i]));
      band.closeSubpath();
      QColor c = s.color;
      c.setAlpha(55);
      p.fillPath(band, c);
    }
    if (s.style == Style::Points) {
      p.setPen(Qt::NoPen);
      p.setBrush(s.color);
      for (int i = 0; i < s.x.size(); ++i)
        if (std::isfinite(s.y[i])) p.drawEllipse(toPixel(r, s.x[i], s.y[i]), 2.2, 2.2);
      continue;
    }
    QPen pen(s.color, s.width);
    if (s.style == Style::Dashed) pen.setStyle(Qt::DashLine);
    p.setPen(pen);
    p.setBrush(Qt::NoBrush);
    QPainterPath path;
    bool started = false;
    for (int i = 0; i < s.x.size(); ++i) {
      if (!std::isfinite(s.y[i])) {
        started = false;
        continue;
      }
      const QPointF q = toPixel(r, s.x[i], s.y[i]);
      if (!started) {
        path.moveTo(q);
        started = true;
      } else if (s.style == Style::Steps) {
        path.lineTo(QPointF(q.x(), path.currentPosition().y()));
        path.lineTo(q);
      } else {
        path.lineTo(q);
      }
    }
    p.drawPath(path);
  }
  p.restore();

  // Рамка, подписи делений
  p.setPen(QPen(fg, 1));
  p.setBrush(Qt::NoBrush);
  p.drawRect(r);
  for (double v : xt) {
    const QPointF q = toPixel(r, v, y0_);
    p.drawLine(q, q + QPointF(0, 4));
    const QString t = tickLabel(v, xs);
    p.drawText(QPointF(q.x() - fm.horizontalAdvance(t) / 2.0, r.bottom() + fm.ascent() + 6), t);
  }
  for (double v : yt) {
    const QPointF q = toPixel(r, x0_, v);
    p.drawLine(q, q - QPointF(4, 0));
    const QString t = tickLabel(v, ys);
    p.drawText(QPointF(r.left() - fm.horizontalAdvance(t) - 7, q.y() + fm.ascent() / 2.5), t);
  }
  if (nmAxis_ && x0_ > 0) {
    // Верхняя ось λ = 10⁷/ν̃ нм
    double ns;
    const auto nt = niceTicks(1e7 / x1_, 1e7 / x0_, std::max(3, int(r.width() / 90)), ns);
    for (double lam : nt) {
      const double nu = 1e7 / lam;
      if (nu < x0_ || nu > x1_) continue;
      const QPointF q = toPixel(r, nu, y1_);
      p.drawLine(q, q - QPointF(0, 4));
      const QString t = tickLabel(lam, ns);
      p.drawText(QPointF(q.x() - fm.horizontalAdvance(t) / 2.0, r.top() - 7), t);
    }
    const QString lt = "λ, nm";
    p.drawText(QPointF(r.center().x() - fm.horizontalAdvance(lt) / 2.0, r.top() - fm.height() - 8), lt);
  }
  if (!xlabel_.isEmpty())
    p.drawText(QPointF(r.center().x() - fm.horizontalAdvance(xlabel_) / 2.0, r.bottom() + 2 * fm.height() + 8), xlabel_);
  if (!ylabel_.isEmpty()) {
    p.save();
    p.translate(area.left() + fm.height() - 2, r.center().y() + fm.horizontalAdvance(ylabel_) / 2.0);
    p.rotate(-90);
    p.drawText(QPointF(0, 0), ylabel_);
    p.restore();
  }
  if (!title_.isEmpty()) {
    QFont f = p.font();
    f.setBold(true);
    p.setFont(f);
    const QString t = QFontMetrics(f).elidedText(title_, Qt::ElideRight, int(area.right() - r.left() - 4));
    p.drawText(QPointF(r.left(), area.top() + QFontMetrics(f).ascent() + 4), t);
    f.setBold(false);
    p.setFont(f);
  }
  // Цветовая шкала тепловой карты
  if (!hz_.isEmpty()) {
    const QRectF cb(r.right() + 14, r.top(), 14, r.height());
    for (int k = 0; k < int(cb.height()); ++k) {
      p.setPen(hdiverging_ ? rdbu(1.0 - k / cb.height()) : viridis(1.0 - k / cb.height()));
      p.drawLine(QPointF(cb.left(), cb.top() + k), QPointF(cb.right(), cb.top() + k));
    }
    p.setPen(fg);
    p.setBrush(Qt::NoBrush);
    p.drawRect(cb);
    double cs;
    for (double v : niceTicks(hzmin_, hzmax_, 6, cs)) {
      const double yy = cb.bottom() - (v - hzmin_) / (hzmax_ - hzmin_) * cb.height();
      p.drawLine(QPointF(cb.right(), yy), QPointF(cb.right() + 3, yy));
      p.drawText(QPointF(cb.right() + 5, yy + fm.ascent() / 2.5), tickLabel(v, cs));
    }
    if (!hlabel_.isEmpty()) {
      p.save();
      p.translate(area.right() - 4, cb.center().y() + fm.horizontalAdvance(hlabel_) / 2.0);
      p.rotate(-90);
      p.drawText(QPointF(0, 0), hlabel_);
      p.restore();
    }
  }
  // Легенда
  if (legend_) {
    int n = 0, w = 0;
    for (const auto& s : series_)
      if (s.visible && !s.name.isEmpty()) {
        ++n;
        w = std::max(w, fm.horizontalAdvance(s.name));
      }
    if (n > 0) {
      const QRectF box(r.right() - w - 40, r.top() + 6, w + 34, n * (fm.height() + 2) + 6);
      QColor b = bg;
      b.setAlpha(215);
      p.setBrush(b);
      QColor e = fg;
      e.setAlpha(60);
      p.setPen(e);
      p.drawRect(box);
      int k = 0;
      for (const auto& s : series_) {
        if (!s.visible || s.name.isEmpty()) continue;
        const double y = box.top() + 4 + k * (fm.height() + 2) + fm.height() / 2.0;
        if (s.style == Style::Points) {
          p.setPen(Qt::NoPen);
          p.setBrush(s.color);
          p.drawEllipse(QPointF(box.left() + 15, y), 3.0, 3.0);
          p.setBrush(Qt::NoBrush);
        } else {
          QPen pen(s.color, 2.2);
          if (s.style == Style::Dashed) pen.setStyle(Qt::DashLine);
          p.setPen(pen);
          p.drawLine(QPointF(box.left() + 6, y), QPointF(box.left() + 24, y));
        }
        p.setPen(fg);
        p.drawText(QPointF(box.left() + 30, y + fm.ascent() / 2.5), s.name);
        ++k;
      }
    }
  }
}

void PlotWidget::mousePressEvent(QMouseEvent* e) {
  if (e->button() == Qt::LeftButton) {
    dragging_ = true;
    dragStart_ = dragNow_ = e->pos();
  } else if (e->button() == Qt::MiddleButton) {
    panning_ = true;
    dragStart_ = e->pos();
    panX0_ = x0_;
    panX1_ = x1_;
  }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* e) {
  const QRectF r = plotRect(rect());
  if (dragging_) {
    dragNow_ = e->pos();
    update();
  } else if (panning_) {
    const double dx = (e->pos().x() - dragStart_.x()) / r.width() * (panX1_ - panX0_);
    setXRange(panX0_ - dx, panX1_ - dx);
  }
  if (r.contains(e->pos())) emit cursorMoved(xFromPixel(r, e->pos().x()), yFromPixel(r, e->pos().y()));
}

void PlotWidget::mouseReleaseEvent(QMouseEvent* e) {
  if (dragging_ && e->button() == Qt::LeftButton) {
    dragging_ = false;
    const QRectF r = plotRect(rect());
    if (std::abs(dragNow_.x() - dragStart_.x()) > 6) {
      const double a = xFromPixel(r, std::min(dragStart_.x(), dragNow_.x()));
      const double b = xFromPixel(r, std::max(dragStart_.x(), dragNow_.x()));
      setXRange(a, b);
    }
    update();
  }
  panning_ = false;
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent*) { autoscale(); }

void PlotWidget::wheelEvent(QWheelEvent* e) {
  const QRectF r = plotRect(rect());
  const double xc = xFromPixel(r, e->position().x());
  const double f = e->angleDelta().y() > 0 ? 0.8 : 1.25;
  setXRange(xc - (xc - x0_) * f, xc + (x1_ - xc) * f);
}

void PlotWidget::contextMenuEvent(QContextMenuEvent* e) {
  QMenu m(this);
  QAction* a1 = m.addAction("Full range");
  QAction* a2 = m.addAction("Save image (PNG/SVG)…");
  QAction* a3 = m.addAction("Save data (CSV)…");
  QAction* a4 = m.addAction(legend_ ? "Hide legend" : "Show legend");
  QAction* r = m.exec(e->globalPos());
  if (r == a1) autoscale();
  else if (r == a2) {
    const QString f = QFileDialog::getSaveFileName(this, "Save plot", "plot.svg", "SVG (*.svg);;PNG (*.png)");
    if (!f.isEmpty()) exportImage(f);
  } else if (r == a3) {
    const QString f = QFileDialog::getSaveFileName(this, "Save data", "plot.csv", "CSV (*.csv)");
    if (!f.isEmpty()) exportCsv(f);
  } else if (r == a4) setLegendVisible(!legend_);
}

bool PlotWidget::exportImage(const QString& path) {
  const QSize sz(std::max(width(), 900), std::max(height(), 520));
  if (path.endsWith(".svg", Qt::CaseInsensitive)) {
    QSvgGenerator gen;
    gen.setFileName(path);
    gen.setSize(sz);
    gen.setViewBox(QRect(QPoint(0, 0), sz));
    gen.setTitle(title_);
    QPainter p(&gen);
    p.setFont(font());
    render(p, QRect(QPoint(0, 0), sz));
    return true;
  }
  QImage img(sz * 2, QImage::Format_ARGB32);
  img.setDevicePixelRatio(2.0);
  QPainter p(&img);
  p.setFont(font());
  render(p, QRect(QPoint(0, 0), sz));
  p.end();
  return img.save(path);
}

bool PlotWidget::exportCsv(const QString& path) const {
  QFile f(path);
  if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
  QTextStream ts(&f);
  ts << "# " << title_ << "\n# x: " << xlabel_ << "; y: " << ylabel_ << "\nseries,x,y,lo,hi\n";
  for (const auto& s : series_)
    for (int i = 0; i < s.x.size(); ++i)
      ts << '"' << s.name << "\"," << s.x[i] << ',' << s.y[i] << ',' << (s.lo.isEmpty() ? QString() : QString::number(s.lo[i]))
         << ',' << (s.hi.isEmpty() ? QString() : QString::number(s.hi[i])) << "\n";
  return true;
}
