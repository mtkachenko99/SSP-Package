#include "geometryview.hpp"

#include <QMatrix4x4>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <numeric>

GeometryView::GeometryView(QWidget* parent) : QWidget(parent) { setMinimumSize(240, 200); }

void GeometryView::setAggregate(const QVector<QVector3D>& r, const QVector<QVector3D>& mu) {
  r_ = r;
  mu_ = mu;
  values_.clear();
  message_.clear();
  update();
}

void GeometryView::setParticle(bool on, double a_xy, double a_z) {
  np_ = on;
  axy_ = a_xy;
  az_ = a_z;
  update();
}

void GeometryView::setValues(const QVector<double>& v, const QString& label) {
  values_ = v;
  valueLabel_ = label;
  update();
}

void GeometryView::setMessage(const QString& m) {
  message_ = m;
  update();
}

void GeometryView::paintEvent(QPaintEvent*) {
  QPainter p(this);
  p.setRenderHint(QPainter::Antialiasing);
  const QColor bg = palette().color(QPalette::Base), fg = palette().color(QPalette::Text);
  p.fillRect(rect(), bg);
  if (!message_.isEmpty()) {
    p.setPen(QColor(200, 60, 60));
    p.drawText(rect().adjusted(8, 8, -8, -8), Qt::AlignCenter | Qt::TextWordWrap, message_);
    return;
  }
  if (r_.isEmpty()) {
    p.setPen(fg);
    p.drawText(rect(), Qt::AlignCenter, "no geometry");
    return;
  }
  QMatrix4x4 R;
  R.rotate(pitch_, 1, 0, 0);
  R.rotate(yaw_, 0, 0, 1);
  // Масштаб: охватить агрегат и частицу
  double ext = 1.0;
  for (const auto& v : r_) ext = std::max(ext, double(v.length()));
  if (np_) ext = std::max(ext, std::max(axy_, az_));
  QVector3D c(0, 0, 0);
  if (!np_) {
    for (const auto& v : r_) c += v;
    c /= r_.size();
    ext = 1.0;
    for (const auto& v : r_) ext = std::max(ext, double((v - c).length()));
  }
  const double s = 0.42 * std::min(width(), height()) / ext * zoom_;
  const QPointF o = rect().center();
  auto proj = [&](const QVector3D& v) {
    const QVector3D q = R.map(v - c);
    return std::make_pair(QPointF(o.x() + s * q.x(), o.y() - s * q.z()), double(q.y()));
  };

  // Оси
  const double al = 0.25 * ext;
  const char* names[3] = {"x", "y", "z"};
  const QVector3D ax[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const QPointF corner(40, height() - 40);
  for (int k = 0; k < 3; ++k) {
    const QVector3D q = R.map(ax[k]);
    const QPointF e = corner + QPointF(25 * q.x(), -25 * q.z());
    p.setPen(QPen(fg, 1.2));
    p.drawLine(corner, e);
    p.drawText(e + QPointF(2, -2), names[k]);
  }
  (void)al;

  // Наночастица: проекция эллипсоида — эллипс, M = (R A⁻¹ Rᵀ) на экранной плоскости (x, z)
  auto drawParticle = [&]() {
    if (!np_ || axy_ <= 0) return;
    const double d[3] = {axy_ * axy_, axy_ * axy_, az_ * az_};
    double m[3][3] = {{0}};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k) m[i][j] += R(i, k) * d[k] * R(j, k);
    // экранные координаты: (x, z) → индексы 0 и 2
    const double a = m[0][0], b = m[0][2], cc = m[2][2];
    const double tr = a + cc, det = a * cc - b * b;
    const double l1 = 0.5 * tr + std::sqrt(std::max(0.0, 0.25 * tr * tr - det)), l2 = tr - l1;
    const double ang = 0.5 * std::atan2(2 * b, a - cc);
    QPainterPath path;
    for (int k = 0; k <= 96; ++k) {
      const double t = 2 * M_PI * k / 96;
      const double u = std::sqrt(l1) * std::cos(t), v = std::sqrt(std::max(l2, 0.0)) * std::sin(t);
      const double X = u * std::cos(ang) - v * std::sin(ang), Z = u * std::sin(ang) + v * std::cos(ang);
      const QPointF q(o.x() + s * X, o.y() - s * Z);
      if (k == 0) path.moveTo(q);
      else path.lineTo(q);
    }
    QRadialGradient g(o, s * std::sqrt(l1));
    g.setColorAt(0, QColor(255, 215, 120, 200));
    g.setColorAt(1, QColor(190, 140, 40, 200));
    p.setPen(QPen(QColor(150, 110, 30), 1.2));
    p.setBrush(g);
    p.drawPath(path);
  };

  // Молекулы по глубине (дальние первыми); частица рисуется между дальними и ближними
  QVector<int> idx(r_.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::vector<double> depth(r_.size());
  for (int i = 0; i < r_.size(); ++i) depth[i] = proj(r_[i]).second;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return depth[a] > depth[b]; });
  double vmin = 0, vmax = 1;
  if (!values_.isEmpty()) {
    vmin = *std::min_element(values_.begin(), values_.end());
    vmax = *std::max_element(values_.begin(), values_.end());
    if (vmax <= vmin) vmax = vmin + 1;
  }
  const double rad = std::clamp(0.18 * s, 2.5, 9.0);
  bool npDrawn = false;
  for (int i : idx) {
    if (!npDrawn && depth[i] < 0) {
      drawParticle();
      npDrawn = true;
    }
    const auto [q, dz] = proj(r_[i]);
    QColor col;
    if (!values_.isEmpty()) {
      const double t = (values_[i] - vmin) / (vmax - vmin);
      col = QColor::fromHsvF(0.66 * (1.0 - t), 0.85, 0.9);
    } else {
      const double t = std::clamp(0.5 + dz / (2 * ext), 0.0, 1.0);
      col = QColor::fromRgbF(0.15 + 0.5 * t, 0.35 + 0.4 * t, 0.75);
    }
    const QVector3D tip = r_[i] + 0.45f * mu_[i];
    const QVector3D tail = r_[i] - 0.45f * mu_[i];
    p.setPen(QPen(fg, 1.4));
    const QPointF a = proj(tail).first, b = proj(tip).first;
    p.drawLine(a, b);
    // наконечник
    const QPointF d = b - a;
    const double L = std::hypot(d.x(), d.y());
    if (L > 3) {
      const QPointF u = d / L, n(-u.y(), u.x());
      QPolygonF head;
      head << b << b - 5 * u + 2.5 * n << b - 5 * u - 2.5 * n;
      p.setBrush(fg);
      p.drawPolygon(head);
    }
    p.setPen(QPen(col.darker(150), 1));
    p.setBrush(col);
    p.drawEllipse(q, rad, rad);
  }
  if (!npDrawn) drawParticle();

  // Масштабная линейка
  {
    const double target = 80.0 / s;  // нм на ~80 px
    const double mag = std::pow(10.0, std::floor(std::log10(target)));
    const double f = target / mag, len = (f < 2 ? 1 : f < 5 ? 2 : 5) * mag;
    const QPointF b0(width() - 20 - len * s, height() - 20), b1(width() - 20, height() - 20);
    p.setPen(QPen(fg, 2));
    p.drawLine(b0, b1);
    p.drawLine(b0 + QPointF(0, -4), b0 + QPointF(0, 4));
    p.drawLine(b1 + QPointF(0, -4), b1 + QPointF(0, 4));
    p.drawText(QPointF(b0.x(), b0.y() - 6), QString("%1 nm").arg(len));
  }
  p.setPen(fg);
  QString info = QString("molecules: %1").arg(r_.size());
  if (!values_.isEmpty()) info += QString("   %1: %2 … %3").arg(valueLabel_).arg(vmin, 0, 'g', 3).arg(vmax, 0, 'g', 3);
  p.drawText(QPointF(8, 16), info);
}

void GeometryView::mousePressEvent(QMouseEvent* e) { last_ = e->pos(); }
void GeometryView::mouseMoveEvent(QMouseEvent* e) {
  if (!(e->buttons() & Qt::LeftButton)) return;
  yaw_ += 0.5 * (e->pos().x() - last_.x());
  pitch_ = std::clamp(pitch_ + 0.5 * (e->pos().y() - last_.y()), -89.0, 89.0);
  last_ = e->pos();
  update();
}
void GeometryView::wheelEvent(QWheelEvent* e) {
  zoom_ *= e->angleDelta().y() > 0 ? 1.15 : 1 / 1.15;
  update();
}
void GeometryView::mouseDoubleClickEvent(QMouseEvent*) {
  yaw_ = -35;
  pitch_ = 25;
  zoom_ = 1.0;
  update();
}
