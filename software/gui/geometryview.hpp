#pragma once
// Ортографический 3D-просмотр агрегата: молекулы, переходные диполи, контур наночастицы (сфера/сфероид).
// Левая кнопка — вращение, колесо — масштаб, двойной щелчок — сброс.

#include <QVector3D>
#include <QVector>
#include <QWidget>

class GeometryView : public QWidget {
  Q_OBJECT
 public:
  explicit GeometryView(QWidget* parent = nullptr);
  void setAggregate(const QVector<QVector3D>& r, const QVector<QVector3D>& mu);
  void setParticle(bool on, double a_xy = 0, double a_z = 0);
  void setValues(const QVector<double>& v, const QString& label);  // раскраска молекул (напр. усиление поля)
  void setMessage(const QString& m);

 protected:
  void paintEvent(QPaintEvent*) override;
  void mousePressEvent(QMouseEvent*) override;
  void mouseMoveEvent(QMouseEvent*) override;
  void wheelEvent(QWheelEvent*) override;
  void mouseDoubleClickEvent(QMouseEvent*) override;

 private:
  QVector<QVector3D> r_, mu_;
  QVector<double> values_;
  QString valueLabel_, message_;
  bool np_ = false;
  double axy_ = 0, az_ = 0;
  double yaw_ = -35, pitch_ = 25, zoom_ = 1.0;
  QPoint last_;
};
