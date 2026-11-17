#pragma once
// Научный 2D-график на QPainter: линии, точки, полосы неопределённости, легенда, масштабирование мышью,
// связывание оси x между графиками (основная панель + панель остатков), экспорт PNG/SVG/CSV.

#include <QColor>
#include <QStringList>
#include <QString>
#include <QVector>
#include <QWidget>

class PlotWidget : public QWidget {
  Q_OBJECT
 public:
  enum class Style { Line, Dashed, Points, Steps };

  struct Series {
    QString name;
    QVector<double> x, y;
    QVector<double> lo, hi;  // полоса (пустые — без полосы)
    QColor color;
    Style style = Style::Line;
    double width = 1.6;
    bool visible = true;
  };

  explicit PlotWidget(QWidget* parent = nullptr);

  void clear();
  int addSeries(const Series& s);
  QVector<Series>& series() { return series_; }
  void setTitle(const QString& t);
  void setLabels(const QString& x, const QString& y);
  void setZeroLine(bool on) { zeroLine_ = on; update(); }
  void setReferenceLines(const QVector<double>& y) { refLines_ = y; update(); }
  void setVerticalLines(const QVector<double>& x, const QStringList& labels) { vLines_ = x; vLabels_ = labels; update(); }
  void setSecondaryNmAxis(bool on) { nmAxis_ = on; update(); }  // верхняя ось λ (нм) для оси x в см⁻¹
  void setLegendVisible(bool on) { legend_ = on; update(); }
  // Тепловая карта под рядами: x[nx], y[ny], z[iy][ix]; подпись цветовой шкалы.
  void setHeatmap(const QVector<double>& x, const QVector<double>& y, const QVector<QVector<double>>& z,
                  const QString& label);
  void clearHeatmap();

  void autoscale();
  void setXRange(double a, double b, bool emitSignal = true);
  double xmin() const { return x0_; }
  double xmax() const { return x1_; }

  bool exportImage(const QString& path);  // .png или .svg
  bool exportCsv(const QString& path) const;

  static QColor paletteColor(int i);

 signals:
  void xRangeChanged(double a, double b);
  void cursorMoved(double x, double y);

 protected:
  void paintEvent(QPaintEvent*) override;
  void mousePressEvent(QMouseEvent*) override;
  void mouseMoveEvent(QMouseEvent*) override;
  void mouseReleaseEvent(QMouseEvent*) override;
  void mouseDoubleClickEvent(QMouseEvent*) override;
  void wheelEvent(QWheelEvent*) override;
  void contextMenuEvent(QContextMenuEvent*) override;

 private:
  void render(QPainter& p, const QRect& area);
  QRectF plotRect(const QRect& area) const;
  QPointF toPixel(const QRectF& r, double x, double y) const;
  double xFromPixel(const QRectF& r, double px) const;
  double yFromPixel(const QRectF& r, double py) const;
  void autoscaleY();

  QVector<Series> series_;
  QVector<double> hx_, hy_;
  QVector<QVector<double>> hz_;
  QString hlabel_;
  double hzmin_ = 0, hzmax_ = 1;
  static QColor viridis(double t);
  QString title_, xlabel_, ylabel_;
  double x0_ = 0, x1_ = 1, y0_ = 0, y1_ = 1;
  bool zeroLine_ = false, nmAxis_ = false, legend_ = true;
  QVector<double> refLines_, vLines_;
  QStringList vLabels_;
  bool hdiverging_ = false;
  static QColor rdbu(double t);
  bool dragging_ = false, panning_ = false;
  QPoint dragStart_, dragNow_;
  double panX0_ = 0, panX1_ = 0;
};
