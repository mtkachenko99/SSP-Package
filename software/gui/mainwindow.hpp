#pragma once

#include "csvtable.hpp"
#include "experiment.hpp"

#include <QFutureWatcher>
#include <QJsonObject>
#include <QMainWindow>
#include <QMap>

#include <atomic>
#include <memory>

namespace ssp::app {
struct CalibrationOutput;
}

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QProgressBar;
class QTableWidget;
class QTabWidget;
class QToolButton;
class QVBoxLayout;
class FormPanel;
class GeometryView;
class PlotWidget;

class MainWindow : public QMainWindow {
  Q_OBJECT
 public:
  explicit MainWindow(QWidget* parent = nullptr);
  bool openConfig(const QString& path);
  void loadResults(const QString& dir);
  void selectResultTab(int i);
  bool loadExperimentFile(const QString& path, const QString& unit);
  void setFitOptions(bool offset, bool shift);
  void startRun() { runCalculation(); }
  void startCalibration(const QString& path);

 signals:
  void runFinished(bool ok);
  void calibrationFinished(bool ok);

 private slots:
  void newConfig();
  void openConfigDialog();
  bool saveConfig();
  bool saveConfigAs();
  void runCalculation();
  void calibrateDialog();
  void applyCalibration();
  void stopCalculation();
  void openResultsDialog();
  void loadExperimentDialog();
  void clearExperiment();

 private:
  struct ConfigResult {
    QString name;
    int n = 0;
    CsvTable spectra, observables, plexciton, dynamics, heomAbs, heomPop, nearField, geometry, states, vibronic, hybridEm;
  };

  void buildUi();
  QWidget* buildConfigurationsTab();
  QWidget* buildResultsArea();
  void normalizeConfig();
  void syncFormsFromConfig();
  void syncConfigFromForms();
  void refreshJsonView();
  void applyJsonText();
  void refreshConfigList(int select);
  void loadEntryEditor(int row);
  void storeEntryEditor();
  void fillParamTable(const QString& type, const QJsonObject& params);
  void updateGeometryPreview();
  void appendLog(const QString& line);
  void setRunning(bool on);
  QString resolvedOutputDir() const;

  void loadRunResults(const QString& dir);
  void loadScanResults(const QString& dir);
  void plotScan();
  void plotScanObservable();
  void refreshResultViews();
  void plotSpectra();
  void plotPlexciton();
  void plotDynamics();
  void plotHeom();
  void plotNearField();
  void fillSummary();
  void runComparison();
  const ConfigResult* resultByName(const QString& n) const;

  // Конфиг
  QJsonObject cfg_;
  QString cfgPath_, baseDir_;
  bool syncing_ = false;
  int currentEntry_ = -1;

  // Левая панель
  QTabWidget* paramTabs_ = nullptr;
  FormPanel *modelForm_ = nullptr, *calcForm_ = nullptr, *npForm_ = nullptr, *heomForm_ = nullptr, *scanForm_ = nullptr;
  QListWidget* entryList_ = nullptr;
  QLineEdit* entryName_ = nullptr;
  QComboBox* entryType_ = nullptr;
  QTableWidget* paramTable_ = nullptr;
  QCheckBox* centerCheck_ = nullptr;
  QDoubleSpinBox *shift_[3] = {}, *euler_[3] = {};
  GeometryView* geomView_ = nullptr;
  QPlainTextEdit* jsonEdit_ = nullptr;

  // Результаты
  QTabWidget* resultTabs_ = nullptr;
  QString resultsDir_;
  QVector<ConfigResult> results_;
  QListWidget* specConfigs_ = nullptr;
  QComboBox *specMode_ = nullptr, *specUnits_ = nullptr;
  QCheckBox *specNormalize_ = nullptr, *specErrors_ = nullptr, *specPerMolecule_ = nullptr, *specSticks_ = nullptr;
  double monomerCm_ = 0.0;
  PlotWidget *specPlot_ = nullptr, *residPlot_ = nullptr;
  QLabel* compareLabel_ = nullptr;
  QComboBox* compareTarget_ = nullptr;
  QCheckBox *fitOffset_ = nullptr, *fitShift_ = nullptr, *residNormalized_ = nullptr;
  QComboBox* plexConfig_ = nullptr;
  PlotWidget* plexPlot_ = nullptr;
  QComboBox* plexView_ = nullptr;  // сечения / вклад агрегата / ФЛ гибрида
  QComboBox* dynConfig_ = nullptr;
  PlotWidget *dynPlot_ = nullptr, *dynEnergyPlot_ = nullptr;
  QComboBox* heomConfig_ = nullptr;
  PlotWidget *heomAbsPlot_ = nullptr, *heomResidPlot_ = nullptr, *heomPopPlot_ = nullptr;
  QLabel* heomLabel_ = nullptr;
  QComboBox* nfConfig_ = nullptr;
  QComboBox* nfQuantity_ = nullptr;  // скорость тушения / Γ_rad / квантовый выход
  PlotWidget *nfEnhPlot_ = nullptr, *nfDecayPlot_ = nullptr;
  GeometryView* nfGeom_ = nullptr;
  QTableWidget* summaryTable_ = nullptr;
  QPlainTextEdit* log_ = nullptr;
  // Скан
  struct ScanShapeData {
    QMap<QString, CsvTable> maps;
    CsvTable branches, coupling, observables;
  };
  QString scanParameter_;
  QVector<double> scanValues_;
  QMap<QString, ScanShapeData> scan_;
  QWidget* scanTab_ = nullptr;
  QComboBox *scanShape_ = nullptr, *scanMap_ = nullptr, *scanObs_ = nullptr;
  QCheckBox *scanNormRows_ = nullptr, *scanOverlay_ = nullptr;
  PlotWidget *scanPlot_ = nullptr, *scanObsPlot_ = nullptr;
  QLabel* scanCoupling_ = nullptr;

  // Эксперимент
  QVector<double> expX_, expY_, expS_;
  QString expName_;

  // Запуск
  QFutureWatcher<void> watcher_;
  // Калибровка по мономеру
  QFutureWatcher<void> calWatcher_;
  std::shared_ptr<ssp::app::CalibrationOutput> calib_;
  QString calibError_;
  QWidget* calTab_ = nullptr;
  PlotWidget *calAbsPlot_ = nullptr, *calEmPlot_ = nullptr, *calResPlot_ = nullptr;
  QTableWidget* calTable_ = nullptr;
  QLabel* calLabel_ = nullptr;
  void showCalibration();
  std::atomic<bool> cancel_{false};
  QString runError_;
  QProgressBar* progress_ = nullptr;
  QAction *runAct_ = nullptr, *stopAct_ = nullptr;
  QLabel* cursorLabel_ = nullptr;
};
