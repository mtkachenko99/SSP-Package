// ssp_gui — графическая оболочка ssp_aggregate (Qt6 Widgets).
// Опции: [config.json] [--results DIR] [--experiment FILE --exp-unit nm|cm|eV] [--tab N]
//        [--screenshot out.png [--size WxH]]   (снимок окна и выход; удобно с QT_QPA_PLATFORM=offscreen)

#include "mainwindow.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFileInfo>
#include <QLocale>
#include <QTimer>

#include <cstdio>

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  QApplication::setApplicationName("ssp_gui");
  QApplication::setOrganizationName("SSP Nano");
  // English-only interface: decimal point and number formats independent of the system locale
  QLocale::setDefault(QLocale(QLocale::English, QLocale::UnitedStates));

  QCommandLineParser p;
  p.addHelpOption();
  p.addPositionalArgument("config", "ssp_aggregate JSON config");
  QCommandLineOption results("results", "Open results folder", "dir");
  QCommandLineOption experiment("experiment", "Experimental spectrum", "file");
  QCommandLineOption expUnit("exp-unit", "Experiment x-axis units (nm|cm|eV)", "unit", "nm");
  QCommandLineOption tab("tab", "Result tab", "n", "0");
  QCommandLineOption paramTab("param-tab", "Parameter tab", "n", "0");
  QCommandLineOption fit("fit", "Fit when comparing: offset,shift", "list");
  QCommandLineOption run("run", "Start the calculation right after opening");
  QCommandLineOption calib("calibrate", "Run a monomer calibration (JSON with monomer_fit)", "file");
  QCommandLineOption shot("screenshot", "Save a window screenshot and exit", "png");
  QCommandLineOption size("size", "Window size for the screenshot", "WxH", "1560x960");
  p.addOptions({results, experiment, expUnit, tab, paramTab, fit, run, calib, shot, size});
  p.process(app);

  MainWindow w;
  const QStringList sz = p.value(size).split('x');
  w.resize(sz.value(0).toInt() > 0 ? sz.value(0).toInt() : 1560, sz.value(1).toInt() > 0 ? sz.value(1).toInt() : 960);
  // Пример по умолчанию: рядом с программой (дистрибутив: examples/), иначе — в исходниках (сборка разработчика)
  QString defCfg = QCoreApplication::applicationDirPath() + "/examples/notebook_demo.json";
#ifdef SSP_EXAMPLES_DIR
  if (!QFileInfo::exists(defCfg)) defCfg = QString(SSP_EXAMPLES_DIR) + "/shapes_with_gold.json";
#endif
  if (!p.positionalArguments().isEmpty()) w.openConfig(QFileInfo(p.positionalArguments().first()).absoluteFilePath());
  else if (!defCfg.isEmpty() && QFileInfo::exists(defCfg)) w.openConfig(defCfg);
  if (p.isSet(results)) w.loadResults(QFileInfo(p.value(results)).absoluteFilePath());
  if (p.isSet(experiment)) w.loadExperimentFile(p.value(experiment), p.value(expUnit));
  if (p.isSet(fit)) w.setFitOptions(p.value(fit).contains("offset"), p.value(fit).contains("shift"));
  w.selectResultTab(p.value(tab).toInt());
  w.show();
  const QString out = p.value(shot);
  const int pt = p.value(paramTab).toInt();
  const int rt = p.value(tab).toInt();
  auto snapshot = [&w, out, pt] {
      const auto all = w.findChildren<QTabWidget*>();
      for (auto* t : all)
        if (t->tabText(0) == "Shapes") t->setCurrentIndex(pt);
      QTimer::singleShot(300, [&w, out] {
        w.grab().save(out);
        QApplication::quit();
      });
  };
  if (p.isSet(calib)) {
    QObject::connect(&w, &MainWindow::calibrationFinished, &w, [snapshot, &p, &shot](bool ok) {
      std::fprintf(stderr, "calibration finished: %s\n", ok ? "successfully" : "with an error");
      if (p.isSet(shot)) QTimer::singleShot(200, snapshot);
    });
    const QString cf = QFileInfo(p.value(calib)).absoluteFilePath();
    QTimer::singleShot(200, &w, [&w, cf] { w.startCalibration(cf); });
  } else if (p.isSet(run)) {
    QObject::connect(&w, &MainWindow::runFinished, &w, [&w, snapshot, rt, &p, &shot](bool ok) {
      std::fprintf(stderr, "calculation finished: %s\n", ok ? "successfully" : "with an error");
      w.selectResultTab(rt);
      if (p.isSet(shot)) QTimer::singleShot(200, snapshot);
    });
    QTimer::singleShot(200, &w, [&w] { w.startRun(); });
  } else if (p.isSet(shot)) {
    QTimer::singleShot(400, snapshot);
  }
  return app.exec();
}
