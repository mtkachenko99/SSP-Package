#include "mainwindow.hpp"

#include "calibrate.hpp"

#include "formpanel.hpp"
#include "geometryview.hpp"
#include "plotwidget.hpp"

#include "config.hpp"
#include "runner.hpp"
#include "scan.hpp"
#include "ssp/excitons/geometry.hpp"

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDesktopServices>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QStatusBar>
#include <QTableWidget>
#include <QTabWidget>
#include <QToolBar>
#include <QUrl>
#include <QVBoxLayout>
#include <QtConcurrent>

#include <cmath>
#include <cstdio>

namespace {

using T = FieldSpec::Type;

FieldSpec F(const QString& path, const QString& label, T type, QVariant def, const QString& tip = {},
            const QString& section = {}) {
  FieldSpec f;
  f.path = path;
  f.label = label;
  f.type = type;
  f.def = def;
  f.tip = tip;
  f.section = section;
  return f;
}
FieldSpec D(const QString& path, const QString& label, double def, double mn, double mx, int dec, const QString& suf,
            const QString& tip = {}, const QString& section = {}) {
  FieldSpec f = F(path, label, T::Double, def, tip, section);
  f.min = mn;
  f.max = mx;
  f.decimals = dec;
  f.suffix = suf;
  return f;
}
FieldSpec I(const QString& path, const QString& label, int def, int mn, int mx, const QString& tip = {},
            const QString& section = {}) {
  FieldSpec f = F(path, label, T::Int, def, tip, section);
  f.min = mn;
  f.max = mx;
  return f;
}
FieldSpec C(const QString& path, const QString& label, const QStringList& vals, const QStringList& labels,
            const QString& def, const QString& tip = {}, const QString& section = {}) {
  FieldSpec f = F(path, label, T::Choice, def, tip, section);
  f.choices = vals;
  f.choiceLabels = labels;
  return f;
}

QList<FieldSpec> modelFields() {
  return {
      D("exciton.monomer_energy_cm", "Monomer energy", 17000, 1000, 60000, 1, "cm⁻¹",
        "Vertical transition energy of an isolated molecule", "Exciton model"),
      D("exciton.transition_dipole_D", "Transition dipole μ", 10, 0.1, 100, 2, "D"),
      D("exciton.disorder_sigma_cm", "Static disorder σ", 60, 0, 5000, 1, "cm⁻¹",
        "Gaussian spread of site energies"),
      D("exciton.n_medium", "Medium refractive index", 1.0, 1.0, 4.0, 3, "",
        "For the absolute absorption cross section and default screening (1/n²)"),
      D("exciton.disorder_corr_length_nm", "Disorder correlation length ξ (0 — none)", 0, 0, 1000, 3, "nm",
        "⟨δ_iδ_j⟩ = σ² exp(−r_ij/ξ)"),
      D("exciton.orientation_jitter_deg", "Orientation jitter (structural)", 0, 0, 90, 2, "°",
        "Random rotation of every dipole, new in each realization"),
      D("exciton.position_jitter_nm", "Position jitter (structural)", 0, 0, 5, 3, "nm"),
      F("exciton.dispersion.enabled", "Site-dependent dispersion shift (van der Waals)", T::Bool, false,
        "ΔE_m = Q·Σ_k V_mk², Q calibrated to the Bayliss continuum shift; site energy = monomer energy in its "
        "solvent + back-to-gas + ΔE_m (Megow et al. 2015)"),
      D("exciton.dispersion.F_eV", "Bayliss factor F", 1.43, 0, 10, 3, "eV"),
      D("exciton.dispersion.n_aggregate", "n of the dye layer", 1.78, 1, 3, 3, ""),
      D("exciton.dispersion.n_solvent", "n of the monomer solvent", 1.33, 1, 3, 3, ""),
      D("exciton.dispersion.dipole_length_nm", "Extended dipole length (dispersion)", 1.0, 0.1, 5, 3, "nm"),
      D("exciton.dispersion.softening_nm", "Charge softening", 0.3, 0, 2, 3, "nm",
        "1/√(r² + a²): transition densities are not point charges"),
      D("exciton.dispersion.cutoff_nm", "Summation radius", 3.0, 0.5, 20, 2, "nm"),
      D("exciton.dispersion.extra_shift_cm", "Extra site shift (solvent contact)", 0, -1e4, 1e4, 1, "cm⁻¹"),
      C("exciton.relaxation", "Relaxation theory", {"redfield", "modified_redfield"},
        {"secular Redfield", "modified Redfield (N ≲ 300)"}, "redfield",
        "Modified Redfield treats the diagonal exciton–phonon coupling non-perturbatively (better for λ ~ J)"),
      C("exciton.coupling.model", "Coupling model", {"point_dipole", "extended_dipole", "legacy_J0", "transition_charges"},
        {"point dipole", "extended dipole", "as in the notebook (J0·(r0/r)³)", "transition charges (TrEsp)"},
        "point_dipole", {}, "Coupling J_ij"),
      F("exciton.coupling.charges_file", "Transition charges file (TrEsp)", T::File, "",
        "Columns x_nm y_nm z_nm q in the molecular frame (x along the transition dipole); rescaled to μ"),
      [&] { auto f = D("exciton.coupling.screening", "Screening f_s", 1.0, 0.0, 2.0, 4, "", "Multiplier of J; usually 1/n²"); f.autoText = "auto (1/n²)"; return f; }(),
      D("exciton.coupling.dipole_length_nm", "Extended dipole length", 0.8, 0.01, 5, 3, "nm"),
      D("exciton.coupling.J0_cm", "J0 (legacy)", 900, -1e5, 1e5, 1, "cm⁻¹"),
      D("exciton.coupling.r0_nm", "r0 (legacy)", 1.0, 0.01, 100, 3, "nm"),
      C("bath.type", "Spectral density", {"drude_lorentz", "ohmic_exponential"}, {"Drude–Lorentz", "ohmic (exp)"},
        "drude_lorentz", {}, "Bath (line shapes, relaxation)"),
      D("bath.lambda_cm", "Reorganization energy λ", 50, 0, 5000, 2, "cm⁻¹"),
      D("bath.cutoff_cm", "Cutoff frequency ω_c", 100, 1, 10000, 1, "cm⁻¹"),
      D("bath.temperature_K", "Temperature", 300, 0, 2000, 1, "K"),
      D("decay.monomer_lifetime_ns", "Monomer τ (0 — from formula)", 0, 0, 1e4, 3, "ns", {}, "Decay and pump"),
      D("decay.nonradiative_ps", "Γ_nr", 0, 0, 1e3, 6, "ps⁻¹"),
      C("pump.type", "Pump", {"oscillator_strength", "laser"}, {"∝ oscillator strength (as in the notebook)", "laser at ν̃_exc"},
        "oscillator_strength"),
      [&] { auto f = D("pump.excitation_cm", "Excitation ν̃", 17500, 1000, 60000, 1, "cm⁻¹"); f.autoText = "auto (= monomer)"; return f; }(),
      F("holstein.enabled", "Vibronic absorption (Holstein, one-particle approx.)", T::Bool, false,
        "Extra spectrum vibronic_absorption.csv; do not put the same mode into bath modes", "Vibronic model"),
      D("holstein.omega_cm", "Mode ω₀", 1300, 10, 5000, 1, "cm⁻¹"),
      D("holstein.huang_rhys", "Huang–Rhys S", 0.6, 0, 5, 3, ""),
      I("holstein.nu_max", "Max. quanta ν", 4, 0, 12),
      I("holstein.realizations", "Realizations", 20, 1, 10000),
  };
}

QList<FieldSpec> calcFields() {
  return {
      F("output_dir", "Results folder", T::String, "results/run", "Relative to the config file", "Output"),
      D("spectrum.wmin_cm", "ν̃ min", 15500, 100, 1e5, 1, "cm⁻¹", {}, "Spectra"),
      D("spectrum.wmax_cm", "ν̃ max", 18000, 100, 1e5, 1, "cm⁻¹"),
      I("spectrum.nw", "Points", 1800, 10, 200000),
      D("spectrum.t_max_ps", "Line-shape time window", 100, 1, 1e5, 1, "ps"),
      I("ensemble.min_realizations", "Min. realizations", 20, 1, 1000000, {}, "Disorder ensemble"),
      I("ensemble.max_realizations", "Max. realizations", 400, 1, 1000000),
      D("ensemble.rel_tolerance", "Tolerance SE/peak", 0.01, 1e-5, 1, 4, ""),
      I("ensemble.seed", "RNG seed", 1234, 0, 2000000000),
      F("dynamics.enabled", "Population kinetics", T::Bool, true, {}, "Kinetics"),
      D("dynamics.t_max_ps", "t max", 20, 0.001, 1e6, 3, "ps"),
      I("dynamics.n_times", "Time points", 200, 2, 100000),
      C("dynamics.initial", "Initial state", {"brightest", "highest", "pump"},
        {"brightest", "highest", "from pump"}, "brightest"),
  };
}

QList<FieldSpec> npFields() {
  return {
      F("nanoparticle.enabled", "Include nanoparticle", T::Bool, false, {}, "Nanoparticle (centred at the origin)"),
      C("nanoparticle.shape", "Shape", {"sphere", "prolate", "oblate"}, {"sphere", "rod (prolate spheroid)",
                                                                       "disk (oblate spheroid)"}, "sphere"),
      D("nanoparticle.semi_axis_xy_nm", "Semi-axis in xy plane", 15, 0.1, 1000, 2, "nm",
        "Sphere: the radius. Rod (prolate): the smaller axis. Disk (oblate): the larger axis."),
      D("nanoparticle.semi_axis_z_nm", "Semi-axis along z", 15, 0.1, 1000, 2, "nm",
        "Symmetry axis. Rod: z > xy. Disk: z < xy. Ignored for a sphere."),
      F("nanoparticle.material_file", "Material (refractiveindex.info YAML)", T::File, "", "Path relative to the config"),
      D("nanoparticle.n_medium", "Medium n", 1.33, 1, 4, 3, ""),
      I("nanoparticle.l_max", "Multipole L_max", 40, 1, 200, "Checked by the truncation estimate in the log"),
      F("nanoparticle.dynamic_dipole", "Retardation (sphere: Mie for l ≤ 8; spheroid: MLWA for l = 1)", T::Bool, true),
      D("nanoparticle.shell_thickness_nm", "Ligand/spacer shell thickness (sphere)", 0, 0, 100, 3, "nm",
        "Dielectric layer on the metal; particle_shell molecules sit above its outer surface"),
      D("nanoparticle.shell_n", "Shell refractive index", 1.45, 1, 4, 3, ""),
      D("nanoparticle.d_perp_nm", "Feibelman d⊥, Re (sphere)", 0, -2, 2, 3, "nm",
        "> 0: electron spill-out (red shift of the plasmon), < 0: spill-in"),
      D("nanoparticle.d_perp_im_nm", "Feibelman d⊥, Im", 0, -2, 2, 3, "nm", "Imaginary part: surface-assisted (Landau) damping"),
      [&] {
        auto f = C("nanoparticle.fast_l", "Solver: exact multipoles up to l", {"auto", "0", "8", "12", "16", "24"}, {},
                   "auto",
                   "auto — fast solver for N ≥ 600; 0 — exact N×N; L > 0 — multipoles l ≤ L exact (Woodbury), "
                   "higher ones diagonal in the exciton basis");
        f.editable = true;
        return f;
      }(),
      D("nanoparticle.size_sigma_rel", "Core size spread σ_R/R (0 — none)", 0, 0, 0.5, 3, "",
        "Spectra averaged over a normal distribution of core radii (Gauss–Hermite nodes)"),
      I("nanoparticle.size_nodes", "Size nodes", 5, 1, 15),
      F("nanoparticle.emission", "Hybrid PL, quantum yield, Purcell", T::Bool, true),
      D("nanoparticle.kreibig_A", "Kreibig A (0 — off)", 0, 0, 5, 3, ""),
      I("nanoparticle.realizations", "Realizations for spectra", 30, 1, 10000),
      [&] { auto f = D("nanoparticle.wmin_cm", "ν̃ min", 14000, 100, 1e5, 1, "cm⁻¹"); f.autoText = "auto (= spectrum)"; return f; }(),
      [&] { auto f = D("nanoparticle.wmax_cm", "ν̃ max", 22000, 100, 1e5, 1, "cm⁻¹"); f.autoText = "auto (= spectrum)"; return f; }(),
      I("nanoparticle.nw", "Points", 600, 10, 100000),
      D("nanoparticle.homogeneous_hwhm_cm", "Homogeneous HWHM (<0 — Redfield)", -1, -1, 1e4, 2, "cm⁻¹"),
      [&] { auto f = D("nanoparticle.near_field_cm", "ν̃ for near-field map", 17000, 100, 1e5, 1, "cm⁻¹"); f.autoText = "auto (= monomer)"; return f; }(),
  };
}

QList<FieldSpec> heomFields() {
  return {
      F("heom.enabled", "Compute HEOM reference", T::Bool, false, "Drude–Lorentz bath without modes only",
        "HEOM (numerically exact reference)"),
      I("heom.max_sites", "Max. molecules", 10, 1, 40, "HEOM is skipped for larger aggregates"),
      I("heom.matsubara", "Matsubara terms K", 1, 0, 10),
      I("heom.depth", "Hierarchy depth L", 4, 1, 20),
      F("heom.terminator", "Terminator", T::Bool, true),
      D("heom.t_max_ps", "t max", 1.0, 0.01, 100, 3, "ps"),
      D("heom.dt_fs", "Step (0 — auto)", 0, 0, 100, 3, "fs"),
      F("heom.absorption", "Absorption", T::Bool, true),
      F("heom.populations", "Populations", T::Bool, true),
      F("heom.initial", "Initial state", T::String, "brightest", "brightest | highest | site:<k>"),
  };
}

QList<FieldSpec> scanFields() {
  FieldSpec par = C("scan.parameter", "Parameter (JSON path)",
                    {"exciton.monomer_energy_cm", "exciton.transition_dipole_D", "exciton.disorder_sigma_cm",
                     "bath.temperature_K", "bath.lambda_cm", "nanoparticle.semi_axis_xy_nm", "nanoparticle.semi_axis_z_nm",
                     "nanoparticle.n_medium", "geometry_params.N", "geometry_params.gap", "geometry_params.theta_max_deg"},
                    {}, "exciton.monomer_energy_cm",
                    "Any numeric key of the config. geometry_params.X — all shapes; configurations[name].geometry_params.X — "
                    "one shape. For a sphere scan semi_axis_xy_nm together with semi_axis_z_nm (also set).");
  par.editable = true;
  return {
      F("scan.enabled", "Run a parameter scan", T::Bool, false,
        "Each value runs the full calculation into output_dir/scan_NN; maps, branches and the coupled-oscillator fit "
        "are written to output_dir/<shape>/",
        "Parameter scan"),
      par,
      D("scan.from", "From", 16500, -1e9, 1e9, 4, ""),
      D("scan.to", "To", 21500, -1e9, 1e9, 4, ""),
      I("scan.steps", "Steps", 11, 2, 1000),
      F("scan.also_set", "Also set (comma-separated paths)", T::String, "",
        "Other paths that receive the same value, e.g. nanoparticle.semi_axis_z_nm"),
  };
}

QJsonObject defaultConfig() {
  QJsonObject c;
  c["output_dir"] = "results/gui_run";
  QJsonObject e{{"name", "ring"}, {"geometry_type", "ring"},
                {"geometry_params", QJsonObject{{"N", 16}, {"radius", 3.5}, {"tilt_deg", 0.0}}}};
  c["configurations"] = QJsonArray{e};
  return c;
}

QString fmt(double v, int p = 4) { return QString::number(v, 'g', p); }

}  // namespace

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
  setWindowTitle("SSP Nano — aggregates and plasmonic nanoparticles");
  buildUi();
  cfg_ = defaultConfig();
  baseDir_ = QDir::currentPath();
  syncFormsFromConfig();
  connect(&watcher_, &QFutureWatcher<void>::finished, this, [this] {
    setRunning(false);
    if (!runError_.isEmpty()) {
      appendLog("ERROR: " + runError_);
      if (QGuiApplication::platformName() == "offscreen") std::fprintf(stderr, "ERROR: %s\n", runError_.toUtf8().constData());
      emit runFinished(false);
      if (QGuiApplication::platformName() != "offscreen") QMessageBox::critical(this, "Calculation error", runError_);
      return;
    }
    loadResults(resolvedOutputDir());
    emit runFinished(true);
  });
  connect(&calWatcher_, &QFutureWatcher<void>::finished, this, [this] {
    runAct_->setEnabled(true);
    progress_->setVisible(false);
    if (!calibError_.isEmpty()) {
      appendLog("ERROR: " + calibError_);
      if (QGuiApplication::platformName() == "offscreen") std::fprintf(stderr, "ERROR: %s\n", calibError_.toUtf8().constData());
      else QMessageBox::critical(this, "Calibration error", calibError_);
      emit calibrationFinished(false);
      return;
    }
    showCalibration();
    emit calibrationFinished(true);
  });
}

// ---------------------------------------------------------------- UI

void MainWindow::buildUi() {
  auto* fileMenu = menuBar()->addMenu("File");
  auto* tb = addToolBar("Main");
  tb->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  auto add = [&](const QString& text, auto slot, QKeySequence key = {}, bool toolbar = true) {
    auto* a = new QAction(text, this);
    a->setShortcut(key);
    connect(a, &QAction::triggered, this, slot);
    fileMenu->addAction(a);
    if (toolbar) tb->addAction(a);
    return a;
  };
  add("New", &MainWindow::newConfig, QKeySequence::New);
  add("Open…", &MainWindow::openConfigDialog, QKeySequence::Open);
  add("Save", &MainWindow::saveConfig, QKeySequence::Save);
  add("Save as…", &MainWindow::saveConfigAs, QKeySequence::SaveAs, false);
  fileMenu->addSeparator();
  tb->addSeparator();
  runAct_ = add("▶ Run", &MainWindow::runCalculation, QKeySequence(Qt::Key_F5));
  stopAct_ = add("■ Stop", &MainWindow::stopCalculation, {});
  stopAct_->setEnabled(false);
  tb->addSeparator();
  add("Open results…", &MainWindow::openResultsDialog);
  add("Experiment…", &MainWindow::loadExperimentDialog);
  add("Calibrate monomer…", &MainWindow::calibrateDialog);
  fileMenu->addSeparator();
  auto* q = fileMenu->addAction("Quit");
  connect(q, &QAction::triggered, this, &QWidget::close);

  progress_ = new QProgressBar;
  progress_->setMaximumWidth(260);
  progress_->setVisible(false);
  cursorLabel_ = new QLabel;
  statusBar()->addPermanentWidget(cursorLabel_);
  statusBar()->addPermanentWidget(progress_);

  // Левая часть — параметры
  paramTabs_ = new QTabWidget;
  paramTabs_->addTab(buildConfigurationsTab(), "Shapes");
  auto wrap = [](QWidget* w) {
    auto* s = new QScrollArea;
    s->setWidget(w);
    s->setWidgetResizable(true);
    return s;
  };
  modelForm_ = new FormPanel(modelFields());
  calcForm_ = new FormPanel(calcFields());
  npForm_ = new FormPanel(npFields());
  heomForm_ = new FormPanel(heomFields());
  scanForm_ = new FormPanel(scanFields());
  paramTabs_->addTab(wrap(modelForm_), "Model");
  paramTabs_->addTab(wrap(calcForm_), "Calculation");
  paramTabs_->addTab(wrap(npForm_), "Nanoparticle");
  paramTabs_->addTab(wrap(heomForm_), "HEOM");
  paramTabs_->addTab(wrap(scanForm_), "Scan");
  for (auto* f : {modelForm_, calcForm_, npForm_, heomForm_, scanForm_})
    connect(f, &FormPanel::changed, this, [this] {
      if (syncing_) return;
      syncConfigFromForms();
      updateGeometryPreview();
    });

  auto* jsonBox = new QWidget;
  auto* jl = new QVBoxLayout(jsonBox);
  jsonEdit_ = new QPlainTextEdit;
  QFont mono("monospace");
  mono.setStyleHint(QFont::Monospace);
  jsonEdit_->setFont(mono);
  auto* applyJson = new QPushButton("Apply JSON → forms");
  jl->addWidget(new QLabel("Full config (the same file ssp_aggregate reads). Unknown keys are errors."));
  jl->addWidget(jsonEdit_);
  jl->addWidget(applyJson);
  connect(applyJson, &QPushButton::clicked, this, &MainWindow::applyJsonText);
  paramTabs_->addTab(jsonBox, "JSON");
  connect(paramTabs_, &QTabWidget::currentChanged, this, [this](int i) {
    if (paramTabs_->tabText(i) == "JSON") refreshJsonView();
  });

  auto* split = new QSplitter;
  split->addWidget(paramTabs_);
  split->addWidget(buildResultsArea());
  split->setStretchFactor(0, 0);
  split->setStretchFactor(1, 1);
  split->setSizes({470, 1100});
  setCentralWidget(split);
}

QWidget* MainWindow::buildConfigurationsTab() {
  auto* w = new QWidget;
  auto* v = new QVBoxLayout(w);
  auto* top = new QHBoxLayout;
  entryList_ = new QListWidget;
  entryList_->setMaximumHeight(130);
  top->addWidget(entryList_);
  auto* btns = new QVBoxLayout;
  auto mk = [&](const QString& t) {
    auto* b = new QPushButton(t);
    btns->addWidget(b);
    return b;
  };
  auto* bAdd = mk("Add");
  auto* bDup = mk("Duplicate");
  auto* bDel = mk("Delete");
  btns->addStretch();
  top->addLayout(btns);
  v->addLayout(top);

  auto* form = new QFormLayout;
  entryName_ = new QLineEdit;
  entryType_ = new QComboBox;
  for (const auto& t : ssp::excitons::geometry_types()) entryType_->addItem(QString::fromStdString(t));
  form->addRow("Name", entryName_);
  form->addRow("Shape type", entryType_);
  v->addLayout(form);
  paramTable_ = new QTableWidget(0, 2);
  paramTable_->setHorizontalHeaderLabels({"Parameter", "Value"});
  paramTable_->horizontalHeader()->setStretchLastSection(true);
  paramTable_->verticalHeader()->setVisible(false);
  paramTable_->setMinimumHeight(150);
  v->addWidget(paramTable_);

  auto* pl = new QGroupBox("Placement (relative to the NP at the origin)");
  auto* g = new QFormLayout(pl);
  centerCheck_ = new QCheckBox("first move the aggregate centre to 0");
  g->addRow(centerCheck_);
  auto row = [&](QDoubleSpinBox** arr, const QString& label, double lo, double hi, const QString& suf) {
    auto* h = new QHBoxLayout;
    for (int k = 0; k < 3; ++k) {
      arr[k] = new QDoubleSpinBox;
      arr[k]->setRange(lo, hi);
      arr[k]->setDecimals(2);
      arr[k]->setSuffix(suf);
      arr[k]->setKeyboardTracking(false);
      h->addWidget(arr[k]);
      connect(arr[k], &QDoubleSpinBox::valueChanged, this, [this] {
        storeEntryEditor();
        updateGeometryPreview();
      });
    }
    g->addRow(label, h);
  };
  row(shift_, "Shift x, y, z", -1000, 1000, " nm");
  row(euler_, "Euler α, β, γ", -360, 360, "°");
  v->addWidget(pl);
  geomView_ = new GeometryView;
  geomView_->setMinimumHeight(260);
  v->addWidget(geomView_, 1);

  connect(entryList_, &QListWidget::currentRowChanged, this, [this](int r) {
    if (syncing_) return;
    storeEntryEditor();
    loadEntryEditor(r);
  });
  connect(entryName_, &QLineEdit::editingFinished, this, [this] {
    storeEntryEditor();
    if (auto* it = entryList_->item(currentEntry_)) it->setText(entryName_->text());
  });
  connect(entryType_, &QComboBox::currentTextChanged, this, [this](const QString& t) {
    if (syncing_) return;
    fillParamTable(t, {});
    storeEntryEditor();
    updateGeometryPreview();
  });
  connect(paramTable_, &QTableWidget::cellChanged, this, [this] {
    if (syncing_) return;
    storeEntryEditor();
    updateGeometryPreview();
  });
  connect(centerCheck_, &QCheckBox::toggled, this, [this] {
    storeEntryEditor();
    updateGeometryPreview();
  });
  connect(bAdd, &QPushButton::clicked, this, [this] {
    storeEntryEditor();
    QJsonArray a = cfg_["configurations"].toArray();
    a.append(QJsonObject{{"name", QString("shape%1").arg(a.size() + 1)}, {"geometry_type", "linear_chain"},
                         {"geometry_params", QJsonObject{}}});
    cfg_["configurations"] = a;
    refreshConfigList(a.size() - 1);
  });
  connect(bDup, &QPushButton::clicked, this, [this] {
    storeEntryEditor();
    QJsonArray a = cfg_["configurations"].toArray();
    if (currentEntry_ < 0 || currentEntry_ >= a.size()) return;
    QJsonObject o = a[currentEntry_].toObject();
    o["name"] = o["name"].toString() + "_copy";
    a.insert(currentEntry_ + 1, o);
    cfg_["configurations"] = a;
    refreshConfigList(currentEntry_ + 1);
  });
  connect(bDel, &QPushButton::clicked, this, [this] {
    QJsonArray a = cfg_["configurations"].toArray();
    if (currentEntry_ < 0 || currentEntry_ >= a.size() || a.size() <= 1) return;
    a.removeAt(currentEntry_);
    cfg_["configurations"] = a;
    currentEntry_ = -1;
    refreshConfigList(0);
  });
  return w;
}

QWidget* MainWindow::buildResultsArea() {
  resultTabs_ = new QTabWidget;
  auto cursor = [this](PlotWidget* p, bool nm) {
    connect(p, &PlotWidget::cursorMoved, this, [this, nm](double x, double y) {
      cursorLabel_->setText(nm && x > 0 ? QString("ν̃ = %1 cm⁻¹ (λ = %2 nm), y = %3").arg(x, 0, 'f', 1).arg(1e7 / x, 0, 'f', 1).arg(y, 0, 'g', 4)
                                         : QString("x = %1, y = %2").arg(x, 0, 'g', 5).arg(y, 0, 'g', 4));
    });
  };

  // --- Спектры + сравнение с экспериментом
  {
    auto* w = new QWidget;
    auto* h = new QHBoxLayout(w);
    auto* side = new QVBoxLayout;
    specConfigs_ = new QListWidget;
    specConfigs_->setMaximumWidth(220);
    side->addWidget(new QLabel("Shapes:"));
    side->addWidget(specConfigs_, 1);
    specMode_ = new QComboBox;
    specMode_->addItems({"Absorption", "Fluorescence", "Absorption + fluorescence"});
    side->addWidget(specMode_);
    specPerMolecule_ = new QCheckBox("Per molecule (σ/N)");
    specPerMolecule_->setChecked(true);
    specPerMolecule_->setToolTip("Shows how aggregation redistributes oscillator strength (superradiance); "
                                 "normalization to the maximum hides it");
    specNormalize_ = new QCheckBox("Normalize to maximum");
    specErrors_ = new QCheckBox("±SE bands (disorder sampling only)");
    specErrors_->setChecked(true);
    specSticks_ = new QCheckBox("Exciton states (one realization)");
    specSticks_->setChecked(true);
    side->addWidget(specPerMolecule_);
    side->addWidget(specNormalize_);
    side->addWidget(specErrors_);
    side->addWidget(specSticks_);
    auto* cmp = new QGroupBox("Comparison with experiment");
    auto* cl = new QVBoxLayout(cmp);
    compareTarget_ = new QComboBox;
    cl->addWidget(new QLabel("Model to compare:"));
    cl->addWidget(compareTarget_);
    fitOffset_ = new QCheckBox("fit baseline");
    fitShift_ = new QCheckBox("fit shift δ");
    residNormalized_ = new QCheckBox("residuals in units of σ");
    residNormalized_->setChecked(true);
    cl->addWidget(fitOffset_);
    cl->addWidget(fitShift_);
    cl->addWidget(residNormalized_);
    auto* bLoad = new QPushButton("Load…");
    auto* bClear = new QPushButton("Remove");
    auto* bh = new QHBoxLayout;
    bh->addWidget(bLoad);
    bh->addWidget(bClear);
    cl->addLayout(bh);
    compareLabel_ = new QLabel("no experiment loaded");
    compareLabel_->setWordWrap(true);
    compareLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    cl->addWidget(compareLabel_);
    side->addWidget(cmp);
    h->addLayout(side);

    auto* plots = new QSplitter(Qt::Vertical);
    specPlot_ = new PlotWidget;
    specPlot_->setSecondaryNmAxis(true);
    residPlot_ = new PlotWidget;
    residPlot_->setZeroLine(true);
    residPlot_->setLegendVisible(false);
    residPlot_->setVisible(false);
    plots->addWidget(specPlot_);
    plots->addWidget(residPlot_);
    plots->setSizes({600, 220});
    h->addWidget(plots, 1);
    connect(specPlot_, &PlotWidget::xRangeChanged, residPlot_, [this](double a, double b) { residPlot_->setXRange(a, b, false); });
    connect(residPlot_, &PlotWidget::xRangeChanged, specPlot_, [this](double a, double b) { specPlot_->setXRange(a, b, false); });
    cursor(specPlot_, true);
    cursor(residPlot_, true);
    connect(specConfigs_, &QListWidget::itemChanged, this, [this] { plotSpectra(); });
    for (auto* c : {specNormalize_, specErrors_, fitOffset_, fitShift_, residNormalized_, specPerMolecule_, specSticks_})
      connect(c, &QCheckBox::toggled, this, [this] { plotSpectra(); });
    connect(specMode_, &QComboBox::currentIndexChanged, this, [this] { plotSpectra(); });
    connect(compareTarget_, &QComboBox::currentIndexChanged, this, [this] { plotSpectra(); });
    connect(bLoad, &QPushButton::clicked, this, &MainWindow::loadExperimentDialog);
    connect(bClear, &QPushButton::clicked, this, &MainWindow::clearExperiment);
    resultTabs_->addTab(w, "Spectra");
  }
  auto withCombo = [&](QComboBox*& combo, QWidget* content, QWidget* extra = nullptr) {
    auto* w = new QWidget;
    auto* v = new QVBoxLayout(w);
    auto* h = new QHBoxLayout;
    combo = new QComboBox;
    combo->setMinimumWidth(200);
    h->addWidget(new QLabel("Shape:"));
    h->addWidget(combo);
    if (extra) h->addWidget(extra);
    h->addStretch();
    v->addLayout(h);
    v->addWidget(content, 1);
    return w;
  };
  // --- Плекситон
  {
    plexPlot_ = new PlotWidget;
    plexPlot_->setSecondaryNmAxis(true);
    cursor(plexPlot_, true);
    plexView_ = new QComboBox;
    plexView_->addItems({"aggregate contribution (system − NP)", "cross sections (system and NP)",
                         "PL: hybrid vs aggregate without NP"});
    resultTabs_->addTab(withCombo(plexConfig_, plexPlot_, plexView_), "Aggregate + NP");
    connect(plexConfig_, &QComboBox::currentIndexChanged, this, [this] { plotPlexciton(); });
    connect(plexView_, &QComboBox::currentIndexChanged, this, [this] { plotPlexciton(); });
  }
  // --- Кинетика
  {
    auto* s = new QSplitter(Qt::Vertical);
    dynPlot_ = new PlotWidget;
    dynEnergyPlot_ = new PlotWidget;
    s->addWidget(dynPlot_);
    s->addWidget(dynEnergyPlot_);
    cursor(dynPlot_, false);
    cursor(dynEnergyPlot_, false);
    connect(dynPlot_, &PlotWidget::xRangeChanged, dynEnergyPlot_, [this](double a, double b) { dynEnergyPlot_->setXRange(a, b, false); });
    connect(dynEnergyPlot_, &PlotWidget::xRangeChanged, dynPlot_, [this](double a, double b) { dynPlot_->setXRange(a, b, false); });
    resultTabs_->addTab(withCombo(dynConfig_, s), "Kinetics");
    connect(dynConfig_, &QComboBox::currentIndexChanged, this, [this] { plotDynamics(); });
  }
  // --- HEOM
  {
    auto* s = new QSplitter(Qt::Horizontal);
    auto* left = new QSplitter(Qt::Vertical);
    heomAbsPlot_ = new PlotWidget;
    heomAbsPlot_->setSecondaryNmAxis(true);
    heomResidPlot_ = new PlotWidget;
    heomResidPlot_->setZeroLine(true);
    left->addWidget(heomAbsPlot_);
    left->addWidget(heomResidPlot_);
    left->setSizes({500, 200});
    heomPopPlot_ = new PlotWidget;
    s->addWidget(left);
    s->addWidget(heomPopPlot_);
    connect(heomAbsPlot_, &PlotWidget::xRangeChanged, heomResidPlot_, [this](double a, double b) { heomResidPlot_->setXRange(a, b, false); });
    connect(heomResidPlot_, &PlotWidget::xRangeChanged, heomAbsPlot_, [this](double a, double b) { heomAbsPlot_->setXRange(a, b, false); });
    cursor(heomAbsPlot_, true);
    cursor(heomPopPlot_, false);
    heomLabel_ = new QLabel;
    resultTabs_->addTab(withCombo(heomConfig_, s, heomLabel_), "HEOM reference");
    connect(heomConfig_, &QComboBox::currentIndexChanged, this, [this] { plotHeom(); });
  }
  // --- Ближнее поле
  {
    auto* s = new QSplitter(Qt::Horizontal);
    auto* left = new QSplitter(Qt::Vertical);
    nfEnhPlot_ = new PlotWidget;
    nfDecayPlot_ = new PlotWidget;
    left->addWidget(nfEnhPlot_);
    left->addWidget(nfDecayPlot_);
    nfGeom_ = new GeometryView;
    s->addWidget(left);
    s->addWidget(nfGeom_);
    cursor(nfEnhPlot_, false);
    cursor(nfDecayPlot_, false);
    nfQuantity_ = new QComboBox;
    nfQuantity_->addItems({"quenching rate −2 Im W_ii", "radiative enhancement Γ_rad/Γ₀ (Purcell)",
                           "monomer quantum yield"});
    resultTabs_->addTab(withCombo(nfConfig_, s, nfQuantity_), "Near field");
    connect(nfQuantity_, &QComboBox::currentIndexChanged, this, [this] { plotNearField(); });
    connect(nfConfig_, &QComboBox::currentIndexChanged, this, [this] { plotNearField(); });
  }
  // --- Скан
  {
    scanTab_ = new QWidget;
    auto* v = new QVBoxLayout(scanTab_);
    auto* h = new QHBoxLayout;
    scanShape_ = new QComboBox;
    scanMap_ = new QComboBox;
    scanMap_->addItems({"System absorption", "System extinction", "System scattering", "System − NP absorption",
                        "Aggregate absorption", "Aggregate fluorescence"});
    scanNormRows_ = new QCheckBox("normalize each row");
    scanNormRows_->setChecked(true);
    scanOverlay_ = new QCheckBox("peaks and fitted branches");
    scanOverlay_->setChecked(true);
    h->addWidget(new QLabel("Shape:"));
    h->addWidget(scanShape_);
    h->addWidget(new QLabel("Map:"));
    h->addWidget(scanMap_);
    h->addWidget(scanNormRows_);
    h->addWidget(scanOverlay_);
    h->addStretch();
    v->addLayout(h);
    auto* split = new QSplitter(Qt::Vertical);
    scanPlot_ = new PlotWidget;
    scanPlot_->setSecondaryNmAxis(true);
    cursor(scanPlot_, true);
    split->addWidget(scanPlot_);
    auto* bottom = new QWidget;
    auto* bl = new QHBoxLayout(bottom);
    scanCoupling_ = new QLabel;
    scanCoupling_->setWordWrap(true);
    scanCoupling_->setTextFormat(Qt::RichText);
    scanCoupling_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    scanCoupling_->setMinimumWidth(420);
    bl->addWidget(scanCoupling_, 1);
    auto* ob = new QVBoxLayout;
    scanObs_ = new QComboBox;
    ob->addWidget(scanObs_);
    scanObsPlot_ = new PlotWidget;
    ob->addWidget(scanObsPlot_, 1);
    bl->addLayout(ob, 1);
    split->addWidget(bottom);
    split->setSizes({560, 260});
    v->addWidget(split, 1);
    scanMap_->setCurrentIndex(3);  // по умолчанию — вклад агрегата (поглощение НЧ его заслоняет)
    resultTabs_->addTab(scanTab_, "Scan");
    connect(scanShape_, &QComboBox::currentIndexChanged, this, [this] {
      plotScan();
      plotScanObservable();
    });
    connect(scanMap_, &QComboBox::currentIndexChanged, this, [this] { plotScan(); });
    connect(scanNormRows_, &QCheckBox::toggled, this, [this] { plotScan(); });
    connect(scanOverlay_, &QCheckBox::toggled, this, [this] { plotScan(); });
    connect(scanObs_, &QComboBox::currentIndexChanged, this, [this] { plotScanObservable(); });
  }
  // --- Калибровка по мономеру
  {
    calTab_ = new QWidget;
    auto* h = new QHBoxLayout(calTab_);
    auto* left = new QSplitter(Qt::Vertical);
    calAbsPlot_ = new PlotWidget;
    calAbsPlot_->setSecondaryNmAxis(true);
    calEmPlot_ = new PlotWidget;
    calResPlot_ = new PlotWidget;
    calResPlot_->setZeroLine(true);
    left->addWidget(calAbsPlot_);
    left->addWidget(calEmPlot_);
    left->addWidget(calResPlot_);
    left->setSizes({400, 300, 200});
    auto* right = new QWidget;
    auto* v = new QVBoxLayout(right);
    calLabel_ = new QLabel("Calibrate the model on the isolated dye: File → Calibrate monomer… (JSON with a "
                           "\"monomer_fit\" section; example: examples/monomer_calibration.json).");
    calLabel_->setWordWrap(true);
    calTable_ = new QTableWidget;
    calTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    auto* apply = new QPushButton("Apply to the aggregate model");
    apply->setToolTip("Sets monomer energy, transition dipole (if absolute ε), refractive index and the bath "
                      "(Drude + modes) in the current configuration");
    connect(apply, &QPushButton::clicked, this, &MainWindow::applyCalibration);
    v->addWidget(calLabel_);
    v->addWidget(calTable_, 1);
    v->addWidget(apply);
    calTable_->horizontalHeader()->setStretchLastSection(true);
    right->setMaximumWidth(470);
    h->addWidget(left, 1);
    h->addWidget(right);
    cursor(calAbsPlot_, true);
    resultTabs_->addTab(calTab_, "Calibration");
  }
  // --- Сводка
  summaryTable_ = new QTableWidget;
  summaryTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  resultTabs_->addTab(summaryTable_, "Shape summary");
  // --- Журнал
  log_ = new QPlainTextEdit;
  log_->setReadOnly(true);
  QFont mono("monospace");
  mono.setStyleHint(QFont::Monospace);
  log_->setFont(mono);
  resultTabs_->addTab(log_, "Log");
  return resultTabs_;
}

// ---------------------------------------------------------------- конфиг ↔ формы

void MainWindow::normalizeConfig() {
  // radius_nm (сокращение для сферы) → явные полуоси, которые редактирует форма
  QJsonObject np = cfg_["nanoparticle"].toObject();
  if (np.contains("radius_nm")) {
    const double r = np["radius_nm"].toDouble();
    if (!np.contains("semi_axis_xy_nm")) np["semi_axis_xy_nm"] = r;
    if (!np.contains("semi_axis_z_nm")) np["semi_axis_z_nm"] = r;
    np.remove("radius_nm");
    cfg_["nanoparticle"] = np;
  }
}

void MainWindow::syncFormsFromConfig() {
  normalizeConfig();
  syncing_ = true;
  for (auto* f : {modelForm_, calcForm_, npForm_, heomForm_, scanForm_}) {
    f->setBaseDir(baseDir_);
    f->load(cfg_);
    f->store(cfg_);  // конфиг дополняется значениями по умолчанию — то, что видно в форме, и будет посчитано
  }
  syncing_ = false;
  currentEntry_ = -1;
  refreshConfigList(0);
  refreshJsonView();
}

void MainWindow::syncConfigFromForms() {
  storeEntryEditor();
  for (auto* f : {modelForm_, calcForm_, npForm_, heomForm_, scanForm_}) f->store(cfg_);
  // Значения для неактивных разделов допустимы, но для сферы z-полуось = xy
  if (jsonGet(cfg_, "nanoparticle.shape").toString() == "sphere")
    jsonSet(cfg_, "nanoparticle.semi_axis_z_nm", jsonGet(cfg_, "nanoparticle.semi_axis_xy_nm"));
}

void MainWindow::refreshJsonView() {
  if (!jsonEdit_) return;
  syncConfigFromForms();
  jsonEdit_->setPlainText(QString::fromUtf8(QJsonDocument(cfg_).toJson(QJsonDocument::Indented)));
}

void MainWindow::applyJsonText() {
  QJsonParseError err;
  const QJsonDocument d = QJsonDocument::fromJson(jsonEdit_->toPlainText().toUtf8(), &err);
  if (err.error != QJsonParseError::NoError || !d.isObject()) {
    QMessageBox::warning(this, "JSON", "Parse error: " + err.errorString() + QString(" (offset %1)").arg(err.offset));
    return;
  }
  cfg_ = d.object();
  syncFormsFromConfig();
}

void MainWindow::refreshConfigList(int select) {
  syncing_ = true;
  entryList_->clear();
  for (const auto& v : cfg_["configurations"].toArray()) entryList_->addItem(v.toObject()["name"].toString());
  syncing_ = false;
  currentEntry_ = -1;
  if (entryList_->count() > 0) {
    const int r = std::clamp(select, 0, entryList_->count() - 1);
    syncing_ = true;
    entryList_->setCurrentRow(r);
    syncing_ = false;
    loadEntryEditor(r);
  }
}

void MainWindow::fillParamTable(const QString& type, const QJsonObject& params) {
  syncing_ = true;
  paramTable_->setRowCount(0);
  try {
    const auto defs = ssp::excitons::default_geometry_params(type.toStdString());
    for (const auto& [k, v] : defs) {
      const int r = paramTable_->rowCount();
      paramTable_->insertRow(r);
      auto* key = new QTableWidgetItem(QString::fromStdString(k));
      key->setFlags(key->flags() & ~Qt::ItemIsEditable);
      paramTable_->setItem(r, 0, key);
      const QString ks = QString::fromStdString(k);
      double val = v;
      if (params.contains(ks)) val = params[ks].isBool() ? (params[ks].toBool() ? 1.0 : 0.0) : params[ks].toDouble();
      const bool autoAxis = type == "particle_shell" && (ks == "a_xy" || ks == "a_z") && !params.contains(ks);
      auto* item = new QTableWidgetItem(autoAxis ? QString("auto (NP)") : QString::number(val, 'g', 10));
      if (autoAxis) item->setToolTip("Taken from the Nanoparticle tab; type a number to override");
      paramTable_->setItem(r, 1, item);
    }
  } catch (const std::exception& e) {
    geomView_->setMessage(e.what());
  }
  syncing_ = false;
}

void MainWindow::loadEntryEditor(int row) {
  const QJsonArray a = cfg_["configurations"].toArray();
  if (row < 0 || row >= a.size()) return;
  currentEntry_ = row;
  const QJsonObject o = a[row].toObject();
  syncing_ = true;
  entryName_->setText(o["name"].toString());
  entryType_->setCurrentText(o["geometry_type"].toString());
  const QJsonObject pl = o["placement"].toObject();
  centerCheck_->setChecked(pl["center"].toBool());
  const QJsonArray sh = pl["shift_nm"].toArray(), eu = pl["euler_deg"].toArray();
  for (int k = 0; k < 3; ++k) {
    shift_[k]->setValue(k < sh.size() ? sh[k].toDouble() : 0.0);
    euler_[k]->setValue(k < eu.size() ? eu[k].toDouble() : 0.0);
  }
  syncing_ = false;
  fillParamTable(entryType_->currentText(), o["geometry_params"].toObject());
  updateGeometryPreview();
}

void MainWindow::storeEntryEditor() {
  if (syncing_ || currentEntry_ < 0) return;
  QJsonArray a = cfg_["configurations"].toArray();
  if (currentEntry_ >= a.size()) return;
  QJsonObject o;
  o["name"] = entryName_->text().isEmpty() ? entryType_->currentText() : entryName_->text();
  o["geometry_type"] = entryType_->currentText();
  QJsonObject params;
  try {
    const auto defs = ssp::excitons::default_geometry_params(entryType_->currentText().toStdString());
    for (int r = 0; r < paramTable_->rowCount(); ++r) {
      const QString k = paramTable_->item(r, 0)->text();
      bool ok = false;
      const double v = paramTable_->item(r, 1) ? paramTable_->item(r, 1)->text().toDouble(&ok) : 0.0;
      if (ok && defs.contains(k.toStdString()) && v != defs.at(k.toStdString())) params[k] = v;
    }
  } catch (...) {
  }
  o["geometry_params"] = params;
  const bool place = centerCheck_->isChecked() || shift_[0]->value() != 0 || shift_[1]->value() != 0 ||
                     shift_[2]->value() != 0 || euler_[0]->value() != 0 || euler_[1]->value() != 0 ||
                     euler_[2]->value() != 0;
  if (place) {
    o["placement"] = QJsonObject{{"center", centerCheck_->isChecked()},
                                 {"shift_nm", QJsonArray{shift_[0]->value(), shift_[1]->value(), shift_[2]->value()}},
                                 {"euler_deg", QJsonArray{euler_[0]->value(), euler_[1]->value(), euler_[2]->value()}}};
  }
  a[currentEntry_] = o;
  cfg_["configurations"] = a;
}

void MainWindow::updateGeometryPreview() {
  const QJsonArray a = cfg_["configurations"].toArray();
  if (currentEntry_ < 0 || currentEntry_ >= a.size()) return;
  const QJsonObject o = a[currentEntry_].toObject();
  try {
    ssp::excitons::ParamMap pm;
    const QJsonObject gp = o["geometry_params"].toObject();
    for (const auto& k : gp.keys()) pm[k.toStdString()] = gp[k].isBool() ? double(gp[k].toBool()) : gp[k].toDouble();
    if (o["geometry_type"].toString() == "particle_shell") {  // как в load_config: полуоси из nanoparticle
      const double sh = jsonGet(cfg_, "nanoparticle.shell_thickness_nm").toDouble(0.0);  // лигандный слой
      if (!pm.contains("a_xy")) pm["a_xy"] = jsonGet(cfg_, "nanoparticle.semi_axis_xy_nm").toDouble(15.0) + sh;
      if (!pm.contains("a_z"))
        pm["a_z"] = jsonGet(cfg_, "nanoparticle.shape").toString() == "sphere"
                        ? pm["a_xy"]
                        : jsonGet(cfg_, "nanoparticle.semi_axis_z_nm").toDouble(15.0) + sh;
    }
    auto agg = ssp::excitons::build_geometry(o["geometry_type"].toString().toStdString(), pm);
    const QJsonObject pl = o["placement"].toObject();
    if (pl["center"].toBool()) ssp::excitons::transform(agg, -ssp::excitons::centroid(agg), 0, 0, 0);
    const QJsonArray sh = pl["shift_nm"].toArray(), eu = pl["euler_deg"].toArray();
    ssp::excitons::Vec3 s = ssp::excitons::Vec3::Zero();
    double e[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) {
      if (k < sh.size()) s(k) = sh[k].toDouble();
      if (k < eu.size()) e[k] = eu[k].toDouble();
    }
    ssp::excitons::transform(agg, s, e[0], e[1], e[2]);
    QVector<QVector3D> r, mu;
    for (std::size_t i = 0; i < agg.size(); ++i) {
      r << QVector3D(agg.positions[i].x(), agg.positions[i].y(), agg.positions[i].z());
      mu << QVector3D(agg.dipoles[i].x(), agg.dipoles[i].y(), agg.dipoles[i].z());
    }
    geomView_->setAggregate(r, mu);
    const bool np = jsonGet(cfg_, "nanoparticle.enabled").toBool();
    geomView_->setParticle(np, jsonGet(cfg_, "nanoparticle.semi_axis_xy_nm").toDouble(),
                           jsonGet(cfg_, "nanoparticle.shape").toString() == "sphere"
                               ? jsonGet(cfg_, "nanoparticle.semi_axis_xy_nm").toDouble()
                               : jsonGet(cfg_, "nanoparticle.semi_axis_z_nm").toDouble());
  } catch (const std::exception& ex) {
    geomView_->setMessage(QString::fromUtf8(ex.what()));
  }
}

// ---------------------------------------------------------------- файлы

void MainWindow::newConfig() {
  cfg_ = defaultConfig();
  cfgPath_.clear();
  syncFormsFromConfig();
  setWindowTitle("SSP Nano — new config");
}

bool MainWindow::openConfig(const QString& path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly)) return false;
  QJsonParseError err;
  const QJsonDocument d = QJsonDocument::fromJson(f.readAll(), &err);
  if (err.error != QJsonParseError::NoError || !d.isObject()) {
    QMessageBox::warning(this, "Config", "JSON error: " + err.errorString());
    return false;
  }
  cfg_ = d.object();
  cfgPath_ = QFileInfo(path).absoluteFilePath();
  baseDir_ = QFileInfo(path).absolutePath();
  syncFormsFromConfig();
  setWindowTitle("SSP Nano — " + QFileInfo(path).fileName());
  const QString out = resolvedOutputDir();
  if (QFile::exists(out + "/summary.csv") || QFile::exists(out + "/scan_info.csv")) loadResults(out);
  return true;
}

void MainWindow::openConfigDialog() {
  const QString p = QFileDialog::getOpenFileName(this, "Open config", baseDir_, "JSON (*.json)");
  if (!p.isEmpty()) openConfig(p);
}

bool MainWindow::saveConfig() {
  if (cfgPath_.isEmpty()) return saveConfigAs();
  syncConfigFromForms();
  QFile f(cfgPath_);
  if (!f.open(QIODevice::WriteOnly)) return false;
  f.write(QJsonDocument(cfg_).toJson(QJsonDocument::Indented));
  statusBar()->showMessage("Saved: " + cfgPath_, 4000);
  return true;
}

bool MainWindow::saveConfigAs() {
  const QString p = QFileDialog::getSaveFileName(this, "Save config", baseDir_ + "/config.json", "JSON (*.json)");
  if (p.isEmpty()) return false;
  // Относительные пути материала пересчитать на новую папку
  const QString mat = jsonGet(cfg_, "nanoparticle.material_file").toString();
  if (!mat.isEmpty() && QFileInfo(mat).isRelative()) {
    const QString abs = QDir(baseDir_).absoluteFilePath(mat);
    jsonSet(cfg_, "nanoparticle.material_file", QFileInfo(p).absoluteDir().relativeFilePath(abs));
  }
  cfgPath_ = QFileInfo(p).absoluteFilePath();
  baseDir_ = QFileInfo(p).absolutePath();
  syncing_ = true;
  npForm_->setBaseDir(baseDir_);
  npForm_->load(cfg_);
  syncing_ = false;
  setWindowTitle("SSP Nano — " + QFileInfo(p).fileName());
  return saveConfig();
}

QString MainWindow::resolvedOutputDir() const {
  const QString o = cfg_["output_dir"].toString("ssp_results");
  return QFileInfo(o).isAbsolute() ? o : QDir(baseDir_).absoluteFilePath(o);
}

// ---------------------------------------------------------------- расчёт

void MainWindow::appendLog(const QString& line) { log_->appendPlainText(line); }

void MainWindow::setRunning(bool on) {
  runAct_->setEnabled(!on);
  stopAct_->setEnabled(on);
  progress_->setVisible(on);
  paramTabs_->setEnabled(!on);
}

void MainWindow::runCalculation() {
  syncConfigFromForms();
  const std::string json = QJsonDocument(cfg_).toJson(QJsonDocument::Indented).toStdString();
  ssp::app::RunConfig rc;
  try {
    rc = ssp::app::load_config_text(json, baseDir_.toStdString(), cfgPath_.toStdString());
  } catch (const std::exception& e) {
    QMessageBox::warning(this, "Config", QString::fromUtf8(e.what()));
    return;
  }
  rc.output_dir = resolvedOutputDir().toStdString();
  log_->clear();
  resultTabs_->setCurrentWidget(log_);
  appendLog("Run: " + QString::fromStdString(rc.output_dir));
  cancel_ = false;
  runError_.clear();
  progress_->setRange(0, int(rc.configurations.size()));
  progress_->setValue(0);
  setRunning(true);
  bool scan = false;
  try {
    scan = ssp::app::parse_scan(json).enabled;
  } catch (const std::exception& e) {
    QMessageBox::warning(this, "Scan", QString::fromUtf8(e.what()));
    return;
  }
  if (scan) {
    const auto spec = ssp::app::parse_scan(json);
    progress_->setRange(0, int(rc.configurations.size() * spec.values.size()));
  }
  const std::string base = baseDir_.toStdString(), cpath = cfgPath_.toStdString();
  watcher_.setFuture(QtConcurrent::run([this, rc, json, scan, base, cpath] {
    try {
      auto logf = [this](const std::string& l) {
        const QString s = QString::fromStdString(l);
        QMetaObject::invokeMethod(this, [this, s] { appendLog(s); }, Qt::QueuedConnection);
      };
      auto prog = [this](int i, int n, const std::string& name) {
        const QString s = QString::fromStdString(name);
        QMetaObject::invokeMethod(
            this,
            [this, i, n, s] {
              progress_->setValue(i);
              progress_->setFormat(QString("%1 (%2/%3)").arg(s).arg(i + 1).arg(n));
            },
            Qt::QueuedConnection);
      };
      if (scan) {
        // output_dir уже абсолютный в rc; в JSON подставляем его же
        QJsonObject root = QJsonDocument::fromJson(QByteArray::fromStdString(json)).object();
        root["output_dir"] = QString::fromStdString(rc.output_dir);
        ssp::app::run_scan(QJsonDocument(root).toJson().toStdString(), base, cpath, logf, &cancel_, prog);
        return;
      }
      ssp::app::run(
          rc, json,
          [this](const std::string& l) {
            const QString s = QString::fromStdString(l);
            QMetaObject::invokeMethod(this, [this, s] { appendLog(s); }, Qt::QueuedConnection);
          },
          &cancel_,
          [this](int i, int n, const std::string& name) {
            const QString s = QString::fromStdString(name);
            QMetaObject::invokeMethod(
                this,
                [this, i, n, s] {
                  progress_->setValue(i);
                  progress_->setFormat(QString("%1 (%2/%3)").arg(s).arg(i + 1).arg(n));
                },
                Qt::QueuedConnection);
          });
    } catch (const std::exception& e) {
      runError_ = QString::fromUtf8(e.what());
    }
  }));
}

void MainWindow::calibrateDialog() {
  const QString f = QFileDialog::getOpenFileName(this, "Monomer calibration config (JSON with \"monomer_fit\")",
                                                 baseDir_, "JSON (*.json)");
  if (!f.isEmpty()) startCalibration(f);
}

void MainWindow::startCalibration(const QString& path) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly)) {
    QMessageBox::warning(this, "Calibration", "Cannot open " + path);
    return;
  }
  const std::string json = f.readAll().toStdString();
  if (!ssp::app::is_calibration_config(json)) {
    QMessageBox::warning(this, "Calibration", "The file has no \"monomer_fit\" section");
    return;
  }
  log_->clear();
  resultTabs_->setCurrentWidget(log_);
  appendLog("Monomer calibration: " + path);
  calibError_.clear();
  runAct_->setEnabled(false);
  progress_->setRange(0, 0);
  progress_->setVisible(true);
  const std::string base = ssp::app::config_dir_of(path.toStdString());
  calWatcher_.setFuture(QtConcurrent::run([this, json, base] {
    try {
      auto out = std::make_shared<ssp::app::CalibrationOutput>(ssp::app::run_calibration(json, base, [this](const std::string& l) {
        const QString s = QString::fromStdString(l);
        QMetaObject::invokeMethod(this, [this, s] { appendLog(s); }, Qt::QueuedConnection);
      }));
      QMetaObject::invokeMethod(this, [this, out] { calib_ = out; }, Qt::BlockingQueuedConnection);
    } catch (const std::exception& e) {
      calibError_ = QString::fromUtf8(e.what());
    }
  }));
}

void MainWindow::showCalibration() {
  if (!calib_) return;
  const auto& C = *calib_;
  const auto& F = C.fit;
  auto vec = [](const std::vector<double>& v) { return QVector<double>(v.begin(), v.end()); };
  calAbsPlot_->clear();
  calEmPlot_->clear();
  calResPlot_->clear();
  calAbsPlot_->addSeries({"absorption (data)", vec(C.abs.w), vec(C.abs.y), {}, {}, PlotWidget::paletteColor(0),
                          PlotWidget::Style::Points, 1, true});
  calAbsPlot_->addSeries({"model", vec(C.abs.w), vec(F.abs_model), {}, {}, PlotWidget::paletteColor(1),
                          PlotWidget::Style::Line, 1.8, true});
  calAbsPlot_->setTitle(QString("Monomer absorption: cumulant fit, χ²_ν = %1").arg(F.chi2_red, 0, 'f', 2));
  calAbsPlot_->setLabels("ν̃, cm⁻¹", "absorption (data units)");
  calAbsPlot_->autoscale();
  auto resid = [](const ssp::excitons::Spectrum& s, const std::vector<double>& m) {
    QVector<double> r;
    const double mx = s.y.empty() ? 1.0 : *std::max_element(s.y.begin(), s.y.end());
    for (std::size_t i = 0; i < s.w.size(); ++i)
      r << (s.y[i] - m[i]) / (s.sigma.size() == s.w.size() ? s.sigma[i] : 0.02 * mx);
    return r;
  };
  calResPlot_->addSeries({"absorption", vec(C.abs.w), resid(C.abs, F.abs_model), {}, {}, PlotWidget::paletteColor(0),
                          PlotWidget::Style::Points, 1, true});
  if (!C.em.w.empty()) {
    calEmPlot_->addSeries({"PL (data)", vec(C.em.w), vec(C.em.y), {}, {}, PlotWidget::paletteColor(2),
                           PlotWidget::Style::Points, 1, true});
    calEmPlot_->addSeries({"model", vec(C.em.w), vec(F.em_model), {}, {}, PlotWidget::paletteColor(1),
                           PlotWidget::Style::Line, 1.8, true});
    calResPlot_->addSeries({"PL", vec(C.em.w), resid(C.em, F.em_model), {}, {}, PlotWidget::paletteColor(2),
                            PlotWidget::Style::Points, 1, true});
  }
  calEmPlot_->setTitle(C.em.w.empty() ? "No PL spectrum given" : "Monomer PL (per wavenumber interval)");
  calEmPlot_->setLabels("ν̃, cm⁻¹", "PL (data units)");
  calEmPlot_->autoscale();
  calResPlot_->setTitle("Residuals (data − model)/σ — structure means missing physics");
  calResPlot_->setLabels("ν̃, cm⁻¹", "residual, σ");
  calResPlot_->autoscale();
  calTable_->clear();
  calTable_->setColumnCount(3);
  calTable_->setHorizontalHeaderLabels({"parameter", "value", "± error"});
  QList<QStringList> rows;
  for (std::size_t k = 0; k < F.names.size(); ++k)
    rows << QStringList{QString::fromStdString(F.names[k]), QString::number(F.values[k], 'g', 6),
                        QString::number(F.errors[k], 'g', 3)};
  double lt = F.p.lambda_cm;
  for (const auto& m : F.p.modes) lt += m.huang_rhys * m.omega_cm;
  rows << QStringList{"λ_tot (Drude + modes), cm⁻¹", QString::number(lt, 'g', 5), ""};
  if (F.mu_D > 0) {
    rows << QStringList{"μ, D (∫ε/ν̃)", QString::number(F.mu_D, 'g', 5), ""};
    rows << QStringList{"τ_rad, ns (model formula)", QString::number(F.tau_rad_ns, 'g', 4), ""};
    rows << QStringList{"τ_rad, ns (Strickler–Berg)", QString::number(F.tau_sb_ns, 'g', 4), ""};
  }
  rows << QStringList{"χ²_ν", QString::number(F.chi2_red, 'g', 4), F.converged ? "converged" : "NOT converged"};
  calTable_->setRowCount(rows.size());
  for (int r = 0; r < rows.size(); ++r)
    for (int c = 0; c < 3; ++c) calTable_->setItem(r, c, new QTableWidgetItem(rows[r].value(c)));
  calTable_->resizeColumnsToContents();
  calLabel_->setText("Results: " + QString::fromStdString(C.output_dir) +
                     "\nσ_inh of the monomer is an upper bound for the site disorder σ of the aggregate.");
  resultTabs_->setCurrentWidget(calTab_);
}

void MainWindow::applyCalibration() {
  if (!calib_) {
    QMessageBox::information(this, "Calibration", "Run File → Calibrate monomer… first");
    return;
  }
  syncConfigFromForms();
  const auto& F = calib_->fit;
  jsonSet(cfg_, "exciton.monomer_energy_cm", F.p.E0_cm);
  if (F.mu_D > 0) jsonSet(cfg_, "exciton.transition_dipole_D", F.mu_D);
  jsonSet(cfg_, "bath.type", "drude_lorentz");
  jsonSet(cfg_, "bath.lambda_cm", F.p.lambda_cm);
  jsonSet(cfg_, "bath.cutoff_cm", F.p.cutoff_cm);
  jsonSet(cfg_, "bath.temperature_K", F.p.temperature_K);
  QJsonArray modes;
  for (const auto& m : F.p.modes)
    modes.append(QJsonObject{{"omega_cm", m.omega_cm}, {"huang_rhys", m.huang_rhys}, {"gamma_cm", m.gamma_cm}});
  jsonSet(cfg_, "bath.modes", modes);
  syncFormsFromConfig();
  appendLog("Calibration applied: monomer energy, μ, bath (Drude + modes). Site disorder σ is not changed "
            "(monomer σ_inh = " + QString::number(F.p.sigma_inh_cm, 'f', 1) + " cm⁻¹ is an upper bound).");
}

void MainWindow::stopCalculation() {
  cancel_ = true;
  appendLog("Stopping after the current shape…");
}

// ---------------------------------------------------------------- результаты

void MainWindow::openResultsDialog() {
  const QString d = QFileDialog::getExistingDirectory(this, "Results folder", resultsDir_.isEmpty() ? baseDir_ : resultsDir_);
  if (!d.isEmpty()) loadResults(d);
}

void MainWindow::loadResults(const QString& dir) {
  scan_.clear();
  scanValues_.clear();
  if (QFile::exists(dir + "/scan_info.csv")) {
    CsvTable info;
    info.load(dir + "/scan_info.csv");
    const QString first = info.rows() ? info.text[0].value(2) : QString();
    if (!first.isEmpty() && QFile::exists(dir + "/" + first + "/summary.csv")) loadRunResults(dir + "/" + first);
    loadScanResults(dir);
    resultTabs_->setCurrentWidget(scanTab_);
    return;
  }
  loadRunResults(dir);
}

void MainWindow::loadScanResults(const QString& dir) {
  CsvTable info;
  info.load(dir + "/scan_info.csv");
  scanParameter_ = info.comments.isEmpty() ? QString("scan value") : info.comments[0].section("scan over ", 1);
  scanValues_ = info.col("value");
  scan_.clear();
  QDir d(dir);
  for (const QString& sub : d.entryList(QDir::Dirs | QDir::NoDotAndDotDot)) {
    if (sub.startsWith("scan_")) continue;
    ScanShapeData sd;
    for (const QString& m : QStringList{"abs_sys", "ext", "sca", "abs_diff", "abs", "em"}) {
      CsvTable t;
      if (t.load(dir + "/" + sub + "/scan_" + m + ".csv")) sd.maps[m] = t;
    }
    sd.branches.load(dir + "/" + sub + "/scan_branches.csv");
    sd.coupling.load(dir + "/" + sub + "/scan_coupling.csv");
    sd.observables.load(dir + "/" + sub + "/scan_observables.csv");
    if (!sd.maps.isEmpty()) scan_[sub] = sd;
  }
  {
    QSignalBlocker b1(scanShape_), b2(scanObs_);
    scanShape_->clear();
    scanShape_->addItems(scan_.keys());
    scanObs_->clear();
    if (!scan_.isEmpty()) {
      const auto& h = scan_.first().observables.header;
      for (int i = 1; i < h.size(); ++i) scanObs_->addItem(h[i]);
      const int k = scanObs_->findText("abs_peak_cm");
      if (k >= 0) scanObs_->setCurrentIndex(k);
    }
  }
  plotScan();
  plotScanObservable();
  statusBar()->showMessage("Scan loaded: " + dir, 5000);
}

void MainWindow::plotScan() {
  scanPlot_->clear();
  scanPlot_->clearHeatmap();
  scanCoupling_->clear();
  if (!scan_.contains(scanShape_->currentText())) return;
  const auto& sd = scan_[scanShape_->currentText()];
  static const QStringList keys = {"abs_sys", "ext", "sca", "abs_diff", "abs", "em"};
  const QString key = keys.value(scanMap_->currentIndex());
  if (!sd.maps.contains(key)) {
    scanCoupling_->setText("This map is not available (no nanoparticle in the calculation?).");
    return;
  }
  const CsvTable& t = sd.maps[key];
  const QVector<double> w = t.cols.value(0);
  QVector<QVector<double>> z;
  for (int c = 1; c < t.cols.size(); ++c) {
    QVector<double> row = t.cols[c];
    if (scanNormRows_->isChecked()) {
      double m = 0;
      for (double v : row) m = std::max(m, std::abs(v));
      if (m > 0)
        for (double& v : row) v /= m;
    }
    z << row;
  }
  QVector<double> vals = scanValues_.mid(0, z.size());
  scanPlot_->setHeatmap(w, vals, z, scanNormRows_->isChecked() ? "normalized per row" : "σ, nm²");
  bool resolved = false;
  for (int r = 0; r < sd.coupling.rows(); ++r)
    if (sd.coupling.text[r].value(0) == "abs") resolved = sd.coupling.text[r].value(sd.coupling.column("resolved")) == "yes";
  if (scanOverlay_->isChecked() && sd.branches.rows() > 0) {
    const auto& b = sd.branches;
    const QVector<double> v = b.col("value");
    scanPlot_->addSeries({"reference: NP dipolar mode (pole fit)", b.col("w_p"), v, {}, {}, QColor(30, 30, 30), PlotWidget::Style::Points, 1, true});
    scanPlot_->addSeries({"reference: aggregate alone", b.col("w_x"), v, {}, {}, QColor(0, 150, 60), PlotWidget::Style::Points, 1, true});
    const QString pk = key == "ext" ? "ext_peak" : key == "sca" ? "sca_peak" : "abs_peak";
    for (int j = 1; j <= 3; ++j)
      scanPlot_->addSeries({j == 1 ? QString("maxima of the system %1").arg(pk.section('_', 0, 0)) : QString(),
                            b.col(pk + QString::number(j)), v, {}, {}, QColor(230, 120, 0), PlotWidget::Style::Points, 1, true});
    if (resolved) {
      scanPlot_->addSeries({"coupled-oscillator branches (fit)", b.col("model_lower"), v, {}, {}, QColor(30, 30, 30), PlotWidget::Style::Dashed, 1.4, true});
      scanPlot_->addSeries({QString(), b.col("model_upper"), v, {}, {}, QColor(30, 30, 30), PlotWidget::Style::Dashed, 1.4, true});
    }
  }
  static const QStringList titles = {"System absorption", "System extinction", "System scattering",
                                     "System − NP absorption (aggregate contribution)", "Aggregate absorption (no NP)",
                                     "Aggregate fluorescence (no NP)"};
  scanPlot_->setTitle(titles.value(scanMap_->currentIndex()) + " vs " + scanParameter_);
  scanPlot_->setLabels("ν̃, cm⁻¹", scanParameter_);
  scanPlot_->autoscale();

  // Таблица связи
  const auto& c = sd.coupling;
  if (c.rows() > 0) {
    QString h = "<b>Coupling analysis</b> (γ = FWHM; γp, ωp from the NP dipolar-mode pole; γx from the aggregate's "
                "ensemble spectrum; values at the smallest detuning)<br>"
                "<table border='1' cellspacing='0' cellpadding='3'><tr><th>method</th><th>g, cm⁻¹</th><th>interval</th>"
                "<th>branches resolved</th><th>γp</th><th>γx</th><th>2g &gt; (γp+γx)/2</th><th>g &gt; |γp−γx|/4</th><th>C</th></tr>";
    for (int r = 0; r < c.rows(); ++r) {
      auto cell = [&](const QString& n, int prec = 4) {
        const QString v = c.text[r].value(c.column(n));
        bool ok;
        const double x = v.toDouble(&ok);
        return ok ? QString::number(x, 'g', prec) : v;
      };
      auto yes = [&](const QString& n) {
        const bool y = cell(n) == "yes";
        return QString("<font color='%1'>%2</font>").arg(y ? "#1a8f3a" : "#b03030", y ? "yes" : "no");
      };
      const QString kind = cell("spectrum");
      const bool micro = kind == "micro";
      const bool res = cell("resolved") == "yes";
      const QString label = micro ? "<b>microscopic</b> (mode × dipoles)" : "peak fit: " + kind;
      const QString g = (!micro && !res) ? QString("<i>undefined</i>") : cell("g_cm");
      const QString iv = micro ? QString("—") : (res ? "[" + cell("g_lo_cm") + ", " + cell("g_hi_cm") + "]" : QString("—"));
      h += "<tr><td>" + label + "</td><td>" + g + "</td><td>" + iv + "</td><td>" +
           (micro ? QString("n/a") : yes("resolved") + " (" + cell("two_peak_points") + " pts)") + "</td><td>" +
           cell("gamma_p_cm") + "</td><td>" + cell("gamma_x_cm") + "</td><td>" +
           ((!micro && !res) ? QString("—") : yes("visible_splitting")) + "</td><td>" +
           ((!micro && !res) ? QString("—") : yes("anticrossing")) + "</td><td>" +
           ((!micro && !res) ? QString("—") : cell("cooperativity", 3)) + "</td></tr>";
    }
    h += "</table>The microscopic g needs no resolved peaks. A peak-based g is reported only if both branches are seen "
         "in ≥ 2 scan points. Strong coupling is claimed only from absorption.";
    scanCoupling_->setText(h);
  }
}

void MainWindow::plotScanObservable() {
  scanObsPlot_->clear();
  if (!scan_.contains(scanShape_->currentText())) return;
  const auto& o = scan_[scanShape_->currentText()].observables;
  if (o.rows() == 0) return;
  PlotWidget::Series s;
  s.name = scanObs_->currentText();
  s.x = o.col("value");
  s.y = o.col(scanObs_->currentText());
  s.color = PlotWidget::paletteColor(0);
  scanObsPlot_->addSeries(s);
  PlotWidget::Series p = s;
  p.name.clear();
  p.style = PlotWidget::Style::Points;
  scanObsPlot_->addSeries(p);
  scanObsPlot_->setTitle(scanObs_->currentText() + " vs " + scanParameter_);
  scanObsPlot_->setLabels(scanParameter_, scanObs_->currentText());
  scanObsPlot_->setLegendVisible(false);
  scanObsPlot_->autoscale();
}

void MainWindow::loadRunResults(const QString& dir) {
  CsvTable summary;
  QString err;
  if (!summary.load(dir + "/summary.csv", &err)) {
    QMessageBox::warning(this, "Results", "No summary.csv in " + dir);
    return;
  }
  resultsDir_ = dir;
  results_.clear();
  monomerCm_ = 0.0;
  {
    QFile f(dir + "/config_used.json");
    if (f.open(QIODevice::ReadOnly))
      monomerCm_ = jsonGet(QJsonDocument::fromJson(f.readAll()).object(), "exciton.monomer_energy_cm").toDouble(0.0);
  }
  const int nc = summary.column("name");
  for (int r = 0; r < summary.rows(); ++r) {
    ConfigResult cr;
    cr.name = summary.text[r][nc];
    const QString d = dir + "/" + cr.name;
    cr.spectra.load(d + "/spectra.csv");
    cr.observables.load(d + "/observables.csv");
    cr.plexciton.load(d + "/plexciton.csv");
    cr.dynamics.load(d + "/dynamics.csv");
    cr.heomAbs.load(d + "/heom_absorption.csv");
    cr.heomPop.load(d + "/heom_populations.csv");
    cr.nearField.load(d + "/molecules_near_field.csv");
    cr.vibronic.load(d + "/vibronic_absorption.csv");
    cr.hybridEm.load(d + "/hybrid_emission.csv");
    cr.geometry.load(d + "/geometry.csv");
    cr.states.load(d + "/states_representative.csv");
    cr.n = summary.text[r].value(summary.column("N")).toInt();
    results_ << cr;
  }
  // Сводка: строки — величины, столбцы — формы
  summaryTable_->clear();
  summaryTable_->setColumnCount(results_.size() + 1);
  QStringList hdr{"quantity"};
  for (const auto& c : results_) hdr << c.name;
  summaryTable_->setHorizontalHeaderLabels(hdr);
  QStringList rows;
  for (int c = 1; c < summary.header.size(); ++c)
    if (!summary.header[c].endsWith("_se") && summary.header[c] != "name") rows << summary.header[c];
  summaryTable_->setRowCount(rows.size());
  QMap<QString, QString> desc;
  if (!results_.isEmpty()) {
    const auto& ob = results_[0].observables;
    for (int r = 0; r < ob.rows(); ++r) desc[ob.text[r][0]] = ob.text[r].value(4);
  }
  for (int i = 0; i < rows.size(); ++i) {
    auto* it = new QTableWidgetItem(rows[i]);
    it->setToolTip(desc.value(rows[i]));
    summaryTable_->setItem(i, 0, it);
    const int c = summary.column(rows[i]), cs = summary.column(rows[i] + "_se");
    for (int r = 0; r < summary.rows(); ++r) {
      QString v = summary.text[r].value(c);
      bool ok;
      const double dv = v.toDouble(&ok);
      if (ok) v = QString::number(dv, 'g', 6);
      if (cs >= 0 && !summary.text[r].value(cs).isEmpty()) v += " ± " + QString::number(summary.text[r].value(cs).toDouble(), 'g', 2);
      summaryTable_->setItem(i, r + 1, new QTableWidgetItem(v));
    }
  }
  summaryTable_->resizeColumnsToContents();
  refreshResultViews();
  statusBar()->showMessage("Results loaded: " + dir, 5000);
}

const MainWindow::ConfigResult* MainWindow::resultByName(const QString& n) const {
  for (const auto& r : results_)
    if (r.name == n) return &r;
  return nullptr;
}

void MainWindow::refreshResultViews() {
  QSignalBlocker b1(specConfigs_), b2(plexConfig_), b3(dynConfig_), b4(heomConfig_), b5(nfConfig_), b6(compareTarget_);
  specConfigs_->clear();
  for (auto* c : {plexConfig_, dynConfig_, heomConfig_, nfConfig_, compareTarget_}) c->clear();
  for (const auto& r : results_) {
    auto* it = new QListWidgetItem(r.name);
    it->setFlags(it->flags() | Qt::ItemIsUserCheckable);
    it->setCheckState(Qt::Checked);
    specConfigs_->addItem(it);
    compareTarget_->addItem(r.name);
    if (r.plexciton.rows()) plexConfig_->addItem(r.name);
    if (r.dynamics.rows()) dynConfig_->addItem(r.name);
    if (r.heomAbs.rows() || r.heomPop.rows()) heomConfig_->addItem(r.name);
    if (r.nearField.rows()) nfConfig_->addItem(r.name);
  }
  plotSpectra();
  plotPlexciton();
  plotDynamics();
  plotHeom();
  plotNearField();
  resultTabs_->setCurrentIndex(0);
}

void MainWindow::plotSpectra() {
  specPlot_->clear();
  const int mode = specMode_->currentIndex();
  const bool norm = specNormalize_->isChecked(), err = specErrors_->isChecked();
  const bool perMol = specPerMolecule_->isChecked() && !norm;
  int ci = 0;
  double ymax = 0.0;
  const ConfigResult* first = nullptr;
  for (int i = 0; i < specConfigs_->count(); ++i) {
    const auto* r = resultByName(specConfigs_->item(i)->text());
    if (!r || r->spectra.rows() == 0) continue;
    const QColor col = PlotWidget::paletteColor(ci++);
    if (specConfigs_->item(i)->checkState() != Qt::Checked) continue;
    if (!first) first = r;
    auto addOne = [&](const QString& ycol, const QString& secol, const QString& suffix, PlotWidget::Style st) {
      PlotWidget::Series s;
      s.name = r->name + suffix;
      s.x = r->spectra.col("w_cm");
      s.y = r->spectra.col(ycol);
      const QVector<double> se = r->spectra.col(secol);
      double mx = 0;
      for (double v : s.y) mx = std::max(mx, v);
      double f = 1.0;
      if (norm && mx > 0) f = 1.0 / mx;
      else if (perMol && r->n > 0 && ycol == "abs_nm2") f = 1.0 / r->n;
      for (double& v : s.y) {
        v *= f;
        ymax = std::max(ymax, v);
      }
      if (err) {
        for (int k = 0; k < s.y.size(); ++k) {
          s.lo << s.y[k] - se[k] * f;
          s.hi << s.y[k] + se[k] * f;
        }
      }
      s.color = col;
      s.style = st;
      specPlot_->addSeries(s);
    };
    if (mode == 0 || mode == 2) addOne("abs_nm2", "abs_se", mode == 2 ? " (abs.)" : "", PlotWidget::Style::Line);
    if ((mode == 0 || mode == 2) && r->vibronic.rows() > 0) {  // модель Хольстейна (OPA)
      PlotWidget::Series v;
      v.name = r->name + " (Holstein OPA)";
      v.x = r->vibronic.col("w_cm");
      v.y = r->vibronic.col("abs_nm2");
      double mx = 0;
      for (double y : v.y) mx = std::max(mx, y);
      const double f = norm && mx > 0 ? 1.0 / mx : (perMol && r->n > 0 ? 1.0 / r->n : 1.0);
      for (double& y : v.y) {
        y *= f;
        ymax = std::max(ymax, y);
      }
      v.color = col;
      v.style = PlotWidget::Style::Dashed;
      v.width = 1.1;
      specPlot_->addSeries(v);
    }
    if (mode == 1 || mode == 2) addOne("em", "em_se", mode == 2 ? " (PL)" : "", mode == 2 ? PlotWidget::Style::Dashed : PlotWidget::Style::Line);
  }
  // «Палочки» экситонных состояний одной реализации (|μ_n|², масштаб: максимум = 35 % высоты графика)
  if (specSticks_->isChecked() && first && first->states.rows() > 0 && ymax > 0) {
    const QVector<double> e = first->states.col("energy_cm"), m2 = first->states.col("mu2_over_mu2mono");
    double mm = 0;
    for (double v : m2) mm = std::max(mm, v);
    PlotWidget::Series st;
    st.name = first->name + ": exciton states |μ|² (1 realization, scaled)";
    for (int k = 0; k < e.size(); ++k) {
      st.x << e[k] << e[k] << e[k];
      st.y << 0.0 << 0.35 * ymax * m2[k] / std::max(mm, 1e-300) << std::numeric_limits<double>::quiet_NaN();
    }
    st.color = QColor(120, 120, 120);
    st.width = 1.0;
    specPlot_->addSeries(st);
  }
  if (monomerCm_ > 0) specPlot_->setVerticalLines({monomerCm_}, {"monomer"});
  specPlot_->setLabels("ν̃, cm⁻¹", norm ? "normalized intensity"
                                       : (mode == 1 ? "PL, photons/(ps·cm⁻¹)"
                                                    : (perMol ? "σ_abs per molecule, nm²" : "σ_abs, nm² per aggregate")));
  specPlot_->setTitle((mode == 0 ? QString("Absorption") : mode == 1 ? QString("Fluorescence") : QString("Absorption and fluorescence")) +
                      " — ensemble average; bands: ±SE of disorder sampling only (not model uncertainty)");
  runComparison();
  specPlot_->autoscale();
}

void MainWindow::runComparison() {
  residPlot_->clear();
  if (expX_.isEmpty()) {
    residPlot_->setVisible(false);
    compareLabel_->setText("no experiment loaded");
    return;
  }
  const auto* r = resultByName(compareTarget_->currentText());
  if (!r || r->spectra.rows() == 0) return;
  const bool em = specMode_->currentIndex() == 1;
  const QVector<double> mx = r->spectra.col("w_cm"), my = r->spectra.col(em ? "em" : "abs_nm2");
  CompareOptions opt;
  opt.fit_offset = fitOffset_->isChecked();
  opt.fit_shift = fitShift_->isChecked();
  const CompareResult c = compareSpectra(mx, my, expX_, expY_, expS_, opt);
  if (!c.ok) {
    compareLabel_->setText(c.message);
    return;
  }
  // Пересчёт в единицы основного графика (нормировка на максимум модели)
  double mmax = 0;
  for (double v : my) mmax = std::max(mmax, v);
  const double f = (specNormalize_->isChecked() && mmax > 0 && c.scale != 0) ? 1.0 / (c.scale * mmax) : 1.0;
  PlotWidget::Series e;
  e.name = "experiment: " + expName_;
  e.x = c.x;
  for (double v : c.y) e.y << (v - c.offset) * f;
  e.color = QColor(30, 30, 30);
  if (palette().color(QPalette::Base).lightness() < 128) e.color = QColor(235, 235, 235);
  e.style = PlotWidget::Style::Points;
  specPlot_->addSeries(e);
  PlotWidget::Series m;
  m.name = QString("model (fit)%1").arg(opt.fit_shift ? QString(", δ = %1 cm⁻¹").arg(c.shift, 0, 'f', 1) : QString());
  m.x = c.x;
  for (double v : c.model) m.y << (v - c.offset) * f;
  m.color = QColor(200, 30, 30);
  m.width = 1.2;
  specPlot_->addSeries(m);

  residPlot_->setVisible(true);
  PlotWidget::Series rs;
  rs.x = c.x;
  rs.y = residNormalized_->isChecked() ? c.normalized_residual : c.residual;
  rs.color = QColor(200, 30, 30);
  rs.style = PlotWidget::Style::Points;
  residPlot_->addSeries(rs);
  residPlot_->setReferenceLines(residNormalized_->isChecked() ? QVector<double>{-2, -1, 1, 2} : QVector<double>{});
  residPlot_->setLabels("ν̃, cm⁻¹", residNormalized_->isChecked() ? "(y − model)/σ" : "y − model");
  residPlot_->setTitle("Residuals");
  residPlot_->autoscale();

  QString t = QString("<b>%1</b> vs %2 (%3)<br>").arg(compareTarget_->currentText(), expName_, em ? "PL" : "absorption");
  t += QString("points: %1, d.o.f.: %2<br>scale s = %3").arg(c.n).arg(c.dof).arg(fmt(c.scale));
  if (opt.fit_offset) t += QString(", baseline b = %1").arg(fmt(c.offset));
  if (opt.fit_shift) t += QString(", shift δ = %1 cm⁻¹").arg(c.shift, 0, 'f', 1);
  t += "<br>";
  if (c.have_sigma) t += QString("χ²_ν = %1<br>").arg(fmt(c.chi2_nu));
  else t += QString("RMSE = %1 (no σ given — χ²_ν undefined)<br>").arg(fmt(c.rmse));
  t += QString("Durbin–Watson = %1, autocorr.(1) = %2").arg(c.durbin_watson, 0, 'f', 2).arg(c.lag1_autocorr, 0, 'f', 2);
  if (c.durbin_watson < 1.0) t += "<br><font color='#c03030'>residuals are structured — the model misses something (M7)</font>";
  compareLabel_->setText(t);
}

bool MainWindow::loadExperimentFile(const QString& path, const QString& unit) {
  QString err;
  QVector<double> x, y, s;
  if (!loadXYTable(path, x, y, s, &err)) {
    QMessageBox::warning(this, "Experiment", err);
    return false;
  }
  convertToWavenumber(x, y, s, unit);
  expX_ = x;
  expY_ = y;
  expS_ = s;
  expName_ = QFileInfo(path).fileName();
  plotSpectra();
  return true;
}

void MainWindow::loadExperimentDialog() {
  const QString p = QFileDialog::getOpenFileName(this, "Experimental spectrum", baseDir_,
                                                 "Tables (*.csv *.txt *.dat *.tsv);;All files (*)");
  if (p.isEmpty()) return;
  bool ok;
  const QString u = QInputDialog::getItem(this, "x-axis units", "x axis in the file:", {"cm⁻¹", "nm", "eV"}, 1, false, &ok);
  if (!ok) return;
  if (loadExperimentFile(p, u == "nm" ? "nm" : u == "eV" ? "eV" : "cm")) resultTabs_->setCurrentIndex(0);
}

void MainWindow::clearExperiment() {
  expX_.clear();
  expY_.clear();
  expS_.clear();
  plotSpectra();
}

void MainWindow::plotPlexciton() {
  plexPlot_->clear();
  const auto* r = resultByName(plexConfig_->currentText());
  if (!r) return;
  if (plexView_->currentIndex() == 2) {
    const auto& h = r->hybridEm;
    if (h.rows() == 0) {
      plexPlot_->setTitle("No hybrid_emission.csv (nanoparticle.emission = false or an older run)");
      return;
    }
    const QVector<double> w = h.col("w_cm");
    PlotWidget::Series a{"PL next to the NP", w, h.col("pl_hybrid"), {}, {}, PlotWidget::paletteColor(0),
                         PlotWidget::Style::Line, 1.6, true};
    PlotWidget::Series b{"PL without the NP", w, h.col("pl_free"), {}, {}, PlotWidget::paletteColor(1),
                         PlotWidget::Style::Dashed, 1.2, true};
    plexPlot_->addSeries(a);
    plexPlot_->addSeries(b);
    plexPlot_->setZeroLine(false);
    plexPlot_->setTitle("Photoluminescence for the same absorbed pump; " + h.comments.value(1));
    plexPlot_->setLabels("ν̃, cm⁻¹", "photons/(ps·cm⁻¹)");
    plexPlot_->autoscale();
    return;
  }
  const auto& t = r->plexciton;
  const QVector<double> w = t.col("w_cm");
  const bool diff = plexView_->currentIndex() == 0;
  const char* names[3][3] = {{"ext", "ext_se", "np_ext"}, {"abs", "abs_se", "np_abs"}, {"sca", "sca_se", "np_sca"}};
  const QString labels[3] = {"extinction", "absorption", "scattering"};
  for (int k = 0; k < 3; ++k) {
    PlotWidget::Series s;
    s.name = labels[k] + (diff ? " (system − NP)" : " (system)");
    s.x = w;
    const auto y = t.col(names[k][0]), se = t.col(names[k][1]), np = t.col(names[k][2]);
    for (int i = 0; i < w.size(); ++i) {
      const double v = diff ? y[i] - np[i] : y[i];
      s.y << v;
      s.lo << v - se[i];
      s.hi << v + se[i];
    }
    s.color = PlotWidget::paletteColor(k);
    plexPlot_->addSeries(s);
    if (!diff) {
      PlotWidget::Series n;
      n.name = labels[k] + " (NP only)";
      n.x = w;
      n.y = np;
      n.color = PlotWidget::paletteColor(k);
      n.style = PlotWidget::Style::Dashed;
      n.width = 1.1;
      plexPlot_->addSeries(n);
    }
  }
  plexPlot_->setZeroLine(diff);
  plexPlot_->setTitle("Aggregate + nanoparticle: cross sections (orientation average), ±SE bands over disorder");
  plexPlot_->setLabels("ν̃, cm⁻¹", "σ, nm²");
  plexPlot_->autoscale();
}

void MainWindow::plotDynamics() {
  dynPlot_->clear();
  dynEnergyPlot_->clear();
  const auto* r = resultByName(dynConfig_->currentText());
  if (!r) return;
  const auto& d = r->dynamics;
  const QVector<double> t = d.col("t_ps");
  const QString cols[3] = {"total", "lowest", "brightest"};
  const QString labs[3] = {"total population", "lowest state", "brightest state"};
  for (int k = 0; k < 3; ++k) {
    PlotWidget::Series s;
    s.name = labs[k];
    s.x = t;
    s.y = d.col(cols[k]);
    const QVector<double> se = d.col(cols[k] + "_se");
    if (se.size() == s.y.size())
      for (int i = 0; i < s.y.size(); ++i) {
        s.lo << s.y[i] - se[i];
        s.hi << s.y[i] + se[i];
      }
    dynPlot_->addSeries(s);
  }
  dynPlot_->setTitle("Population kinetics — disorder-ensemble average ±SE (Redfield + decay)");
  dynPlot_->setLabels("t, ps", "population");
  dynPlot_->autoscale();
  PlotWidget::Series e;
  e.name = "mean energy";
  e.x = t;
  e.y = d.col("mean_energy_cm");
  {
    const QVector<double> se = d.col("mean_energy_se");
    if (se.size() == e.y.size())
      for (int i = 0; i < e.y.size(); ++i) {
        e.lo << e.y[i] - se[i];
        e.hi << e.y[i] + se[i];
      }
  }
  e.color = PlotWidget::paletteColor(3);
  dynEnergyPlot_->addSeries(e);
  dynEnergyPlot_->setTitle("Intraband relaxation");
  dynEnergyPlot_->setLabels("t, ps", "⟨E⟩, cm⁻¹");
  dynEnergyPlot_->autoscale();
}

void MainWindow::plotHeom() {
  heomAbsPlot_->clear();
  heomResidPlot_->clear();
  heomPopPlot_->clear();
  heomLabel_->clear();
  const auto* r = resultByName(heomConfig_->currentText());
  if (!r) return;
  if (r->heomAbs.rows()) {
    const QVector<double> w = r->heomAbs.col("w_cm"), h = r->heomAbs.col("heom"), c = r->heomAbs.col("cumulant");
    PlotWidget::Series sh{"HEOM (reference)", w, h, {}, {}, PlotWidget::paletteColor(6), PlotWidget::Style::Line, 2.2, true};
    PlotWidget::Series sc{"cumulant + Redfield", w, c, {}, {}, PlotWidget::paletteColor(1), PlotWidget::Style::Dashed, 1.6, true};
    heomAbsPlot_->addSeries(sh);
    heomAbsPlot_->addSeries(sc);
    heomAbsPlot_->setTitle("Absorption of one realization: fast theory vs HEOM");
    heomAbsPlot_->setLabels("ν̃, cm⁻¹", "normalized line shape, cm");
    heomAbsPlot_->autoscale();
    PlotWidget::Series d;
    d.name = "cumulant − HEOM";
    d.x = w;
    double pk = 0, md = 0;
    for (int i = 0; i < w.size(); ++i) {
      d.y << c[i] - h[i];
      pk = std::max(pk, h[i]);
      md = std::max(md, std::abs(c[i] - h[i]));
    }
    d.color = PlotWidget::paletteColor(1);
    heomResidPlot_->addSeries(d);
    heomResidPlot_->setTitle("Residual");
    heomResidPlot_->setLabels("ν̃, cm⁻¹", "Δ");
    heomResidPlot_->autoscale();
    heomLabel_->setText(QString("max|Δ|/peak = %1 %").arg(100.0 * md / std::max(pk, 1e-300), 0, 'f', 1));
  }
  if (r->heomPop.rows()) {
    const QVector<double> t = r->heomPop.col("t_ps");
    int k = 0;
    while (r->heomPop.column(QString("heom_%1").arg(k)) >= 0) {
      const QColor col = PlotWidget::paletteColor(k);
      heomPopPlot_->addSeries({QString("HEOM %1").arg(k), t, r->heomPop.col(QString("heom_%1").arg(k)), {}, {}, col,
                               PlotWidget::Style::Line, 2.0, true});
      heomPopPlot_->addSeries({QString("Redfield %1").arg(k), t, r->heomPop.col(QString("redfield_%1").arg(k)), {}, {},
                               col, PlotWidget::Style::Dashed, 1.3, true});
      ++k;
    }
    heomPopPlot_->setTitle("Exciton state populations (solid — HEOM, dashed — Redfield)");
    heomPopPlot_->setLabels("t, ps", "population");
    heomPopPlot_->setLegendVisible(k <= 4);
    heomPopPlot_->autoscale();
  }
}

void MainWindow::plotNearField() {
  nfEnhPlot_->clear();
  nfDecayPlot_->clear();
  const auto* r = resultByName(nfConfig_->currentText());
  if (!r) {
    nfGeom_->setAggregate({}, {});
    return;
  }
  const auto& t = r->nearField;
  PlotWidget::Series e{"⟨|E/E₀|²⟩", t.col("gap_nm"), t.col("field_enhancement"), {}, {}, PlotWidget::paletteColor(0),
                       PlotWidget::Style::Points, 1, true};
  nfEnhPlot_->addSeries(e);
  nfEnhPlot_->setTitle("Field enhancement at the molecules");
  nfEnhPlot_->setLabels("gap to the outer NP surface (incl. ligand shell), nm", "⟨|E/E₀|²⟩");
  nfEnhPlot_->autoscale();
  const int qk = nfQuantity_->currentIndex();
  const bool hasQy = t.column("quantum_yield") >= 0;
  const QString col = qk == 1 && hasQy ? "rad_enhancement" : qk == 2 && hasQy ? "quantum_yield" : "decay_to_np_ps";
  PlotWidget::Series d{qk == 0 ? "decay into NP" : col, t.col("gap_nm"), t.col(col), {}, {}, PlotWidget::paletteColor(1),
                       PlotWidget::Style::Points, 1, true};
  nfDecayPlot_->addSeries(d);
  if (qk == 2 && hasQy && t.rows() > 0) {
    const double q0 = t.col("quantum_yield_free").value(0);
    nfDecayPlot_->addSeries({"without NP", {t.col("gap_nm").first(), t.col("gap_nm").last()}, {q0, q0}, {}, {},
                             PlotWidget::paletteColor(2), PlotWidget::Style::Dashed, 1, true});
  }
  if (qk > 0 && !hasQy) nfDecayPlot_->setTitle("Older run: no Purcell/QY columns — showing the quenching rate");
  else if (qk == 1) nfDecayPlot_->setTitle("Radiative rate enhancement |L_iᵀμ̂_i|² (emission through the NP)");
  else if (qk == 2) nfDecayPlot_->setTitle("Quantum yield of an isolated molecule at each position");
  else nfDecayPlot_->setTitle("Molecular quenching rate (−2 Im W_ii)");
  nfDecayPlot_->setLabels("gap to the outer NP surface (incl. ligand shell), nm",
                          col == "rad_enhancement" ? "Γ_rad/Γ₀" : col == "quantum_yield" ? "QY" : "Γ_NP, ps⁻¹");
  nfDecayPlot_->autoscale();
  // 3D с раскраской по усилению поля
  const auto& g = r->geometry;
  QVector<QVector3D> pos, mu;
  for (int i = 0; i < g.rows(); ++i) {
    pos << QVector3D(g.cols[0][i], g.cols[1][i], g.cols[2][i]);
    mu << QVector3D(g.cols[3][i], g.cols[4][i], g.cols[5][i]);
  }
  nfGeom_->setAggregate(pos, mu);
  const double axy = jsonGet(cfg_, "nanoparticle.semi_axis_xy_nm").toDouble(15);
  const double az = jsonGet(cfg_, "nanoparticle.shape").toString() == "sphere" ? axy : jsonGet(cfg_, "nanoparticle.semi_axis_z_nm").toDouble(axy);
  nfGeom_->setParticle(true, axy, az);
  nfGeom_->setValues(t.col("field_enhancement"), "|E/E₀|²");
}

void MainWindow::selectResultTab(int i) { resultTabs_->setCurrentIndex(i); }

void MainWindow::setFitOptions(bool offset, bool shift) {
  fitOffset_->setChecked(offset);
  fitShift_->setChecked(shift);
}
