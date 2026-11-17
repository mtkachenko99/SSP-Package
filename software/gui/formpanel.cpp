#include "formpanel.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>

QJsonValue jsonGet(const QJsonObject& root, const QString& path) {
  const QStringList p = path.split('.');
  QJsonObject o = root;
  for (int i = 0; i + 1 < p.size(); ++i) {
    if (!o.contains(p[i]) || !o[p[i]].isObject()) return QJsonValue(QJsonValue::Undefined);
    o = o[p[i]].toObject();
  }
  return o.contains(p.back()) ? o[p.back()] : QJsonValue(QJsonValue::Undefined);
}

static void setRec(QJsonObject& o, const QStringList& p, int i, const QJsonValue& v) {
  if (i + 1 == p.size()) {
    o[p[i]] = v;
    return;
  }
  QJsonObject child = o[p[i]].toObject();
  setRec(child, p, i + 1, v);
  o[p[i]] = child;
}

void jsonSet(QJsonObject& root, const QString& path, const QJsonValue& v) { setRec(root, path.split('.'), 0, v); }

static void removeRec(QJsonObject& o, const QStringList& p, int i) {
  if (i + 1 == p.size()) {
    o.remove(p[i]);
    return;
  }
  if (!o.contains(p[i])) return;
  QJsonObject child = o[p[i]].toObject();
  removeRec(child, p, i + 1);
  o[p[i]] = child;
}

void jsonRemove(QJsonObject& root, const QString& path) { removeRec(root, path.split('.'), 0); }

FormPanel::FormPanel(const QList<FieldSpec>& fields, QWidget* parent) : QWidget(parent), fields_(fields) {
  auto* form = new QFormLayout(this);
  form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  for (const auto& f : fields_) {
    if (!f.section.isEmpty()) {
      auto* h = new QLabel("<b>" + f.section + "</b>");
      h->setContentsMargins(0, 8, 0, 2);
      form->addRow(h);
    }
    QWidget* w = nullptr;
    switch (f.type) {
      case FieldSpec::Type::Double: {
        auto* s = new QDoubleSpinBox;
        s->setRange(f.autoText.isEmpty() ? f.min : f.min - 1.0, f.max);
        if (!f.autoText.isEmpty()) s->setSpecialValueText(f.autoText);
        s->setDecimals(f.decimals);
        s->setSuffix(f.suffix.isEmpty() ? QString() : " " + f.suffix);
        s->setKeyboardTracking(false);
        connect(s, &QDoubleSpinBox::valueChanged, this, [this] { if (!loading_) emit changed(); });
        w = s;
        break;
      }
      case FieldSpec::Type::Int: {
        auto* s = new QSpinBox;
        s->setRange(int(f.min), int(std::min(f.max, 2e9)));
        s->setSuffix(f.suffix.isEmpty() ? QString() : " " + f.suffix);
        s->setKeyboardTracking(false);
        connect(s, &QSpinBox::valueChanged, this, [this] { if (!loading_) emit changed(); });
        w = s;
        break;
      }
      case FieldSpec::Type::Bool: {
        auto* c = new QCheckBox;
        connect(c, &QCheckBox::toggled, this, [this] { if (!loading_) emit changed(); });
        w = c;
        break;
      }
      case FieldSpec::Type::Choice: {
        auto* c = new QComboBox;
        c->setEditable(f.editable);
        for (int i = 0; i < f.choices.size(); ++i)
          c->addItem(i < f.choiceLabels.size() ? f.choiceLabels[i] : f.choices[i], f.choices[i]);
        connect(c, &QComboBox::currentIndexChanged, this, [this] { if (!loading_) emit changed(); });
        if (f.editable)
          connect(c->lineEdit(), &QLineEdit::editingFinished, this, [this] { if (!loading_) emit changed(); });
        w = c;
        break;
      }
      case FieldSpec::Type::String:
      case FieldSpec::Type::File: {
        auto* box = new QWidget;
        auto* h = new QHBoxLayout(box);
        h->setContentsMargins(0, 0, 0, 0);
        auto* e = new QLineEdit;
        e->setObjectName("edit");
        h->addWidget(e);
        connect(e, &QLineEdit::editingFinished, this, [this] { if (!loading_) emit changed(); });
        if (f.type == FieldSpec::Type::File) {
          auto* b = new QPushButton("…");
          b->setFixedWidth(28);
          h->addWidget(b);
          connect(b, &QPushButton::clicked, this, [this, e] {
            const QString path = QFileDialog::getOpenFileName(this, "File", baseDir_);
            if (path.isEmpty()) return;
            e->setText(baseDir_.isEmpty() ? path : QDir(baseDir_).relativeFilePath(path));
            emit changed();
          });
        }
        w = box;
        break;
      }
    }
    w->setToolTip(f.tip);
    editors_ << w;
    auto* lab = new QLabel(f.label);
    lab->setToolTip(f.tip);
    form->addRow(lab, w);
  }
}

void FormPanel::load(const QJsonObject& root) {
  loading_ = true;
  for (int i = 0; i < fields_.size(); ++i) {
    const auto& f = fields_[i];
    QJsonValue v = jsonGet(root, f.path);
    if (v.isUndefined()) v = QJsonValue::fromVariant(f.def);
    QWidget* w = editors_[i];
    switch (f.type) {
      case FieldSpec::Type::Double: {
        auto* sb = static_cast<QDoubleSpinBox*>(w);
        const bool missing = jsonGet(root, f.path).isUndefined();
        sb->setValue(!f.autoText.isEmpty() && missing ? sb->minimum() : v.toDouble());
        break;
      }
      case FieldSpec::Type::Int: static_cast<QSpinBox*>(w)->setValue(int(v.toDouble())); break;
      case FieldSpec::Type::Bool: static_cast<QCheckBox*>(w)->setChecked(v.toBool()); break;
      case FieldSpec::Type::Choice: {
        auto* c = static_cast<QComboBox*>(w);
        const QString sv = v.isDouble() ? QString::number(v.toDouble(), 'g', 12) : v.toString();
        const int k = c->findData(sv);
        if (k >= 0) c->setCurrentIndex(k);
        else if (f.editable) c->setEditText(sv);
        else c->setCurrentIndex(0);
        break;
      }
      case FieldSpec::Type::String:
      case FieldSpec::Type::File: w->findChild<QLineEdit*>("edit")->setText(v.toString()); break;
    }
  }
  loading_ = false;
}

void FormPanel::store(QJsonObject& root) const {
  for (int i = 0; i < fields_.size(); ++i) {
    const auto& f = fields_[i];
    QWidget* w = editors_[i];
    QJsonValue v;
    switch (f.type) {
      case FieldSpec::Type::Double: {
        auto* sb = static_cast<QDoubleSpinBox*>(w);
        if (!f.autoText.isEmpty() && sb->value() == sb->minimum()) {
          jsonRemove(root, f.path);  // «авто»: значение по умолчанию вычисляет ssp_aggregate
          continue;
        }
        v = sb->value();
        break;
      }
      case FieldSpec::Type::Int: v = static_cast<QSpinBox*>(w)->value(); break;
      case FieldSpec::Type::Bool: v = static_cast<QCheckBox*>(w)->isChecked(); break;
      case FieldSpec::Type::Choice: {
        auto* c = static_cast<QComboBox*>(w);
        v = f.editable ? c->currentText().trimmed() : c->currentData().toString();
        bool num = false;  // числовой текст в редактируемом списке — число в JSON
        const double d = v.toString().toDouble(&num);
        if (f.editable && num) v = d;
        break;
      }
      case FieldSpec::Type::String:
      case FieldSpec::Type::File: {
        const QString t = w->findChild<QLineEdit*>("edit")->text();
        if (t.isEmpty()) continue;  // пустые строки не пишем
        v = t;
        break;
      }
    }
    jsonSet(root, f.path, v);
  }
}
