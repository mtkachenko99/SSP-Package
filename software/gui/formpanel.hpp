#pragma once
// Форма параметров, построенная по схеме: каждое поле привязано к пути в JSON-конфиге ("exciton.coupling.model").

#include <QJsonObject>
#include <QStringList>
#include <QVariant>
#include <QWidget>

#include <functional>

struct FieldSpec {
  enum class Type { Double, Int, Bool, Choice, String, File };
  QString path, label;
  Type type = Type::Double;
  QVariant def;
  QStringList choices;       // для Choice: значения в JSON
  QStringList choiceLabels;  // подписи (если пусто — сами значения)
  double min = -1e12, max = 1e12;
  int decimals = 3;
  QString suffix, tip;
  QString section;           // заголовок группы (отображается перед полем)
  bool editable = false;     // Choice: можно ввести своё значение
  QString autoText;          // Double: если задан — минимум спинбокса означает «авто» (ключ не пишется в конфиг)
};

class FormPanel : public QWidget {
  Q_OBJECT
 public:
  FormPanel(const QList<FieldSpec>& fields, QWidget* parent = nullptr);
  void load(const QJsonObject& root);  // значения из конфига (или значения по умолчанию)
  void store(QJsonObject& root) const;
  // База для относительных путей файлов
  void setBaseDir(const QString& d) { baseDir_ = d; }

 signals:
  void changed();

 private:
  QList<FieldSpec> fields_;
  QList<QWidget*> editors_;
  QString baseDir_;
  bool loading_ = false;
};

// Чтение/запись значения по пути "a.b.c" во вложенных объектах.
QJsonValue jsonGet(const QJsonObject& root, const QString& path);
void jsonSet(QJsonObject& root, const QString& path, const QJsonValue& v);
void jsonRemove(QJsonObject& root, const QString& path);
