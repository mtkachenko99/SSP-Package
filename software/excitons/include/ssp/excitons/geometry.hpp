#pragma once
// Геометрия молекулярного агрегата: положения молекул (нм) и единичные направления переходных диполей.
// Формы и их параметры совпадают с testcalc.ipynb (build_*), длины трактуются в нм.

#include <Eigen/Dense>

#include <map>
#include <string>
#include <vector>

namespace ssp::excitons {

using Vec3 = Eigen::Vector3d;

struct Aggregate {
  std::vector<Vec3> positions;  // нм
  std::vector<Vec3> dipoles;    // единичные векторы
  std::vector<Vec3> normals;    // необязательно: нормаль к плоскости молекулы (для переходных зарядов)
  std::size_t size() const { return positions.size(); }
  // Система координат молекулы i: x — вдоль диполя, z — нормаль (если не задана — перпендикуляр, выбранный по ẑ/x̂).
  Eigen::Matrix3d frame(std::size_t i) const;
};

// Параметры формы: имя → число (логические параметры: 0/1). Неизвестный ключ — ошибка.
using ParamMap = std::map<std::string, double>;

// Список поддерживаемых типов в порядке ноутбука.
const std::vector<std::string>& geometry_types();

// Построить агрегат заданного типа. Отсутствующие параметры берутся по умолчанию (как в ноутбуке).
// Для форм со случайностью (tube_radial_disorder) используется параметр seed.
Aggregate build_geometry(const std::string& type, const ParamMap& params);

// Параметры по умолчанию для типа (для справки и вывода в отчёт).
ParamMap default_geometry_params(const std::string& type);

// Перенос и поворот агрегата: r' = R·r + shift, μ' = R·μ. Поворот задаётся углами Эйлера z-y-z (градусы).
void transform(Aggregate& agg, const Vec3& shift_nm, double alpha_deg, double beta_deg, double gamma_deg);

// Центр масс положений.
Vec3 centroid(const Aggregate& agg);

// Минимальное расстояние между молекулами (нм).
double min_pair_distance(const Aggregate& agg);

}  // namespace ssp::excitons
