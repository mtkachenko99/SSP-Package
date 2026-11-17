#include "ssp/excitons/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

namespace ssp::excitons {
namespace {

constexpr double kPi = std::numbers::pi;
double deg(double d) { return d * kPi / 180.0; }

Vec3 normalized(const Vec3& v) {
  const double n = v.norm();
  return n > 0.0 ? Vec3(v / n) : v;
}

// Доступ к параметрам с проверкой: каждый прочитанный ключ отмечается, лишние ключи — ошибка.
class Params {
 public:
  Params(const std::string& type, const ParamMap& user, ParamMap defaults)
      : type_(type), values_(std::move(defaults)) {
    for (const auto& [k, v] : user) {
      if (!values_.contains(k)) throw std::invalid_argument("geometry '" + type + "': unknown parameter '" + k + "'");
      values_[k] = v;
    }
  }
  double operator[](const std::string& k) const { return values_.at(k); }
  int integer(const std::string& k) const {
    const double v = values_.at(k);
    if (v < 0 || std::abs(v - std::round(v)) > 1e-9)
      throw std::invalid_argument("geometry '" + type_ + "': parameter '" + k + "' must be a non-negative integer");
    return static_cast<int>(std::lround(v));
  }
  bool flag(const std::string& k) const { return values_.at(k) != 0.0; }

 private:
  std::string type_;
  ParamMap values_;
};

const std::map<std::string, ParamMap>& defaults_table() {
  static const std::map<std::string, ParamMap> t = {
      {"linear_chain", {{"N", 20}, {"spacing", 1.0}, {"dipole_angle_deg", 20.0}}},
      {"ring", {{"N", 20}, {"radius", 4.0}, {"tangential", 1}, {"tilt_deg", 20.0}}},
      {"elliptic_ring", {{"N", 20}, {"a", 5.0}, {"b", 3.0}, {"tangential", 1}, {"tilt_deg", 20.0}}},
      {"cylinder",
       {{"n_rings", 8}, {"molecules_per_ring", 8}, {"radius", 3.5}, {"pitch", 1.1},
        {"helical_twist_deg", 18.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}}},
      {"double_wall_tube",
       {{"n_rings", 8}, {"molecules_per_ring_inner", 8}, {"molecules_per_ring_outer", 10},
        {"radius_inner", 3.0}, {"radius_outer", 4.0}, {"pitch", 1.1}, {"helical_twist_deg", 18.0},
        {"wall_offset_deg", 10.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}}},
      {"elliptic_tube",
       {{"n_rings", 8}, {"molecules_per_ring", 10}, {"a", 5.0}, {"b", 3.0}, {"pitch", 1.1},
        {"helical_twist_deg", 18.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}}},
      {"helix",
       {{"N", 40}, {"radius", 3.0}, {"pitch_per_turn", 5.0}, {"n_turns", 4.0}, {"dipole_tilt_deg", 20.0},
        {"dipole_follow_tangent", 1}}},
      {"ribbon_2d", {{"n_x", 10}, {"n_y", 4}, {"spacing_x", 1.0}, {"spacing_y", 1.2}, {"dipole_angle_deg", 20.0}}},
      {"brick_layer_sheet",
       {{"n_x", 10}, {"n_y", 4}, {"spacing_x", 1.0}, {"spacing_y", 1.2}, {"slip_x", 0.5},
        {"dipole_angle_deg", 20.0}}},
      {"stacked_rings",
       {{"n_rings", 6}, {"molecules_per_ring", 12}, {"radius", 4.0}, {"ring_spacing", 1.0},
        {"ring_rotation_deg", 10.0}, {"tangential", 1}, {"tilt_deg", 20.0}}},
      {"tube_radial_disorder",
       {{"n_rings", 8}, {"molecules_per_ring", 8}, {"radius", 3.5}, {"radial_sigma", 0.2}, {"pitch", 1.1},
        {"helical_twist_deg", 18.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}, {"seed", 0}}},
      {"tube_variable_pitch",
       {{"n_rings", 8}, {"molecules_per_ring", 8}, {"radius", 3.5}, {"pitch0", 0.9}, {"pitch1", 1.4},
        {"helical_twist_deg", 18.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}}},
      {"particle_shell",
       {{"a_xy", 15.0}, {"a_z", 15.0}, {"gap", 1.5}, {"N", 300}, {"density_per_nm2", 0.0}, {"theta_max_deg", 180.0},
        {"dipole_mode", 0}, {"tilt_deg", 0.0}, {"domains", 0}, {"domain_jitter_deg", 10.0}, {"seed", 0}}},
      {"brickwork_tube",
       {{"walls", 2}, {"radius_inner", 3.55}, {"radius_outer", 4.75}, {"length", 25.0}, {"stack_spacing", 0.5},
        {"row_spacing", 2.08}, {"beta_inner_deg", 44.3}, {"beta_outer_deg", 34.5}, {"slip", 0.9}, {"seamless", 1}}},
      {"tapered_tube",
       {{"n_rings", 8}, {"molecules_per_ring", 8}, {"radius0", 3.0}, {"radius1", 5.0}, {"pitch", 1.1},
        {"helical_twist_deg", 18.0}, {"dipole_tilt_deg", 20.0}, {"dipole_azimuthal", 1}}},
  };
  return t;
}

// Диполь: cos(tilt)·in_plane + sin(tilt)·ẑ, нормированный.
Vec3 tilted(const Vec3& in_plane, double tilt) {
  return normalized(std::cos(tilt) * in_plane + std::sin(tilt) * Vec3::UnitZ());
}

Aggregate linear_chain(const Params& p) {
  const int N = p.integer("N");
  const double th = deg(p["dipole_angle_deg"]);
  const Vec3 mu = normalized(Vec3(std::cos(th), 0.0, std::sin(th)));
  Aggregate a;
  for (int i = 0; i < N; ++i) {
    a.positions.emplace_back(i * p["spacing"], 0.0, 0.0);
    a.dipoles.push_back(mu);
  }
  return a;
}

Aggregate ring_impl(int N, double radius, bool tangential, double tilt_deg) {
  Aggregate a;
  const double tilt = deg(tilt_deg);
  for (int n = 0; n < N; ++n) {
    const double phi = 2.0 * kPi * n / N;
    a.positions.emplace_back(radius * std::cos(phi), radius * std::sin(phi), 0.0);
    const Vec3 radial(std::cos(phi), std::sin(phi), 0.0);
    const Vec3 tangent(-std::sin(phi), std::cos(phi), 0.0);
    a.dipoles.push_back(tilted(tangential ? tangent : radial, tilt));
  }
  return a;
}

Aggregate ring(const Params& p) {
  return ring_impl(p.integer("N"), p["radius"], p.flag("tangential"), p["tilt_deg"]);
}

Aggregate elliptic_ring(const Params& p) {
  const int N = p.integer("N");
  const double A = p["a"], B = p["b"], tilt = deg(p["tilt_deg"]);
  Aggregate a;
  for (int n = 0; n < N; ++n) {
    const double phi = 2.0 * kPi * n / N;
    const double x = A * std::cos(phi), y = B * std::sin(phi);
    a.positions.emplace_back(x, y, 0.0);
    const Vec3 radial = normalized(Vec3(x / std::max(A, 1e-12), y / std::max(B, 1e-12), 0.0));
    const Vec3 tangent = normalized(Vec3(-A * std::sin(phi), B * std::cos(phi), 0.0));
    a.dipoles.push_back(tilted(p.flag("tangential") ? tangent : radial, tilt));
  }
  return a;
}

// Общая трубка: кольца на высотах z_m, радиус r_m (функция номера кольца), поворот m·twist.
template <class RadiusFn, class ZFn>
Aggregate tube(int n_rings, int per_ring, RadiusFn radius_of, ZFn z_of, double twist_deg, double tilt_deg,
               bool azimuthal) {
  Aggregate a;
  const double twist = deg(twist_deg), tilt = deg(tilt_deg);
  for (int m = 0; m < n_rings; ++m) {
    for (int n = 0; n < per_ring; ++n) {
      const double phi = 2.0 * kPi * n / per_ring + m * twist;
      const double r = radius_of(m, n);
      a.positions.emplace_back(r * std::cos(phi), r * std::sin(phi), z_of(m));
      const Vec3 radial(std::cos(phi), std::sin(phi), 0.0);
      const Vec3 tangent(-std::sin(phi), std::cos(phi), 0.0);
      a.dipoles.push_back(tilted(azimuthal ? tangent : radial, tilt));
    }
  }
  return a;
}

Aggregate cylinder_impl(int n_rings, int per_ring, double radius, double pitch, double twist_deg, double tilt_deg,
                        bool azimuthal) {
  return tube(
      n_rings, per_ring, [&](int, int) { return radius; }, [&](int m) { return m * pitch; }, twist_deg, tilt_deg,
      azimuthal);
}

Aggregate cylinder(const Params& p) {
  return cylinder_impl(p.integer("n_rings"), p.integer("molecules_per_ring"), p["radius"], p["pitch"],
                       p["helical_twist_deg"], p["dipole_tilt_deg"], p.flag("dipole_azimuthal"));
}

Aggregate double_wall_tube(const Params& p) {
  Aggregate inner = cylinder_impl(p.integer("n_rings"), p.integer("molecules_per_ring_inner"), p["radius_inner"],
                                  p["pitch"], p["helical_twist_deg"], p["dipole_tilt_deg"],
                                  p.flag("dipole_azimuthal"));
  const Aggregate outer = cylinder_impl(p.integer("n_rings"), p.integer("molecules_per_ring_outer"),
                                        p["radius_outer"], p["pitch"], p["helical_twist_deg"] + p["wall_offset_deg"],
                                        p["dipole_tilt_deg"], p.flag("dipole_azimuthal"));
  inner.positions.insert(inner.positions.end(), outer.positions.begin(), outer.positions.end());
  inner.dipoles.insert(inner.dipoles.end(), outer.dipoles.begin(), outer.dipoles.end());
  return inner;
}

Aggregate elliptic_tube(const Params& p) {
  const int n_rings = p.integer("n_rings"), per_ring = p.integer("molecules_per_ring");
  const double A = p["a"], B = p["b"], twist = deg(p["helical_twist_deg"]), tilt = deg(p["dipole_tilt_deg"]);
  Aggregate a;
  for (int m = 0; m < n_rings; ++m) {
    for (int n = 0; n < per_ring; ++n) {
      const double phi = 2.0 * kPi * n / per_ring + m * twist;
      const double x = A * std::cos(phi), y = B * std::sin(phi);
      a.positions.emplace_back(x, y, m * p["pitch"]);
      const Vec3 radial = normalized(Vec3(x / std::max(A, 1e-12), y / std::max(B, 1e-12), 0.0));
      const Vec3 tangent = normalized(Vec3(-A * std::sin(phi), B * std::cos(phi), 0.0));
      a.dipoles.push_back(tilted(p.flag("dipole_azimuthal") ? tangent : radial, tilt));
    }
  }
  return a;
}

Aggregate helix(const Params& p) {
  const int N = p.integer("N");
  const double R = p["radius"], pitch = p["pitch_per_turn"], turns = p["n_turns"], tilt = deg(p["dipole_tilt_deg"]);
  Aggregate a;
  for (int i = 0; i < N; ++i) {
    const double t = 2.0 * kPi * turns * i / std::max(N - 1, 1);
    a.positions.emplace_back(R * std::cos(t), R * std::sin(t), pitch * t / (2.0 * kPi));
    const Vec3 tangent = normalized(Vec3(-R * std::sin(t), R * std::cos(t), pitch / (2.0 * kPi)));
    const Vec3 radial = normalized(Vec3(std::cos(t), std::sin(t), 0.0));
    a.dipoles.push_back(tilted(p.flag("dipole_follow_tangent") ? tangent : radial, tilt));
  }
  return a;
}

Aggregate sheet(const Params& p, double slip) {
  const int nx = p.integer("n_x"), ny = p.integer("n_y");
  const double th = deg(p["dipole_angle_deg"]);
  const Vec3 mu = normalized(Vec3(std::cos(th), 0.0, std::sin(th)));
  Aggregate a;
  for (int iy = 0; iy < ny; ++iy) {
    const double shift = (iy % 2 == 1) ? slip : 0.0;
    for (int ix = 0; ix < nx; ++ix) {
      a.positions.emplace_back(ix * p["spacing_x"] + shift, iy * p["spacing_y"], 0.0);
      a.dipoles.push_back(mu);
    }
  }
  return a;
}

Aggregate stacked_rings(const Params& p) {
  Aggregate a;
  const double rot = deg(p["ring_rotation_deg"]);
  for (int m = 0; m < p.integer("n_rings"); ++m) {
    const Aggregate base = ring_impl(p.integer("molecules_per_ring"), p["radius"], p.flag("tangential"), p["tilt_deg"]);
    const Eigen::Matrix3d R = Eigen::AngleAxisd(m * rot, Vec3::UnitZ()).toRotationMatrix();
    for (std::size_t k = 0; k < base.size(); ++k) {
      Vec3 r = R * base.positions[k];
      r.z() += m * p["ring_spacing"];
      a.positions.push_back(r);
      a.dipoles.push_back(R * base.dipoles[k]);
    }
  }
  return a;
}

Aggregate tube_radial_disorder(const Params& p) {
  std::mt19937_64 rng(static_cast<std::uint64_t>(p.integer("seed")));
  std::normal_distribution<double> gauss(0.0, p["radial_sigma"]);
  const double R = p["radius"], pitch = p["pitch"];
  return tube(
      p.integer("n_rings"), p.integer("molecules_per_ring"), [&](int, int) { return R + gauss(rng); },
      [&](int m) { return m * pitch; }, p["helical_twist_deg"], p["dipole_tilt_deg"], p.flag("dipole_azimuthal"));
}

Aggregate tube_variable_pitch(const Params& p) {
  const int n_rings = p.integer("n_rings");
  // Как в ноутбуке: шаги linspace(pitch0, pitch1, max(n_rings-1, 1)), z_m — накопленная сумма.
  const int n_steps = std::max(n_rings - 1, 1);
  std::vector<double> z(n_rings, 0.0);
  for (int m = 1; m < n_rings; ++m) {
    const double t = n_steps > 1 ? double(m - 1) / (n_steps - 1) : 0.0;
    z[m] = z[m - 1] + p["pitch0"] + t * (p["pitch1"] - p["pitch0"]);
  }
  const double R = p["radius"];
  return tube(
      n_rings, p.integer("molecules_per_ring"), [&](int, int) { return R; }, [&](int m) { return z[m]; },
      p["helical_twist_deg"], p["dipole_tilt_deg"], p.flag("dipole_azimuthal"));
}

// Слой молекул на поверхности сфероида (оси a_xy, a_xy, a_z) на расстоянии gap по нормали.
// Точки распределены равномерно по площади (последовательность Фибоначчи + обращение функции распределения
// площади по полярному углу), покрывают полярную «шапочку» θ ≤ θ_max (180° — вся поверхность).
// Диполи в касательной плоскости: dipole_mode 0 — случайная ориентация, 1 — по азимуту φ̂, 2 — по меридиану θ̂,
// 3 — по нормали; tilt_deg — наклон к нормали.
Aggregate particle_shell(const Params& p) {
  const double a = p["a_xy"], c = p["a_z"], gap = p["gap"];
  const double tmax = std::clamp(deg(p["theta_max_deg"]), 1e-6, kPi);
  const int mode = p.integer("dipole_mode");
  const int n_dom = p.integer("domains");
  // Число молекул: N или плотность × площадь покрываемой части смещённой поверхности (численно)
  int N = p.integer("N");
  if (p["density_per_nm2"] > 0.0) {
    const double A = a + gap, C = c + gap;
    double area = 0.0;
    const int M0 = 4000;
    for (int k = 0; k < M0; ++k) {
      const double t = tmax * (k + 0.5) / M0;
      area += 2.0 * kPi * A * std::sin(t) * std::sqrt(C * C * std::sin(t) * std::sin(t) + A * A * std::cos(t) * std::cos(t)) * tmax / M0;
    }
    N = std::max(1, static_cast<int>(std::lround(p["density_per_nm2"] * area)));
  }
  if (a <= 0 || c <= 0) throw std::invalid_argument("particle_shell: a_xy and a_z must be > 0");
  if (gap <= 0) throw std::invalid_argument("particle_shell: gap must be > 0");
  if (mode > 3) throw std::invalid_argument("particle_shell: dipole_mode must be 0..3");
  // Функция распределения площади по θ: dA ∝ sinθ·sqrt(c² sin²θ + a² cos²θ) dθ
  const int M = 4000;
  std::vector<double> th(M + 1), cdf(M + 1, 0.0);
  for (int k = 0; k <= M; ++k) th[k] = tmax * k / M;
  for (int k = 1; k <= M; ++k) {
    auto dens = [&](double t) { return std::sin(t) * std::sqrt(c * c * std::sin(t) * std::sin(t) + a * a * std::cos(t) * std::cos(t)); };
    cdf[k] = cdf[k - 1] + 0.5 * (dens(th[k - 1]) + dens(th[k])) * (th[k] - th[k - 1]);
  }
  for (auto& v : cdf) v /= cdf.back();
  std::mt19937_64 rng(static_cast<std::uint64_t>(p.integer("seed")));
  std::uniform_real_distribution<double> uni(0.0, 2.0 * kPi);
  std::normal_distribution<double> jit(0.0, deg(p["domain_jitter_deg"]));
  const double golden = kPi * (3.0 - std::sqrt(5.0));
  const double tilt = deg(p["tilt_deg"]);
  // Домены: центры — точки Фибоначчи на единичной сфере (направления), директор — случайный касательный вектор
  std::vector<Vec3> dom_c, dom_d;
  for (int k = 0; k < n_dom; ++k) {
    const double z = 1.0 - 2.0 * (k + 0.5) / n_dom, rr = std::sqrt(std::max(0.0, 1.0 - z * z)), ph = k * golden;
    const Vec3 u(rr * std::cos(ph), rr * std::sin(ph), z);
    const Vec3 e1 = (std::abs(u.z()) < 0.9 ? Vec3::UnitZ() : Vec3::UnitX()).cross(u).normalized(), e2 = u.cross(e1);
    const double psi = uni(rng);
    dom_c.push_back(u);
    dom_d.push_back(std::cos(psi) * e1 + std::sin(psi) * e2);
  }
  Aggregate agg;
  for (int k = 0; k < N; ++k) {
    const double u = (k + 0.5) / N;
    const auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
    const std::size_t j = std::clamp<std::size_t>(static_cast<std::size_t>(it - cdf.begin()), 1, M);
    const double f = (u - cdf[j - 1]) / std::max(cdf[j] - cdf[j - 1], 1e-300);
    const double theta = th[j - 1] + f * (th[j] - th[j - 1]);
    const double phi = k * golden;
    const Vec3 s(a * std::sin(theta) * std::cos(phi), a * std::sin(theta) * std::sin(phi), c * std::cos(theta));
    const Vec3 n = Vec3(s.x() / (a * a), s.y() / (a * a), s.z() / (c * c)).normalized();
    const Vec3 ephi(-std::sin(phi), std::cos(phi), 0.0);
    const Vec3 etheta = ephi.cross(n).normalized();
    Vec3 d;
    if (n_dom > 0 && mode != 3) {
      // ближайший домен (по направлению), его директор — в локальную касательную плоскость, плюс разброс
      const Vec3 u = s.normalized();
      int best = 0;
      for (int q = 1; q < n_dom; ++q)
        if (u.dot(dom_c[q]) > u.dot(dom_c[best])) best = q;
      Vec3 t = dom_d[best] - dom_d[best].dot(n) * n;
      if (t.norm() < 1e-6) t = etheta;
      t.normalize();
      const double ang = jit(rng);
      d = std::cos(ang) * t + std::sin(ang) * n.cross(t);
    } else if (mode == 0) {
      const double psi = uni(rng);
      d = std::cos(psi) * etheta + std::sin(psi) * ephi;
    } else if (mode == 1) d = ephi;
    else if (mode == 2) d = etheta;
    else d = n;
    if (mode != 3) d = std::cos(tilt) * d + std::sin(tilt) * n;
    agg.positions.push_back(s + gap * n);
    agg.dipoles.push_back(normalized(d));
    agg.normals.push_back(n);
  }
  return agg;
}

Aggregate tapered_tube(const Params& p) {
  const int n_rings = p.integer("n_rings");
  const double r0 = p["radius0"], r1 = p["radius1"], pitch = p["pitch"];
  return tube(
      n_rings, p.integer("molecules_per_ring"),
      [&](int m, int) { return n_rings > 1 ? r0 + (r1 - r0) * m / (n_rings - 1) : r0; },
      [&](int m) { return m * pitch; }, p["helical_twist_deg"], p["dipole_tilt_deg"], p.flag("dipole_azimuthal"));
}

// Трубчатый агрегат из «кирпичной» решётки, свёрнутой на цилиндр (модель C8S3: Megow et al., Didraga et al.):
// молекулы длиной row_spacing стоят стопками с шагом stack_spacing; длинная ось (диполь) — под углом β к оси
// трубки, касательно к поверхности; соседние в стопке сдвинуты вдоль длинной оси на ±slip/2 (две молекулы в
// ячейке; slip = 0 — H-тип). walls = 1 — только внутренняя стенка. seamless = 1: вектор свёртки — вектор
// решётки n₁·(2d·ê_стопка) + n₂·(w·ê_ось), поэтому β и R подстраиваются (обычно на доли градуса и нм) — без шва;
// иначе на шве неизбежны слишком близкие пары (они дают ложные ловушки), и такие молекулы (< 0.35 нм) отбрасываются.
Aggregate brickwork_tube(const Params& p) {
  const int walls = p.integer("walls");
  if (walls < 1 || walls > 2) throw std::invalid_argument("brickwork_tube: walls must be 1 or 2");
  const double L = p["length"], d = p["stack_spacing"], w = p["row_spacing"], slip = p["slip"];
  if (L <= 0 || d <= 0 || w <= 0) throw std::invalid_argument("brickwork_tube: length and spacings must be > 0");
  Aggregate a;
  for (int wall = 0; wall < walls; ++wall) {
    double R = wall == 0 ? p["radius_inner"] : p["radius_outer"];
    double beta = deg(wall == 0 ? p["beta_inner_deg"] : p["beta_outer_deg"]);
    if (R <= 0) throw std::invalid_argument("brickwork_tube: radii must be > 0");
    if (p.flag("seamless")) {
      // ближайший к (2πR, β) вектор решётки, лежащий вдоль окружности: tan β' = n₂w/(2n₁d)
      double best = 1e300, Rb = R, bb = beta;
      for (int n1 = 1; n1 <= 400; ++n1) {
        const int n2 = static_cast<int>(std::lround(2.0 * d * n1 * std::tan(beta) / w));
        const double b2 = std::atan2(n2 * w, 2.0 * d * n1);
        const double C2 = 2.0 * d * n1 * std::cos(b2) + n2 * w * std::sin(b2);
        const double cost = std::abs(C2 / (2.0 * kPi) - R) / R + std::abs(b2 - beta) / deg(3.0);
        if (C2 > 0 && cost < best) {
          best = cost;
          Rb = C2 / (2.0 * kPi);
          bb = b2;
        }
      }
      R = Rb;
      beta = bb;
    }
    const double C = 2.0 * kPi * R;
    // единичные векторы в развёртке (s — вдоль окружности, z — вдоль оси)
    const double eds = std::sin(beta), edz = std::cos(beta);  // длинная ось
    const double ess = std::cos(beta), esz = -std::sin(beta);  // направление стопки
    const int ni = static_cast<int>(std::ceil((C + L) / d)) + 2, nj = static_cast<int>(std::ceil((C + L) / w)) + 2;
    const std::size_t first = a.size();
    for (int i = -ni; i <= ni; ++i)
      for (int j = -nj; j <= nj; ++j) {
        const double off = (i % 2 == 0 ? 0.5 : -0.5) * slip;
        const double s = i * d * ess + (j * w + off) * eds;
        const double z = i * d * esz + (j * w + off) * edz;
        if (s < 0 || s >= C || z < 0 || z >= L) continue;
        const double ph = s / R;
        const Vec3 pos(R * std::cos(ph), R * std::sin(ph), z);
        const Vec3 es(-std::sin(ph), std::cos(ph), 0.0), n(std::cos(ph), std::sin(ph), 0.0);
        bool clash = false;
        for (std::size_t k = first; k < a.size() && !clash; ++k) clash = (a.positions[k] - pos).squaredNorm() < 0.35 * 0.35;
        if (clash) continue;
        a.positions.push_back(pos);
        a.dipoles.push_back(normalized(eds * es + edz * Vec3::UnitZ()));
        a.normals.push_back(n);
      }
  }
  return a;
}

}  // namespace

const std::vector<std::string>& geometry_types() {
  static const std::vector<std::string> v = {
      "linear_chain",  "ring",          "elliptic_ring",        "cylinder",           "double_wall_tube",
      "elliptic_tube", "helix",         "ribbon_2d",            "brick_layer_sheet",  "stacked_rings",
      "tube_radial_disorder", "tube_variable_pitch", "tapered_tube", "particle_shell", "brickwork_tube"};
  return v;
}

ParamMap default_geometry_params(const std::string& type) {
  const auto& t = defaults_table();
  const auto it = t.find(type);
  if (it == t.end()) throw std::invalid_argument("unknown geometry type: " + type);
  return it->second;
}

Aggregate build_geometry(const std::string& type, const ParamMap& params) {
  const Params p(type, params, default_geometry_params(type));
  Aggregate a;
  if (type == "linear_chain") a = linear_chain(p);
  else if (type == "ring") a = ring(p);
  else if (type == "elliptic_ring") a = elliptic_ring(p);
  else if (type == "cylinder") a = cylinder(p);
  else if (type == "double_wall_tube") a = double_wall_tube(p);
  else if (type == "elliptic_tube") a = elliptic_tube(p);
  else if (type == "helix") a = helix(p);
  else if (type == "ribbon_2d") a = sheet(p, 0.0);
  else if (type == "brick_layer_sheet") a = sheet(p, p["slip_x"]);
  else if (type == "stacked_rings") a = stacked_rings(p);
  else if (type == "tube_radial_disorder") a = tube_radial_disorder(p);
  else if (type == "tube_variable_pitch") a = tube_variable_pitch(p);
  else if (type == "tapered_tube") a = tapered_tube(p);
  else if (type == "particle_shell") a = particle_shell(p);
  else if (type == "brickwork_tube") a = brickwork_tube(p);
  if (a.size() == 0) throw std::invalid_argument("geometry '" + type + "': empty aggregate");
  return a;
}

void transform(Aggregate& agg, const Vec3& shift_nm, double alpha_deg, double beta_deg, double gamma_deg) {
  const Eigen::Matrix3d R = (Eigen::AngleAxisd(deg(alpha_deg), Vec3::UnitZ()) *
                             Eigen::AngleAxisd(deg(beta_deg), Vec3::UnitY()) *
                             Eigen::AngleAxisd(deg(gamma_deg), Vec3::UnitZ()))
                                .toRotationMatrix();
  for (std::size_t i = 0; i < agg.size(); ++i) {
    agg.positions[i] = R * agg.positions[i] + shift_nm;
    agg.dipoles[i] = R * agg.dipoles[i];
    if (i < agg.normals.size()) agg.normals[i] = R * agg.normals[i];
  }
}

Eigen::Matrix3d Aggregate::frame(std::size_t i) const {
  const Vec3 x = dipoles[i].normalized();
  Vec3 z = i < normals.size() ? normals[i] : (std::abs(x.z()) < 0.9 ? Vec3::UnitZ() : Vec3::UnitX());
  z = (z - z.dot(x) * x).normalized();
  Eigen::Matrix3d F;
  F.col(0) = x;
  F.col(1) = z.cross(x);
  F.col(2) = z;
  return F;
}

Vec3 centroid(const Aggregate& agg) {
  Vec3 c = Vec3::Zero();
  for (const auto& r : agg.positions) c += r;
  return c / static_cast<double>(agg.size());
}

double min_pair_distance(const Aggregate& agg) {
  double dmin = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < agg.size(); ++i)
    for (std::size_t j = i + 1; j < agg.size(); ++j) dmin = std::min(dmin, (agg.positions[i] - agg.positions[j]).norm());
  return dmin;
}

}  // namespace ssp::excitons
