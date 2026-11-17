#pragma once
// Конфигурация запуска ssp_aggregate (JSON). Неизвестные ключи — ошибка (защита от опечаток).

#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/geometry.hpp"
#include "ssp/excitons/heom.hpp"
#include "ssp/excitons/plexciton.hpp"
#include "ssp/excitons/vibronic.hpp"
#include "ssp/plasmonics/material.hpp"

#include <string>
#include <vector>

namespace ssp::app {

struct Placement {
  excitons::Vec3 shift_nm = excitons::Vec3::Zero();
  double euler_deg[3] = {0.0, 0.0, 0.0};
  bool center = false;  // сначала перенести центр агрегата в начало координат
};

struct Configuration {
  std::string name, type;
  excitons::ParamMap params;
  Placement placement;
};

struct DynamicsConfig {
  bool enabled = true;
  double t_max_ps = 20.0;
  int n_times = 200;
  std::string initial = "brightest";  // brightest | highest | pump
};

struct NanoparticleConfig {
  bool enabled = false;
  excitons::NanoparticleSpec spec;
  std::string material_file;          // путь к YAML refractiveindex.info (относительно файла конфига)
  bool use_drude = false;
  plasmonics::DrudeParams drude;
  double kreibig_A = 0.0;             // 0 — без размерного затухания
  double fermi_velocity_m_s = 1.40e6;
  double kreibig_omega_p_eV = 9.0;
  double kreibig_gamma_bulk_eV = 0.07;
  int realizations = 30;
  double wmin_cm = 0.0, wmax_cm = 0.0;  // 0 — как у экситонных спектров
  std::size_t nw = 600;
  double homogeneous_hwhm_cm = -1.0;  // < 0 — из Редфилда
  double near_field_cm = 0.0;         // 0 — энергия мономера
  bool emission = true;               // ФЛ гибрида, квантовый выход, Парселл
};

struct HolsteinConfig {
  bool enabled = false;
  excitons::HolsteinParams params;
  int realizations = 20;
};

struct HeomConfig {
  bool enabled = false;
  excitons::HeomParams params;
  int max_sites = 10;
  bool absorption = true;
  bool populations = true;
  std::string initial = "brightest";  // brightest | highest | site:<k>
  int n_out = 200;
};

struct RunConfig {
  std::string config_path, config_dir, output_dir;
  std::vector<Configuration> configurations;
  excitons::ExcitonModel model;
  excitons::SpectrumGrid grid;
  excitons::EnsembleSettings ensemble;
  DynamicsConfig dynamics;
  NanoparticleConfig nanoparticle;
  HeomConfig heom;
  HolsteinConfig holstein;
};

RunConfig load_config(const std::string& path);
// Разбор JSON-текста; относительные пути (material_file) — относительно base_dir.
RunConfig load_config_text(const std::string& json, const std::string& base_dir, const std::string& path = {});

}  // namespace ssp::app
