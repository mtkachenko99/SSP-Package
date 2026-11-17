// ssp_aggregate — спектры и кинетика молекулярных агрегатов произвольной формы, опционально рядом с плазмонной НЧ.
// Использование: ssp_aggregate config.json         (расчёт; при разделе "monomer_fit" — калибровка по мономеру)
//                ssp_aggregate --list-geometries

#include "calibrate.hpp"
#include "config.hpp"
#include "runner.hpp"
#include "scan.hpp"

#include "ssp/excitons/geometry.hpp"
#include "ssp/excitons/huckel.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
  using namespace ssp;
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json | --list-geometries\n";
    return 2;
  }
  const std::string arg = argv[1];
  if (arg == "--huckel-tresp") {  // переходные заряды по геометрии: --huckel-tresp in.xyz out.txt [заряд π-системы]
    if (argc < 4) {
      std::cerr << "Usage: " << argv[0] << " --huckel-tresp molecule.xyz charges.txt [pi_charge=1]\n";
      return 2;
    }
    try {
      excitons::HuckelOptions o;
      if (argc > 4) o.pi_charge = std::stoi(argv[4]);
      const auto atoms = excitons::read_xyz(argv[2]);
      const auto h = excitons::huckel_transition_charges(atoms, o);
      excitons::write_transition_charges(argv[3], h, std::string("from ") + argv[2]);
      std::cout << "π-system: " << h.pi_atoms.size() << " atoms, HOMO " << h.homo + 1 << ", HOMO–LUMO gap "
                << h.energies(h.homo + 1) - h.energies(h.homo) << " |β|, |Σqr| = " << h.dipole_e_nm << " e·nm\n";
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << "\n";
      return 1;
    }
    return 0;
  }
  if (arg == "--list-geometries") {
    for (const auto& t : excitons::geometry_types()) {
      std::cout << t << ":";
      for (const auto& [k, v] : excitons::default_geometry_params(t)) std::cout << " " << k << "=" << v;
      std::cout << "\n";
    }
    return 0;
  }
  try {
    std::ifstream f(arg);
    std::stringstream ss;
    ss << f.rdbuf();
    const auto log = [](const std::string& line) { std::cout << line << std::endl; };
    if (app::is_calibration_config(ss.str())) {  // калибровка по мономеру
      app::run_calibration(ss.str(), app::config_dir_of(arg), log);
      return 0;
    }
    const app::RunConfig rc = app::load_config(arg);
    if (app::parse_scan(ss.str()).enabled)
      app::run_scan(ss.str(), rc.config_dir, arg, log);
    else
      app::run(rc, ss.str(), log);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
