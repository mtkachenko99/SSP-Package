SSP Nano — modelling of exciton spectra of molecular aggregates and plasmon–exciton hybrids
==========================================================================================

Programs
  ssp_gui        graphical interface (start here): parameters, run, plots, comparison with experiment
  ssp_aggregate  command line:  ssp_aggregate examples/notebook_demo.json
                 monomer calibration: ssp_aggregate examples/monomer_calibration.json
                 transition charges from geometry: ssp_aggregate --huckel-tresp molecule.xyz charges.txt
  ssp_benchmark  validation against exact Mie theory / datasheet values

Examples (folder examples/): notebook_demo, heom_benchmark, shapes_with_gold, shell_scan,
  dense_ligand_shell, monomer_calibration (synthetic data). Paths inside a config are relative to the config file;
  results are written to output_dir (relative to the config).

Gold optical constants: examples that use "material_file" expect the refractiveindex.info database at the path
  given in the config; replace it with "drude": {"eps_inf": 9.5, "omega_p_eV": 9.0, "gamma_eV": 0.07}
  or point it to a downloaded YAML file (https://refractiveindex.info, e.g. Au/McPeak).

Documentation: see the accompanying Word document (model, methods, validation, examples).
