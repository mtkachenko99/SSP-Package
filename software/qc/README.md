# Quantum-chemistry helpers

Optional tools used to produce transition charges (TrEsp) for the exciton model. They are **not** needed to build or
run `ssp_aggregate` / `ssp_gui`.

| File | Purpose |
|---|---|
| `c8s3_tresp.py` | PySCF TD-DFT (CAM-B3LYP/6-31G*, TDA) for the C8S3 chromophore → ESP-fitted transition charges |
| `nwchem.sh` | wrapper that starts NWChem from a local conda-forge environment with the basis-library paths set |

Environments are not stored in git. To recreate them:

```bash
python3 -m venv software/qc/venv && software/qc/venv/bin/pip install pyscf matplotlib
```

```bash
# NWChem without root: micromamba + conda-forge
curl -sSfL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C software/qc bin/micromamba
MAMBA_ROOT_PREFIX=software/qc/mamba software/qc/bin/micromamba create -y -p software/qc/nwchem-env -c conda-forge nwchem
```
