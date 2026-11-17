#!/usr/bin/env bash
# NWChem из conda-forge (software/qc/nwchem-env): задаёт пути к библиотекам базисов и запускает расчёт.
# Использование: software/qc/nwchem.sh input.nw [число MPI-процессов, по умолчанию 1] > output.out
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env="$here/nwchem-env"
export NWCHEM_BASIS_LIBRARY="$env/share/nwchem/libraries/"
export NWCHEM_NWPW_LIBRARY="$env/share/nwchem/libraryps/"
export OMP_NUM_THREADS=1
np="${2:-1}"
if [ "$np" -gt 1 ]; then
  exec "$env/bin/mpirun" --oversubscribe -np "$np" "$env/bin/nwchem" "$1"
else
  exec "$env/bin/nwchem" "$1"
fi
