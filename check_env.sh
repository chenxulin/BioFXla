#!/usr/bin/env bash
set -e

MAIN_ENV="aidd"

# Ensure conda/mamba is available even in non-interactive shells.
if ! command -v mamba >/dev/null 2>&1 || ! command -v conda >/dev/null 2>&1; then
	for CONDA_ROOT in \
		"${CONDA_PREFIX:-}" \
		"/data/cxl/miniforge3" \
		"$HOME/miniconda3" \
		"$HOME/anaconda3"
	do
		if [[ -n "$CONDA_ROOT" && -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
			# shellcheck disable=SC1090
			source "$CONDA_ROOT/etc/profile.d/conda.sh"
			break
		fi
	done
fi

PKG_MGR=""
if command -v mamba >/dev/null 2>&1; then
	PKG_MGR="mamba"
elif command -v conda >/dev/null 2>&1; then
	PKG_MGR="conda"
else
	echo "ERROR: neither mamba nor conda is available in PATH."
	echo "Hint: initialize conda in your shell or install mamba/conda first."
	exit 127
fi

echo "== mamba / uv =="
"$PKG_MGR" --version
uv --version

echo "== ${MAIN_ENV} =="
"$PKG_MGR" run -n "${MAIN_ENV}" python -c "from rdkit import Chem; print('RDKit ok')"
"$PKG_MGR" run -n "${MAIN_ENV}" python -c "import openmm; print('OpenMM ok')"
"$PKG_MGR" run -n "${MAIN_ENV}" python -m openmm.testInstallation
"$PKG_MGR" run -n "${MAIN_ENV}" python -c "import MDAnalysis, mdtraj; print('trajectory tools ok')"
"$PKG_MGR" run -n "${MAIN_ENV}" python -c "import pandas, numpy, scipy, sklearn; print('data stack ok')"
"$PKG_MGR" run -n "${MAIN_ENV}" python -c "import openff.toolkit; print('OpenFF ok')"
"$PKG_MGR" run -n "${MAIN_ENV}" obabel -V
"$PKG_MGR" run -n "${MAIN_ENV}" admet_predict --help >/dev/null

echo "== gnina =="
which gnina
gnina --help >/dev/null

echo "== openfe =="
"$PKG_MGR" run -n openfe openfe --help >/dev/null
"$PKG_MGR" run -n openfe python -c "print('OpenFe ok')"

echo "== reinvent4 =="
cd /data/cxl/BioFXla/tools/REINVENT4
uv run reinvent --help >/dev/null
uv run python -c "import torch; print('torch:', torch.__version__, 'cuda:', torch.cuda.is_available())"

echo "ALL CHECKS PASSED"