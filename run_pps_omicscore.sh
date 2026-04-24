#!/bin/bash
#$ -N pps_omicscore
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=15:00:00
#$ -t 1-30
#$ -l h_vmem=15G
#$ -o logs/
#$ -e logs/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PPS_IDX="${SGE_TASK_ID:-${1:-}}"

if [ -z "$PPS_IDX" ]; then
  echo "Usage: bash run_pps_omicscore.sh <idx>" >&2
  echo "For SGE arrays, submit with: qsub run_pps_omicscore.sh" >&2
  exit 2
fi

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  # Set PPS_CONDA_ENV="" to skip conda activation.
  if [ "${PPS_CONDA_ENV+r}" ]; then
    if [ -n "$PPS_CONDA_ENV" ]; then
      source "$HOME/miniconda3/etc/profile.d/conda.sh"
      conda activate "$PPS_CONDA_ENV"
    fi
  else
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate r-gcc
  fi
fi

DATA_DIR="${PPS_DATA_DIR:-data}"
NCORES="${PPS_NCORES:-${NSLOTS:-4}}"
SEED="${PPS_SEED:-1427}"

Rscript "$SCRIPT_DIR/pps_omicscore.R" \
  --idx "$PPS_IDX" \
  --data-dir "$DATA_DIR" \
  --prot-file "${PPS_PROT_FILE:-$DATA_DIR/proteins/proteins_visit_0.fst}" \
  --metab-file "${PPS_METAB_FILE:-$DATA_DIR/metabolites/nmr_threephases.csv}" \
  --clusters-file "${PPS_CLUSTERS_FILE:-$DATA_DIR/clusters.csv}" \
  --ids-file "${PPS_IDS_FILE:-$DATA_DIR/id_lists.csv}" \
  --out-dir "${PPS_OUT_DIR:-results/omics_score/pps/disease_free}" \
  --ncores "$NCORES" \
  --seed "$SEED"
