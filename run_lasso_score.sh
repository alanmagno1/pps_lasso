#!/bin/bash
#$ -N lasso_score
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=15:00:00
#$ -l h_vmem=15G
#$ -o logs/
#$ -e logs/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  # Set LASSO_CONDA_ENV="" to skip conda activation.
  if [ "${LASSO_CONDA_ENV+r}" ]; then
    if [ -n "$LASSO_CONDA_ENV" ]; then
      source "$HOME/miniconda3/etc/profile.d/conda.sh"
      conda activate "$LASSO_CONDA_ENV"
    fi
  else
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate r-gcc
  fi
fi

Rscript "$SCRIPT_DIR/lasso_score.R" "$@"

