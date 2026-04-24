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
exec bash "$SCRIPT_DIR/run_pps_omicscore.sh" "$@"
