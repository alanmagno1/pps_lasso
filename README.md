# PPS LASSO Omics Score Pipeline

This repository contains an R pipeline for training PPS prediction models from
protein and metabolite features. The main script is `pps_omicscore.R`.

For each PPS column, the pipeline:

- loads PPS clusters, protein features, metabolite features, and ID lists;
- filters features and rows with more than 10% missingness;
- restricts the merged dataset to disease-free IDs;
- splits train/test using `id_lists.csv` (`protein` vs `diet` IDs);
- applies training-derived preprocessing only: imputation, inverse-rank normal
  transform, and z-score scaling;
- trains LASSO (`alpha = 1`) and elastic net (`alpha = 0.5`) models with nested
  cross-validation;
- saves fold IDs, preprocessing parameters, model weights, predictions, and
  performance metrics.

## Repository Structure

```text
.
  pps_omicscore.R          # Main proteins + metabolites pipeline
  run_pps_omicscore.sh     # Local/SGE launcher
  run_pps_omiscore.sh      # Backward-compatible launcher with old spelling
  metscore.R               # Companion metabolite-only pipeline
  collect.R                # Aggregate PPS score files
  collate.R                # Collate scores, metrics, and summaries
  requirements.R           # Installs required R packages
  DESCRIPTION              # Project dependency metadata
  data/README.md           # Expected private data layout
```

## Data

By default, `pps_omicscore.R` expects:

```text
data/clusters.csv
data/id_lists.csv
data/proteins/proteins_visit_0.fst
data/metabolites/nmr_threephases.csv
```

Required columns:

- `clusters.csv`: `f.eid` plus one or more PPS columns.
- `id_lists.csv`: at least `disease_free`, `protein`, and `diet` columns.
- protein file: `f.eid` plus protein feature columns.
- metabolite file: `f.eid`, `visit`, NMR feature columns, and metadata columns.

Do not commit individual-level data, UK Biobank data, model outputs, logs, or
large derived files. The `.gitignore` is configured for that.

## Install

Use R 4.4 or newer. The reproducible environment is recorded in `renv.lock`.
Recommended setup:

```bash
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore()'
```

For a simpler non-locked install, use:

```bash
Rscript requirements.R
```

## Run One PPS Locally

```bash
Rscript pps_omicscore.R \
  --idx 1 \
  --data-dir data \
  --out-dir results/omics_score/pps/disease_free \
  --ncores 4 \
  --seed 1427
```

You can override individual files instead of using `--data-dir`:

```bash
Rscript pps_omicscore.R \
  --idx 1 \
  --prot-file /path/to/proteins_visit_0.fst \
  --metab-file /path/to/nmr_threephases.csv \
  --clusters-file /path/to/clusters.csv \
  --ids-file /path/to/id_lists.csv \
  --out-dir results/omics_score/pps/disease_free
```

## Run With the Shell Launcher

Local run:

```bash
bash run_pps_omicscore.sh 1
```

SGE array run:

```bash
mkdir -p logs
qsub run_pps_omicscore.sh
```

The launcher accepts these environment variables:

- `PPS_DATA_DIR`
- `PPS_PROT_FILE`
- `PPS_METAB_FILE`
- `PPS_CLUSTERS_FILE`
- `PPS_IDS_FILE`
- `PPS_OUT_DIR`
- `PPS_NCORES`
- `PPS_SEED`
- `PPS_CONDA_ENV`

Set `PPS_CONDA_ENV=""` to skip conda activation in the shell launcher.

## Aggregate Outputs

After all PPS jobs finish, aggregate prediction columns and collate summaries:

```bash
Rscript collect.R --pps-root results/omics_score/pps/disease_free
Rscript collate.R --pps-root results/omics_score/pps/disease_free
```

## Outputs

For each PPS, outputs are written under:

```text
results/omics_score/pps/disease_free/<PPS_NAME>/
  ids/
  lasso/
    folds/
    metrics/
    scores/
    weights/
  elastic_net/
    folds/
    metrics/
    scores/
    weights/
  <PPS_NAME>_preproc_params_V3.rds
  performance_summary_V3.csv
```

## Notes Before Publishing

- Confirm whether you want a public license, for example MIT, Apache-2.0, or no
  reuse license.
- Keep real input data outside git.
- Add a short example dataset only if it is synthetic or fully shareable.
