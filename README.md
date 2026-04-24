# Omics LASSO Score Pipeline

This repository contains a general R pipeline for training LASSO-based omics
scores from protein and/or metabolite predictors.

The main script is `lasso_score.R`. It can run three predictor modes:

- `proteins`: protein-only score, named `protscore` by default.
- `metabolites`: metabolite-only score, named `metscore` by default.
- `combined`: proteins + metabolites score, named `omicscore` by default.

The outcome is not hard-coded. You provide any numeric outcome column in an
outcome dataset, and the script creates a reproducible 70/30 train/test split.

The pipeline:

- loads an outcome dataset plus protein and/or metabolite predictor datasets;
- merges samples by an ID column;
- filters predictors and samples by missingness using the training set;
- applies predictor preprocessing:
  - proteins: KNN imputation, inverse-rank normal transform, z-score scaling;
  - metabolites: half-minimum imputation, inverse-rank normal transform,
    z-score scaling;
- scales the outcome within the training set;
- fits a LASSO model with nested cross-validation to summarize lambda values;
- evaluates predictions on the held-out 30% test set;
- saves IDs, preprocessing parameters, feature map, weights, predictions,
  metrics, and fitted model objects.

## Repository Structure

```text
.
  lasso_score.R            # Main generalized score pipeline
  run_lasso_score.sh       # Local/SGE launcher
  requirements.R           # Installs required R packages
  DESCRIPTION              # Project dependency metadata
  data/README.md           # Expected private data layout
```

## Data

Example private data layout:

```text
data/outcome.csv
data/proteins/proteins_visit_0.fst
data/metabolites/nmr_threephases.csv
```

Required columns:

- outcome file: sample ID column plus the numeric outcome column to predict.
- protein file: sample ID column plus numeric protein feature columns.
- metabolite file: sample ID column plus numeric metabolite feature columns.

The default sample ID column is `f.eid`, but you can change it with `--id-col`,
`--protein-id-col`, and `--metabolite-id-col`.

Do not commit individual-level data, restricted biobank data, model outputs,
logs, or large derived files. The `.gitignore` is configured for that.

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

## Run A Score

Protein-only score:

```bash
Rscript lasso_score.R \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --id-col f.eid \
  --mode proteins \
  --protein-file data/proteins/proteins_visit_0.fst \
  --out-dir results
```

Metabolite-only score:

```bash
Rscript lasso_score.R \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode metabolites \
  --metabolite-file data/metabolites/nmr_threephases.csv \
  --metabolite-exclude-cols phase,sample_id,plate_id,plate_position,visit,spectrometer \
  --metabolite-visit-col visit \
  --metabolite-visit-value "Main Phase" \
  --drop-pct-cols \
  --out-dir results
```

Combined protein + metabolite score:

```bash
Rscript lasso_score.R \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode combined \
  --protein-file data/proteins/proteins_visit_0.fst \
  --metabolite-file data/metabolites/nmr_threephases.csv \
  --out-dir results
```

## Run With the Shell Launcher

The launcher simply activates the optional conda environment and forwards all
arguments to `lasso_score.R`.

```bash
bash run_lasso_score.sh \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode combined \
  --protein-file data/proteins/proteins_visit_0.fst \
  --metabolite-file data/metabolites/nmr_threephases.csv
```

SGE run:

```bash
mkdir -p logs
qsub run_lasso_score.sh --outcome-file data/outcome.csv --outcome-col bmi --mode proteins --protein-file data/proteins/proteins_visit_0.fst
```

Set `LASSO_CONDA_ENV=""` to skip conda activation in the shell launcher.

## Useful Options

- `--train-prop`: training proportion. Default is `0.70`.
- `--alpha`: glmnet alpha. Default is `1` for LASSO; use `0.5` for elastic net.
- `--lambda-choice`: `lambda_min`, `lambda_1se`, or `both`. Default is `both`.
- `--feature-missing-threshold`: maximum training-set missingness per feature.
  Default is `0.10`.
- `--sample-missing-threshold`: maximum missingness per sample. Default is
  `0.10`.
- `--outer-folds` and `--inner-folds`: nested cross-validation fold counts.
- `--score-name`: output folder/name. If omitted, defaults to
  `<outcome>_protscore`, `<outcome>_metscore`, or `<outcome>_omicscore`.

## Outputs

For one run, outputs are written under:

```text
results/<score_name>/
  ids/
  folds/
  metrics/
  models/
  scores/
  weights/
  run_config.csv
  <score_name>_feature_map.csv
  <score_name>_preproc_params.rds
  performance_summary.csv
```

Feature names are prefixed internally as `prot__` and `metab__` to avoid name
collisions. The feature map keeps the original names.

## Notes Before Publishing

- Confirm whether you want a public license, for example MIT, Apache-2.0, or no
  reuse license.
- Keep real input data outside git.
- Add a short example dataset only if it is synthetic or fully shareable.
