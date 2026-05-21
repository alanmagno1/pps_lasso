# Omics LASSO Score Pipeline

This repository contains a general R pipeline for training omics-based
prediction scores with penalized linear regression. It can use protein
features, metabolite features, or both to predict any user-specified numeric
outcome.

The main script is `lasso_score.R`.

## What The Pipeline Does

Given an outcome table and one or more omics predictor tables, the pipeline:

- merges samples by a shared ID column;
- creates a reproducible train/test split, 70/30 by default;
- filters predictors and samples by missingness using the training set;
- preprocesses predictors using parameters learned only from training data;
- fits a penalized Gaussian regression model with `glmnet`;
- chooses lambda values from nested cross-validation summaries;
- evaluates predictions in the held-out test set;
- writes IDs, fold assignments, weights, predictions, metrics, and
  preprocessing parameters.

The default model is LASSO (`alpha = 1`). You can also run elastic net by
setting `--alpha 0.5`, or any `glmnet` alpha value between 0 and 1.

## Score Modes

Use `--mode` to choose which predictor block to train on.

| Mode | Required file(s) | Default score suffix | Meaning |
| --- | --- | --- | --- |
| `proteins` | `--protein-file` | `protscore` | Protein-only score |
| `metabolites` | `--metabolite-file` | `metscore` | Metabolite-only score |
| `combined` | `--protein-file` and `--metabolite-file` | `omicscore` | Protein + metabolite score |

If `--score-name` is not supplied, output folders are named from the outcome and
mode, for example `bmi_protscore`, `bmi_metscore`, or `bmi_omicscore`.

## Repository Structure

```text
.
  lasso_score.R              # Main generalized score pipeline
  generate_synthetic_data.R  # Generates synthetic test datasets
  run_lasso_score.sh         # Optional local/SGE shell launcher
  requirements.R             # Installs required R packages
  data/README.md             # Notes about private input data
  logs/.gitkeep              # Empty logs directory placeholder
  results/.gitkeep           # Empty results directory placeholder
```

## Installation

Only tested using R 4.3.

Install required packages with:

```bash
Rscript requirements.R
```

The required R packages are:

- `caret`
- `data.table`
- `doParallel`
- `fst`
- `getopt`
- `glmnet`
- `impute`

## Input Data

The pipeline expects one outcome dataset and at least one predictor dataset.
Files can be `.fst`, `.csv`, `.tsv`, or `.txt`. Files are read with `fst` when
the extension is `.fst`; otherwise they are read with `data.table::fread()`.

Example private layout:

```text
data/
  outcome.csv
  proteins/
    proteins_visit_0.fst
  metabolites/
    nmr_threephases.csv
```

These files are intentionally ignored by git. Do not commit individual-level
data, restricted cohort data, generated model objects, or large result files.

### Outcome File

The outcome file must contain:

- a sample ID column;
- one numeric outcome column to predict.

Minimal shape:

```text
f.eid,bmi
1001,26.4
1002,31.2
1003,24.8
```

Use `--outcome-col` to select the outcome. The outcome must be numeric or
coercible to numeric. Samples with missing outcome values are removed before
modeling.

### Protein File

The protein file must contain:

- a sample ID column;
- numeric protein feature columns.

Minimal shape:

```text
f.eid,protein_a,protein_b,protein_c
1001,0.23,1.42,0.91
1002,0.18,,1.10
1003,0.31,1.20,0.77
```

Protein preprocessing uses KNN imputation, inverse-rank normal transformation,
and z-score scaling.

### Metabolite File

The metabolite file must contain:

- a sample ID column;
- numeric metabolite feature columns;
- optional metadata columns that can be excluded.

Minimal shape:

```text
f.eid,visit,metabolite_a,metabolite_b
1001,Main Phase,14.2,0.81
1002,Main Phase,,0.67
1003,Main Phase,12.9,0.73
```

Metabolite preprocessing uses half-minimum imputation, inverse-rank normal
transformation, and z-score scaling.

Use `--metabolite-exclude-cols` to exclude metadata columns from predictors.
Use `--metabolite-visit-col` and `--metabolite-visit-value` to keep a specific
visit. Use `--drop-pct-cols` to remove metabolite columns ending in `_pct`.

## Basic Usage

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
  --id-col f.eid \
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
  --id-col f.eid \
  --mode combined \
  --protein-file data/proteins/proteins_visit_0.fst \
  --metabolite-file data/metabolites/nmr_threephases.csv \
  --metabolite-exclude-cols phase,sample_id,plate_id,plate_position,visit,spectrometer \
  --metabolite-visit-col visit \
  --metabolite-visit-value "Main Phase" \
  --drop-pct-cols \
  --out-dir results
```

Elastic net example:

```bash
Rscript lasso_score.R \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode combined \
  --protein-file data/proteins/proteins_visit_0.fst \
  --metabolite-file data/metabolites/nmr_threephases.csv \
  --alpha 0.5
```

Custom train/test split:

```bash
Rscript lasso_score.R \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode proteins \
  --protein-file data/proteins/proteins_visit_0.fst \
  --train-ids my_train_ids.csv \
  --test-ids my_test_ids.csv
```

The ID files must be CSVs with a column matching `--id-col` (default `f.eid`).
When `--train-ids` and `--test-ids` are both provided, `--train-prop` is
ignored.

## Shell Launcher

`run_lasso_score.sh` is a thin wrapper around `lasso_score.R`. It optionally
activates a conda environment and forwards all arguments.

```bash
bash run_lasso_score.sh \
  --outcome-file data/outcome.csv \
  --outcome-col bmi \
  --mode combined \
  --protein-file data/proteins/proteins_visit_0.fst \
  --metabolite-file data/metabolites/nmr_threephases.csv
```

By default, the launcher tries to activate a conda environment named `lasso` if
Miniconda is available. To use a different environment:

```bash
LASSO_CONDA_ENV=my_r_env bash run_lasso_score.sh ...
```

To skip conda activation:

```bash
LASSO_CONDA_ENV="" bash run_lasso_score.sh ...
```

For SGE:

```bash
mkdir -p logs
qsub run_lasso_score.sh --outcome-file data/outcome.csv --outcome-col bmi --mode proteins --protein-file data/proteins/proteins_visit_0.fst
```

## Command-Line Arguments

Required arguments:

| Argument | Description |
| --- | --- |
| `--outcome-file`, `-y` | File containing sample IDs and the outcome column |
| `--outcome-col`, `-O` | Numeric outcome column to predict |
| `--mode`, `-M` | One of `proteins`, `metabolites`, or `combined` |

Predictor file arguments:

| Argument | Description |
| --- | --- |
| `--protein-file`, `-p` | Protein predictor file; required for `proteins` and `combined` |
| `--metabolite-file`, `-m` | Metabolite predictor file; required for `metabolites` and `combined` |

Column and filtering arguments:

| Argument | Default | Description |
| --- | --- | --- |
| `--id-col`, `-I` | `f.eid` | Shared sample ID column in the outcome file |
| `--protein-id-col`, `-A` | value of `--id-col` | Protein ID column if named differently |
| `--metabolite-id-col`, `-B` | value of `--id-col` | Metabolite ID column if named differently |
| `--protein-exclude-cols`, `-X` | none | Comma-separated protein columns to exclude |
| `--metabolite-exclude-cols`, `-Y` | none | Comma-separated metabolite columns to exclude |
| `--metabolite-visit-col`, `-v` | none | Optional visit column in the metabolite file |
| `--metabolite-visit-value`, `-V` | none | Visit value to keep |
| `--drop-pct-cols`, `-Z` | off | Drop metabolite columns ending in `_pct` |

Model arguments:

| Argument | Default | Description |
| --- | --- | --- |
| `--train-prop`, `-t` | `0.70` | Proportion of samples assigned to training |
| `--train-ids`, `-J` | none | CSV with pre-defined training sample IDs (must contain `--id-col`) |
| `--test-ids`, `-U` | none | CSV with pre-defined test sample IDs (must contain `--id-col`) |
| `--alpha`, `-a` | `1` | `glmnet` alpha; `1` is LASSO, `0.5` is elastic net |
| `--lambda-choice`, `-L` | `both` | `lambda_min`, `lambda_1se`, or `both` |
| `--feature-missing-threshold`, `-F` | `0.10` | Max allowed missingness per feature in training |
| `--sample-missing-threshold`, `-R` | `0.10` | Max allowed missingness per sample |
| `--outer-folds`, `-k` | `10` | Outer CV folds used to summarize lambda values |
| `--inner-folds`, `-K` | `10` | Inner folds passed to `cv.glmnet()` |
| `--knn-k`, `-q` | `10` | K for protein KNN imputation |

Runtime and output arguments:

| Argument | Default | Description |
| --- | --- | --- |
| `--score-name`, `-N` | `<outcome>_<mode suffix>` | Output score name |
| `--out-dir`, `-o` | `results` | Parent output directory |
| `--ncores`, `-n` | `NSLOTS` or `4` | Parallel workers for `cv.glmnet()` |
| `--seed`, `-S` | `1427` | Random seed for reproducible splitting/CV |
| `--help`, `-h` | off | Print CLI help |

## Preprocessing And Leakage Control

All preprocessing parameters are learned from the training set only. The same
parameters are then applied to the held-out test set.

This matters because the test set should mimic new unseen samples. Estimating
imputation values, rank-normal references, means, or standard deviations from
the full dataset would leak test-set information into training.

Training-set preprocessing:

- predictors with missingness above `--feature-missing-threshold` are removed;
- samples with missingness above `--sample-missing-threshold` are removed;
- proteins are KNN-imputed, inverse-rank normalized, and z-scored;
- metabolites are half-minimum imputed, inverse-rank normalized, and z-scored;
- zero-variance predictors are removed;
- the outcome is z-scored for model fitting.

Test-set preprocessing:

- only retained training-set features are used;
- training-derived imputation, rank-normal, and scaling parameters are applied;
- predictions are transformed back to the original outcome scale.

## Model Fitting

The script fits penalized Gaussian regression models with `glmnet`.

For nested cross-validation:

- outer folds split the training set into train/validation folds;
- each outer fold runs `cv.glmnet()` on the fold-specific training data;
- `lambda.min` and `lambda.1se` are collected across outer folds;
- final lambda values are the median `lambda.min` and median `lambda.1se`;
- final models are fit on the full training set.

`lambda_min` usually gives a less penalized model with more predictors.
`lambda_1se` usually gives a more regularized model with fewer predictors.

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

Output files:

| Path | Description |
| --- | --- |
| `run_config.csv` | Runtime configuration and sample/feature counts |
| `<score_name>_feature_map.csv` | Maps internal feature names to original names and feature type |
| `<score_name>_preproc_params.rds` | Saved preprocessing parameters for reproducibility/scoring |
| `ids/*_train_ids.csv` | Training sample IDs |
| `ids/*_test_ids.csv` | Held-out test sample IDs |
| `ids/*_dropped_*_ids.csv` | Samples dropped by missingness filters, if any |
| `folds/*_fold*_train_ids.csv` | Outer-CV fold training IDs |
| `folds/*_fold*_valid_ids.csv` | Outer-CV fold validation IDs |
| `metrics/*_nested_cv_metrics.csv` | Outer-CV metrics and fold-specific lambda values |
| `models/*_model_lambda_*.rds` | Final fitted `glmnet` model objects |
| `weights/*_weights_lambda_*.csv` | Non-zero model coefficients |
| `scores/*_scores_lambda_*.csv` | Test-set observed values, predictions, and residuals |
| `performance_summary.csv` | Test-set RMSE, MAE, correlation, R2, and number of predictors |

Internal feature names are prefixed as `prot__` and `metab__` to avoid name
collisions between predictor files. Use `<score_name>_feature_map.csv` to map
them back to original feature names.

## Interpreting Results

Start with `performance_summary.csv`.

Important columns:

- `Lambda`: whether results use `lambda_min` or `lambda_1se`;
- `Predictors`: number of non-zero predictors selected by the model;
- `RMSE`: root mean squared error on the original outcome scale;
- `MAE`: mean absolute error on the original outcome scale;
- `Corr`: Pearson correlation between observed and predicted test-set values;
- `R2`: squared Pearson correlation.

Then inspect `weights/` to see which predictors were selected and their
coefficients. Positive weights increase the predicted outcome; negative weights
decrease it, conditional on the preprocessed feature scale.

## Limitations

- The current pipeline supports numeric continuous outcomes only.
- It does not currently fit logistic, Cox, or other non-Gaussian models.
- It does not currently include covariate adjustment.
- It assumes rows are independent samples.

## Troubleshooting

`Outcome file is missing ID column`

Check `--id-col`. If your outcome ID column is named `eid`, use:

```bash
--id-col eid
```

`protein file is missing ID column` or `metabolite file is missing ID column`

If predictor files use different ID names, set:

```bash
--protein-id-col protein_id
--metabolite-id-col metabolite_id
```

`No numeric predictor columns found`

Most likely metadata columns were not excluded, or features were read as
character columns. Check the file format and use `--protein-exclude-cols` or
`--metabolite-exclude-cols`.

`No predictors passed the training-set missingness filter`

Too many values are missing. Increase `--feature-missing-threshold` cautiously,
or inspect missingness in the predictor file.

`Too few samples remain after sample missingness filtering`

Increase `--sample-missing-threshold` cautiously, or check whether the outcome
and predictor files overlap on the ID column.

`Outcome has zero variance in training set`

The selected outcome does not vary after filtering/splitting, or the wrong
column was selected with `--outcome-col`.

## Synthetic Data

To generate synthetic datasets for testing:

```bash
Rscript generate_synthetic_data.R
```

This creates `data/outcome.csv`, `data/proteins/proteins_visit_0.csv`, and
`data/metabolites/nmr_threephases.csv` with 500 samples by default. Options:
`--n-samples`, `--n-proteins`, `--n-metabolites`, `--seed`, `--out-dir`.

## TODOs

- Add covariate adjustment support
- Add logistic/Cox regression modes
