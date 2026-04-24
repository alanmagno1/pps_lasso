#!/usr/bin/env Rscript
# =============================================================================
# General omics score pipeline
#   Modes:
#     proteins     -> protein-only score (protscore)
#     metabolites  -> metabolite-only score (metscore)
#     combined     -> protein + metabolite score (omicscore)
# =============================================================================

suppressPackageStartupMessages({
  library(getopt)
  library(data.table)
  library(fst)
  library(glmnet)
  library(doParallel)
  library(caret)
  library(impute)
})

setDTthreads(1)

# ---------- 0. CLI -----------------------------------------------------------
spec <- matrix(c(
  "outcome-file",              "y", 1, "character",
  "outcome-col",               "O", 1, "character",
  "id-col",                    "I", 1, "character",
  "mode",                      "M", 1, "character",
  "protein-file",              "p", 1, "character",
  "metabolite-file",           "m", 1, "character",
  "protein-id-col",            "A", 1, "character",
  "metabolite-id-col",         "B", 1, "character",
  "protein-exclude-cols",      "X", 1, "character",
  "metabolite-exclude-cols",   "Y", 1, "character",
  "metabolite-visit-col",      "v", 1, "character",
  "metabolite-visit-value",    "V", 1, "character",
  "drop-pct-cols",             "Z", 0, "logical",
  "score-name",                "N", 1, "character",
  "out-dir",                   "o", 1, "character",
  "train-prop",                "t", 1, "double",
  "alpha",                     "a", 1, "double",
  "lambda-choice",             "L", 1, "character",
  "feature-missing-threshold", "F", 1, "double",
  "sample-missing-threshold",  "R", 1, "double",
  "outer-folds",               "k", 1, "integer",
  "inner-folds",               "K", 1, "integer",
  "knn-k",                     "q", 1, "integer",
  "ncores",                    "n", 1, "integer",
  "seed",                      "S", 1, "integer",
  "help",                      "h", 0, "logical"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

usage <- paste(
  "Usage: Rscript lasso_score.R --outcome-file FILE --outcome-col COL --mode MODE [options]",
  "",
  "Required:",
  "  -y, --outcome-file FILE          Dataset with the outcome column",
  "  -O, --outcome-col COL            Numeric outcome to predict",
  "  -M, --mode MODE                  proteins, metabolites, or combined",
  "",
  "Predictor files:",
  "  -p, --protein-file FILE          Protein predictor dataset (.fst, .csv, .tsv, .txt)",
  "  -m, --metabolite-file FILE       Metabolite predictor dataset (.fst, .csv, .tsv, .txt)",
  "",
  "Column options:",
  "  -I, --id-col COL                 Shared sample ID column (default: f.eid)",
  "  -A, --protein-id-col COL         Protein sample ID column if different",
  "  -B, --metabolite-id-col COL      Metabolite sample ID column if different",
  "  -X, --protein-exclude-cols CSV   Protein columns to exclude from predictors",
  "  -Y, --metabolite-exclude-cols CSV Metabolite columns to exclude from predictors",
  "  -v, --metabolite-visit-col COL   Optional metabolite visit column",
  "  -V, --metabolite-visit-value VAL Optional metabolite visit value to keep",
  "  -Z, --drop-pct-cols              Drop metabolite columns ending in '_pct'",
  "",
  "Model options:",
  "  -t, --train-prop NUM             Training split proportion (default: 0.70)",
  "  -a, --alpha NUM                  glmnet alpha; 1=LASSO, 0.5=elastic net (default: 1)",
  "  -L, --lambda-choice VALUE        lambda_min, lambda_1se, or both (default: both)",
  "  -F, --feature-missing-threshold NUM Max train missingness per feature (default: 0.10)",
  "  -R, --sample-missing-threshold NUM  Max missingness per sample (default: 0.10)",
  "  -k, --outer-folds INTEGER        Outer CV folds for lambda summary (default: 10)",
  "  -K, --inner-folds INTEGER        Inner cv.glmnet folds (default: 10)",
  "  -q, --knn-k INTEGER              KNN imputation k for proteins (default: 10)",
  "",
  "Runtime/output:",
  "  -N, --score-name NAME            Output score name",
  "  -o, --out-dir DIR                Output directory (default: results)",
  "  -n, --ncores INTEGER             Parallel workers (default: NSLOTS or 4)",
  "  -S, --seed INTEGER               Random seed (default: 1427)",
  "  -h, --help                       Show this help message",
  sep = "\n"
)

if (!is.null(opt$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

arg_value <- function(name, default = NULL) {
  value <- opt[[name]]
  if (!is.null(value) && !is.na(value) && nzchar(as.character(value))) {
    return(value)
  }
  default
}

require_arg <- function(name, label) {
  value <- arg_value(name)
  if (is.null(value)) {
    stop(sprintf("Missing required argument %s\n\n%s", label, usage), call. = FALSE)
  }
  value
}

parse_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return(character())
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

check_probability <- function(x, name) {
  if (is.na(x) || x <= 0 || x >= 1) {
    stop(sprintf("%s must be > 0 and < 1", name), call. = FALSE)
  }
  x
}

check_threshold <- function(x, name) {
  if (is.na(x) || x < 0 || x > 1) {
    stop(sprintf("%s must be between 0 and 1", name), call. = FALSE)
  }
  x
}

read_table_auto <- function(path) {
  path <- path.expand(path)
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  if (grepl("\\.fst$", path, ignore.case = TRUE)) {
    return(as.data.table(read.fst(path)))
  }
  as.data.table(fread(path))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

# ---------- 1. Runtime config ------------------------------------------------
outcome_file <- path.expand(require_arg("outcome-file", "--outcome-file"))
outcome_col <- require_arg("outcome-col", "--outcome-col")
id_col <- arg_value("id-col", "f.eid")
mode <- tolower(require_arg("mode", "--mode"))

valid_modes <- c("proteins", "metabolites", "combined")
if (!mode %in% valid_modes) {
  stop("--mode must be one of: proteins, metabolites, combined", call. = FALSE)
}

protein_file <- arg_value("protein-file")
metabolite_file <- arg_value("metabolite-file")
if (mode %in% c("proteins", "combined") && is.null(protein_file)) {
  stop("--protein-file is required for proteins and combined modes", call. = FALSE)
}
if (mode %in% c("metabolites", "combined") && is.null(metabolite_file)) {
  stop("--metabolite-file is required for metabolites and combined modes", call. = FALSE)
}

protein_id_col <- arg_value("protein-id-col", id_col)
metabolite_id_col <- arg_value("metabolite-id-col", id_col)
protein_exclude_cols <- parse_csv(arg_value("protein-exclude-cols"))
metabolite_exclude_cols <- parse_csv(arg_value("metabolite-exclude-cols"))
metabolite_visit_col <- arg_value("metabolite-visit-col")
metabolite_visit_value <- arg_value("metabolite-visit-value")
drop_pct_cols <- isTRUE(arg_value("drop-pct-cols", FALSE))

train_prop <- check_probability(as.numeric(arg_value("train-prop", 0.70)), "--train-prop")
alpha <- as.numeric(arg_value("alpha", 1))
if (is.na(alpha) || alpha < 0 || alpha > 1) {
  stop("--alpha must be between 0 and 1", call. = FALSE)
}

lambda_choice <- tolower(arg_value("lambda-choice", "both"))
valid_lambda_choices <- c("lambda_min", "lambda_1se", "both")
if (!lambda_choice %in% valid_lambda_choices) {
  stop("--lambda-choice must be lambda_min, lambda_1se, or both", call. = FALSE)
}

feature_missing_threshold <- check_threshold(
  as.numeric(arg_value("feature-missing-threshold", 0.10)),
  "--feature-missing-threshold"
)
sample_missing_threshold <- check_threshold(
  as.numeric(arg_value("sample-missing-threshold", 0.10)),
  "--sample-missing-threshold"
)

outer_folds <- as.integer(arg_value("outer-folds", 10))
inner_folds <- as.integer(arg_value("inner-folds", 10))
knn_k <- as.integer(arg_value("knn-k", 10))
ncores <- as.integer(arg_value("ncores", Sys.getenv("NSLOTS", "4")))
seed <- as.integer(arg_value("seed", 1427))

if (any(is.na(c(outer_folds, inner_folds, knn_k, ncores, seed)))) {
  stop("--outer-folds, --inner-folds, --knn-k, --ncores, and --seed must be integers", call. = FALSE)
}
if (outer_folds < 2 || inner_folds < 2 || knn_k < 1 || ncores < 1) {
  stop("Fold counts, --knn-k, and --ncores must be positive; folds must be >= 2", call. = FALSE)
}

default_suffix <- switch(
  mode,
  proteins = "protscore",
  metabolites = "metscore",
  combined = "omicscore"
)
score_name <- arg_value("score-name", paste(safe_name(outcome_col), default_suffix, sep = "_"))
out_dir <- path.expand(arg_value("out-dir", "results"))
run_dir <- file.path(out_dir, score_name)

registerDoParallel(ncores)
set.seed(seed)

# ---------- 2. Data loading --------------------------------------------------
load_outcome <- function(path, id_col, outcome_col) {
  dt <- read_table_auto(path)
  if (!id_col %in% names(dt)) {
    stop("Outcome file is missing ID column: ", id_col, call. = FALSE)
  }
  if (!outcome_col %in% names(dt)) {
    stop("Outcome file is missing outcome column: ", outcome_col, call. = FALSE)
  }
  dt <- dt[, c(id_col, outcome_col), with = FALSE]
  dt <- dt[!duplicated(get(id_col))]
  if (!is.numeric(dt[[outcome_col]])) {
    converted <- suppressWarnings(as.numeric(dt[[outcome_col]]))
    if (sum(!is.na(converted)) == sum(!is.na(dt[[outcome_col]]))) {
      dt[, (outcome_col) := converted]
    } else {
      stop("--outcome-col must be numeric or coercible to numeric", call. = FALSE)
    }
  }
  dt <- dt[!is.na(get(outcome_col))]
  if (nrow(dt) < 10) {
    stop("Need at least 10 samples with non-missing outcome", call. = FALSE)
  }
  dt
}

load_predictor_block <- function(path, source, id_col, source_id_col, exclude_cols,
                                 visit_col = NULL, visit_value = NULL,
                                 drop_pct_cols = FALSE) {
  dt <- read_table_auto(path)
  if (!source_id_col %in% names(dt)) {
    stop(sprintf("%s file is missing ID column: %s", source, source_id_col), call. = FALSE)
  }
  if (source_id_col != id_col) {
    setnames(dt, source_id_col, id_col)
  }

  if (!is.null(visit_col)) {
    if (!visit_col %in% names(dt)) {
      stop(sprintf("%s file is missing visit column: %s", source, visit_col), call. = FALSE)
    }
    if (!is.null(visit_value)) {
      dt <- dt[get(visit_col) == visit_value]
    }
    exclude_cols <- unique(c(exclude_cols, visit_col))
  }

  exclude_cols <- unique(c(id_col, outcome_col, exclude_cols))
  candidate_cols <- setdiff(names(dt), exclude_cols)
  if (drop_pct_cols) {
    candidate_cols <- candidate_cols[!grepl("_pct$", candidate_cols)]
  }

  numeric_cols <- candidate_cols[vapply(dt[, ..candidate_cols], is.numeric, logical(1))]
  dropped_non_numeric <- setdiff(candidate_cols, numeric_cols)
  if (!length(numeric_cols)) {
    stop(sprintf("No numeric predictor columns found in %s file", source), call. = FALSE)
  }

  dt <- dt[, c(id_col, numeric_cols), with = FALSE]
  dt <- dt[!duplicated(get(id_col))]

  prefix <- if (source == "protein") "prot__" else "metab__"
  feature_cols <- paste0(prefix, make.unique(numeric_cols))
  setnames(dt, numeric_cols, feature_cols)

  feature_map <- data.table(
    Feature = feature_cols,
    OriginalFeature = numeric_cols,
    Type = source
  )

  list(
    data = dt,
    features = feature_cols,
    feature_map = feature_map,
    dropped_non_numeric = dropped_non_numeric
  )
}

message("Loading outcome and predictor data")
outcome_dt <- load_outcome(outcome_file, id_col, outcome_col)

blocks <- list()
if (mode %in% c("proteins", "combined")) {
  blocks$protein <- load_predictor_block(
    protein_file, "protein", id_col, protein_id_col, protein_exclude_cols
  )
}
if (mode %in% c("metabolites", "combined")) {
  blocks$metabolite <- load_predictor_block(
    metabolite_file, "metabolite", id_col, metabolite_id_col,
    metabolite_exclude_cols, metabolite_visit_col, metabolite_visit_value,
    drop_pct_cols
  )
}

merged <- Reduce(
  function(x, y) merge(x, y, by = id_col, all = FALSE, sort = FALSE),
  c(list(outcome_dt), lapply(blocks, `[[`, "data"))
)
merged <- merged[!duplicated(get(id_col))]

protein_features <- if (!is.null(blocks$protein)) blocks$protein$features else character()
metabolite_features <- if (!is.null(blocks$metabolite)) blocks$metabolite$features else character()
feature_map <- rbindlist(lapply(blocks, `[[`, "feature_map"), use.names = TRUE)

if (nrow(merged) < 10) {
  stop("Need at least 10 samples after merging outcome and predictor files", call. = FALSE)
}

message(sprintf("Merged %d samples and %d predictor columns", nrow(merged), nrow(feature_map)))

# ---------- 3. Train/test split and feature filters --------------------------
train_idx <- createDataPartition(merged[[outcome_col]], p = train_prop, list = FALSE)
train <- copy(merged[train_idx])
test <- copy(merged[-train_idx])

if (nrow(test) < 2) {
  stop("Test set has fewer than 2 samples; reduce --train-prop or add samples", call. = FALSE)
}

filter_features_by_train_missingness <- function(train_dt, features, threshold) {
  if (!length(features)) {
    return(character())
  }
  miss <- colMeans(is.na(train_dt[, ..features]))
  names(miss)[miss <= threshold]
}

protein_features <- filter_features_by_train_missingness(
  train, protein_features, feature_missing_threshold
)
metabolite_features <- filter_features_by_train_missingness(
  train, metabolite_features, feature_missing_threshold
)
all_features <- c(protein_features, metabolite_features)
if (!length(all_features)) {
  stop("No predictors passed the training-set missingness filter", call. = FALSE)
}
feature_map <- feature_map[Feature %in% all_features]

drop_sparse_samples <- function(dt, features, threshold) {
  prop_missing <- rowMeans(is.na(dt[, ..features]))
  list(
    kept = dt[prop_missing <= threshold],
    dropped = dt[prop_missing > threshold, .(sample_id = get(id_col), missing_fraction = prop_missing[prop_missing > threshold])]
  )
}

train_filter <- drop_sparse_samples(train, all_features, sample_missing_threshold)
test_filter <- drop_sparse_samples(test, all_features, sample_missing_threshold)
train <- train_filter$kept
test <- test_filter$kept

if (nrow(train) < 10 || nrow(test) < 2) {
  stop("Too few samples remain after sample missingness filtering", call. = FALSE)
}

message(sprintf("Training samples: %d; test samples: %d", nrow(train), nrow(test)))
message(sprintf("Protein predictors: %d; metabolite predictors: %d", length(protein_features), length(metabolite_features)))

# ---------- 4. Preprocessing -------------------------------------------------
build_rin_ref <- function(v) {
  x_sorted <- sort(v)
  p <- (seq_along(x_sorted) - 0.5) / length(x_sorted)
  up <- tapply(p, x_sorted, mean)
  list(vals = as.numeric(names(up)), p = as.numeric(up))
}

apply_rin_from_ref <- function(ref, x) {
  if (length(ref$vals) == 1L) {
    return(rep(0, length(x)))
  }
  p <- approx(x = ref$vals, y = ref$p, xout = x, method = "linear", rule = 2)$y
  qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

median_values <- function(dt, features) {
  vapply(
    features,
    function(feature) median(dt[[feature]], na.rm = TRUE),
    numeric(1)
  )
}

median_impute <- function(dt, features, medians) {
  if (!length(features)) {
    return(dt)
  }
  for (feature in features) {
    replacement <- medians[[feature]]
    if (is.na(replacement) || is.infinite(replacement)) {
      replacement <- 0
    }
    values <- dt[[feature]]
    values[is.na(values)] <- replacement
    set(dt, j = feature, value = values)
  }
  dt
}

protein_knn_impute <- function(dt, features, k, medians = NULL) {
  if (!length(features)) {
    return(dt)
  }
  if (is.null(medians)) {
    medians <- median_values(dt, features)
  }
  if (!anyNA(dt[, ..features])) {
    return(dt)
  }
  if (length(features) < 2 || nrow(dt) < 2) {
    return(median_impute(dt, features, medians))
  }

  mat <- as.matrix(dt[, ..features])
  k_use <- min(k, max(1L, length(features) - 1L))
  imputed <- tryCatch(
    impute.knn(t(mat), k = k_use)$data,
    error = function(e) NULL
  )
  if (is.null(imputed)) {
    return(median_impute(dt, features, medians))
  }
  dt[, (features) := as.data.table(t(imputed))]
  median_impute(dt, features, medians)
}

preprocess_proteins_train <- function(train_dt, features, k) {
  if (!length(features)) {
    return(list(train = train_dt, features = character(), params = list()))
  }
  medians <- median_values(train_dt, features)
  train_dt <- protein_knn_impute(train_dt, features, k, medians)

  rin_ref <- lapply(train_dt[, ..features], build_rin_ref)
  train_dt[, (features) := Map(apply_rin_from_ref, rin_ref, .SD), .SDcols = features]

  mu <- sapply(train_dt[, ..features], mean)
  sd <- sapply(train_dt[, ..features], sd)
  zero_var <- names(sd)[is.na(sd) | sd == 0]
  if (length(zero_var)) {
    features <- setdiff(features, zero_var)
    train_dt[, (zero_var) := NULL]
    mu <- mu[features]
    sd <- sd[features]
    rin_ref <- rin_ref[features]
    medians <- medians[features]
  }
  if (length(features)) {
    train_dt[, (features) := Map(function(x, m, s) (x - m) / s, .SD, mu, sd), .SDcols = features]
  }

  list(
    train = train_dt,
    features = features,
    params = list(mu = mu, sd = sd, rin_ref = rin_ref, medians = medians, k = k)
  )
}

preprocess_metabolites_train <- function(train_dt, features) {
  if (!length(features)) {
    return(list(train = train_dt, features = character(), params = list()))
  }
  mins <- sapply(train_dt[, ..features], min, na.rm = TRUE) / 2
  mins[is.infinite(mins) | is.na(mins)] <- 0
  train_dt[, (features) := Map(function(x, m) {
    x[is.na(x)] <- m
    x
  }, .SD, mins), .SDcols = features]

  rin_ref <- lapply(train_dt[, ..features], build_rin_ref)
  train_dt[, (features) := Map(apply_rin_from_ref, rin_ref, .SD), .SDcols = features]

  mu <- sapply(train_dt[, ..features], mean)
  sd <- sapply(train_dt[, ..features], sd)
  zero_var <- names(sd)[is.na(sd) | sd == 0]
  if (length(zero_var)) {
    features <- setdiff(features, zero_var)
    train_dt[, (zero_var) := NULL]
    mu <- mu[features]
    sd <- sd[features]
    mins <- mins[features]
    rin_ref <- rin_ref[features]
  }
  if (length(features)) {
    train_dt[, (features) := Map(function(x, m, s) (x - m) / s, .SD, mu, sd), .SDcols = features]
  }

  list(
    train = train_dt,
    features = features,
    params = list(mu = mu, sd = sd, mins = mins, rin_ref = rin_ref)
  )
}

protein_prep <- preprocess_proteins_train(train, protein_features, knn_k)
train <- protein_prep$train
protein_features <- protein_prep$features

metabolite_prep <- preprocess_metabolites_train(train, metabolite_features)
train <- metabolite_prep$train
metabolite_features <- metabolite_prep$features

all_features <- c(protein_features, metabolite_features)
feature_map <- feature_map[Feature %in% all_features]
if (!length(all_features)) {
  stop("No predictors remain after preprocessing zero-variance filters", call. = FALSE)
}

outcome_mu <- mean(train[[outcome_col]])
outcome_sd <- sd(train[[outcome_col]])
if (is.na(outcome_sd) || outcome_sd == 0) {
  stop("Outcome has zero variance in training set", call. = FALSE)
}
train[, .outcome_scaled := (get(outcome_col) - outcome_mu) / outcome_sd]

apply_preprocessing_to_test <- function(test_dt, protein_features, metabolite_features,
                                        protein_params, metabolite_params) {
  if (length(protein_features)) {
    test_dt <- protein_knn_impute(
      test_dt, protein_features, protein_params$k, protein_params$medians
    )
    test_dt[, (protein_features) := Map(
      apply_rin_from_ref, protein_params$rin_ref, .SD
    ), .SDcols = protein_features]
    test_dt[, (protein_features) := Map(
      function(x, m, s) (x - m) / s, .SD, protein_params$mu, protein_params$sd
    ), .SDcols = protein_features]
  }
  if (length(metabolite_features)) {
    test_dt[, (metabolite_features) := Map(function(x, m) {
      x[is.na(x)] <- m
      x
    }, .SD, metabolite_params$mins), .SDcols = metabolite_features]
    test_dt[, (metabolite_features) := Map(
      apply_rin_from_ref, metabolite_params$rin_ref, .SD
    ), .SDcols = metabolite_features]
    test_dt[, (metabolite_features) := Map(
      function(x, m, s) (x - m) / s, .SD, metabolite_params$mu, metabolite_params$sd
    ), .SDcols = metabolite_features]
  }
  test_dt
}

test <- apply_preprocessing_to_test(
  test, protein_features, metabolite_features,
  protein_prep$params, metabolite_prep$params
)
test[, .outcome_scaled := (get(outcome_col) - outcome_mu) / outcome_sd]

# ---------- 5. Output setup --------------------------------------------------
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "ids"), showWarnings = FALSE)
dir.create(file.path(run_dir, "folds"), showWarnings = FALSE)
dir.create(file.path(run_dir, "metrics"), showWarnings = FALSE)
dir.create(file.path(run_dir, "models"), showWarnings = FALSE)
dir.create(file.path(run_dir, "scores"), showWarnings = FALSE)
dir.create(file.path(run_dir, "weights"), showWarnings = FALSE)

fwrite(
  data.table(sample_id = train[[id_col]]),
  file.path(run_dir, "ids", sprintf("%s_train_ids.csv", score_name))
)
fwrite(
  data.table(sample_id = test[[id_col]]),
  file.path(run_dir, "ids", sprintf("%s_test_ids.csv", score_name))
)
if (nrow(train_filter$dropped)) {
  fwrite(train_filter$dropped, file.path(run_dir, "ids", sprintf("%s_dropped_train_ids.csv", score_name)))
}
if (nrow(test_filter$dropped)) {
  fwrite(test_filter$dropped, file.path(run_dir, "ids", sprintf("%s_dropped_test_ids.csv", score_name)))
}

run_config <- data.table(
  Parameter = c(
    "score_name", "mode", "outcome_file", "outcome_col", "id_col",
    "protein_file", "metabolite_file", "train_prop", "alpha",
    "lambda_choice", "feature_missing_threshold", "sample_missing_threshold",
    "outer_folds", "inner_folds", "knn_k", "ncores", "seed",
    "n_merged", "n_train", "n_test", "n_protein_features", "n_metabolite_features"
  ),
  Value = as.character(c(
    score_name, mode, outcome_file, outcome_col, id_col,
    ifelse(is.null(protein_file), "", protein_file),
    ifelse(is.null(metabolite_file), "", metabolite_file),
    train_prop, alpha, lambda_choice, feature_missing_threshold,
    sample_missing_threshold, outer_folds, inner_folds, knn_k, ncores, seed,
    nrow(merged), nrow(train), nrow(test), length(protein_features), length(metabolite_features)
  ))
)
fwrite(run_config, file.path(run_dir, "run_config.csv"))
fwrite(feature_map, file.path(run_dir, sprintf("%s_feature_map.csv", score_name)))

preproc_params <- list(
  score_name = score_name,
  mode = mode,
  outcome = list(column = outcome_col, mu = outcome_mu, sd = outcome_sd),
  id_col = id_col,
  protein = c(protein_prep$params, list(features = protein_features)),
  metabolite = c(metabolite_prep$params, list(features = metabolite_features)),
  feature_map = feature_map
)
saveRDS(preproc_params, file.path(run_dir, sprintf("%s_preproc_params.rds", score_name)))

# ---------- 6. Model fitting -------------------------------------------------
metrics_fun <- function(pred, obs) {
  rmse <- sqrt(mean((pred - obs)^2))
  mae <- mean(abs(pred - obs))
  pred_sd <- sd(pred)
  obs_sd <- sd(obs)
  r <- if (length(pred) < 2 || is.na(pred_sd) || is.na(obs_sd) || pred_sd == 0 || obs_sd == 0) {
    NA_real_
  } else {
    cor(pred, obs)
  }
  data.table(RMSE = rmse, MAE = mae, Corr = r, R2 = r^2)
}

fit_cv <- function(x, y, alpha, nfolds, parallel) {
  cv.glmnet(
    x, y,
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    parallel = parallel,
    standardize = FALSE
  )
}

outer_folds <- min(outer_folds, nrow(train))
inner_folds <- min(inner_folds, max(2L, nrow(train) - 1L))
folds_out <- createFolds(train$.outcome_scaled, k = outer_folds, list = TRUE, returnTrain = TRUE)
lam_min <- lam_1se <- numeric(length(folds_out))
fold_stats <- data.table()

for (i in seq_along(folds_out)) {
  tr <- folds_out[[i]]
  va <- setdiff(seq_len(nrow(train)), tr)

  xtr <- as.matrix(train[tr, ..all_features])
  ytr <- train$.outcome_scaled[tr]
  xva <- as.matrix(train[va, ..all_features])
  yva <- train$.outcome_scaled[va]

  cvfit <- fit_cv(xtr, ytr, alpha, min(inner_folds, nrow(xtr)), TRUE)
  lam_min[i] <- cvfit$lambda.min
  lam_1se[i] <- cvfit$lambda.1se

  pred_scaled <- predict(cvfit, xva, s = cvfit$lambda.1se)[, 1]
  pred_orig <- pred_scaled * outcome_sd + outcome_mu
  obs_orig <- yva * outcome_sd + outcome_mu
  mets <- metrics_fun(pred_orig, obs_orig)

  fold_stats <- rbind(
    fold_stats,
    data.table(
      fold = i,
      LambdaMin = cvfit$lambda.min,
      Lambda1se = cvfit$lambda.1se,
      RMSE = mets$RMSE,
      MAE = mets$MAE,
      Corr = mets$Corr,
      R2 = mets$R2
    )
  )

  fwrite(
    data.table(sample_id = train[[id_col]][tr]),
    file.path(run_dir, "folds", sprintf("%s_fold%02d_train_ids.csv", score_name, i))
  )
  fwrite(
    data.table(sample_id = train[[id_col]][va]),
    file.path(run_dir, "folds", sprintf("%s_fold%02d_valid_ids.csv", score_name, i))
  )
}

fwrite(fold_stats, file.path(run_dir, "metrics", sprintf("%s_nested_cv_metrics.csv", score_name)))

lambda_values <- list(
  lambda_min = median(lam_min),
  lambda_1se = median(lam_1se)
)
if (lambda_choice != "both") {
  lambda_values <- lambda_values[lambda_choice]
}

x_train <- as.matrix(train[, ..all_features])
y_train <- train$.outcome_scaled
x_test <- as.matrix(test[, ..all_features])
y_test_scaled <- test$.outcome_scaled

results <- data.table()

save_coefs <- function(fit, tag, lambda_value) {
  raw_coefs <- coef(fit)
  weights <- as.numeric(raw_coefs)
  features <- rownames(raw_coefs)
  nz <- which(weights != 0)
  coef_dt <- data.table(
    Score = score_name,
    Outcome = outcome_col,
    Mode = mode,
    Feature = features[nz],
    Weight = weights[nz],
    Lambda = lambda_value,
    Alpha = alpha
  )
  coef_dt <- merge(
    coef_dt,
    rbind(
      data.table(Feature = "(Intercept)", OriginalFeature = "(Intercept)", Type = "intercept"),
      feature_map
    ),
    by = "Feature",
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(coef_dt, c(
    "Score", "Outcome", "Mode", "Feature", "OriginalFeature", "Type",
    "Weight", "Lambda", "Alpha"
  ))
  fwrite(coef_dt, file.path(run_dir, "weights", sprintf("%s_weights_%s.csv", score_name, tag)))
  coef_dt
}

evaluate_model <- function(coefs, tag) {
  intercept <- if ("(Intercept)" %in% coefs$Feature) {
    coefs[Feature == "(Intercept)", Weight][1]
  } else {
    0
  }
  w <- coefs[Feature != "(Intercept)"]
  if (nrow(w)) {
    pred_scaled <- as.vector(x_test[, w$Feature, drop = FALSE] %*% w$Weight + intercept)
  } else {
    pred_scaled <- rep(intercept, nrow(x_test))
  }
  pred_orig <- pred_scaled * outcome_sd + outcome_mu
  obs_orig <- y_test_scaled * outcome_sd + outcome_mu
  mets <- metrics_fun(pred_orig, obs_orig)

  scores <- data.table(
    sample_id = test[[id_col]],
    Observed = obs_orig,
    Predicted = pred_orig,
    Residual = obs_orig - pred_orig
  )
  fwrite(scores, file.path(run_dir, "scores", sprintf("%s_scores_%s.csv", score_name, tag)))

  data.table(
    Score = score_name,
    Outcome = outcome_col,
    Mode = mode,
    Lambda = tag,
    Alpha = alpha,
    Predictors = nrow(w),
    NTrain = nrow(train),
    NTest = nrow(test),
    RMSE = mets$RMSE,
    MAE = mets$MAE,
    Corr = mets$Corr,
    R2 = mets$R2
  )
}

for (tag in names(lambda_values)) {
  lambda_value <- lambda_values[[tag]]
  fit <- glmnet(
    x_train, y_train,
    family = "gaussian",
    alpha = alpha,
    lambda = lambda_value,
    standardize = FALSE
  )
  saveRDS(fit, file.path(run_dir, "models", sprintf("%s_model_%s.rds", score_name, tag)))
  coefs <- save_coefs(fit, tag, lambda_value)
  results <- rbind(results, evaluate_model(coefs, tag))
}

fwrite(results, file.path(run_dir, "performance_summary.csv"))

message(sprintf(
  "Finished %s: %d train, %d test, %d predictors",
  score_name, nrow(train), nrow(test), length(all_features)
))
message("Outputs written to: ", run_dir)
