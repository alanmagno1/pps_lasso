#!/usr/bin/env Rscript
# =============================================================================
# PPS Metabolite Models (V3)
# =============================================================================

suppressPackageStartupMessages({
  library(getopt)
  library(data.table); setDTthreads(1)
  library(readr)
  library(glmnet)
  library(doParallel)
  library(caret)
})

# ---------- 0 · Command-line interface --------------------------------------
spec <- matrix(c(
  "idx",           "i", 1, "integer",
  "data-dir",      "d", 1, "character",
  "metab-file",    "m", 1, "character",
  "clusters-file", "c", 1, "character",
  "ids-file",      "s", 1, "character",
  "out-dir",       "o", 1, "character",
  "ncores",        "n", 1, "integer",
  "seed",          "S", 1, "integer",
  "help",          "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

usage <- paste(
  "Usage: Rscript metscore.R --idx <1-n_PPS> [options]",
  "",
  "Required:",
  "  -i, --idx INTEGER             1-based PPS column index from clusters.csv",
  "",
  "Input options:",
  "  -d, --data-dir PATH           Base input directory (default: data or PPS_DATA_DIR)",
  "  -m, --metab-file PATH         NMR/metabolite CSV file",
  "  -c, --clusters-file PATH      PPS clusters CSV file",
  "  -s, --ids-file PATH           ID lists CSV file",
  "",
  "Output/runtime options:",
  "  -o, --out-dir PATH            Output directory (default: results/metabolite_score/pps)",
  "  -n, --ncores INTEGER          Parallel workers (default: NSLOTS, PPS_NCORES, or 4)",
  "  -S, --seed INTEGER            Random seed (default: PPS_SEED or 1427)",
  "  -h, --help                    Show this help message",
  sep = "\n"
)

if (!is.null(opt$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}
if (is.null(opt$idx)) {
  stop(usage, call. = FALSE)
}
pps_idx <- as.integer(opt$idx)

arg_or_env <- function(opt_name, env_name, default) {
  value <- opt[[opt_name]]
  if (!is.null(value) && !is.na(value) && nzchar(as.character(value))) {
    return(value)
  }
  env_value <- Sys.getenv(env_name, unset = NA_character_)
  if (!is.na(env_value) && nzchar(env_value)) {
    return(env_value)
  }
  default
}

# ---------- 1 · Parallel resources ------------------------------------------
ncores <- as.integer(arg_or_env("ncores", "PPS_NCORES", Sys.getenv("NSLOTS", "4")))
if (is.na(ncores) || ncores < 1) {
  stop("--ncores must be a positive integer", call. = FALSE)
}
registerDoParallel(ncores)

seed <- as.integer(arg_or_env("seed", "PPS_SEED", "1427"))
if (is.na(seed)) {
  stop("--seed must be an integer", call. = FALSE)
}
set.seed(seed)

# ---------- 2 · Paths and constants -----------------------------------------
DATA_DIR <- path.expand(arg_or_env("data-dir", "PPS_DATA_DIR", "data"))
METAB_FILE <- path.expand(arg_or_env(
  "metab-file", "PPS_METAB_FILE",
  file.path(DATA_DIR, "metabolites", "nmr_threephases.csv")
))
CLUSTERS_FILE <- path.expand(arg_or_env(
  "clusters-file", "PPS_CLUSTERS_FILE",
  file.path(DATA_DIR, "clusters.csv")
))
IDS_FILE <- path.expand(arg_or_env(
  "ids-file", "PPS_IDS_FILE",
  file.path(DATA_DIR, "id_lists.csv")
))
BASE_DIR <- path.expand(arg_or_env(
  "out-dir", "PPS_OUT_DIR",
  file.path("results", "metabolite_score", "pps")
))

input_files <- c(
  METAB_FILE = METAB_FILE,
  CLUSTERS_FILE = CLUSTERS_FILE,
  IDS_FILE = IDS_FILE
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop(
    "Missing input file(s):\n",
    paste(sprintf("  %s: %s", names(missing_files), missing_files), collapse = "\n"),
    call. = FALSE
  )
}

pps_all  <- fread(CLUSTERS_FILE)
pps_cols <- setdiff(names(pps_all), "f.eid")
if (is.na(pps_idx) || pps_idx < 1 || pps_idx > length(pps_cols)) {
  stop(
    sprintf("--idx must be between 1 and %d for %s", length(pps_cols), CLUSTERS_FILE),
    call. = FALSE
  )
}
pps_name <- pps_cols[pps_idx]
message(sprintf("Running PPS %d/%d: %s", pps_idx, length(pps_cols), pps_name))

PPS_DIR <- file.path(BASE_DIR, pps_name)
IDS_DIR <- file.path(PPS_DIR, "ids")
dir.create(PPS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(IDS_DIR,            showWarnings = FALSE)

# ---------- 3 · Load PPS and metabolite data --------------------------------
nmr      <- as.data.table(read_csv(METAB_FILE))
nmr_main <- nmr[visit == "Main Phase"]

exclude_meta <- c(
  "phase","eid_30418","V2","V3","V4","V5","V6",
  "f.eid","V8","V9","V10","V11","V12",
  "sample_id","plate_id","plate_position","visit",
  "spectrometer","Glucose"
)
pct_cols   <- grep("_pct$", names(nmr_main), value = TRUE)
features   <- setdiff(names(nmr_main), unique(c(exclude_meta, pct_cols)))

na_col     <- colMeans(is.na(nmr_main[, ..features]))
keep_metab <- names(na_col)[na_col <= 0.10]
na_row     <- rowMeans(is.na(nmr_main[, ..keep_metab]))
metab      <- nmr_main[na_row <= 0.10, c("f.eid", keep_metab), with = FALSE]

merged <- merge(pps_all, metab, by = "f.eid")[!duplicated(f.eid)]
merged <- merged[!is.na( merged[[pps_name]] )]

ids <- fread(IDS_FILE)
merged <- merged[f.eid %in% ids$disease_free]

# ---------- 4 · 70/30 split -------------------------------------------------
train_idx <- createDataPartition(merged[[pps_name]], p = 0.7, list = FALSE)
train70   <- merged[train_idx]
test30    <- merged[-train_idx]

fwrite(data.table(f.eid = test30$f.eid),
       file = file.path(IDS_DIR, sprintf("%s_test30_ids_V3.csv", pps_name)))

# ---------- 5 · Pre‑processing (training only) ------------------------------
# 5.1 Half‑minimum imputation
mins_train <- sapply(train70[, ..keep_metab], min, na.rm = TRUE) / 2
train70[, (keep_metab) := Map(\(x,m){x[is.na(x)]<-m;x}, .SD, mins_train),
        .SDcols = keep_metab]

# 5.2 Inverse‑rank normal transformation
build_rin_ref <- function(v) {
  ord <- order(v)
  x_sorted <- v[ord]
  p <- (seq_along(x_sorted) - 0.5) / length(x_sorted)
  # Colapsar duplicados para evitar problemas en approx()
  ux <- unique(x_sorted)
  up <- tapply(p, x_sorted, mean) # Promedio de p para valores idénticos
  list(vals = as.numeric(names(up)), p = as.numeric(up))
}
rin_ref_list <- lapply(train70[, ..keep_metab], build_rin_ref)

apply_rin_from_ref <- function(ref, x) {
  p <- approx(x = ref$vals, y = ref$p, xout = x, method = "linear", rule = 2)$y
  qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
}
train70[, (keep_metab) := Map(apply_rin_from_ref, rin_ref_list, .SD), .SDcols = keep_metab]

# 5.3 Z‑score scaling of metabolites
mu_train <- sapply(train70[, ..keep_metab], mean)
sd_train <- sapply(train70[, ..keep_metab], sd)
zero_var <- names(sd_train)[sd_train == 0]
if (length(zero_var)) {
  keep_metab <- setdiff(keep_metab, zero_var)
  train70[, (zero_var) := NULL]
  mu_train <- mu_train[keep_metab]
  sd_train <- sd_train[keep_metab]
  ins_train <- mins_train[keep_metab]
}
train70[, (keep_metab) := Map(\(x,m,s) (x - m) / s, .SD, mu_train, sd_train),
        .SDcols = keep_metab]

# 5.4 Z‑score scaling of the PPS (response)
pps_mu <- mean(train70[[pps_name]])
pps_sd <- sd(train70[[pps_name]])
train70[, (pps_name) := (get(pps_name) - pps_mu) / pps_sd]

# 5.5 Save preprocessing parameters
saveRDS(list(mu = mu_train,
             sd = sd_train,
             mins = mins_train,
             pps_mu = pps_mu,
             pps_sd = pps_sd,
             rin_ref = rin_ref_list,
             feats  = keep_metab),
        file = file.path(PPS_DIR, sprintf("%s_preproc_params_V3.rds", pps_name)))

# ---------- 6 · Utility: performance metrics -------------------------------
metrics_fun <- function(pred, obs) {
  rmse <- sqrt(mean((pred - obs)^2))
  mae  <- mean(abs(pred - obs))
  r    <- ifelse(sd(pred) == 0 || sd(obs) == 0, NA, cor(pred, obs))
  data.table(RMSE = rmse, MAE = mae, Corr = r, R2 = r^2)
}

# ---------- 7 · Model variants: LASSO and Elastic Net ----------------------
folds_out_global <- createFolds(train70[[pps_name]], k = 10,
                                list = TRUE, returnTrain = TRUE)

variants <- list(lasso = 1, elastic_net = 0.5)   # name -> alpha
results  <- data.table()

for (var in names(variants)) {
  alpha_val <- variants[[var]]
  folds_out <- folds_out_global

  # --- Create variant‑specific directories
  VAR_DIR   <- file.path(PPS_DIR, var); dir.create(VAR_DIR, showWarnings = FALSE)
  FOLDS_DIR <- file.path(VAR_DIR, "folds");   dir.create(FOLDS_DIR,   showWarnings = FALSE)
  WEIG_DIR  <- file.path(VAR_DIR, "weights"); dir.create(WEIG_DIR,    showWarnings = FALSE)
  METR_DIR  <- file.path(VAR_DIR, "metrics"); dir.create(METR_DIR,    showWarnings = FALSE)
  SCOR_DIR  <- file.path(VAR_DIR, "scores");  dir.create(SCOR_DIR,    showWarnings = FALSE)

  # --- Nested CV (outer folds = 10)
  lam_min <- lam_1se <- numeric(length(folds_out))
  fold_stats <- data.table(fold = integer(), RMSE = double(), Corr = double())

  for (i in seq_along(folds_out)) {
    tr <- folds_out[[i]]
    va <- setdiff(seq_len(nrow(train70)), tr)

    Xtr <- as.matrix(train70[tr, ..keep_metab])
    ytr <- train70[[pps_name]][tr]
    Xva <- as.matrix(train70[va, ..keep_metab])
    yva <- train70[[pps_name]][va]

    cvfit <- cv.glmnet(
      Xtr, ytr,
      family      = "gaussian",
      alpha       = alpha_val,
      nfolds      = 10,
      parallel    = TRUE,
      standardize = FALSE          
    )

    lam_min[i] <- cvfit$lambda.min
    lam_1se[i] <- cvfit$lambda.1se

    preds <- predict(cvfit, Xva, s = cvfit$lambda.1se)[, 1]
    fold_stats <- rbind(
      fold_stats,
      data.table(fold = i,
                 RMSE = sqrt(mean((preds - yva)^2)),
                 Corr = cor(preds, yva))
    )

    # Save fold‑specific IDs
    fwrite(data.table(f.eid = train70$f.eid[tr]),
           file = file.path(FOLDS_DIR,
                            sprintf("%s_%s_fold%02d_train_ids_V3.csv", pps_name, var, i)))
    fwrite(data.table(f.eid = train70$f.eid[va]),
           file = file.path(FOLDS_DIR,
                            sprintf("%s_%s_fold%02d_valid_ids_V3.csv", pps_name, var, i)))
  }

  fwrite(fold_stats,
         file = file.path(METR_DIR,
                          sprintf("%s_%s_nested_metrics_V3.csv", pps_name, var)))

  # --- Final models on full 70 % training data
  lam_final_min  <- median(lam_min)
  lam_final_1se  <- median(lam_1se)
  X70 <- as.matrix(train70[, ..keep_metab])
  y70 <- train70[[pps_name]]

  fit_min  <- glmnet(X70, y70, family = "gaussian",
                     alpha = alpha_val, lambda = lam_final_min,
                     standardize = FALSE)
  fit_1se  <- glmnet(X70, y70, family = "gaussian",
                     alpha = alpha_val, lambda = lam_final_1se,
                     standardize = FALSE)

  # --- Save weights
  save_coefs <- function(fit, tag, lam) {
    b  <- as.numeric(coef(fit))
    nz <- which(b != 0)
    dt <- data.table(PPS        = pps_name,
                     Metabolite = rownames(coef(fit))[nz],
                     Weight     = b[nz],
                     Lambda     = lam)
    fwrite(dt,
           file = file.path(WEIG_DIR,
                            sprintf("%s_%s_weights_%s_V3.csv",
                                    pps_name, var, tag)))
    dt
  }
  coefs_min <- save_coefs(fit_min,  "lambda_min",  lam_final_min)
  coefs_1se <- save_coefs(fit_1se,  "lambda_1se", lam_final_1se)

  # ---------- 8 · Apply identical preprocessing to TEST30 ---------------
  params <- readRDS(file.path(PPS_DIR,
                              sprintf("%s_preproc_params_V3.rds", pps_name)))
  test_dt <- copy(test30)
  rin_ref <- params$rin_ref

  # 8.1 Half‑minimum imputation
  test_dt[, (keep_metab) := Map(\(x,m){x[is.na(x)] <- m; x},
                                .SD, params$mins),
          .SDcols = keep_metab]

  # 8.2 Inverse‑rank normal
  apply_rin_from_ref <- function(ref, x) {
    p <- approx(x = ref$vals, y = ref$p, xout = x, method = "linear", rule = 2)$y
    qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
  }
  test_dt[, (keep_metab) := Map(apply_rin_from_ref, rin_ref, .SD), .SDcols = keep_metab]

  # 8.3 Z‑score metabolites
  test_dt[, (keep_metab) := Map(\(x,m,s) (x - m) / s,
                                .SD, params$mu, params$sd),
          .SDcols = keep_metab]

  # 8.4 Z‑score PPS
  test_dt[, (pps_name) := (get(pps_name) - params$pps_mu) / params$pps_sd]

  Xtest        <- as.matrix(test_dt[, ..keep_metab])
  ytest_scaled <- test_dt[[pps_name]]  # scaled

  # ---------- 9 · Evaluate on TEST30 and save scores -------------------
  eval_model <- function(coefs, tag) {
    intercept <- if ("(Intercept)" %in% coefs$Metabolite)
                   coefs[Metabolite == "(Intercept)", Weight] else 0
    w         <- coefs[Metabolite != "(Intercept)"]

    pred_scaled <- as.vector(
      Xtest[, w$Metabolite, drop = FALSE] %*% w$Weight + intercept
    )

    # Back‑transform to original PPS scale
    pred_orig <- pred_scaled * params$pps_sd + params$pps_mu
    obs_orig  <- ytest_scaled  * params$pps_sd + params$pps_mu

    mets <- metrics_fun(pred_orig, obs_orig)

    # Save per‑individual scores
    scores_dt <- data.table(f.eid     = test_dt$f.eid,
                            Observed  = obs_orig,
                            Predicted = pred_orig)
    fwrite(scores_dt,
           file = file.path(SCOR_DIR,
                            sprintf("%s_%s_%s_scores_V3.csv",
                                    pps_name, var, tag)))

    # Append to global results
    results <<- rbind(
      results,
      data.table(PPS        = pps_name,
                 Variant    = var,
                 Lambda     = tag,
                 Predictors = nrow(w),
                 RMSE       = mets$RMSE,
                 MAE        = mets$MAE,
                 Corr       = mets$Corr,
                 R2         = mets$R2)
    )
  }

  eval_model(coefs_min,  "lambda_min")
  eval_model(coefs_1se,  "lambda_1se")

  message(sprintf("[%s] %s: lambda_min=%.4g (%d predictors)  lambda_1se=%.4g (%d predictors)",
                  pps_name, var,
                  lam_final_min,  nrow(coefs_min),
                  lam_final_1se,  nrow(coefs_1se)))
}

# ---------- 10 · Save performance summary -----------------------------------
summary_file <- file.path(PPS_DIR, "performance_summary_V3.csv")
fwrite(results, summary_file)
cat("Performance summary saved to", summary_file, "\n")
