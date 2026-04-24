#!/usr/bin/env Rscript
# =============================================================================
# PPS Omics Score (V3): proteins + metabolites as predictors
# Variants: LASSO (alpha = 1) and Elastic Net (alpha = 0.5)
# =============================================================================

suppressPackageStartupMessages({
  library(getopt)
  library(data.table); setDTthreads(1)
  library(readr)
  library(fst)
  library(glmnet)
  library(doParallel)
  library(caret)
  library(impute)
})

# ---------- 0 · CLI ----------------------------------------------------------
spec <- matrix(c(
  "idx",           "i", 1, "integer",
  "data-dir",      "d", 1, "character",
  "prot-file",     "p", 1, "character",
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
  "Usage: Rscript pps_omicscore.R --idx <1-n_PPS> [options]",
  "",
  "Required:",
  "  -i, --idx INTEGER             1-based PPS column index from clusters.csv",
  "",
  "Input options:",
  "  -d, --data-dir PATH           Base input directory (default: data or PPS_DATA_DIR)",
  "  -p, --prot-file PATH          Protein .fst file",
  "  -m, --metab-file PATH         NMR/metabolite CSV file",
  "  -c, --clusters-file PATH      PPS clusters CSV file",
  "  -s, --ids-file PATH           ID lists CSV file",
  "",
  "Output/runtime options:",
  "  -o, --out-dir PATH            Output directory (default: results/omics_score/pps/disease_free)",
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
PROT_FILE <- path.expand(arg_or_env(
  "prot-file", "PPS_PROT_FILE",
  file.path(DATA_DIR, "proteins", "proteins_visit_0.fst")
))
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
  file.path("results", "omics_score", "pps", "disease_free")
))

input_files <- c(
  PROT_FILE = PROT_FILE,
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

# ---------- 3 · Load PPS and data -------------------------------------------
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

# --- Proteins
prot <- as.data.table(read.fst(PROT_FILE))
prot_features <- setdiff(names(prot), "f.eid")
na_col_prot   <- colMeans(is.na(prot[, ..prot_features]))
keep_prots    <- names(na_col_prot)[na_col_prot <= 0.10]
na_row_prot   <- rowMeans(is.na(prot[, ..keep_prots]))
prot          <- prot[na_row_prot <= 0.10, c("f.eid", keep_prots), with = FALSE]

# --- Metabolites 
nmr      <- as.data.table(read_csv(METAB_FILE, show_col_types = FALSE))
nmr_main <- nmr[visit == "Main Phase"]
exclude_meta <- c("phase","eid_30418","V2","V3","V4","V5","V6",
                  "f.eid","V8","V9","V10","V11","V12",
                  "sample_id","plate_id","plate_position","visit",
                  "spectrometer","Glucose")
pct_cols       <- grep("_pct$", names(nmr_main), value = TRUE)
metab_features <- setdiff(names(nmr_main), unique(c(exclude_meta, pct_cols)))
na_col_metab   <- colMeans(is.na(nmr_main[, ..metab_features]))
keep_metab     <- names(na_col_metab)[na_col_metab <= 0.10]
na_row_metab   <- rowMeans(is.na(nmr_main[, ..keep_metab]))
metab          <- nmr_main[na_row_metab <= 0.10, c("f.eid", keep_metab), with = FALSE]

# --- Merge: PPS + proteins + metabolites
merged <- Reduce(function(x, y) merge(x, y, by = "f.eid", all = FALSE),
                 list(pps_all, prot, metab))
merged <- merged[!duplicated(f.eid)]
merged <- merged[!is.na(get(pps_name))]

# --- Disease-free filter
ids <- fread(IDS_FILE)
merged <- merged[f.eid %in% ids$disease_free]

# ---------- 4 Split-------------------------
train_ids <- setdiff(ids$protein, ids$diet)
test_ids  <- intersect(ids$protein, ids$diet)

train <- merged[f.eid %in% train_ids]
test  <- merged[f.eid %in% test_ids]

# --- Output dirs + save ID lists --------------------------
PPS_DIR <- file.path(BASE_DIR, pps_name)
IDS_DIR <- file.path(PPS_DIR, "ids")
dir.create(PPS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(IDS_DIR,            showWarnings = FALSE)

fwrite(data.table(f.eid = train$f.eid),
       file = file.path(IDS_DIR, sprintf("%s_train_ids_V3.csv", pps_name)))
fwrite(data.table(f.eid = test$f.eid),
       file = file.path(IDS_DIR, sprintf("%s_test_ids_V3.csv", pps_name)))

# ---------- 5 · Preprocessing (TRAIN only) ----------------------------------
# Rank-normal helpers
build_rin_ref <- function(v) {
  ord <- order(v); x_sorted <- v[ord]
  p <- (seq_along(x_sorted) - 0.5) / length(x_sorted)
  up <- tapply(p, x_sorted, mean)
  list(vals = as.numeric(names(up)), p = as.numeric(up))
}
apply_rin_from_ref <- function(ref, x) {
  p <- approx(x = ref$vals, y = ref$p, xout = x, method = "linear", rule = 2)$y
  qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

# ----- 5A · Proteins 
if (length(keep_prots)) {
  matP_tr <- as.matrix(train[, ..keep_prots])
  impP_tr_t <- impute.knn(t(matP_tr), k = 10)$data
  train[, (keep_prots) := as.data.table(t(impP_tr_t))]
  rinP_ref <- lapply(train[, ..keep_prots], build_rin_ref)
  train[, (keep_prots) := Map(apply_rin_from_ref, rinP_ref, .SD), .SDcols = keep_prots]
  muP <- sapply(train[, ..keep_prots], mean)
  sdP <- sapply(train[, ..keep_prots], sd)
  zeroP <- names(sdP)[sdP == 0]
  if (length(zeroP)) {
    keep_prots <- setdiff(keep_prots, zeroP)
    train[, (zeroP) := NULL]
    muP <- muP[keep_prots]; sdP <- sdP[keep_prots]; rinP_ref <- rinP_ref[keep_prots]
  }
  if (length(keep_prots)) {
    train[, (keep_prots) := Map(function(x,m,s) (x - m)/s, .SD, muP, sdP), .SDcols = keep_prots]
  }
} else { muP <- sdP <- rinP_ref <- list() }

# ----- 5B · Metabolites
if (length(keep_metab)) {
  minsM <- sapply(train[, ..keep_metab], min, na.rm = TRUE) / 2
  train[, (keep_metab) := Map(function(x,m){x[is.na(x)] <- m; x}, .SD, minsM), .SDcols = keep_metab]
  rinM_ref <- lapply(train[, ..keep_metab], build_rin_ref)
  train[, (keep_metab) := Map(apply_rin_from_ref, rinM_ref, .SD), .SDcols = keep_metab]
  muM <- sapply(train[, ..keep_metab], mean)
  sdM <- sapply(train[, ..keep_metab], sd)
  zeroM <- names(sdM)[sdM == 0]
  if (length(zeroM)) {
    keep_metab <- setdiff(keep_metab, zeroM)
    train[, (zeroM) := NULL]
    muM <- muM[keep_metab]; sdM <- sdM[keep_metab]; minsM <- minsM[keep_metab]; rinM_ref <- rinM_ref[keep_metab]
  }
  if (length(keep_metab)) {
    train[, (keep_metab) := Map(function(x,m,s) (x - m)/s, .SD, muM, sdM), .SDcols = keep_metab]
  }
} else { muM <- sdM <- minsM <- rinM_ref <- list() }

# ----- 5C · Scale PPS  --------------------------------------------
pps_mu <- mean(train[[pps_name]])
pps_sd <- sd(train[[pps_name]])
train[, (pps_name) := (get(pps_name) - pps_mu) / pps_sd]

# ----- 5D · Save preprocessing params ---------------------------------------
saveRDS(
  list(
    prot  = list(mu = muP, sd = sdP, rin_ref = rinP_ref, feats = keep_prots),
    metab = list(mu = muM, sd = sdM, mins = minsM, rin_ref = rinM_ref, feats = keep_metab),
    pps_mu = pps_mu, pps_sd = pps_sd
  ),
  file = file.path(PPS_DIR, sprintf("%s_preproc_params_V3.rds", pps_name))
)

# ---------- 6 · Metrics ------------------------------------------------------
metrics_fun <- function(pred, obs) {
  rmse <- sqrt(mean((pred - obs)^2))
  mae  <- mean(abs(pred - obs))
  r    <- ifelse(sd(pred) == 0 || sd(obs) == 0, NA, cor(pred, obs))
  data.table(RMSE = rmse, MAE = mae, Corr = r, R2 = r^2)
}

# ---------- 7 · Nested CV + Models -------------------------
folds_out_global <- createFolds(train[[pps_name]], k = 10, list = TRUE, returnTrain = TRUE)
variants <- list(lasso = 1, elastic_net = 0.5)
results  <- data.table()

for (var in names(variants)) {
  alpha_val <- variants[[var]]
  folds_out <- folds_out_global
  
  VAR_DIR   <- file.path(PPS_DIR, var); dir.create(VAR_DIR, showWarnings = FALSE)
  FOLDS_DIR <- file.path(VAR_DIR, "folds");   dir.create(FOLDS_DIR,   showWarnings = FALSE)
  WEIG_DIR  <- file.path(VAR_DIR, "weights"); dir.create(WEIG_DIR,    showWarnings = FALSE)
  METR_DIR  <- file.path(VAR_DIR, "metrics"); dir.create(METR_DIR,    showWarnings = FALSE)
  SCOR_DIR  <- file.path(VAR_DIR, "scores");  dir.create(SCOR_DIR,    showWarnings = FALSE)
  
  lam_min <- lam_1se <- numeric(length(folds_out))
  fold_stats <- data.table(fold = integer(), RMSE = double(), Corr = double())
  
  feats_all <- c(keep_prots, keep_metab)
  for (i in seq_along(folds_out)) {
    tr <- folds_out[[i]]
    va <- setdiff(seq_len(nrow(train)), tr)
    
    Xtr <- as.matrix(train[tr, ..feats_all])
    ytr <- train[[pps_name]][tr]
    Xva <- as.matrix(train[va, ..feats_all])
    yva <- train[[pps_name]][va]
    
    cvfit <- cv.glmnet(
      Xtr, ytr,
      family = "gaussian",
      alpha = alpha_val,
      nfolds = 10,
      parallel = TRUE,
      standardize = FALSE
    )
    
    lam_min[i] <- cvfit$lambda.min
    lam_1se[i] <- cvfit$lambda.1se
    
    preds <- predict(cvfit, Xva, s = cvfit$lambda.1se)[, 1]
    fold_stats <- rbind(fold_stats,
                        data.table(fold = i,
                                   RMSE = sqrt(mean((preds - yva)^2)),
                                   Corr = cor(preds, yva)))
    
    fwrite(data.table(f.eid = train$f.eid[tr]),
           file = file.path(FOLDS_DIR, sprintf("%s_%s_fold%02d_train_ids_V3.csv", pps_name, var, i)))
    fwrite(data.table(f.eid = train$f.eid[va]),
           file = file.path(FOLDS_DIR, sprintf("%s_%s_fold%02d_valid_ids_V3.csv", pps_name, var, i)))
  }
  
  fwrite(fold_stats,
         file = file.path(METR_DIR, sprintf("%s_%s_nested_metrics_V3.csv", pps_name, var)))
  
  # Final models on full TRAIN
  lam_final_min <- median(lam_min)
  lam_final_1se <- median(lam_1se)
  Xtr_full <- as.matrix(train[, ..feats_all])
  ytr_full <- train[[pps_name]]
  
  fit_min <- glmnet(Xtr_full, ytr_full, family = "gaussian",
                    alpha = alpha_val, lambda = lam_final_min, standardize = FALSE)
  fit_1se <- glmnet(Xtr_full, ytr_full, family = "gaussian",
                    alpha = alpha_val, lambda = lam_final_1se, standardize = FALSE)
  
  # Save weights
  save_coefs <- function(fit, tag, lam) {
    b  <- as.numeric(coef(fit))
    rn <- rownames(coef(fit))
    nz <- which(b != 0)
    feats_nz <- rn[nz]
    type <- ifelse(feats_nz == "(Intercept)", "intercept",
                   ifelse(feats_nz %in% keep_prots, "protein",
                          ifelse(feats_nz %in% keep_metab, "metabolite", "unknown")))
    dt <- data.table(PPS = pps_name, Feature = feats_nz, type = type,
                     Weight = b[nz], Lambda = lam)
    fwrite(dt, file = file.path(WEIG_DIR,
                                sprintf("%s_%s_weights_%s_V3.csv", pps_name, var, tag)))
    dt
  }
  coefs_min <- save_coefs(fit_min, "lambda_min",  lam_final_min)
  coefs_1se <- save_coefs(fit_1se, "lambda_1se",  lam_final_1se)
  
  # ---------- 8 · Apply preprocessing to TEST  -------------
  params <- readRDS(file.path(PPS_DIR, sprintf("%s_preproc_params_V3.rds", pps_name)))
  test_dt <- copy(test)
  
  if (length(params$prot$feats)) {
    fp <- params$prot$feats
    matP_te <- as.matrix(test_dt[, ..fp])
    impP_te_t <- impute.knn(t(matP_te), k = 10)$data
    test_dt[, (fp) := as.data.table(t(impP_te_t))]
    test_dt[, (fp) := Map(apply_rin_from_ref, params$prot$rin_ref, .SD), .SDcols = fp]
    test_dt[, (fp) := Map(function(x,m,s) (x - m)/s, .SD, params$prot$mu, params$prot$sd), .SDcols = fp]
  }
  if (length(params$metab$feats)) {
    fm <- params$metab$feats
    test_dt[, (fm) := Map(function(x,m){x[is.na(x)] <- m; x}, .SD, params$metab$mins), .SDcols = fm]
    test_dt[, (fm) := Map(apply_rin_from_ref, params$metab$rin_ref, .SD), .SDcols = fm]
    test_dt[, (fm) := Map(function(x,m,s) (x - m)/s, .SD, params$metab$mu, params$metab$sd), .SDcols = fm]
  }
  
  test_dt[, (pps_name) := (get(pps_name) - params$pps_mu) / params$pps_sd]
  
  feats_use <- c(params$prot$feats, params$metab$feats)
  Xtest        <- as.matrix(test_dt[, ..feats_use])
  ytest_scaled <- test_dt[[pps_name]]
  
  # ---------- 9 · Evaluate & save scores ------------------------------------
  eval_model <- function(coefs, tag) {
    intercept <- if ("(Intercept)" %in% coefs$Feature)
      coefs[Feature == "(Intercept)", Weight] else 0
    w <- coefs[Feature != "(Intercept)"]
    
    pred_scaled <- as.vector(Xtest[, w$Feature, drop = FALSE] %*% w$Weight + intercept)
    
    pred_orig <- pred_scaled * params$pps_sd + params$pps_mu
    obs_orig  <- ytest_scaled * params$pps_sd + params$pps_mu
    mets <- metrics_fun(pred_orig, obs_orig)
    
    fwrite(data.table(f.eid = test_dt$f.eid, Observed = obs_orig, Predicted = pred_orig),
           file = file.path(SCOR_DIR, sprintf("%s_%s_%s_scores_V3.csv", pps_name, var, tag)))
    
    results <<- rbind(
      results,
      data.table(PPS = pps_name, Variant = var, Lambda = tag,
                 Predictors = nrow(w),
                 RMSE = mets$RMSE, MAE = mets$MAE, Corr = mets$Corr, R2 = mets$R2)
    )
  }
  
  eval_model(coefs_min,  "lambda_min")
  eval_model(coefs_1se,  "lambda_1se")
  
  message(sprintf("[%s] %s: lambda_min=%.4g (%d preds) | lambda_1se=%.4g (%d preds)",
                  pps_name, var,
                  lam_final_min,  sum(coefs_min$Feature != "(Intercept)"),
                  lam_final_1se,  sum(coefs_1se$Feature != "(Intercept)")))
}

# ---------- 10 · Save performance summary -----------------------------------
summary_file <- file.path(PPS_DIR, "performance_summary_V3.csv")
fwrite(results, summary_file)
cat("Performance summary saved to", summary_file, "\n")
