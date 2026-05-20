#!/usr/bin/env Rscript
# =============================================================================
# Generate synthetic data for testing the LASSO score pipeline
#
# Creates:
#   data/outcome.csv                    — outcome table (f.eid + bmi)
#   data/proteins/proteins_visit_0.csv  — protein predictors
#   data/metabolites/nmr_threephases.csv — metabolite predictors
#
# Usage:
#   Rscript generate_synthetic_data.R [--n-samples 500] [--n-proteins 50]
#          [--n-metabolites 30] [--seed 42] [--out-dir data]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Parse arguments --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1])
  }
  default
}

n_samples    <- as.integer(get_arg("--n-samples", 10000))
n_proteins   <- as.integer(get_arg("--n-proteins", 2000))
n_metabolites <- as.integer(get_arg("--n-metabolites", 100))
seed_val     <- as.integer(get_arg("--seed", 42))
out_dir      <- get_arg("--out-dir", "data")

set.seed(seed_val)

message(sprintf(
  "Generating synthetic data: %d samples, %d proteins, %d metabolites",
  n_samples, n_proteins, n_metabolites
))

# ---- Helper: inject missingness ---------------------------------------------
inject_na <- function(mat, prob = 0.03) {
  mask <- matrix(runif(length(mat)) < prob, nrow = nrow(mat))
  mat[mask] <- NA
  mat
}

# ---- 1. Latent factors that drive outcome and some predictors ---------------
n_latent <- 5
latent <- matrix(rnorm(n_samples * n_latent), nrow = n_samples)

# ---- 2. Protein predictors --------------------------------------------------
# Some proteins are correlated with latent factors (signal), rest are noise
n_signal_prot <- min(10, n_proteins)
prot_mat <- matrix(rnorm(n_samples * n_proteins), nrow = n_samples)
for (j in seq_len(n_signal_prot)) {
  factor_idx <- ((j - 1) %% n_latent) + 1
  loading <- runif(1, 0.3, 0.8) * sample(c(-1, 1), 1)
  prot_mat[, j] <- loading * latent[, factor_idx] + rnorm(n_samples, sd = 0.5)
}
prot_mat <- inject_na(prot_mat, prob = 0.03)

prot_names <- paste0("protein_", sprintf("%03d", seq_len(n_proteins)))
prot_dt <- data.table(f.eid = 1000L + seq_len(n_samples))
prot_dt[, (prot_names) := as.data.table(prot_mat)]

# ---- 3. Metabolite predictors -----------------------------------------------
n_signal_met <- min(8, n_metabolites)
met_mat <- matrix(abs(rnorm(n_samples * n_metabolites, mean = 5, sd = 2)),
                  nrow = n_samples)
for (j in seq_len(n_signal_met)) {
  factor_idx <- ((j - 1) %% n_latent) + 1
  loading <- runif(1, 0.2, 0.6)
  met_mat[, j] <- abs(loading * latent[, factor_idx] + rnorm(n_samples, mean = 5, sd = 1))
}
met_mat <- inject_na(met_mat, prob = 0.02)

met_names <- paste0("metabolite_", sprintf("%03d", seq_len(n_metabolites)))
# Add a few _pct columns so --drop-pct-cols can be tested
pct_names <- paste0("metabolite_", sprintf("%03d", seq_len(min(3, n_metabolites))), "_pct")
pct_mat <- matrix(runif(n_samples * length(pct_names), 0, 100), nrow = n_samples)

met_dt <- data.table(f.eid = 1000L + seq_len(n_samples))
met_dt[, visit := "Main Phase"]
met_dt[, phase := "phase_1"]
met_dt[, sample_id := paste0("S", seq_len(n_samples))]
met_dt[, plate_id := sample(paste0("PLT", 1:5), n_samples, replace = TRUE)]
met_dt[, plate_position := sample(paste0(LETTERS[1:8], sprintf("%02d", 1:12)), n_samples, replace = TRUE)]
met_dt[, spectrometer := sample(c("spec_A", "spec_B"), n_samples, replace = TRUE)]
met_dt[, (met_names) := as.data.table(met_mat)]
met_dt[, (pct_names) := as.data.table(pct_mat)]

# Add a second visit so --metabolite-visit-col filtering can be tested
met_dt2 <- copy(met_dt)
met_dt2[, visit := "Repeat Phase"]
# Perturb values slightly
for (col in met_names) {
  met_dt2[, (col) := get(col) + rnorm(.N, sd = 0.3)]
}
met_dt <- rbindlist(list(met_dt, met_dt2))

# ---- 4. Outcome (BMI) -------------------------------------------------------
# BMI = linear combination of latent factors + noise
bmi_weights <- c(1.5, -0.8, 0.6, -0.3, 0.4)
bmi <- 25 + latent %*% bmi_weights + rnorm(n_samples, sd = 2)

outcome_dt <- data.table(
  f.eid = 1000L + seq_len(n_samples),
  bmi = round(as.numeric(bmi), 1)
)
# Inject a few missing outcomes
na_idx <- sample(seq_len(n_samples), size = max(1, round(n_samples * 0.01)))
outcome_dt[na_idx, bmi := NA]

# ---- 5. Write files ---------------------------------------------------------
dir.create(file.path(out_dir, "proteins"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "metabolites"), recursive = TRUE, showWarnings = FALSE)

outcome_path <- file.path(out_dir, "outcome.csv")
protein_path <- file.path(out_dir, "proteins", "proteins_visit_0.csv")
metabolite_path <- file.path(out_dir, "metabolites", "nmr_threephases.csv")

fwrite(outcome_dt, outcome_path)
fwrite(prot_dt, protein_path)
fwrite(met_dt, metabolite_path)

message(sprintf("Outcome:      %s  (%d rows, %d cols)", outcome_path, nrow(outcome_dt), ncol(outcome_dt)))
message(sprintf("Proteins:     %s  (%d rows, %d cols)", protein_path, nrow(prot_dt), ncol(prot_dt)))
message(sprintf("Metabolites:  %s  (%d rows, %d cols)", metabolite_path, nrow(met_dt), ncol(met_dt)))
message("Done. You can now run the pipeline with:")
message("")
message("  # Protein-only")
message(sprintf("  Rscript lasso_score.R --outcome-file %s --outcome-col bmi --mode proteins --protein-file %s", outcome_path, protein_path))
message("")
message("  # Metabolite-only")
message(sprintf("  Rscript lasso_score.R --outcome-file %s --outcome-col bmi --mode metabolites --metabolite-file %s --metabolite-exclude-cols phase,sample_id,plate_id,plate_position,visit,spectrometer --metabolite-visit-col visit --metabolite-visit-value 'Main Phase' --drop-pct-cols", outcome_path, metabolite_path))
message("")
message("  # Combined")
message(sprintf("  Rscript lasso_score.R --outcome-file %s --outcome-col bmi --mode combined --protein-file %s --metabolite-file %s --metabolite-exclude-cols phase,sample_id,plate_id,plate_position,visit,spectrometer --metabolite-visit-col visit --metabolite-visit-value 'Main Phase' --drop-pct-cols", outcome_path, protein_path, metabolite_path))
