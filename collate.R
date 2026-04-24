#!/usr/bin/env Rscript
# =============================================================================
# Collate ALL PPS outputs (V3)
#   • Combines the 4×25 score files into ONE wide table
#   • Stacks all performance‑summary files into ONE long table
#   • Generates descriptive statistics for every PPS column
# =============================================================================

suppressPackageStartupMessages({
  library(getopt)
  library(data.table)
  setDTthreads(0)
  library(stringr)
})

## ---- 0 · CLI and core paths -----------------------------------------------
spec <- matrix(c(
  "pps-root", "p", 1, "character",
  "help",     "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

usage <- paste(
  "Usage: Rscript collate.R [options]",
  "",
  "Options:",
  "  -p, --pps-root PATH    Directory containing PPS outputs and aggregate score files",
  "                         (default: results/omics_score/pps/disease_free or PPS_ROOT)",
  "  -h, --help             Show this help message",
  sep = "\n"
)

if (!is.null(opt$help)) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

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

PPS_ROOT <- path.expand(arg_or_env(
  "pps-root", "PPS_ROOT",
  file.path("results", "omics_score", "pps", "disease_free")
))
AGG_ROOT <- PPS_ROOT # where the 4 score files were written
OUT_SCORES <- file.path(PPS_ROOT, "pps_scores_ALL_V3.csv")
OUT_METR <- file.path(PPS_ROOT, "metrics_ALL_PPS_V3.csv")
OUT_STAT <- file.path(PPS_ROOT, "pps_scores_summary_V3.csv")

## ---- 1 ·  Merge the 4 aggregated score files ------------------------------
agg_files <- c(
  lasso_min  = "pps_score_lasso_lambda_min_V3.csv",
  lasso_1se  = "pps_score_lasso_lambda_1se_V3.csv",
  enet_min   = "pps_score_enet_lambda_min_V3.csv",
  enet_1se   = "pps_score_enet_lambda_1se_V3.csv"
)
agg_paths <- file.path(AGG_ROOT, agg_files)
names(agg_paths) <- names(agg_files)

if (!all(file.exists(agg_paths))) {
  stop(
    "Missing aggregated score file(s):\n",
    paste(agg_paths[!file.exists(agg_paths)], collapse = "\n")
  )
}

merge_scores <- NULL
for (tag in names(agg_paths)) {
  dt <- fread(agg_paths[tag])
  variant <- str_extract(tag, "lasso|enet")
  lambda <- str_extract(tag, "min|1se")
  pps_cols <- setdiff(names(dt), "f.eid")
  # Rename columns: <pps>__<variant>_<lambda>
  newnames <- paste0(pps_cols, "__", variant, "_", lambda)
  setnames(dt, pps_cols, newnames)
  merge_scores <-
    if (is.null(merge_scores)) {
      dt
    } else {
      merge(merge_scores, dt, by = "f.eid", all = TRUE, sort = FALSE)
    }
}
setcolorder(merge_scores, c("f.eid", sort(setdiff(names(merge_scores), "f.eid"))))
fwrite(merge_scores, OUT_SCORES)
cat("Combined score table: ", OUT_SCORES, "\n")

## ---- 2 ·  Stack all performance‑summary files ------------------------------
perf_files <- list.files(
  PPS_ROOT,
  pattern = "^performance_summary_V3\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
if (!length(perf_files)) {
  stop("No performance_summary_V3.csv files found under ", PPS_ROOT)
}

metrics_all <- rbindlist(lapply(perf_files, fread), use.names = TRUE, fill = TRUE)

# Ensure Variant / Lambda columns are present
if (!"Variant" %in% names(metrics_all)) {
  metrics_all[, Variant := str_extract(perf_files, "(lasso|elastic_net)")[.I]]
}
if (!"Lambda" %in% names(metrics_all)) {
  metrics_all[, Lambda := ifelse(str_detect(perf_files, "lambda_1se"), "lambda_1se", "lambda_min")[.I]]
}

setcolorder(
  metrics_all,
  c(
    "PPS", "Variant", "Lambda",
    setdiff(names(metrics_all), c("PPS", "Variant", "Lambda"))
  )
)
fwrite(metrics_all, OUT_METR)
cat("Combined metrics table: ", OUT_METR, "\n")

## ---- 3 ·  Descriptive statistics per PPS column ---------------------------
num_cols <- setdiff(names(merge_scores), "f.eid")

summary_dt <- rbindlist(lapply(num_cols, function(col) {
  x <- merge_scores[[col]]
  data.table(
    Column = col,
    N      = sum(!is.na(x)),
    Mean   = mean(x, na.rm = TRUE),
    SD     = sd(x, na.rm = TRUE),
    Min    = min(x, na.rm = TRUE),
    Q1     = quantile(x, 0.25, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Q3     = quantile(x, 0.75, na.rm = TRUE),
    Max    = max(x, na.rm = TRUE)
  )
}))
fwrite(summary_dt, OUT_STAT)
cat("Descriptive statistics: ", OUT_STAT, "\n")
