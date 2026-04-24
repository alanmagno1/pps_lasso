#!/usr/bin/env Rscript
# =============================================================================
# Aggregate predicted PPS scores (test set) into one file per variant + lambda
# Output files match the names expected by the Cox‑associations script:
#   pps_score_<lasso|enet>_lambda_<min|1se>_V3.csv
# Each file has columns: f.eid + one column per PPS (Predicted values)
# =============================================================================

suppressPackageStartupMessages({
    library(getopt)
    library(data.table)
    setDTthreads(0)
})

## ---- 0 · CLI and paths -----------------------------------------------------
spec <- matrix(c(
    "pps-root", "p", 1, "character",
    "help",     "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

usage <- paste(
    "Usage: Rscript collect.R [options]",
    "",
    "Options:",
    "  -p, --pps-root PATH    Parent directory containing one subdirectory per PPS",
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
OUT_FMT <- "pps_score_%s_lambda_%s_V3.csv" # output name template

## ---- 1 · Variant / lambda mapping ------------------------------------------
variant_dirs <- list( # label in file  → subfolder name inside each PPS
    lasso = "lasso",
    enet  = "elastic_net"
)
lambda_tags <- c( # label → suffix used in filenames
    min = "lambda_min",
    "1se" = "lambda_1se"
)

## ---- 2 · List all PPS directories ------------------------------------------
pps_dirs <- list.dirs(PPS_ROOT, full.names = TRUE, recursive = FALSE)
if (!length(pps_dirs)) {
    stop("No PPS directories found under ", PPS_ROOT)
}

## ---- 3 · Aggregation helper -------------------------------------------------
aggregate_variant <- function(var_label, var_subdir, lam_label, lam_tag) {
    merged <- NULL

    for (pps_dir in pps_dirs) {
        pps_name <- basename(pps_dir) # e.g. fi_7_adiposity_driven_hyperinsulinem
        score_file <- file.path(
            pps_dir, var_subdir, "scores",
            sprintf("%s_%s_%s_scores_V3.csv", pps_name, var_subdir, lam_tag)
        )
        if (!file.exists(score_file)) next

        dt <- fread(score_file, select = c("f.eid", "Predicted"))
        setnames(dt, "Predicted", pps_name)
        merged <- if (is.null(merged)) {
            dt
        } else {
            merge(merged, dt, by = "f.eid", all = TRUE, sort = FALSE)
        }
    }

    if (is.null(merged)) {
        stop("No score files found for variant ", var_label, " / lambda ", lam_label)
    }

    out_path <- file.path(PPS_ROOT, sprintf(OUT_FMT, var_label, lam_label))
    fwrite(merged, out_path)
    message(sprintf(
        "Wrote %s (%d PPS, %d individuals)",
        basename(out_path), ncol(merged) - 1L, nrow(merged)
    ))
}

## ---- 4 · Run for all 4 combinations ----------------------------------------
for (var_label in names(variant_dirs)) {
    for (lam_label in names(lambda_tags)) {
        aggregate_variant(
            var_label, variant_dirs[[var_label]],
            lam_label, lambda_tags[[lam_label]]
        )
    }
}
