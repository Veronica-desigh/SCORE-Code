# Run all S04 stages in order, with a fresh R process for each stage.
run_s04 <- function() {
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  args <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- if (length(frames)) dirname(tail(frames, 1L)[[1L]]) else
    if (length(args)) dirname(sub("^--file=", "", args[1L])) else "."
  root <- normalizePath(file.path(script_dir, "../../.."), winslash = "/", mustWork = TRUE)
  previous <- setwd(root)
  on.exit(setwd(previous), add = TRUE)
  code <- file.path("LJW_code", "Sensitivity_code", "S04_ELSA")
  stages <- c("01_Data_quality_S04.R", "02_Guideline_strategies_S04.R",
              "03_Risk_transitions_S04.R", "04_RM_imputation_distribution_S04.R",
              "05_RM_strategies_S04.R", "06_RM_tables_matrices_contributions_S04.R",
              "07_RM_vector_matrices_S04.R")
  paths <- file.path(code, stages)
  if (!all(file.exists(paths))) stop("Missing S04 scripts: ", paste(paths[!file.exists(paths)], collapse = ", "))
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # Avoid local socket requirements by default. Honor explicit worker settings.
  worker_vars <- c("S04_WORKERS", "S04_RM_THREADS")
  missing_vars <- worker_vars[!nzchar(Sys.getenv(worker_vars))]
  if (length(missing_vars)) {
    previous_workers <- Sys.getenv(missing_vars, unset = NA_character_)
    do.call(Sys.setenv, as.list(setNames(rep("1", length(missing_vars)), missing_vars)))
    on.exit({
      Sys.unsetenv(missing_vars)
      existing <- !is.na(previous_workers)
      if (any(existing)) do.call(Sys.setenv, as.list(previous_workers[existing]))
    }, add = TRUE)
  }
  for (i in seq_along(paths)) {
    message(sprintf("[S04 %d/%d] %s", i, length(paths), stages[i]))
    status <- system2(rscript, shQuote(paths[i]))
    if (status != 0L) stop(sprintf("%s failed (exit %d).", stages[i], status), call. = FALSE)
  }
  message("S04 completed: all seven stages passed.")
}

run_s04()
