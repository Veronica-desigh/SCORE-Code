# Shared performance and timing utilities for M01 (01--06).
# These functions do not alter any study-specific records or random streams.
# Intended to be sourced by main scripts after creating data_clean_dir.

m01_runtime_init <- function(module, data_clean_dir) {
  dir.create(file.path(data_clean_dir, "quality_control"), recursive = TRUE,
             showWarnings = FALSE)
  options(m01.runtime.module = module,
          m01.runtime.id = paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()),
          m01.runtime.log = file.path(data_clean_dir, "quality_control", "runtime_timing_M01.csv"),
          m01.runtime.started = Sys.time(),
          m01.runtime.progress = !identical(Sys.getenv("M01_PROGRESS", "1"), "0"))
  m01_record_timing("module_started", 0, "START")
  invisible(NULL)
}

m01_record_timing <- function(step, elapsed_sec, status = "OK", note = "") {
  path <- getOption("m01.runtime.log", NULL)
  if (is.null(path)) return(invisible(NULL))
  entry <- data.table::data.table(
    run_id = getOption("m01.runtime.id", NA_character_),
    module = getOption("m01.runtime.module", NA_character_),
    step = as.character(step), status = as.character(status),
    started_or_logged_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    elapsed_seconds = round(as.numeric(elapsed_sec), 3L),
    note = as.character(note))
  # Parent process writes; parallel workers write only their own progress logs.
  data.table::fwrite(entry, path, append = file.exists(path),
                     col.names = !file.exists(path), bom = FALSE)
  invisible(NULL)
}

m01_timed <- function(step, expr) {
  started <- Sys.time()
  cat(sprintf("[%s] [%s] START %s\n", format(started, "%H:%M:%S"),
              getOption("m01.runtime.module", "M01"), step))
  flush.console()
  status <- "ERROR"
  on.exit({
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    m01_record_timing(step, elapsed, status)
    cat(sprintf("[%s] [%s] %s %s | %.2f s (%.2f min)\n",
                format(Sys.time(), "%H:%M:%S"),
                getOption("m01.runtime.module", "M01"), status,
                step, elapsed, elapsed / 60))
    flush.console()
  }, add = TRUE)
  result <- force(expr)
  status <- "OK"
  result
}

m01_checkpoint <- function(step, current, total, detail = "") {
  if (!isTRUE(getOption("m01.runtime.progress", TRUE))) return(invisible(NULL))
  if (!is.numeric(current) || !is.numeric(total) || total <= 0L) return(invisible(NULL))
  cat(sprintf("[%s] [%s] PROGRESS %s: %d/%d (%.1f%%)%s\n",
      format(Sys.time(), "%H:%M:%S"),
      getOption("m01.runtime.module", "M01"), step,
      as.integer(current), as.integer(total), 100 * current / total,
      if (nzchar(detail)) paste0(" | ", detail) else ""))
  flush.console()
  invisible(NULL)
}

m01_configure_threads <- function() {
  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(physical) || physical < 1L)
    physical <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(physical) || physical < 1L) physical <- 2L
  requested <- Sys.getenv("M01_DT_THREADS", "auto")
  if (tolower(requested) == "auto") {
    n <- max(1L, min(8L, as.integer(physical) - 1L))
  } else {
    n <- suppressWarnings(as.integer(requested))
    if (!is.finite(n) || n < 1L) stop("M01_DT_THREADS must be a positive integer or 'auto'")
    n <- min(n, as.integer(physical))
  }
  data.table::setDTthreads(n)
  cat(sprintf("[M01] data.table=%d threads (physical cores=%d; requested=%s)\n",
              data.table::getDTthreads(), as.integer(physical), requested))
  invisible(data.table::getDTthreads())
}

m01_runtime_finish <- function() {
  began <- getOption("m01.runtime.started", Sys.time())
  elapsed <- as.numeric(difftime(Sys.time(), began, units = "secs"))
  m01_record_timing("module_finished", elapsed, "OK")
  cat(sprintf("[%s] [%s] COMPLETED | %.2f min\n",
              format(Sys.time(), "%H:%M:%S"),
              getOption("m01.runtime.module", "M01"), elapsed / 60))
  invisible(NULL)
}
