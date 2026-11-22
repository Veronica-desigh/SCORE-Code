# ============================================================================
# S04 / ELSA | 04 | Risk-modifier missing-data multiple imputation & diagnostics
# One analysis, ONE self-contained R script. Run after 02_Guideline_strategies_S04.R
# Existing 01-03 inputs, scores, sample membership and treatment transitions
# are READ-ONLY. No cleaning, re-screening or re-scoring occurs in Step 04.
# Only Wave 6 hs-CRP is used (no Wave 4 / persistence information).
# ============================================================================

# -------------------------- 01. USER PARAMETERS ------------------------------
THREADS <- 5L                     # independent parallel MI chains (PSOCK)
M <- 5L                           # number of imputed data sets
MAXIT <- 5L                       # MICE iterations within each chain
SEED <- 20261113L
MIN_DONORS <- 15L                 # minimum observed values for an imputation target
MIN_BINARY_CLASS <- 5L            # minimum observed 0 and 1 for binary logreg
MIN_CONTINUOUS_UNIQUE <- 5L       # minimum distinct observations for PMM
P_MM_DONORS <- 5L                 # predictive mean matching donors
# Run only a sample for a smoke test (0 = original full fixed analysis cohort).
TEST_N <- 0L
# ----------------------------------------------------------------------------

THREADS <- as.integer(Sys.getenv("S04_RM_THREADS", as.character(THREADS)))
M <- as.integer(Sys.getenv("S04_RM_M", as.character(M)))
MAXIT <- as.integer(Sys.getenv("S04_RM_MAXIT", as.character(MAXIT)))
TEST_N <- as.integer(Sys.getenv("S04_RM_TEST_N", as.character(TEST_N)))
stopifnot(!anyNA(c(THREADS, M, MAXIT, TEST_N)), THREADS >= 1L,
          M >= 1L, MAXIT >= 1L, TEST_N >= 0L,
          MIN_DONORS >= 2L, MIN_BINARY_CLASS >= 2L,
          MIN_CONTINUOUS_UNIQUE >= 2L)

FIXED <- c("age", "sex", "current_smoking", "tc", "hdl", "ldl", "sbp", "dbp")
RM_BINARY_COMMON <- c("rm_social_deprivation_proxy", "rm_major_psychiatric_proxy",
                      "rm_physical_inactivity_proxy", "rm_rheumatoid_proxy")
RM_MENOPAUSE <- "rm_natural_menopause_lt40_proxy"
RM_DERIVED <- c("rm_obesity_central", "rm_hscrp_gt2")
RM_ALL <- c(RM_BINARY_COMMON, RM_DERIVED, RM_MENOPAUSE)
RM_SOURCE_CONT <- c("bmi", "waist", "hscrp_w6")
TARGET_MAIN <- c(RM_SOURCE_CONT, RM_BINARY_COMMON)
TARGET_ALL <- c(TARGET_MAIN, RM_MENOPAUSE)
RM_2019 <- c("rm_social_deprivation_proxy", "rm_major_psychiatric_proxy",
             "rm_obesity_central", "rm_physical_inactivity_proxy", "rm_rheumatoid_proxy")
RM_2025 <- c(RM_2019, "rm_hscrp_gt2", RM_MENOPAUSE)
# The lists above represent ONLY the ELSA-assessable RM subset, not every
# modifier specified in either clinical guideline. Unavailable RMs are not 0.

script_location <- function() {
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(args)) return(normalizePath(sub("^--file=", "", args[1L]), mustWork = FALSE))
  fr <- Filter(Negate(is.null), lapply(sys.frames(), function(e) e$ofile))
  if (length(fr)) return(normalizePath(tail(fr, 1L)[[1L]], mustWork = FALSE))
  normalizePath(getwd(), mustWork = FALSE)
}
path0 <- script_location()
script_dir <- if (dir.exists(path0)) path0 else dirname(path0)
root_override <- Sys.getenv("S04_PROJECT_ROOT", "")
root <- if (nzchar(root_override)) normalizePath(root_override, mustWork = FALSE) else
  normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
input_dir <- Sys.getenv("S04_INPUT_DIR", file.path(root, "DATA", "ELSA"))
analysis_dir <- Sys.getenv("S04_ANALYSIS_DIR", file.path(root, "LJW_outcome", "data_intervention", "S04_ELSA", "complete_case"))
out_dir <- Sys.getenv("S04_RM_OUT_DIR", file.path(root, "LJW_outcome", "data_intervention", "S04_ELSA", "rm_imputation"))
if (TEST_N > 0L) out_dir <- file.path(out_dir, sprintf("test_%d", TEST_N))
for (folder in c(out_dir, file.path(out_dir, "quality_control"),
                 file.path(out_dir, "per_imputation"), file.path(out_dir, "logs"))) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
log_path <- file.path(out_dir, "logs", "04_RM_runtime_S04.log")
start_clock <- Sys.time()

log_msg <- function(...) {
  secs <- as.numeric(difftime(Sys.time(), start_clock, units = "secs"))
  txt <- sprintf("[%s +%.1fs] %s", format(Sys.time(), "%H:%M:%S"), secs,
                 paste0(..., collapse = ""))
  message(txt)
  cat(txt, "\n", file = log_path, append = TRUE)
  invisible(txt)
}
write_out <- function(data, filename, folder = out_dir) {
  path <- file.path(folder, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # UTF-8 with BOM, to preserve Chinese labels in Excel on Windows.
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  con1 <- file(tmp, open = "rb")
  on.exit(close(con1), add = TRUE)
  bytes <- readBin(con1, what = "raw", n = file.info(tmp)$size)
  con2 <- file(path, open = "wb")
  on.exit(close(con2), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con2)
  writeBin(bytes, con2)
  log_msg("Saved ", filename, ": ", nrow(data), " rows")
  invisible(path)
}
read_elsa <- function() {
  dict_path <- file.path(input_dir, "ELSA_variable.csv")
  rds_path <- file.path(input_dir, "ELSA_data.rds")
  if (!file.exists(dict_path) || !file.exists(rds_path)) {
    stop("Required ELSA source missing. S04_INPUT_DIR = ", input_dir,
         "; files: ELSA_data.rds and ELSA_variable.csv", call. = FALSE)
  }
  dict <- read.csv(dict_path, fileEncoding = "UTF-8-BOM", check.names = FALSE,
                   stringsAsFactors = FALSE)
  required <- c("output_column_name", "source_field", "ELSA_support")
  if (!all(required %in% names(dict))) stop("Dictionary missing: ",
    paste(setdiff(required, names(dict)), collapse = ", "))
  if (anyNA(dict[c("output_column_name", "source_field")]) ||
      anyDuplicated(dict$output_column_name) || anyDuplicated(dict$source_field)) {
    stop("ELSA dictionary has missing/duplicated mappings")
  }
  data <- as.data.frame(readRDS(rds_path))
  if (!all(names(data) %in% dict$output_column_name)) {
    stop("ELSA data have unmapped fields: ",
         paste(setdiff(names(data), dict$output_column_name), collapse = ", "))
  }
  names(data) <- dict$source_field[match(names(data), dict$output_column_name)]
  names(data)[names(data) == "__derived_high_bp__"] <- "high_bp"
  if (!"idauniq" %in% names(data)) stop("Missing idauniq")
  data$idauniq <- as.character(data$idauniq)
  if (anyNA(data$idauniq) || any(!nzchar(data$idauniq)) ||
      anyDuplicated(data$idauniq)) stop("Invalid/duplicate ELSA participant IDs")
  list(data = data, dict = dict)
}
as_binary <- function(x, label) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  if (is.character(x)) {
    good <- is.na(x) | x %in% c("0", "1")
    if (!all(good)) stop("Non-0/1 binary value in ", label)
    x <- as.integer(x)
  }
  if (!is.numeric(x) || any(!is.na(x) & !(x %in% c(0, 1)))) {
    stop("Expected 0/1/NA in ", label)
  }
  as.integer(x)
}
check_continuous <- function(x, v) {
  if (!is.numeric(x)) stop("Expected a numeric source: ", v)
  if (any(!is.finite(x) & !is.na(x))) stop("Infinite value in ", v)
  if (v == "hscrp_w6" && any(!is.na(x) & x < 0)) stop("Negative hs-CRP in source")
  if (v %in% c("bmi", "waist") && any(!is.na(x) & x <= 0)) {
    stop("Non-positive observed ", v, "; input should already be cleaned")
  }
  x
}
resolve_obesity <- function(d, observed_rm) {
  cut <- ifelse(d$sex == "male", 102, 88)
  yes <- (!is.na(d$bmi) & d$bmi >= 30) |
         (!is.na(d$waist) & d$waist >= cut)
  no <- !is.na(d$bmi) & !is.na(d$waist) & d$bmi < 30 & d$waist < cut
  inferred <- ifelse(yes, 1L, ifelse(no, 0L, NA_integer_))
  ifelse(!is.na(observed_rm), observed_rm, inferred)
}
resolve_hscrp <- function(d, observed_rm) {
  inferred <- ifelse(is.na(d$hscrp_w6), NA_integer_,
                     as.integer(d$hscrp_w6 > 2))
  ifelse(!is.na(observed_rm), observed_rm, inferred)
}
applicable_rows <- function(d, v) {
  if (identical(v, RM_MENOPAUSE)) d$sex == "female" else rep(TRUE, nrow(d))
}
# For each variable, publish observed/imputed counts and distribution on
# its applicable population. MENOPAUSE has female denominator only.
summarize_one <- function(d, v, stage, imp = 0L) {
  applicable <- applicable_rows(d, v)
  x <- d[[v]][applicable]
  x <- if (v %in% c(RM_ALL, "current_smoking")) as_binary(x, v) else x
  observed <- x[!is.na(x)]
  continuous <- v %in% c(RM_SOURCE_CONT, "age", "tc", "hdl", "ldl", "sbp", "dbp")
  numerical <- if (continuous) as.numeric(observed) else numeric(0)
  npos <- if (v %in% c(RM_ALL, "current_smoking")) sum(observed == 1L) else NA_integer_
  nneg <- if (v %in% c(RM_ALL, "current_smoking")) sum(observed == 0L) else NA_integer_
  q <- if (continuous && length(numerical))
    unname(stats::quantile(numerical, c(.25, .5, .75))) else rep(NA_real_, 3L)
  data.frame(imputation_id = imp, stage = stage, variable = v,
             n_total = nrow(d), n_applicable = length(x),
             n_not_applicable = sum(!applicable), n_observed = length(observed),
             n_missing = sum(is.na(x)),
             missing_pct = if (length(x)) 100 * mean(is.na(x)) else NA_real_,
             n_positive = npos, n_negative = nneg,
             positive_pct_observed = if (!is.na(npos) && length(observed))
               100 * npos / length(observed) else NA_real_,
             mean = if (length(numerical)) mean(numerical) else NA_real_,
             sd = if (length(numerical) > 1L) stats::sd(numerical) else NA_real_,
             q25 = q[1L], median = q[2L], q75 = q[3L],
             min = if (length(numerical)) min(numerical) else NA_real_,
             max = if (length(numerical)) max(numerical) else NA_real_)
}
summarize_set <- function(d, variables, stage, imp = 0L) {
  do.call(rbind, lapply(variables, function(v) summarize_one(d, v, stage, imp)))
}
category_set <- function(d, stage, imp = 0L) {
  vars <- c("sex", "current_smoking", RM_ALL)
  do.call(rbind, lapply(vars, function(v) {
    applicable <- applicable_rows(d, v)
    x <- as.character(d[[v]][applicable])
    levels <- if (v == "sex") c("female", "male") else c("0", "1")
    tab <- table(factor(x[!is.na(x)], levels = levels))
    data.frame(imputation_id = imp, stage = stage, variable = v,
               category = levels, n_applicable = length(x),
               n_observed = sum(!is.na(x)), n_missing = sum(is.na(x)),
               n = as.integer(tab),
               pct_observed = if (sum(!is.na(x))) 100 * as.integer(tab) / sum(!is.na(x)) else NA_real_)
  }))
}
# In observed-data support decisions, TRUE means the variable has usable
# donors; √ / △ in a dictionary alone is not enough evidence for estimation.
support_one <- function(d, v, available, dict_support, female_only = FALSE) {
  x <- d[[v]][if (female_only) d$sex == "female" else rep(TRUE, nrow(d))]
  obs <- x[!is.na(x)]
  binary <- v %in% c(RM_BINARY_COMMON, RM_MENOPAUSE)
  n0 <- if (binary) sum(obs == 0L) else NA_integer_
  n1 <- if (binary) sum(obs == 1L) else NA_integer_
  distinct <- length(unique(obs))
  ok <- available && length(obs) >= MIN_DONORS &&
    (if (binary) n0 >= MIN_BINARY_CLASS && n1 >= MIN_BINARY_CLASS else
       distinct >= MIN_CONTINUOUS_UNIQUE)
  reason <- if (!available) "not_in_source" else if (!length(obs))
    "no_observed_donors" else if (anyNA(x) && !ok)
    "insufficient_donors_or_binary_class" else if (!anyNA(x))
    "complete_no_imputation_needed" else "model_supported"
  data.frame(variable = v, elsa_feasibility = dict_support,
             field_in_source = available, applicable_n = length(x),
             observed_n = length(obs), missing_n = sum(is.na(x)),
             observed_zeros = n0, observed_ones = n1,
             distinct_observed = distinct,
             imputation_supported = ok, active = ok && anyNA(x),
             reason = reason, stringsAsFactors = FALSE)
}
# Captures per-chain, per-variable logged statistics; not a formal proof of
# convergence, particularly with M=5 and MAXIT=5.
trace_table <- function(arr, var_arr, imp_id, model_stage) {
  if (is.null(arr) || length(dim(arr)) != 3L) return(NULL)
  ix <- which(!is.na(arr), arr.ind = TRUE)
  if (!nrow(ix)) return(NULL)
  d <- dimnames(arr)
  vn <- if (!is.null(d[[1L]])) d[[1L]][ix[, 1L]] else as.character(ix[, 1L])
  data.frame(imputation_id = imp_id, stage = model_stage,
             variable = vn, iteration = ix[, 2L], chain = ix[, 3L],
             chain_mean = as.numeric(arr[ix]),
             chain_variance = if (!is.null(var_arr)) as.numeric(var_arr[ix]) else NA_real_)
}
# Observed and imputed-only comparisons should prompt review, not be treated
# as a mechanical reason to change the analysis model or drop observations.
evaluate_target <- function(original, current, v, imp) {
  applicable <- applicable_rows(original, v)
  before <- original[[v]][applicable]
  after <- current[[v]][applicable]
  obs <- before[!is.na(before)]
  filled <- after[is.na(before) & !is.na(after)]
  binary <- v %in% c(RM_ALL)
  if (binary) {
    m0 <- if (length(obs)) mean(obs == 1L) else NA_real_
    m1 <- if (length(filled)) mean(filled == 1L) else NA_real_
    delta <- if (is.finite(m0) && is.finite(m1)) 100 * (m1 - m0) else NA_real_
    flag <- is.finite(delta) && abs(delta) > 20
    kind <- "binary_positive_rate"
  } else {
    m0 <- if (length(obs)) mean(obs) else NA_real_
    m1 <- if (length(filled)) mean(filled) else NA_real_
    s0 <- if (length(obs) > 1L) stats::sd(obs) else NA_real_
    delta <- if (is.finite(s0) && s0 > 0 && is.finite(m1)) (m1 - m0) / s0 else NA_real_
    flag <- is.finite(delta) && abs(delta) > 1
    kind <- "continuous_standardized_mean_difference"
  }
  data.frame(imputation_id = imp, variable = v, measure = kind,
             observed_n = length(obs), newly_completed_n = length(filled),
             remaining_missing_n = sum(is.na(after)),
             observed_value = m0, completed_only_value = m1,
             difference = delta, review_large_shift = flag)
}
rm_joint_status <- function(d, vars, label, imp) {
  # "No known positive" is NOT equivalent to "all guideline RMs negative".
  xx <- as.matrix(d[, vars, drop = FALSE])
  # Do not count male menopause N/A as unresolved, or convert the stored
  # male RM values to zeros. Applicability is handled only in this summary.
  app <- matrix(TRUE, nrow = nrow(xx), ncol = ncol(xx))
  if (RM_MENOPAUSE %in% vars) app[d$sex == "male", match(RM_MENOPAUSE, vars)] <- FALSE
  known_pos <- rowSums((xx == 1) & app, na.rm = TRUE)
  unresolved <- rowSums(is.na(xx) & app)
  outcome <- ifelse(known_pos > 0L, "confirmed_positive",
                    ifelse(unresolved == 0L, "all_assessable_negative", "unresolved"))
  tab <- table(factor(outcome,
    levels = c("confirmed_positive", "all_assessable_negative", "unresolved")))
  data.frame(imputation_id = imp, guideline = label,
             status = names(tab), n = as.integer(tab),
             pct = 100 * as.integer(tab) / nrow(d),
             mean_known_positive = mean(known_pos),
             mean_unresolved = mean(unresolved))
}
correlation_rows <- function(d, imp) {
  combos <- utils::combn(RM_ALL, 2L, simplify = FALSE)
  do.call(rbind, lapply(combos, function(v) {
    keep <- applicable_rows(d, v[1L]) & applicable_rows(d, v[2L]) &
      !is.na(d[[v[1L]]]) & !is.na(d[[v[2L]]])
    a <- d[[v[1L]]][keep]; b <- d[[v[2L]]][keep]
    r <- if (length(a) >= 10L && length(unique(a)) > 1L &&
             length(unique(b)) > 1L) stats::cor(a, b) else NA_real_
    data.frame(imputation_id = imp, rm1 = v[1L], rm2 = v[2L],
               n_pair_observed = sum(keep), phi = r)
  }))
}
# Parallel worker. Two phases are used for the sex-specific menopause RM:
# (1) joint MICE for the seven common targets, using other observed/imputed
#     ELSA RM sources as mutual predictors;
# (2) MICE among women only for menopause, conditioned on the completed
#     common RM sources. Men are never coded as missing women / filled with 0.
# Menopause is deliberately not a predictor in phase (1), because values
# for men are structurally not applicable and could not enter an all-sex FCS.
run_chain <- function(task, spec) {
  chain_t0 <- Sys.time()
  worker_log <- function(stage) {
    sec <- round(as.numeric(difftime(Sys.time(), chain_t0, units = "secs")), 1)
    msg <- sprintf("[%s] Chain %d/%d: %s (%.1f s)",
                   format(Sys.time(), "%H:%M:%S"), task$id, spec$m, stage, sec)
    cat(msg, "\n", file = file.path(spec$progress_dir,
                  sprintf("04_RM_chain_%02d_progress_S04.log", task$id)), append = TRUE)
    message(msg)
  }
  worker_log("started")
  if (length(spec$active_main)) {
    # Per-target PMM settings via blots; never forward PMM-only arguments to
    # Logistic imputation methods.
    pmmt <- intersect(spec$active_main, c("bmi", "waist", "hscrp_w6"))
    blots <- stats::setNames(lapply(pmmt, function(v)
      list(donors = spec$donors)), pmmt)
    fit <- mice::mice(spec$model, m = 1L, maxit = spec$maxit,
                      method = spec$methods, predictorMatrix = spec$predmat,
                      where = spec$where, seed = task$seed, blots = blots,
                      printFlag = FALSE,
                      remove.constant = TRUE, remove.collinear = TRUE)
    vals <- mice::complete(fit, 1L)
    main_trace <- fit$chainMean
    main_trace_var <- fit$chainVar
    main_events <- fit$loggedEvents
  } else {
    vals <- spec$model
    main_trace <- NULL; main_trace_var <- NULL; main_events <- NULL
  }
  # Only return imputed targets; the master process protects all source data.
  main_filled <- list()
  for (v in spec$active_main) {
    z <- vals[[v]]
    if (is.factor(z)) z <- as.integer(as.character(z))
    main_filled[[v]] <- z
  }
  worker_log("common RM imputation completed")
  menopause_filled <- NULL
  menopause_trace <- NULL; menopause_trace_var <- NULL
  menopause_events <- NULL; menopause_predictors <- character()
  if (spec$menopause_active) {
    women <- spec$women
    df <- spec$female_base
    for (v in spec$active_main) df[[v]] <- main_filled[[v]][women]
    # Remaining unsupported/incomplete variables cannot serve as predictors.
    candidate <- setdiff(names(df), spec$menopause_var)
    candidate <- candidate[vapply(df[candidate], function(z)
      !anyNA(z) && length(unique(z)) > 1L, logical(1))]
    if (!length(candidate)) stop("No usable predictors for female menopause")
    menopause_predictors <- candidate
    fm <- df[, c(candidate, spec$menopause_var), drop = FALSE]
    for (v in intersect(candidate, spec$binary_vars)) {
      fm[[v]] <- factor(fm[[v]], levels = 0:1)
    }
    fm[[spec$menopause_var]] <- factor(fm[[spec$menopause_var]], levels = 0:1)
    for (v in intersect(c("sex", "current_smoking"), names(fm)))
      fm[[v]] <- factor(fm[[v]])
    meth <- stats::setNames(rep("", ncol(fm)), names(fm))
    meth[spec$menopause_var] <- "logreg"
    pm <- matrix(0L, ncol(fm), ncol(fm), dimnames = list(names(fm), names(fm)))
    pm[spec$menopause_var, candidate] <- 1L
    wh <- is.na(fm)
    wh[, candidate] <- FALSE
    fitf <- mice::mice(fm, m = 1L, maxit = spec$maxit, method = meth,
                       predictorMatrix = pm, where = wh,
                       seed = task$seed + 50000L, printFlag = FALSE,
                       remove.constant = TRUE, remove.collinear = TRUE)
    zz <- mice::complete(fitf, 1L)[[spec$menopause_var]]
    menopause_filled <- if (is.factor(zz)) as.integer(as.character(zz)) else zz
    menopause_trace <- fitf$chainMean
    menopause_trace_var <- fitf$chainVar
    menopause_events <- fitf$loggedEvents
  }
  worker_log("female-specific model / all calculations completed")
  list(id = task$id, seed = task$seed, main = main_filled,
       menopause = menopause_filled, female_predictors = menopause_predictors,
       trace_main = main_trace, tracevar_main = main_trace_var,
       events_main = main_events, trace_female = menopause_trace,
       tracevar_female = menopause_trace_var, events_female = menopause_events)
}

main <- function() {
  on.exit(log_msg("S04 Step 04 exit; elapsed = ",
                  round(as.numeric(difftime(Sys.time(), start_clock, units = "secs")), 1),
                  " seconds"), add = TRUE)
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required; install.packages('mice')")
  }
  log_msg("Start S04 RM: M=", M, ", maxit=", MAXIT, ", threads=", THREADS,
          ", seed=", SEED, ", test_n=", TEST_N)
  loaded <- read_elsa()
  raw <- loaded$data
  dict <- loaded$dict
  score_path <- file.path(analysis_dir, "02_score_strategy_participants_S04.csv")
  if (!file.exists(score_path)) stop("Run Step 02 first: missing ", score_path)
  score <- read.csv(score_path, colClasses = c(idauniq = "character"),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!"idauniq" %in% names(score) || !nrow(score) ||
      anyNA(score$idauniq) || anyDuplicated(score$idauniq)) {
    stop("Step 02 output has invalid IDs")
  }
  ix <- match(score$idauniq, raw$idauniq)
  if (anyNA(ix)) stop("Step 02 participant ID does not exist in ELSA_data.rds")
  d <- raw[ix, , drop = FALSE]
  rownames(d) <- NULL
  stopifnot(identical(d$idauniq, score$idauniq), !anyDuplicated(d$idauniq))
  if (TEST_N > 0L) {
    set.seed(SEED)
    chosen <- sort(sample.int(nrow(d), min(TEST_N, nrow(d))))
    d <- d[chosen, , drop = FALSE]
    score <- score[chosen, , drop = FALSE]
    rownames(d) <- NULL; rownames(score) <- NULL
  }
  if (!all(FIXED %in% names(d))) stop("Missing fixed predictors: ",
     paste(setdiff(FIXED, names(d)), collapse = ", "))
  if (anyNA(d[FIXED])) stop("Step 02 analysis sample has missing fixed predictors")
  if (!all(d$sex %in% c("female", "male"))) stop("Unexpected sex coding")
  if (!all(d$current_smoking %in% 0:1)) stop("Unexpected smoking values")
  for (v in setdiff(FIXED, "sex")) {
    if (!is.numeric(d[[v]]) || any(!is.finite(d[[v]]))) {
      stop("Fixed predictor not finite numeric: ", v)
    }
  }
  if ("ldl" %in% names(score)) {
    if (anyNA(score$ldl) || any(abs(score$ldl - d$ldl) > 1e-8)) {
      stop("Step 02 LDL values differ from source; do not proceed")
    }
  }
  log_msg("Fixed analysis cohort matched: n=", nrow(d),
          "; female=", sum(d$sex == "female"))
  # Dictionary feasibility √ or △ determines RM candidacy, NOT donor sufficiency.
  dict_var <- dict$source_field
  dict_var[dict_var == "__derived_high_bp__"] <- "high_bp"
  selected <- unique(c(RM_ALL, RM_SOURCE_CONT))
  if (!all(selected %in% dict_var)) stop("ELSA dictionary missing specified RM: ",
       paste(setdiff(selected, dict_var), collapse = ", "))
  support_map <- stats::setNames(as.character(dict$ELSA_support), dict_var)
  if (!all(support_map[selected] %in% c("√", "△"))) {
    stop("Some selected RM are not √ or △ in the feasibility dictionary")
  }
  provided <- stats::setNames(selected %in% names(d), selected)
  missing_source_fields <- setdiff(selected, names(d))
  for (v in missing_source_fields) {
    d[[v]] <- if (v %in% RM_ALL) NA_integer_ else NA_real_
  }
  for (v in RM_ALL) d[[v]] <- as_binary(d[[v]], v)
  for (v in RM_SOURCE_CONT) d[[v]] <- check_continuous(d[[v]], v)
  source_original <- d
  before <- summarize_set(d, unique(c(FIXED, RM_SOURCE_CONT, RM_ALL)), "before_imputation")
  write_out(before, "04_variable_distribution_missingness_before_S04.csv")
  write_out(category_set(d, "before_imputation"),
            "04_categorical_distribution_before_S04.csv")

  # Deterministic logic is applied only where original RM is missing.
  d$rm_obesity_central <- as.integer(resolve_obesity(d, source_original$rm_obesity_central))
  d$rm_hscrp_gt2 <- as.integer(resolve_hscrp(d, source_original$rm_hscrp_gt2))
  after_deterministic <- summarize_set(d, unique(c(FIXED, RM_SOURCE_CONT, RM_ALL)),
                                       "after_deterministic")
  write_out(after_deterministic, "04_variable_distribution_after_deterministic_S04.csv")
  write_out(category_set(d, "after_deterministic"),
            "04_categorical_distribution_after_deterministic_S04.csv")
  mismatch <- data.frame(variable = RM_DERIVED, mismatch_observed = c(
    sum(!is.na(source_original$rm_obesity_central) &
        !is.na(resolve_obesity(source_original, rep(NA_integer_, nrow(d)))) &
        source_original$rm_obesity_central !=
          resolve_obesity(source_original, rep(NA_integer_, nrow(d))), na.rm = TRUE),
    sum(!is.na(source_original$rm_hscrp_gt2) &
        !is.na(source_original$hscrp_w6) &
        source_original$rm_hscrp_gt2 != as.integer(source_original$hscrp_w6 > 2),
        na.rm = TRUE)))
  write_out(mismatch, "04_derived_RM_consistency_S04.csv",
            file.path(out_dir, "quality_control"))
  if (any(mismatch$mismatch_observed > 0L)) {
    log_msg("WARNING: Original derived RM disagree with observed source values: ",
            paste(mismatch$variable[mismatch$mismatch_observed > 0L], collapse = ", "),
            ". Original observed RM values are preserved; investigate QC file.")
  }

  # All common target variables are evaluated for inclusion. No 2019/2025
  # unavailable RM will be imputed or used as a predictor.
  sup <- do.call(rbind, lapply(TARGET_ALL, function(v) {
    support_one(d, v, available = unname(provided[v]),
                dict_support = unname(support_map[v]),
                female_only = identical(v, RM_MENOPAUSE))
  }))
  write_out(sup, "04_RM_imputation_support_S04.csv",
            file.path(out_dir, "quality_control"))
  observable <- do.call(rbind, lapply(RM_ALL, function(v) {
    applicable <- applicable_rows(d, v)
    x <- d[[v]][applicable]
    data.frame(variable = v, elsa_feasibility = unname(support_map[v]),
               source_field_present = unname(provided[v]),
               n_applicable = length(x), n_positive = sum(x == 1L, na.rm = TRUE),
               n_negative = sum(x == 0L, na.rm = TRUE), n_missing = sum(is.na(x)),
               observed_support = sum(!is.na(x)) > 0L,
               imputation_role = if (v %in% RM_DERIVED) "derived_from_continuous_source" else
                 if (v == RM_MENOPAUSE) "female_only_target" else "common_binary_target")
  }))
  write_out(observable, "04_RM_observability_S04.csv",
            file.path(out_dir, "quality_control"))
  active_main <- intersect(TARGET_MAIN, sup$variable[sup$active])
  menopause_active <- isTRUE(sup$active[sup$variable == RM_MENOPAUSE])
  log_msg("Common RM targets to impute: ", if (length(active_main))
          paste(active_main, collapse = ", ") else "none")
  log_msg("Female menopause target: ", if (menopause_active) "active" else "not supported / no missing")

  # Matrix of fixed observed clinical variables + all usable ELSA common RM
  # source/indicator variables; composite obesity/CRP are not included beside
  # their deterministic parents (BMI/waist/CRP), to avoid circular equations.
  eligible_common <- TARGET_MAIN[vapply(TARGET_MAIN, function(v) {
    row <- sup[sup$variable == v, , drop = FALSE]
    (nrow(row) == 1L) && (isTRUE(row$imputation_supported) ||
      (!anyNA(d[[v]]) && length(unique(d[[v]])) > 1L))
  }, logical(1))]
  model_vars <- unique(c(FIXED, eligible_common))
  model <- d[, model_vars, drop = FALSE]
  for (v in intersect(RM_BINARY_COMMON, model_vars))
    model[[v]] <- factor(model[[v]], levels = 0:1)
  model$sex <- factor(model$sex, levels = c("female", "male"))
  model$current_smoking <- factor(model$current_smoking, levels = 0:1)
  methods <- stats::setNames(rep("", length(model_vars)), model_vars)
  pm <- matrix(0L, length(model_vars), length(model_vars),
               dimnames = list(model_vars, model_vars))
  where <- is.na(model)
  where[, FIXED] <- FALSE
  for (v in active_main) {
    methods[v] <- if (v %in% RM_BINARY_COMMON) "logreg" else "pmm"
    pred <- setdiff(model_vars, v)
    pred <- pred[vapply(model[pred], function(x)
      length(unique(stats::na.omit(x))) > 1L, logical(1))]
    # All remaining RM sources are either fully observed or themselves
    # imputed in the same MICE chain; incomplete unsupported columns are absent.
    pm[v, pred] <- 1L
  }
  # All targets not scheduled for MI and all fixed clinical fields stay fixed.
  write_out(data.frame(variable = rownames(pm), method = unname(methods),
                       pm, check.names = FALSE),
            "04_RM_predictor_matrix_common_S04.csv",
            file.path(out_dir, "quality_control"))
  write_out(data.frame(variable = model_vars,
                       type = ifelse(model_vars %in% FIXED, "fixed_predictor", "ELSA_RM_source_or_proxy"),
                       included = TRUE), "04_RM_model_variables_S04.csv",
            file.path(out_dir, "quality_control"))

  # The menopause model is female-only. It uses other RM values completed in
  # the same chain; male records keep their original N/A meaning.
  women <- which(d$sex == "female")
  female_candidates <- unique(c(FIXED, TARGET_MAIN))
  female_base <- d[women, unique(c(female_candidates, RM_MENOPAUSE)), drop = FALSE]
  spec <- list(model = model, methods = methods, predmat = pm, where = where,
               active_main = active_main, female_base = female_base,
               women = women, menopause_var = RM_MENOPAUSE,
               menopause_active = menopause_active, maxit = MAXIT,
               donors = P_MM_DONORS, binary_vars = RM_BINARY_COMMON,
               m = M, progress_dir = file.path(out_dir, "logs"))
  tasks <- lapply(seq_len(M), function(i) list(id = i, seed = SEED + i))
  write_out(data.frame(imputation_id = seq_len(M), seed = SEED + seq_len(M),
                       iterations = MAXIT, threads = THREADS, donors = P_MM_DONORS,
                       cohort_n = nrow(d), test_n = TEST_N),
            "04_RM_run_settings_S04.csv", file.path(out_dir, "quality_control"))

  workers <- min(THREADS, M)
  cl <- NULL
  if (workers > 1L && (length(active_main) || menopause_active)) {
    # outfile='' keeps worker messages visible on platforms supporting PSOCK
    # console forwarding; independent per-chain files always record progress.
    cl <- parallel::makePSOCKcluster(workers, outfile = "")
    on.exit({ if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
    parallel::clusterEvalQ(cl, requireNamespace("mice", quietly = TRUE))
    parallel::clusterExport(cl, "run_chain", envir = environment(run_chain))
  }
  chains <- vector("list", M)
  groups <- split(seq_len(M), ceiling(seq_len(M) / workers))
  for (ids in groups) {
    if (!is.null(cl)) {
      chains[ids] <- parallel::parLapplyLB(cl, tasks[ids], run_chain, spec = spec)
    } else {
      chains[ids] <- lapply(tasks[ids], run_chain, spec = spec)
    }
    log_msg("MICE chain progress: ", max(ids), "/", M,
            " (", round(100 * max(ids) / M, 1), "%).")
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
    # on.exit stopCluster is harmlessly redundant for a stopped cluster.
  }

  completed <- vector("list", M)
  evaluation <- list(); dist_after <- list(); cats_after <- list()
  checks <- list(); joint <- list(); correlations <- list()
  traces <- list(); events <- list(); female_matrices <- list()
  derived_after_qc <- list()
  all_describe <- unique(c(FIXED, RM_SOURCE_CONT, RM_ALL))
  for (i in seq_len(M)) {
    ch <- chains[[i]]
    x <- d
    for (v in active_main) {
      fill <- ch$main[[v]]
      if (length(fill) != nrow(x)) stop("Chain returned wrong length for ", v)
      cells <- is.na(x[[v]])
      x[[v]][cells] <- fill[cells]
      if (anyNA(x[[v]])) stop("Active target remained missing: ", v, "; chain ", i)
    }
    if (menopause_active) {
      fill <- ch$menopause
      if (length(fill) != length(women)) stop("Incorrect menopause chain length")
      cells <- is.na(x[[RM_MENOPAUSE]][women])
      x[[RM_MENOPAUSE]][women[cells]] <- fill[cells]
      if (anyNA(x[[RM_MENOPAUSE]][women])) stop("Female menopause MI incomplete")
    }
    x$rm_obesity_central <- as.integer(resolve_obesity(x, source_original$rm_obesity_central))
    x$rm_hscrp_gt2 <- as.integer(resolve_hscrp(x, source_original$rm_hscrp_gt2))
    for (v in RM_ALL) x[[v]] <- as_binary(x[[v]], v)
    for (v in intersect(active_main, RM_SOURCE_CONT)) {
      obs_range <- range(source_original[[v]], na.rm = TRUE)
      if (any(!is.na(x[[v]]) & (x[[v]] < obs_range[1L] | x[[v]] > obs_range[2L]))) {
        stop("PMM output exceeded source donor range: ", v, "; chain ", i)
      }
    }
    # Distinguish an imputed continuous measurement disagreeing with a
    # previously observed binary derived RM; keep the observed RM unchanged.
    unob <- resolve_obesity(x, rep(NA_integer_, nrow(x)))
    uncrp <- resolve_hscrp(x, rep(NA_integer_, nrow(x)))
    derived_after_qc[[i]] <- data.frame(
      imputation_id = i, variable = RM_DERIVED,
      n_original_RM_inconsistent_with_completed_source = c(
        sum(!is.na(source_original$rm_obesity_central) & !is.na(unob) &
            source_original$rm_obesity_central != unob, na.rm = TRUE),
        sum(!is.na(source_original$rm_hscrp_gt2) & !is.na(uncrp) &
            source_original$rm_hscrp_gt2 != uncrp, na.rm = TRUE)))
    # Every observed original value must remain exactly the same, including
    # female RM, SCORE clinical predictors and unrelated source columns.
    if (!identical(x$idauniq, source_original$idauniq) ||
        !identical(x$sex, source_original$sex)) stop("ID/sex changed in MI")
    for (v in names(source_original)) {
      a <- source_original[[v]]; b <- x[[v]]
      known <- !is.na(a)
      if (!identical(b[known], a[known])) {
        stop("Observed value was changed in chain ", i, ": ", v)
      }
      if (!v %in% c(TARGET_ALL, RM_DERIVED) && !identical(a, b)) {
        stop("A fixed/non-target source variable changed: ", v)
      }
    }
    if (!identical(x[[RM_MENOPAUSE]][d$sex == "male"],
                   source_original[[RM_MENOPAUSE]][d$sex == "male"])) {
      stop("Non-applicable male menopause statuses were changed")
    }
    # Attach the frozen Step 02 participant-level SCORE/treatment values.
    # Never recompute or overwrite the source clinical variables.
    joined_cols <- setdiff(names(score), c("idauniq", names(x)))
    for (v in joined_cols) x[[v]] <- score[[v]]
    if (!identical(x$idauniq, score$idauniq)) stop("Step 02 scores misaligned")
    completed[[i]] <- x
    dist_after[[i]] <- summarize_set(x, all_describe, "post_imputation", i)
    cats_after[[i]] <- category_set(x, "post_imputation", i)
    for (v in unique(c(TARGET_ALL, RM_DERIVED))) {
      evaluation[[length(evaluation) + 1L]] <- evaluate_target(source_original, x, v, i)
      valid <- applicable_rows(d, v)
      checks[[length(checks) + 1L]] <- data.frame(
        imputation_id = i, variable = v,
        n_observed_protected = sum(!is.na(source_original[[v]]) & valid),
        n_completed_deterministically = sum(is.na(source_original[[v]]) &
                                            !is.na(d[[v]]) & valid),
        n_completed_by_MI_or_MI_source = sum(is.na(d[[v]]) &
                                            !is.na(x[[v]]) & valid),
        n_completed_total = sum(is.na(source_original[[v]]) & !is.na(x[[v]]) & valid),
        n_missing_after = sum(is.na(x[[v]]) & valid),
        n_not_applicable = sum(!valid),
        n_observed_changed = 0L)
    }
    joint[[length(joint) + 1L]] <- rbind(
      rm_joint_status(x, RM_2019, "2019_ELSA_assessable_subset", i),
      rm_joint_status(x, RM_2025, "2025_ELSA_assessable_subset", i))
    correlations[[i]] <- correlation_rows(x, i)
    if (!is.null(ch$trace_main)) traces[[length(traces) + 1L]] <-
      trace_table(ch$trace_main, ch$tracevar_main, i, "common_RMs")
    if (!is.null(ch$trace_female)) traces[[length(traces) + 1L]] <-
      trace_table(ch$trace_female, ch$tracevar_female, i, "female_menopause")
    if (!is.null(ch$events_main) && nrow(ch$events_main)) {
      events[[length(events) + 1L]] <- data.frame(imputation_id = i,
        stage = "common_RMs", ch$events_main, check.names = FALSE)
    }
    if (!is.null(ch$events_female) && nrow(ch$events_female)) {
      events[[length(events) + 1L]] <- data.frame(imputation_id = i,
        stage = "female_menopause", ch$events_female, check.names = FALSE)
    }
    if (menopause_active) {
      female_matrices[[i]] <- data.frame(imputation_id = i,
        target = RM_MENOPAUSE, predictor = ch$female_predictors)
    }
    log_msg("Completed dataset / QC ", i, "/", M)
  }
  names(completed) <- sprintf("imp%02d", seq_len(M))
  saveRDS(completed, file.path(out_dir, "04_RM_imputed_cohorts_S04.rds"), compress = "xz")
  write_out(do.call(rbind, dist_after),
            "04_variable_distribution_missingness_after_per_imputation_S04.csv",
            file.path(out_dir, "per_imputation"))
  write_out(do.call(rbind, cats_after),
            "04_categorical_distribution_after_per_imputation_S04.csv",
            file.path(out_dir, "per_imputation"))
  pooled <- stats::aggregate(
    do.call(rbind, dist_after)[c("n_observed", "n_missing", "missing_pct",
                                 "n_positive", "n_negative", "positive_pct_observed",
                                 "mean", "sd", "median", "q25", "q75")],
    by = do.call(rbind, dist_after)["variable"],
    FUN = function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE))
  pooled$m <- M
  pooled$n_total <- nrow(d)
  post_all <- do.call(rbind, dist_after)
  between <- stats::aggregate(post_all$positive_pct_observed,
    by = list(variable = post_all$variable),
    FUN = function(z) if (sum(!is.na(z)) < 2L) NA_real_ else stats::sd(z, na.rm = TRUE))
  names(between)[2L] <- "between_imputation_sd_positive_pct"
  pooled <- merge(pooled, between, by = "variable", all.x = TRUE, sort = FALSE)
  write_out(pooled, "04_variable_distribution_after_mean_S04.csv")
  write_out(do.call(rbind, evaluation),
            "04_imputed_value_distribution_evaluation_S04.csv",
            file.path(out_dir, "quality_control"))
  # Compact, publication-QA oriented overview: original vs MI missingness,
  # positive prevalence, and between-imputation variation (descriptive).
  rm_pre <- before[before$variable %in% RM_ALL,
                   c("variable", "n_applicable", "n_missing", "missing_pct",
                     "positive_pct_observed"), drop = FALSE]
  names(rm_pre)[names(rm_pre) == "n_missing"] <- "missing_before_n"
  names(rm_pre)[names(rm_pre) == "missing_pct"] <- "missing_before_pct"
  names(rm_pre)[names(rm_pre) == "positive_pct_observed"] <- "positive_before_pct"
  rm_post <- pooled[pooled$variable %in% RM_ALL,
                    c("variable", "n_missing", "missing_pct", "positive_pct_observed",
                      "between_imputation_sd_positive_pct"), drop = FALSE]
  names(rm_post)[names(rm_post) == "n_missing"] <- "missing_after_mean_n"
  names(rm_post)[names(rm_post) == "missing_pct"] <- "missing_after_mean_pct"
  names(rm_post)[names(rm_post) == "positive_pct_observed"] <- "positive_after_mean_pct"
  eval_summary <- merge(rm_pre, rm_post, by = "variable", sort = FALSE)
  ev <- do.call(rbind, evaluation)
  ev <- ev[ev$variable %in% RM_ALL, , drop = FALSE]
  ev_flag <- stats::aggregate(as.integer(ev$review_large_shift),
           by = list(variable = ev$variable), FUN = sum)
  names(ev_flag)[2L] <- "n_chains_with_large_observed_imputed_shift"
  eval_summary <- merge(eval_summary, ev_flag, by = "variable", all.x = TRUE, sort = FALSE)
  eval_summary <- merge(eval_summary, observable[c("variable", "observed_support")],
                        by = "variable", all.x = TRUE, sort = FALSE)
  write_out(eval_summary, "04_RM_imputation_evaluation_summary_S04.csv")
  if (any(eval_summary$n_chains_with_large_observed_imputed_shift > 0L, na.rm = TRUE)) {
    log_msg("Review recommended: at least one RM exhibits a large observed-vs-imputed shift.")
  }
  write_out(do.call(rbind, checks), "04_RM_imputation_checks_S04.csv",
            file.path(out_dir, "quality_control"))
  write_out(do.call(rbind, derived_after_qc),
            "04_derived_RM_consistency_after_S04.csv",
            file.path(out_dir, "quality_control"))
  if (any(do.call(rbind, derived_after_qc)$n_original_RM_inconsistent_with_completed_source > 0L)) {
    log_msg("WARNING: Some completed continuous values disagree with original observed binary RMs; ",
            "see 04_derived_RM_consistency_after_S04.csv")
  }
  joint_pre <- rbind(rm_joint_status(source_original, RM_2019,
                                   "2019_ELSA_assessable_subset", 0L),
                     rm_joint_status(source_original, RM_2025,
                                   "2025_ELSA_assessable_subset", 0L))
  write_out(rbind(joint_pre, do.call(rbind, joint)),
            "04_RM_any_positive_and_unresolved_S04.csv")
  write_out(rbind(correlation_rows(source_original, 0L),
                  do.call(rbind, correlations)),
            "04_RM_pairwise_correlations_S04.csv",
            file.path(out_dir, "quality_control"))
  if (length(traces)) write_out(do.call(rbind, traces),
                               "04_RM_chain_trace_S04.csv",
                               file.path(out_dir, "quality_control"))
  if (length(events)) {
    write_out(do.call(rbind, events), "04_RM_logged_events_S04.csv",
              file.path(out_dir, "quality_control"))
  } else {
    write_out(data.frame(status = "no_logged_events"),
              "04_RM_logged_events_S04.csv", file.path(out_dir, "quality_control"))
  }
  if (length(female_matrices)) write_out(do.call(rbind, female_matrices),
                                        "04_RM_predictor_matrix_female_menopause_S04.csv",
                                        file.path(out_dir, "quality_control"))
  # Detailed metadata: test mode, model options, original input and versions.
  metadata <- data.frame(
    item = c("analysis", "cohort_n", "m", "maxit", "threads", "seed",
             "source_rds_md5", "source_dictionary_md5", "fixed_score_csv_md5",
             "mice_version", "R_version", "hscrp_wave",
             "female_menopause_model", "clinical_scores_recomputed",
             "status"),
    value = c("ELSA_S04_RM_MICE", nrow(d), M, MAXIT, workers, SEED,
              unname(tools::md5sum(file.path(input_dir, "ELSA_data.rds"))),
              unname(tools::md5sum(file.path(input_dir, "ELSA_variable.csv"))),
              unname(tools::md5sum(score_path)),
              as.character(utils::packageVersion("mice")), R.version.string,
              "wave6_only", "female_only_conditional_on_common_RM_imputations",
              "FALSE", if (TEST_N > 0L) "TEST_ONLY" else "FULL_COHORT"))
  metadata <- rbind(metadata, data.frame(item='imputed_cohorts_md5',
      value=unname(tools::md5sum(file.path(out_dir, '04_RM_imputed_cohorts_S04.rds')))))
  write_out(metadata, "04_RM_run_metadata_S04.csv",
            file.path(out_dir, "quality_control"))
  log_msg("SUCCESS. All observed fields and frozen SCORE data protected; ",
          "M=", M, " completed cohorts. Total runtime = ",
          round(as.numeric(difftime(Sys.time(), start_clock, units = "secs")), 1), " s")
}

main()
