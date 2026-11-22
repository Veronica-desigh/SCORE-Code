# S04 | 02 | SCORE / SCORE2 and guideline-based treatment demand
# Run: Rscript 02_Guideline_strategies_S04.R
# Includes setup/helpers/analysis, with the S03 adaptor embedded; the canonical M01 scoring engine is the only external algorithm dependency.
# ============================= COMMON FUNCTIONS (LOCAL TO THIS SCRIPT) =============================
S04_WORKERS <- suppressWarnings(as.integer(Sys.getenv("S04_WORKERS", "5")))
if (is.na(S04_WORKERS) || S04_WORKERS < 1L) stop("S04_WORKERS must be a positive integer")
S04_BATCH_SIZE <- suppressWarnings(as.integer(Sys.getenv("S04_BATCH_SIZE", "1000")))
if (is.na(S04_BATCH_SIZE) || S04_BATCH_SIZE < 1L) stop("S04_BATCH_SIZE must be positive")

script_path <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE))
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(e) e$ofile))
  if (length(frames)) return(normalizePath(tail(frames, 1L)[[1L]], mustWork=FALSE))
  normalizePath(getwd(), mustWork=FALSE)
}
script_dir <- if (dir.exists(script_path())) script_path() else dirname(script_path())
root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
input_dir <- Sys.getenv("S04_INPUT_DIR", file.path(root, "DATA", "ELSA"))
analysis_dir <- Sys.getenv("S04_ANALYSIS_DIR", file.path(root, "LJW_outcome", "data_intervention", "S04_ELSA", "complete_case"))
qc_dir <- Sys.getenv("S04_QC_DIR", file.path(root, "LJW_outcome", "data_clean", "S04_ELSA", "quality_control"))

elapsed_text <- function(seconds) sprintf("%02d:%02d:%05.2f", floor(seconds/3600), floor(seconds%%3600/60), seconds%%60)
start_clock <- Sys.time()
log_message <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep="")
progress <- function(done,total,label="Completed") {
  pct <- if (total) 100*done/total else 100
  log_message(sprintf("%s: %d/%d (%.1f%%); elapsed %s",label,done,total,pct,elapsed_text(as.numeric(difftime(Sys.time(),start_clock,units="secs")))))
}
finish <- function(script_label) log_message(script_label," finished; elapsed ",elapsed_text(as.numeric(difftime(Sys.time(),start_clock,units="secs"))))
write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  write.csv(x,path,row.names=FALSE,na="")
  log_message("Wrote ",path," [",nrow(x)," rows]")
}
read_elsa <- function() {
  dictionary_file <- file.path(input_dir,"ELSA_variable.csv")
  data_file <- file.path(input_dir,"ELSA_data.rds")
  for (p in c(dictionary_file,data_file)) if (!file.exists(p)) stop("Missing ELSA input: ",p)
  dictionary <- read.csv(dictionary_file,fileEncoding="UTF-8-BOM",check.names=FALSE,stringsAsFactors=FALSE)
  needed <- c("output_column_name","source_field")
  if (!all(needed %in% names(dictionary))) stop("Dictionary missing: ",paste(setdiff(needed,names(dictionary)),collapse=", "))
  if (anyNA(dictionary[needed]) || anyDuplicated(dictionary$output_column_name) || anyDuplicated(dictionary$source_field)) stop("Dictionary has missing or duplicate mappings")
  d <- as.data.frame(readRDS(data_file))
  if (!all(names(d) %in% dictionary$output_column_name)) stop("Unmapped ELSA data fields: ",paste(setdiff(names(d),dictionary$output_column_name),collapse=", "))
  names(d) <- dictionary$source_field[match(names(d),dictionary$output_column_name)]
  names(d)[names(d)=="__derived_high_bp__"] <- "high_bp"
  if (!"idauniq" %in% names(d)) stop("Missing participant ID idauniq")
  d$idauniq <- as.character(d$idauniq)
  if (anyNA(d$idauniq) || any(!nzchar(d$idauniq)) || anyDuplicated(d$idauniq)) stop("Missing or duplicate participant identifiers")
  list(data=d, dictionary=dictionary)
}
# Locate the original validated function without assuming a fixed S03 folder layout.
resolve_external_function <- function(env_name, file_name, preferred) {
  override <- Sys.getenv(env_name, unset = "")
  if (nzchar(override)) {
    path <- normalizePath(override, winslash = "/", mustWork = FALSE)
    if (!file.exists(path)) stop(env_name, " points to a nonexistent file: ", path, call. = FALSE)
    log_message("Function path supplied by ", env_name, ": ", path)
    return(path)
  }
  candidates <- unique(c(
    preferred,
    file.path(script_dir, file_name),
    file.path(dirname(script_dir), "S03_UKHLS", file_name),
    file.path(root, "LJW_code", "Sensitivity_code", "S03_UKHLS", file_name),
    file.path(root, "LJW_code", "Main_code", "M01_ASCVD", file_name)
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) {
    path <- normalizePath(existing[1L], winslash = "/")
    log_message("Located function: ", path)
    return(path)
  }
  roots <- unique(c(script_dir, dirname(script_dir), file.path(root, "LJW_code")))
  roots <- roots[dir.exists(roots)]
  hits <- unique(unlist(lapply(roots, function(folder)
    list.files(folder, pattern = paste0("^", file_name, "$"), full.names = TRUE,
               recursive = TRUE, ignore.case = TRUE)), use.names = FALSE))
  if (length(hits) == 1L) {
    path <- normalizePath(hits, winslash = "/")
    log_message("Located function by project search: ", path)
    return(path)
  }
  if (length(hits) > 1L) stop(
    "Multiple possible ", file_name, " files found. Set ", env_name,
    " explicitly to avoid choosing a wrong version:\n", paste(hits, collapse = "\n"), call. = FALSE)
  stop("Cannot locate validated function ", file_name,
       ". No scoring or transition algorithm has been substituted.\n",
       "Check that the S03 function really exists or set ", env_name,
       " to its full filename.\nSearched folders:\n",
       paste(roots, collapse = "\n"), call. = FALSE)
}

# SCORE and treatment-demand adaptor embedded from validated S03 logic;
# only the canonical M01 risk engine is sourced externally (scientific dependency).
# 只适配字段，评分与降脂策略从 M01 原函数计算。
score_route_s04 <- function(d, m01, workers) {
  if (!nrow(d)) stop("队列为空，无法评分")
  score_env <- new.env(parent = .GlobalEnv)
  source(m01, local = score_env, encoding = "UTF-8")
  data.table::setDTthreads(1L)
  x <- data.table::data.table(
    participant_id = d$pidp, age = d$age, sex = ifelse(d$sex == 2, "female", "male"),
    # M01 的风险引擎只用 current_smoker，0 统一传入非当前类别。
    smoking_status = ifelse(d$current_smoking == 1, "current", "never"),
    systolic_blood_pressure = d$sbp, diastolic_blood_pressure = d$dbp,
    total_cholesterol = d$tc, hdl_cholesterol = d$hdl, ldl_cholesterol = d$ldl)
  score_env$validate_score_inputs(x)
  score_env$validate_score_eligibility_boundaries(x)
  x <- score_env$derive_current_smoker(x)
  x <- score_env$calculate_score_values(x, risk_region = "Low", n_workers = workers)
  x <- score_env$add_score_categories(x)
  cfg <- list(ldl_delta = .4, ldl_min = .6, ldl_targets = c(L = 3, M = 2.6, H = 1.8, VH = 1.4))
  # S0 不涉及随机干预后 LDL；D 只取 M01 的最小必要降幅。
  a <- score_env$m07_one(x$ldl_cholesterol, x$score_2019_category, rep(.5, nrow(x)), cfg)
  b <- score_env$m07_one(x$ldl_cholesterol, x$score2_2025_category, rep(.5, nrow(x)), cfg)
  transition <- score_env$m07_derive_transition(a$A, b$A, a$D, b$D)
  stopifnot(all(a$D >= 0), all(b$D >= 0), identical(x$participant_id, d$pidp))
  data.frame(pidp = d$pidp, wave = d$wave, ldl = x$ldl_cholesterol,
             score_2019_percent = x$score_2019_percent, score2_2025_percent = x$score2_2025_percent,
             score_2019_category = x$score_2019_category, score2_2025_category = x$score2_2025_category,
             target_2019 = a$G, target_2025 = b$G, trigger_2019 = a$A, trigger_2025 = b$A,
             minimum_reduction_2019 = a$D, minimum_reduction_2025 = b$D, transition = transition)
}


# ============================= ANALYSIS =============================
log_message("Step 02 starting; configured scoring workers=",S04_WORKERS,", batch size=",S04_BATCH_SIZE)
for (pkg in c("data.table","RiskScorescvd")) if (!requireNamespace(pkg,quietly=TRUE)) stop("Missing R package: ",pkg)
d <- read_elsa()$data
fields <- c("idauniq","age","sex","current_smoking","sbp","dbp","tc","hdl","ldl")
if (!all(fields %in% names(d))) stop("Required fields unavailable: ",paste(setdiff(fields,names(d)),collapse=", "))
cc <- complete.cases(d[,fields,drop=FALSE])
flow <- data.frame(stage=c("Already cleaned source","Complete SCORE/treatment fields"),n=c(nrow(d),sum(cc)),excluded=c(0L,sum(!cc)))
write_csv(flow,file.path(analysis_dir,"02_complete_case_flow_S04.csv"))
d <- d[cc,,drop=FALSE]
if (!nrow(d)) stop("No complete-case observations for scoring")
if (any(!is.finite(d$age)) || any(d$age<50|d$age>69)) stop("Age 50–69 assumption violated; inspect supplied extract")
if (!all(d$sex %in% c("male","female"))) stop("Unexpected sex coding")
if (!all(d$current_smoking %in% c(0,1))) stop("Unexpected smoking coding")
for (v in c("sbp","dbp","tc","hdl","ldl")) if (any(!is.finite(d[[v]]))) stop("Non-finite scoring input: ",v)
# Compatibility fields are not clinical eligibility re-screening.
d$wave <- 6L; d$eligible <- TRUE
d$pidp <- as.character(d$idauniq)
d$sex <- match(d$sex,c("male","female"))
m01 <- resolve_external_function("S04_M01_SCORING_FUNCTION", "Function_07_Guideline_strategies_M01.R", file.path(root,"LJW_code","Main_code","M01_ASCVD","Function_07_Guideline_strategies_M01.R"))
if (!file.exists(m01)) stop("Missing canonical M01 scoring engine: ",m01)
# Batch-level progress reporting. The established S03 scoring function handles
# the canonical risk engine parallel processing using workers=5 by default.
# Smaller batches may add overhead; S04_BATCH_SIZE tunes this tradeoff.
ids <- split(seq_len(nrow(d)),ceiling(seq_len(nrow(d))/S04_BATCH_SIZE))
results <- vector("list",length(ids))
for (i in seq_along(ids)) {
  index <- ids[[i]]
  log_message("Scoring batch ",i,"/",length(ids)," (",length(index)," participants; workers=",S04_WORKERS,")")
  chunk <- score_route_s04(d[index,,drop=FALSE],m01,workers=S04_WORKERS)
  if (!"pidp" %in% names(chunk)) stop("Scoring function did not return pidp")
  if (nrow(chunk)!=length(index) || !identical(as.character(chunk$pidp),as.character(d$pidp[index]))) stop("Batch returned mismatched IDs/order; no output written")
  results[[i]] <- chunk
  progress(sum(lengths(ids[seq_len(i)])),nrow(d),"Scoring")
}
result <- do.call(rbind,results)
row.names(result) <- NULL
names(result)[names(result)=="pidp"] <- "idauniq"
if (!identical(as.character(result$idauniq),as.character(d$idauniq)) || any(!is.finite(result$score_2019_percent)) || any(!is.finite(result$score2_2025_percent))) stop("Scoring output validation failed")
write_csv(result,file.path(analysis_dir,"02_score_strategy_participants_S04.csv"))
settings <- data.frame(n_source=flow$n[1],n_complete=nrow(d),strategy="S0",risk_region="Low",weighted=FALSE,
                      workers=S04_WORKERS,batch_size=S04_BATCH_SIZE,batches=length(ids),
                      m01_md5=unname(tools::md5sum(m01)),
                      risk_scores_cvd_version=as.character(utils::packageVersion("RiskScorescvd")))
write_csv(settings,file.path(analysis_dir,"02_scoring_settings_S04.csv"))
# Manifest to permit exact downstream sample and source verification.
source_manifest <- data.frame(item=c('source_rds_md5','dictionary_md5','scoring_output_md5','n_scored'),
   value=c(unname(tools::md5sum(file.path(input_dir,'ELSA_data.rds'))),
           unname(tools::md5sum(file.path(input_dir,'ELSA_variable.csv'))),
           unname(tools::md5sum(file.path(analysis_dir,'02_score_strategy_participants_S04.csv'))),
           as.character(nrow(result))))
write_csv(source_manifest,file.path(analysis_dir,'02_input_manifest_S04.csv'))
finish("02 guideline strategies")
