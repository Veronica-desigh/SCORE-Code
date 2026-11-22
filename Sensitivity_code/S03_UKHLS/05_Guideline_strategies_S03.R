# S03 05：SCORE/SCORE2评分及风险分层（独立脚本）
# ======================= 前端可调参数 =======================
THREADS <- 5L                 # 默认并行上限（评分引擎内部使用）
TEST_N <- 0L                  # 0=全量；>0=抽样测试
RISK_REGION <- "Low"         # 保持既有M01风险区域
LDL_DELTA <- 0.4             # 与原S03一致
LDL_MIN <- 0.6               # 与原S03一致
# ===========================================================
threads <- as.integer(Sys.getenv("S03_THREADS", as.character(THREADS)))
test_n <- as.integer(Sys.getenv("S03_ANALYSIS_TEST_N", as.character(TEST_N)))
stopifnot(is.finite(threads), threads >= 1L, is.finite(test_n), test_n >= 0L)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
frames <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
start <- if (length(frames)) dirname(tail(frames, 1)[[1]]) else
  if (length(script_arg)) dirname(sub("^--file=", "", script_arg[1])) else getwd()
root <- normalizePath(start, winslash = "/", mustWork = TRUE)
root_override <- Sys.getenv("S03_PROJECT_ROOT", "")
if (nzchar(root_override)) {
  root <- normalizePath(root_override, winslash = "/", mustWork = TRUE)
}

while (!dir.exists(file.path(root, "LJW_code", "Main_code", "M01_ASCVD"))) {
  parent <- dirname(root)
  if (identical(root, parent)) stop("无法定位M01代码所在项目根目录")
  root <- parent
}
clean <- file.path(root, "LJW_outcome", "data_clean", "S03_UKHLS")
out <- file.path(root, "LJW_outcome", "data_intervention", "S03_UKHLS")
if (test_n > 0) out <- file.path(out, paste0("test_", test_n))
if (Sys.getenv("S03_SMOKING_MODE", "same_wave") == "wave2_prior") out <- file.path(out,"wave2_smoking_sensitivity")
for (p in c("complete_case", "quality_control", "logs")) dir.create(file.path(out, p), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out, "logs", "05_Guideline_strategies_S03_runtime.log")
t0 <- proc.time()[["elapsed"]]
log <- function(...) {
  msg <- sprintf("[%s +%.1fs] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 proc.time()[["elapsed"]] - t0, paste0(..., collapse = ""))
  message(msg); cat(msg, "\n", file = log_file, append = TRUE)
}

# Excel-compatible UTF-8 CSV with BOM. Never rewrite source data in GBK.
s03_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  bytes <- readBin(tmp, what = "raw", n = file.info(tmp)$size)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(bytes, con)
  invisible(path)
}

score_route_s03 <- function(d, m01, workers, risk_region, ldl_delta, ldl_min, log) {
  if (!nrow(d)) stop("04队列为空")
  score_env <- new.env(parent = .GlobalEnv)
  source(m01, local = score_env, encoding = "UTF-8")
  data.table::setDTthreads(1L)
  x <- data.table::data.table(
    participant_id = d$pidp, age = d$age,
    sex = ifelse(d$sex == 2, "female", "male"),
    smoking_status = ifelse(d$current_smoking == 1, "current", "never"),
    systolic_blood_pressure = d$sbp, diastolic_blood_pressure = d$dbp,
    total_cholesterol = d$tc, hdl_cholesterol = d$hdl, ldl_cholesterol = d$ldl)
  log("校验M01评分输入及资格边界")
  score_env$validate_score_inputs(x)
  score_env$validate_score_eligibility_boundaries(x)
  x <- score_env$derive_current_smoker(x)
  log("开始计算SCORE/SCORE2；线程上限=", workers)
  x <- score_env$calculate_score_values(x, risk_region = risk_region, n_workers = workers)
  x <- score_env$add_score_categories(x)
  log("完成风险评分；计算风险分层和LDL-C必要降幅")
  cfg <- list(ldl_delta = ldl_delta, ldl_min = ldl_min,
              ldl_targets = c(L = 3, M = 2.6, H = 1.8, VH = 1.4))
  a <- score_env$m07_one(x$ldl_cholesterol, x$score_2019_category, rep(.5, nrow(x)), cfg)
  b <- score_env$m07_one(x$ldl_cholesterol, x$score2_2025_category, rep(.5, nrow(x)), cfg)
  tr <- score_env$m07_derive_transition(a$A, b$A, a$D, b$D)
  stopifnot(all(a$D >= 0), all(b$D >= 0), identical(as.numeric(x$participant_id), as.numeric(d$pidp)))
  data.frame(pidp = d$pidp, wave = d$wave, ldl = x$ldl_cholesterol,
             score_2019_percent = x$score_2019_percent,
             score2_2025_percent = x$score2_2025_percent,
             score_2019_category = x$score_2019_category,
             score2_2025_category = x$score2_2025_category,
             target_2019 = a$G, target_2025 = b$G,
             trigger_2019 = a$A, trigger_2025 = b$A,
             minimum_reduction_2019 = a$D, minimum_reduction_2025 = b$D,
             transition = tr)
}

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("RiskScorescvd", quietly = TRUE))
  stop("需要安装 data.table 和 RiskScorescvd")
sample_n <- as.integer(Sys.getenv("S03_SAMPLE_N", "0"))
if (is.na(sample_n) || sample_n < 0L) stop("S03_SAMPLE_N invalid")
if (sample_n > 0L) clean <- file.path(clean, paste0("smoke_test_", sample_n))
smoking_mode <- Sys.getenv("S03_SMOKING_MODE", "same_wave")
if (!smoking_mode %in% c("same_wave","wave2_prior")) stop("Unsupported smoking mode")
if (smoking_mode == "wave2_prior") clean <- file.path(clean, "wave2_smoking_sensitivity")
input <- file.path(clean, "04_eligible_cohort_S03.csv")
if (!file.exists(input)) stop("缺少04完整病例队列：", input)
d <- read.csv(input, na.strings = "", check.names = FALSE)
needed <- c("pidp", "wave", "age", "sex", "current_smoking", "sbp", "dbp", "tc", "hdl", "ldl")
if (!all(needed %in% names(d)) || anyNA(d[needed]) || anyDuplicated(d$pidp)) stop("04完整病例队列字段完整性校验失败")
if (test_n > 0L && sample_n == 0L) {
  set.seed(20261112)
  d <- d[sample.int(nrow(d), min(test_n, nrow(d))), , drop = FALSE]
}
m01 <- file.path(root, "LJW_code", "Main_code", "M01_ASCVD", "Function_07_Guideline_strategies_M01.R")
if (!file.exists(m01)) stop("未找到M01评分引擎：", m01)
log("开始05；参与者=", nrow(d), "；线程=", threads, "；测试上限=", test_n)
a <- score_route_s03(d, m01, threads, RISK_REGION, LDL_DELTA, LDL_MIN, log)
p <- file.path(out, "complete_case", "05_score_strategy_participants_S03.csv")
s03_write_csv(a, p)
s03_write_csv(data.frame(engine = m01, md5 = unname(tools::md5sum(m01)),
                     risk_region = RISK_REGION, strategy = "S0",
                     package_version = as.character(utils::packageVersion("RiskScorescvd"))),
          file.path(out, "quality_control", "05_M01_scoring_engine_S03.csv"))
log("完成05；输出=", p, "；总耗时=", round(proc.time()[["elapsed"]] - t0, 1), "秒")
