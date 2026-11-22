# S03 06：风险分层及六类干预变化统计（单队列，独立脚本）
# ======================= 前端可调参数 =======================
THREADS <- 5L                 # 默认线程（分组统计可并行）
TEST_N <- 0L                  # 与05一致，0=全量
# ===========================================================
threads <- as.integer(Sys.getenv("S03_THREADS", as.character(THREADS)))
test_n <- as.integer(Sys.getenv("S03_ANALYSIS_TEST_N", as.character(TEST_N)))
stopifnot(threads >= 1L, test_n >= 0L)
arg <- grep("^--file=", commandArgs(), value = TRUE)
fr <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
start <- if (length(fr)) dirname(tail(fr, 1)[[1]]) else if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd()
root <- normalizePath(start, winslash = "/", mustWork = TRUE)
root_override <- Sys.getenv("S03_PROJECT_ROOT", "")
if (nzchar(root_override)) {
  root <- normalizePath(root_override, winslash = "/", mustWork = TRUE)
}

while (!dir.exists(file.path(root, "LJW_code", "Main_code", "M01_ASCVD"))) {
  up <- dirname(root); if (identical(up, root)) stop("找不到项目根目录"); root <- up
}
out <- file.path(root, "LJW_outcome", "data_intervention", "S03_UKHLS")
if (test_n > 0L) out <- file.path(out, paste0("test_", test_n))
if (Sys.getenv("S03_SMOKING_MODE", "same_wave") == "wave2_prior") out <- file.path(out,"wave2_smoking_sensitivity")
for (p in c("complete_case", "logs", "quality_control")) dir.create(file.path(out, p), recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]
logfile <- file.path(out, "logs", "06_Treatment_transitions_S03_runtime.log")
log <- function(...) {msg <- sprintf("[%s +%.1fs] %s", format(Sys.time(), "%H:%M:%S"), proc.time()[["elapsed"]]-t0, paste0(...,collapse=""));message(msg);cat(msg,"\n",file=logfile,append=TRUE)}

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

summarize_level <- function(d, v) {
  labels <- c("L", "M", "H", "VH", "H+VH")
  n <- nrow(d)
  counts <- vapply(labels, function(level) sum(d[[v]] %in% if (level == "H+VH") c("H", "VH") else level), integer(1))
  if (sum(counts[1:4]) != n) stop("存在无效风险等级：", v)
  data.frame(model = v, category = labels, n_total = n, n = unname(counts), percent = 100*unname(counts)/n)
}
summarize_transition <- function(d, level) {
  keep <- d$transition %in% level
  n <- sum(keep)
  v19 <- d$minimum_reduction_2019[keep]
  v25 <- d$minimum_reduction_2025[keep]
  data.frame(transition = level, n_total = nrow(d), n = n, percent = 100*n/nrow(d),
             mean_minimum_reduction_2019 = if(n) mean(v19) else NA_real_,
             mean_minimum_reduction_2025 = if(n) mean(v25) else NA_real_,
             mean_difference_2025_minus_2019 = if(n) mean(v25-v19) else NA_real_)
}
parallel_map <- function(tasks, fn, nthreads, log) {
  if (nthreads <= 1L || length(tasks) <= 1L) return(lapply(seq_along(tasks), function(i) {a <- fn(tasks[[i]]); log("完成",i,"/",length(tasks));a}))
  cl <- parallel::makePSOCKcluster(min(nthreads, length(tasks)))
  on.exit(parallel::stopCluster(cl), add=TRUE)
  parallel::clusterExport(cl, c("d", "summarize_level", "summarize_transition"), envir = .GlobalEnv)
  ans <- vector("list", length(tasks))
  for (inds in split(seq_along(tasks), ceiling(seq_along(tasks)/min(nthreads,length(tasks))))) {
    ans[inds] <- parallel::parLapplyLB(cl, tasks[inds], fn)
    log("完成", max(inds), "/", length(tasks))
  }
  ans
}
path <- file.path(out, "complete_case", "05_score_strategy_participants_S03.csv")
if (!file.exists(path)) stop("缺少05结果：", path)
d <- read.csv(path, na.strings = "")
stopifnot(nrow(d) > 0, !anyDuplicated(d$pidp), !anyNA(d[c("score_2019_category", "score2_2025_category", "transition", "minimum_reduction_2019", "minimum_reduction_2025")]))
log("开始06；样本=",nrow(d),"；线程=",threads)
risk <- do.call(rbind,parallel_map(list("score_2019_category","score2_2025_category"), function(v) summarize_level(d,v), threads,log))
labels <- c("unchanged_untreated", "newly_intervention", "intervention_removed", "treatment_intensified", "unchanged_intensity", "treatment_deintensified")
trans <- do.call(rbind,parallel_map(as.list(labels), function(v) summarize_transition(d,v), threads,log))
if (sum(trans$n) != nrow(d)) stop("六类干预变化计数与总人数不一致")
p1 <- file.path(out,"complete_case","06_risk_level_H_VH_S03.csv")
p2 <- file.path(out,"complete_case","06_six_transitions_LDL_reduction_S03.csv")
s03_write_csv(risk,p1)
s03_write_csv(trans,p2)
log("完成06；总耗时=",round(proc.time()[["elapsed"]]-t0,1),"秒；输出=",p1,"、",p2)
