# S03 第02步 Clean_nurse：单文件独立运行（功能函数与主程序合并）
# 默认5线程；可在下方 THREADS 参数或环境变量 S03_THREADS 调整。
# 测试参数沿用 S03_SAMPLE_N。
THREADS <- 5L  # 修改此处设置默认并行线程数

# 只读取所需字段。保留原始数值供缺失原因统计。
read_s03 <- function(path, columns) {
  header <- names(read.csv(path, nrows = 0L, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  absent <- setdiff(columns, header)
  if (length(absent)) stop("缺少字段 ", paste(absent, collapse = ", "), "：", path)
  read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM",
           colClasses = ifelse(header %in% columns, "numeric", "NULL"))
}

clean_number <- function(x, minimum = 0, allowed = NULL) {
  valid <- is.finite(x) & x >= minimum
  if (!is.null(allowed)) valid <- valid & x %in% allowed
  replace(x, !valid, NA_real_)
}

clean_date <- function(year, month, day) {
  year <- clean_number(year, allowed = 1900:2100)
  month <- clean_number(month, allowed = 1:12)
  day <- clean_number(day, allowed = 1:31)
  value <- sprintf("%04d-%02d-%02d", year, month, day)
  date <- as.Date(value, format = "%Y-%m-%d")
  replace(date, is.na(date) | format(date, "%Y-%m-%d") != value, as.Date(NA))
}

# 汇总源文件负值编码及空白；保留每种编码的数量以便回查值标签。
missing_audit <- function(raw, source) {
  do.call(rbind, lapply(names(raw), function(name) {
    value <- raw[[name]]
    code <- ifelse(is.na(value), "blank", ifelse(value < 0, as.character(value), "nonnegative"))
    counts <- table(code)
    data.frame(source = rep(source, length(counts)), variable = rep(name, length(counts)),
               raw_code = names(counts), n = as.integer(counts))
  }))
}

# 冒烟测试以护士记录为抽样框，两个波次均抽取，所有脚本共享同一批 pidp。
sample_ids <- function(path, n) {
  if (n == 0L) return(NULL)
  d <- read_s03(path, c("pidp", "wave"))
  set.seed(20261111)
  unique(unlist(lapply(split(d$pidp[d$pidp > 0], d$wave[d$pidp > 0]), function(ids) {
    ids <- unique(ids)
    ids[sample.int(length(ids), min(n, length(ids)))]
  }), use.names = FALSE))
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

# Windows/Linux 通用：仅向PSOCK子进程导出函数，避免重新执行主脚本。
s03_parallel <- function(jobs, worker, threads, function_path = NULL, progress = NULL) {
  if (!length(jobs)) return(list())
  nw <- min(as.integer(threads), length(jobs))
  run_one <- function(j) worker(j)
  if (nw == 1L) {
    ans <- vector("list", length(jobs))
    for (i in seq_along(jobs)) {
      ans[[i]] <- run_one(jobs[[i]])
      if (!is.null(progress)) progress(i, length(jobs))
    }
    return(ans)
  }
  cl <- parallel::makePSOCKcluster(nw, setup_timeout = 30)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  funcs <- ls(.GlobalEnv, all.names = TRUE)
  funcs <- funcs[vapply(funcs, function(nm) is.function(get(nm, .GlobalEnv)), logical(1))]
  parallel::clusterExport(cl, funcs, envir = .GlobalEnv)
  ans <- vector("list", length(jobs))
  # 每批完成后在主进程报告进度；Windows 下同样可以观察。
  for (indices in split(seq_along(jobs), ceiling(seq_along(jobs) / nw))) {
    done <- parallel::parLapplyLB(cl, jobs[indices], worker)
    ans[indices] <- done
    if (!is.null(progress)) progress(max(indices), length(jobs))
  }
  ans
}

write_s03 <- function(data, path) s03_write_csv(data, path)
# 按行切分独立清理任务，合并时恢复输入顺序。
split_rows <- function(d, threads) {
  n <- nrow(d)
  if (!n) return(list(d))
  group <- ceiling(seq_len(n) / ceiling(n / threads))
  split(d, factor(group, levels = unique(group)))
}

nurse_columns <- function() {
  c("pidp", "wave", "nurdayy", "nurdaym", "nurdayd", "omsysval", "omdiaval",
    "bnf7_statins", "medcnjd", "statins", "statina_med",
    "bmival", "wstval", "bpmedc", "psu", "strata", "indnsub_xw")
}

clean_nurse <- function(d) {
  # statina_med 是非处方他汀追问。其阴性无法单独排除处方他汀。
  bnf <- clean_number(d$bnf7_statins, allowed = 0:1)
  reported <- clean_number(d$statina_med, allowed = 1:2)
  # 7251 护士问卷：无用药且非处方他汀为否时，BNF 不适用属于已核实的跳题。
  # 非处方筛问为是时，以其过去 7 天的回答判定这一分支。
  otc_negative <- d$statins %in% 2 | (d$statins %in% 1 & reported %in% 2)
  confirmed_skip <- d$bnf7_statins %in% -8 & d$medcnjd %in% 2 & otc_negative
  statin <- ifelse(bnf %in% 1 | reported %in% 1, 1,
                   ifelse(bnf %in% 0 | confirmed_skip, 0, NA_real_))
  data.frame(
    pidp = clean_number(d$pidp, 1), wave = clean_number(d$wave, allowed = 2:3),
    baseline_date = clean_date(d$nurdayy, d$nurdaym, d$nurdayd),
    sbp = clean_number(d$omsysval, 1), dbp = clean_number(d$omdiaval, 1),
    statin_7d = statin, bnf_statin_7d = bnf, self_report_statin_7d = reported,
    bmi = clean_number(d$bmival, 1), waist_cm = clean_number(d$wstval, 1),
    bp_medication = clean_number(d$bpmedc, allowed = 0:1),
    psu = clean_number(d$psu), strata = clean_number(d$strata),
    nurse_weight = clean_number(d$indnsub_xw)
  )
}

# ------------------------------ 主程序 ------------------------------
# 删除当前环境中的所有对象
# 从项目根目录或任意子目录启动，也支持 Rscript 传入绝对脚本路径。
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
source_frames <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
start <- if (length(source_frames)) dirname(tail(source_frames, 1L)[[1L]]) else
  if (length(script_arg)) dirname(sub("^--file=", "", script_arg[1L])) else getwd()
root <- normalizePath(start, winslash = "/", mustWork = TRUE)
root_override <- Sys.getenv("S03_PROJECT_ROOT", "")
if (nzchar(root_override)) {
  root <- normalizePath(root_override, winslash = "/", mustWork = TRUE)
}

while (!dir.exists(file.path(root, "DATA", "UKHLS", "UKHLS", "csv"))) {
  parent <- dirname(root)
  if (identical(parent, root)) stop("找不到项目中的 UKHLS CSV 目录")
  root <- parent
}

threads <- strtoi(Sys.getenv("S03_THREADS", as.character(THREADS)), base = 10L)
sample_n <- strtoi(Sys.getenv("S03_SAMPLE_N", "0"), base = 10L)
if (is.na(threads) || threads < 1L || is.na(sample_n) || sample_n < 0L) {
  stop("S03_THREADS 必须为正整数，S03_SAMPLE_N 必须为非负整数")
}

input <- file.path(root, "DATA", "UKHLS", "UKHLS", "csv")
stata <- file.path(input, "UKDA-6614-stata", "stata", "stata14_se", "ukhls")
spss <- file.path(input, "UKDA-7251-spss", "spss", "spss28")
output <- file.path(root, "LJW_outcome", "data_clean", "S03_UKHLS")
if (sample_n > 0L) output <- file.path(output, paste0("smoke_test_", sample_n))
dir.create(file.path(output, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "quality_control"), showWarnings = FALSE)
log_path <- file.path(output, "logs", "02_Clean_nurse_S03_runtime.log")
log <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_path, append = TRUE)
}
started <- proc.time()[["elapsed"]]
log("开始；并行上限=", threads, "；每波次测试抽样=", sample_n, "；输入根目录=", input)
# 读 SPSS28 xindresp_ns。护士年月日为基线，有效均值用于血压判定。
# 近7天他汀按 BNF 指标和专项自报判定；BMI、腰围、降压药和设计字段保留。
# 无用药且非处方他汀阴性时，BNF 的不适用编码按已确认的阴性跳题处理。
path <- file.path(spss, "xindresp_ns.csv")
dest <- file.path(output, "02_nurse_baseline_clean_S03.csv")
audit <- file.path(output, "quality_control", "02_nurse_cleaning_audit_S03.csv")
log("读取护士数据：", path)
ids <- sample_ids(path, sample_n)
raw <- read_s03(path, nurse_columns())
if (!is.null(ids)) raw <- raw[raw$pidp %in% ids, , drop = FALSE]
log("清理 ", nrow(raw), " 条护士记录")
parts <- s03_parallel(split_rows(raw, threads), clean_nurse, threads, progress = function(i, n) log("并行任务进度：", i, "/", n))
write_s03(do.call(rbind, parts), dest)
write_s03(missing_audit(raw, basename(path)), audit)
log("完成护士记录清理；输出 ", dest, "；审计 ", audit)
log("总耗时 ", round(proc.time()[["elapsed"]] - started, 2), " 秒")
