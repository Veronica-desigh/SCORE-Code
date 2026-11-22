# S03 第01步 Clean_interview：单文件独立运行（功能函数与主程序合并）
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
# 完整的多选记录支持阴性；只有正向选择的字段不把“不适用”当成阴性。
condition_group <- function(d, stem, codes, full_list) {
  fields <- paste0(stem, codes)
  selected <- rowSums(d[fields] == 1, na.rm = TRUE) > 0
  none <- d[[paste0(stem, "96")]] %in% 1
  answered <- if (full_list) rowSums(d[fields] == 0, na.rm = TRUE) == length(fields) else
    rowSums(d[paste0(stem, c(1:17, 96))] == 1, na.rm = TRUE) > 0
  ifelse(selected, 1, ifelse(none | answered, 0, NA_real_))
}

clean_interview <- function(job) {
  wave <- job$wave
  prefix <- paste0(letters[wave], "_")
  date_cols <- c("intdaty_dv", "intdatm_dv", "intdatd_dv")
  full <- paste0("hcond", c(4:7, 14, 96))
  new <- paste0("hcondn", c(1:17, 96))
  core <- c("age_dv", "sex_dv", "gor_dv", "ethn_dv", "hiqual_dv", "jbstat")
  imd_fields <- c("imd2015qe_dv", "imd2017qni_dv", "imd2016qs_dv", "imd2014qw_dv")
  smoking <- if (wave == 2L) c("smnow", "smever") else c("evrsmo", "smofrq")
  # IMD 是按住址匹配的地区剥夺五分位。Wave 2 是唯一含成人中等强度运动频率的波次。
  fields <- c(date_cols, if (wave == 1L) "sex_dv", if (wave != 2L) c(full, "hcond2"),
              if (wave != 1L) c(new, core, smoking),
              if (wave == 2L) "sportsfreq",
              if (wave == 3L) c("hconde6", "hconde7"))
  raw <- read_s03(job$path, c("pidp", paste0(prefix, fields)))
  area_raw <- if (wave == 1L) NULL else
    read_s03(job$area_path, c("pidp", paste0(prefix, imd_fields)))
  if (!is.null(job$ids)) raw <- raw[raw$pidp %in% job$ids, , drop = FALSE]
  if (!is.null(area_raw) && !is.null(job$ids)) area_raw <- area_raw[area_raw$pidp %in% job$ids, , drop = FALSE]
  d <- raw
  names(d) <- sub(paste0("^", prefix), "", names(d))
  if (!is.null(area_raw)) {
    names(area_raw) <- sub(paste0("^", prefix), "", names(area_raw))
    if (anyDuplicated(area_raw$pidp[!is.na(area_raw$pidp)])) stop("indall 存在重复 pidp，Wave ", wave)
    row <- match(d$pidp, area_raw$pidp)
    for (name in imd_fields) d[[name]] <- area_raw[[name]][row]
  }
  na <- rep(NA_real_, nrow(d))
  smoke <- na
  if (wave == 2L) {
    smoke[d$smnow %in% 1] <- 1
    smoke[d$smnow %in% 2 | (d$smever %in% 2 & !d$smnow %in% 1)] <- 0
  }
  if (wave == 3L) {
    # 6614 Wave 3 问卷 Scyoungadults_w3 模块仅向 16–21 岁自填者提问。
    in_universe <- d$age_dv %in% 16:21
    smoke[in_universe & d$smofrq %in% 3:5] <- 1
    smoke[in_universe & (d$smofrq %in% 1:2 |
      (d$evrsmo %in% 2 & !d$smofrq %in% 3:5))] <- 0
  }
  get <- function(name, allowed = NULL) if (name %in% names(d)) clean_number(d[[name]], allowed = allowed) else na
  imd <- na
  imd_source <- rep(NA_character_, nrow(d))
  if (wave != 1L) {
    for (name in imd_fields) {
      value <- get(name, 1:5)
      take <- is.na(imd) & !is.na(value)
      imd[take] <- value[take]
      imd_source[take] <- name
    }
  }
  sport_frequency <- if (wave == 2L) get("sportsfreq", 1:7) else na
  # 代码 1–2 表示每周至少一次，代码 3–7 表示少于每周一次。
  physical_inactivity <- ifelse(sport_frequency %in% 1:2, 0,
                                ifelse(sport_frequency %in% 3:7, 1, NA_real_))
  asc_history <- if (wave != 2L) condition_group(d, "hcond", 4:7, TRUE) else na
  if (wave == 3L) asc_history[d$hconde6 %in% 1 | d$hconde7 %in% 1] <- 1
  result <- data.frame(
    pidp = clean_number(d$pidp, 1), wave = rep(wave, nrow(d)),
    interview_date = clean_date(d$intdaty_dv, d$intdatm_dv, d$intdatd_dv),
    age = get("age_dv"), sex = get("sex_dv", 1:2), current_smoking = smoke,
    region = get("gor_dv"), ethnicity = get("ethn_dv"), education = get("hiqual_dv"),
    employment = get("jbstat"),
    imd_quintile = imd, imd_source = imd_source,
    moderate_sport_frequency = sport_frequency, physical_inactivity_sport_proxy = physical_inactivity,
    ascvd_history = asc_history,
    diabetes_history = if (wave != 2L) condition_group(d, "hcond", 14, TRUE) else na,
    ascvd_new = if (wave != 1L) condition_group(d, "hcondn", 4:7, FALSE) else na,
    diabetes_new = if (wave != 1L) condition_group(d, "hcondn", 14, FALSE) else na,
    arthritis_history = if (wave != 2L) condition_group(d, "hcond", 2, TRUE) else na,
    arthritis_new = if (wave != 1L) condition_group(d, "hcondn", 2, FALSE) else na
  )
  audit <- missing_audit(raw, basename(job$path))
  if (!is.null(area_raw)) audit <- rbind(audit, missing_audit(area_raw, basename(job$area_path)))
  list(data = result, audit = audit)
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
log_path <- file.path(output, "logs", "01_Clean_interview_S03_runtime.log")
log <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_path, append = TRUE)
}
started <- proc.time()[["elapsed"]]
log("开始；并行上限=", threads, "；每波次测试抽样=", sample_n, "；输入根目录=", input)
# 读 a/b/c_indresp；b/c_indall 补充同波次地区 IMD 五分位。
# 日期按年月日构造；疾病多选按值标签解释，原始缺失原因另存。
# Wave 3 的 Evrsmo/Smofrq 仅适用于 16–21 岁自填问卷；目标年龄的吸烟状态保留缺失。
paths <- file.path(stata, paste0(c("a", "b", "c"), "_indresp.csv"))
area_paths <- file.path(stata, paste0(c("a", "b", "c"), "_indall.csv"))
dest <- file.path(output, "01_interview_clean_S03.csv")
audit <- file.path(output, "quality_control", "01_interview_cleaning_audit_S03.csv")
log("读取并行清理 Wave 1/2/3 主访：", paste(paths, collapse = " | "))
ids <- sample_ids(file.path(spss, "xindresp_ns.csv"), sample_n)
jobs <- Map(function(path, area_path, wave) list(path = path, area_path = area_path, wave = wave, ids = ids),
            paths, area_paths, 1:3)
parts <- s03_parallel(jobs, clean_interview, threads, progress = function(i, n) log("并行任务进度：", i, "/", n))
write_s03(do.call(rbind, lapply(parts, "[[", "data")), dest)
write_s03(do.call(rbind, lapply(parts, "[[", "audit")), audit)
log("完成 3/3 波次，清理记录 ", sum(vapply(parts, function(x) nrow(x$data), integer(1))), "；输出 ", dest, "；审计 ", audit)
log("总耗时 ", round(proc.time()[["elapsed"]] - started, 2), " 秒")
