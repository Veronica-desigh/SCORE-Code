# S03 第04步 Eligibility_criteria：单文件独立运行（功能函数与主程序合并）
# 默认5线程；可在下方 THREADS 参数或环境变量 S03_THREADS 调整。
# 测试参数沿用 S03_SAMPLE_N。
THREADS <- 5L  # 修改此处设置默认并行线程数


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

# 关联键完整的源表必须唯一，防止笛卡尔展开和任意去重。
key_s03 <- function(d) {
  ifelse(is.na(d$pidp) | is.na(d$wave), NA_character_, paste(d$pidp, d$wave, sep = ":"))
}

attach_s03 <- function(base, source, label) {
  key <- key_s03(source)
  if (anyDuplicated(key[!is.na(key)])) stop(label, " 存在重复 pidp + wave")
  match_row <- match(key_s03(base), key[!is.na(key)])
  matched <- which(!is.na(key))[match_row]
  data <- source[matched, setdiff(names(source), c("pidp", "wave")), drop = FALSE]
  rownames(data) <- NULL
  list(data = data, audit = data.frame(source = label, source_rows = nrow(source),
       invalid_key_rows = sum(is.na(key)), base_rows = nrow(base),
       matched_rows = sum(!is.na(matched)), unmatched_rows = sum(is.na(matched))))
}

# 既往阳性可延续；阴性须有可定时的病史问项覆盖。
baseline_history <- function(d, condition) {
  history <- paste0(condition, "_history")
  recent <- paste0(condition, "_new")
  before <- !is.na(d$baseline_date) & !is.na(d$interview_date) & d$interview_date <= d$baseline_date
  positive <- before & (d[[history]] %in% 1 | d[[recent]] %in% 1)
  negative <- before & d[[history]] %in% 0
  for (wave in 1:2) {
    date <- d[[paste0("history_date_", wave)]]
    earlier <- d$wave %in% wave:3 & !is.na(date) &
      !is.na(d$baseline_date) & date <= d$baseline_date
    positive <- positive | (earlier & d[[paste0(history, "_", wave)]] %in% 1)
    if (wave == 2L) positive <- positive | (earlier & d[[paste0(recent, "_2")]] %in% 1)
  }
  # b 波次只有新增病史时，用 a 波次完整病史补齐；c 可沿用完整基线史加两次更新。
  a_date <- d$history_date_1
  b_date <- d$history_date_2
  a_neg <- !is.na(a_date) & !is.na(d$interview_date) & a_date <= d$interview_date &
    d[[paste0(history, "_1")]] %in% 0
  b_neg <- !is.na(b_date) & !is.na(a_date) & !is.na(d$interview_date) &
    b_date >= a_date & b_date <= d$interview_date & d[[paste0(recent, "_2")]] %in% 0
  negative <- negative | (before & a_neg & d[[recent]] %in% 0 &
                             (d$wave %in% 2 | (d$wave %in% 3 & b_neg)))
  ifelse(positive, 1, ifelse(negative, 0, NA_real_))
}

derive_eligibility <- function(d) {
  # 2021 CKD-EPI，无种族系数；性别编码 1=男，2=女。
  female <- d$sex == 2
  scr <- d$creatinine_umol_l / 88.4
  ratio <- scr / ifelse(female, 0.7, 0.9)
  d$egfr <- 142 * pmin(ratio, 1)^ifelse(female, -0.241, -0.302) *
    pmax(ratio, 1)^(-1.200) * 0.9938^d$age * ifelse(female, 1.012, 1)
  ldl <- d$tc - d$hdl - d$tg / 2.2
  d$ldl_reason <- ifelse(is.na(ldl), "missing_lipid_input",
                        ifelse(d$tg >= 4.5, "tg_outside_range",
                               ifelse(ldl <= 0, "nonpositive_estimate", "estimated")))
  d$ldl <- ifelse(d$ldl_reason == "estimated", ldl, NA_real_)
  d$ascvd <- baseline_history(d, "ascvd") == 1
  d$self_report_diabetes <- baseline_history(d, "diabetes") == 1
  d$arthritis <- baseline_history(d, "arthritis")
  # 糖尿病：阳性优先；有任一阴性且无阳性时判阴性；两项都缺才未知。
  hba1c_positive <- d$hba1c >= 48
  positive <- d$self_report_diabetes %in% TRUE | hba1c_positive %in% TRUE
  negative <- d$self_report_diabetes %in% FALSE | hba1c_positive %in% FALSE
  d$diabetes <- ifelse(positive, TRUE, ifelse(negative, FALSE, NA))
  d$low_egfr <- d$egfr < 60
  d$recent_statin <- d$statin_7d == 1
  d$high_tc <- d$tc > 8
  d$high_ldl <- d$ldl > 4.9
  d$high_bp <- d$sbp >= 180 | d$dbp >= 110
  d$extreme_risk <- d$high_tc | d$high_ldl | d$high_bp
  d$rm_obesity_central <- ifelse(d$bmi >= 30 | (d$sex == 1 & d$waist_cm >= 102) |
                                   (d$sex == 2 & d$waist_cm >= 88), 1,
                                  ifelse(!is.na(d$bmi) & !is.na(d$waist_cm) & d$bmi < 30 &
                                           !((d$sex == 1 & d$waist_cm >= 102) |
                                             (d$sex == 2 & d$waist_cm >= 88)), 0, NA_real_))
  d$rm_south_asian <- ifelse(d$ethnicity %in% 9:11, 1,
                             ifelse(d$ethnicity %in% 1:17, 0, NA_real_))
  d$rm_hscrp_gt2 <- ifelse(!is.na(d$hscrp), as.numeric(d$hscrp > 2), NA_real_)
  d$rm_social_deprivation <- ifelse(d$imd_quintile %in% 1, 1,
                                    ifelse(d$imd_quintile %in% 2:5, 0, NA_real_))
  d$rm_physical_inactivity_sport_proxy <- d$physical_inactivity_sport_proxy
  d$rm_chronic_inflammatory_arthritis_proxy <- d$arthritis
  d
}

summarize_missing <- function(d, fields) {
  missing <- vapply(d[fields], function(x) sum(is.na(x)), integer(1))
  data.frame(variable = fields, n_total = nrow(d), n_observed = nrow(d) - missing,
             n_missing = missing, missing_fraction = if (nrow(d)) missing / nrow(d) else NA_real_,
             row.names = NULL)
}

# 这些 RM 仅作描述，不参与本脚本的纳入或排除。
summarize_risk_modifiers <- function(d) {
  metadata <- data.frame(
    variable = c("rm_obesity_central", "rm_south_asian", "rm_hscrp_gt2",
                 "rm_social_deprivation", "rm_physical_inactivity_sport_proxy",
                 "rm_chronic_inflammatory_arthritis_proxy"),
    risk_modifier = c("肥胖或中心性肥胖", "高风险南亚族裔", "hs-CRP 升高",
                      "社会剥夺", "体力活动不足", "慢性炎症疾病（关节炎代理）"),
    guideline = c("2019/2025", "2025", "2025", "2019/2025", "2019/2025", "2019/2025"),
    construction_status = c("已构建", "已构建", "已构建", "已构建", "已构建，Wave 2 限定",
                            "已构建，关节炎代理"),
    definition = c(
      "BMI≥30 kg/m² 或腰围男性≥102 cm、女性≥88 cm。任一阳性即阳性。",
      "同波次主访 ethnicity 代码 9、10、11，即印度、巴基斯坦、孟加拉。",
      "同波次血液 hs-CRP >2 mg/L。",
      "同波次住址匹配的国家 IMD 五分位为 1，即该国最贫困的 20%。英国各国指数不与 UKB Townsend 指数混用。",
      "Wave 2 b_sportsfreq 的中等强度运动少于每周一次，代码 3–7。每周至少一次，代码 1–2，判为阴性。该问项不记录总活动时长，不能等同 UKB Field 22035。",
      "基线前或基线时自报关节炎。UKHLS 未区分类风湿性关节炎、骨关节炎及其他关节炎，也未覆盖银屑病、炎症性肠病等，因此仅作关节炎代理。"
    ),
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(seq_len(nrow(metadata)), function(i) {
    value <- d[[metadata$variable[i]]]
    observed <- !is.na(value)
    n_observed <- sum(observed)
    n_positive <- sum(value %in% 1, na.rm = TRUE)
    n_negative <- sum(value %in% 0, na.rm = TRUE)
    cbind(metadata[i, ], applicable_n = nrow(d), positive_n = n_positive,
          negative_n = n_negative, observed_n = n_observed,
          missing_n = nrow(d) - n_observed,
          positive_percent_observed = if (n_observed) 100 * n_positive / n_observed else NA_real_,
          missing_percent_applicable = if (nrow(d)) 100 * (nrow(d) - n_observed) / nrow(d) else NA_real_)
  }))
}

run_eligibility <- function(inputs, paths, threads, log,
                            use_wave2_smoking = FALSE) {
  tables <- lapply(inputs, function(p) read.csv(p, na.strings = "", check.names = FALSE))
  interview <- tables[[1]]
  nurse <- tables[[2]]
  blood <- tables[[3]]
  valid_pid <- nurse$pidp[!is.na(nurse$pidp)]
  if (anyDuplicated(valid_pid)) stop("护士文件中存在重复 pidp，请先审查基线记录")
  nurse$baseline_date <- as.Date(nurse$baseline_date, format = "%Y-%m-%d")
  interview$interview_date <- as.Date(interview$interview_date, format = "%Y-%m-%d")
  core <- attach_s03(nurse, interview[interview$wave %in% 2:3, ], "same_wave_interview")
  labs <- attach_s03(nurse, blood, "same_wave_blood")
  d <- cbind(nurse, core$data, labs$data)
  for (wave in 1:2) {
    prior <- interview[interview$wave %in% wave, ]
    if (anyDuplicated(prior$pidp[!is.na(prior$pidp)])) stop("历史主访 pidp 重复，Wave ", wave)
    match_row <- match(d$pidp, prior$pidp)
    match_row[is.na(d$pidp)] <- NA_integer_
    d[[paste0("history_date_", wave)]] <- prior$interview_date[match_row]
    for (field in c("ascvd_history", "diabetes_history", "arthritis_history",
                    "ascvd_new", "diabetes_new", "arthritis_new")) {
      d[[paste0(field, "_", wave)]] <- prior[[field]][match_row]
    }
    d[[paste0("sex_", wave)]] <- prior$sex[match_row]
    if (wave == 2L) d$smoking_wave2 <- prior$current_smoking[match_row]
  }
  # 性别缺失时，仅接受基线前已有的有效记录；历史记录冲突时保留未知。
  sex1 <- ifelse(d$history_date_1 <= d$baseline_date, d$sex_1, NA)
  sex2 <- ifelse(d$history_date_2 <= d$baseline_date, d$sex_2, NA)
  conflict <- !is.na(sex1) & !is.na(sex2) & sex1 != sex2
  use_sex2 <- is.na(d$sex) & !conflict & sex2 %in% 1:2
  use_sex1 <- is.na(d$sex) & !conflict & !use_sex2 & sex1 %in% 1:2
  d$sex_source_wave <- ifelse(d$sex %in% 1:2, d$wave, NA_integer_)
  d$sex[use_sex2] <- sex2[use_sex2]
  d$sex[use_sex1] <- sex1[use_sex1]
  d$sex_source_wave[use_sex2] <- 2L
  d$sex_source_wave[use_sex1] <- 1L

  # 吸烟会随时间变化。主分析用同波次基线前记录，敏感性分析允许 Wave 2 前值。
  observed <- d$current_smoking %in% 0:1 & !is.na(d$interview_date) &
    !is.na(d$baseline_date) & d$interview_date <= d$baseline_date
  d$current_smoking[!observed] <- NA
  d$smoking_source_wave <- ifelse(observed, d$wave, NA_integer_)
  d$smoking_source_date <- d$interview_date
  d$smoking_source_date[!observed] <- as.Date(NA)
  if (use_wave2_smoking) {
    use_prior <- d$wave %in% 3 & !observed & d$smoking_wave2 %in% 0:1 &
      !is.na(d$history_date_2) & !is.na(d$baseline_date) & d$history_date_2 <= d$baseline_date
    d$current_smoking[use_prior] <- d$smoking_wave2[use_prior]
    d$smoking_source_wave[use_prior] <- 2L
    d$smoking_source_date[use_prior] <- d$history_date_2[use_prior]
    log("Wave 3 使用基线前 Wave 2 吸烟代理：", sum(use_prior), " 人")
  }
  d$smoking_days_before_baseline <- as.integer(d$baseline_date - d$smoking_source_date)
  log("已关联 ", nrow(d), " 条护士记录，开始并行计算纳排指标")
  d <- do.call(rbind, s03_parallel(split_rows(d, threads), derive_eligibility, threads, progress = function(i, n) log("纳排变量计算：", i, "/", n)))
  required <- c("pidp", "baseline_date", "wave", "age", "sex", "current_smoking",
                "tc", "hdl", "ldl", "sbp", "dbp", "creatinine_umol_l", "egfr",
                "ascvd", "diabetes", "low_egfr", "recent_statin", "extreme_risk")
  d$missing_required_fields <- vapply(seq_len(nrow(d)), function(i) {
    paste(required[is.na(d[i, required, drop = FALSE])], collapse = "|")
  }, character(1))
  rule <- list(
    invalid_id_or_baseline = is.na(d$pidp) | is.na(d$baseline_date) |
      d$baseline_date < as.Date("2010-01-01") | d$baseline_date > as.Date("2012-12-31"),
    age_outside_40_69 = d$age < 40 | d$age > 69,
    prior_ascvd = d$ascvd, diabetes = d$diabetes, low_egfr = d$low_egfr,
    recent_statin = d$recent_statin, extreme_risk = d$extreme_risk,
    required_information_missing = nzchar(d$missing_required_fields)
  )
  for (name in names(rule)) d[[paste0("rule_", name)]] <- rule[[name]]
  d$first_exclusion_reason <- rep(NA_character_, nrow(d))
  remaining <- rep(TRUE, nrow(d))
  flow <- vector("list", length(rule))
  for (i in seq_along(rule)) {
    name <- names(rule)[i]
    if (name == "required_information_missing") {
      required_missing <- summarize_missing(d[remaining, , drop = FALSE], required)
    }
    exclude <- remaining & rule[[i]] %in% TRUE
    before <- sum(remaining)
    d$first_exclusion_reason[exclude] <- name
    remaining[exclude] <- FALSE
    flow[[i]] <- data.frame(step = i, rule = name, n_before = before,
                           n_excluded = sum(exclude), n_remaining = sum(remaining))
    log("纳排 ", i, "/", length(rule), " ", name, "：排除 ", sum(exclude), "，剩余 ", sum(remaining))
  }
  d$eligible <- remaining
  cohort <- d[remaining, , drop = FALSE]
  flow <- do.call(rbind, flow)
  stopifnot(!anyNA(cohort[required]), !anyDuplicated(cohort$pidp),
            all(flow$n_before - flow$n_excluded == flow$n_remaining),
            sum(flow$n_excluded) + nrow(cohort) == nrow(d))
  wave_summary <- do.call(rbind, lapply(c(2L, 3L, NA_integer_), function(wave) {
    records <- d[if (is.na(wave)) is.na(d$wave) else d$wave %in% wave, , drop = FALSE]
    data.frame(wave = wave, n_nurse = nrow(records), n_eligible = sum(records$eligible),
               n_missing_smoking = sum(is.na(records$current_smoking)),
               n_required_info_excluded = sum(records$first_exclusion_reason %in% "required_information_missing"))
  }))
  coverage <- data.frame(
    rule = c("ASCVD", "diabetes", "CKD", "lipid_therapy", "LDL-C", "current_smoking"),
    operationalization = c("pre-baseline self-reported CHD, angina, MI or stroke",
                           "positive self-report or HbA1c >=48 overrides negative; one negative suffices",
                           "CKD-EPI 2021 eGFR <60",
                           "statin use in previous seven days", "Friedewald; TG <4.5 mmol/L",
                           if (use_wave2_smoking) "Wave 3: pre-baseline Wave 2 smoking proxy" else
                             "same-wave interview on or before baseline"),
    limitation = c("TIA, PAD and revascularization not ascertainable",
                   "glucose-lowering treatment not separately ascertainable",
                   "historical N18 not ascertainable", "365-day and non-statin coverage unavailable",
                   "estimated from nonfasting lipid panel",
                   "Wave 3 EVRSMO/SMOFRQ are restricted to ages 16-21")
  )
  results <- list(cohort = cohort, flow = flow, participant = d, required = required_missing,
                  merge = rbind(core$audit, labs$audit), missing = summarize_missing(cohort, names(cohort)),
                  coverage = coverage, wave = wave_summary,
                  risk_modifiers = summarize_risk_modifiers(cohort))
  for (name in names(paths)) {
    write_s03(results[[name]], paths[[name]])
    log("写入 ", paths[[name]])
  }
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
clean_output <- output
smoking_mode <- Sys.getenv("S03_SMOKING_MODE", "same_wave")
if (!smoking_mode %in% c("same_wave", "wave2_prior")) {
  stop("S03_SMOKING_MODE 必须为 same_wave 或 wave2_prior")
}
# 敏感性分析复用 01–03 清理产物，另存其 04 队列和审计表。
if (smoking_mode == "wave2_prior") output <- file.path(output, "wave2_smoking_sensitivity")
dir.create(file.path(output, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "quality_control"), showWarnings = FALSE)
log_path <- file.path(output, "logs", "04_Eligibility_criteria_S03_runtime.log")
log <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
  message(line)
  cat(line, "\n", file = log_path, append = TRUE)
}
started <- proc.time()[["elapsed"]]
log("开始；并行上限=", threads, "；每波次测试抽样=", sample_n, "；输入根目录=", input)
# 以护士为起始样本，按 pidp + wave 关联主访和血液。
# 年龄 40–69；排除基线前自报 ASCVD、糖尿病（含 HbA1c ≥48）、
# eGFR <60、近7天他汀、TC >8、估算 LDL-C >4.9、SBP ≥180 或 DBP ≥110。
# Friedewald LDL-C 用 mmol/L，TG ≥4.5 或 LDL-C ≤0 时设为不可估算。
# 最后仅筛必需字段及可实施的排除状态，其他协变量缺失原样保留。
# 糖尿病按自报和 HbA1c 综合判断；任一阳性优先，任一阴性且无阳性即可判阴性。
# 最后一步要求糖尿病状态可判定，HbA1c 本身允许缺失。
# TIA、PAD、血运重建、降糖药、历史 N18 无可用同定义字段，记录为范围限制。
inputs <- file.path(clean_output, c("01_interview_clean_S03.csv", "02_nurse_baseline_clean_S03.csv",
                             "03_blood_biomarkers_clean_S03.csv"))
flow_dir <- file.path(output, "人群纳排记录_S03")
qc_dir <- file.path(output, "quality_control")
dir.create(flow_dir, showWarnings = FALSE)
paths <- c(
  cohort = file.path(output, "04_eligible_cohort_S03.csv"),
  flow = file.path(flow_dir, "04_eligibility_flow_S03.csv"),
  participant = file.path(flow_dir, "04_eligibility_participant_status_S03.csv"),
  required = file.path(flow_dir, "04_required_information_missingness_S03.csv"),
  wave = file.path(flow_dir, "04_eligibility_by_wave_S03.csv"),
  merge = file.path(qc_dir, "04_merge_audit_S03.csv"),
  missing = file.path(qc_dir, "04_eligible_cohort_missingness_S03.csv"),
  coverage = file.path(qc_dir, "04_eligibility_definition_coverage_S03.csv"),
  risk_modifiers = file.path(qc_dir, "04_RM_distribution_missingness_S03.csv")
)
log("读取清理产物：", paste(inputs, collapse = " | "))
run_eligibility(inputs, paths, threads, log,
                use_wave2_smoking = smoking_mode == "wave2_prior")
log("总耗时 ", round(proc.time()[["elapsed"]] - started, 2), " 秒")
