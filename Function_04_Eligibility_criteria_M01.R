# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_04_Eligibility_criteria_M01.R
# 任务：供04调用的纳排规则、流程审计
# 输入：上游01–03合并记录
# 输出：eligibility标记与审计表
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function_04_Eligibility_criteria_M01_v1.2.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   为 04_Eligibility_criteria_M01.R 提供纳入排除、疾病状态和用药状态判定函数。
#
# 设计原则：
#   1. 主脚本只保留路径、输入输出和纳排流程；
#   2. documented ASCVD、糖尿病、CKD 和显著升高单一危险因素的判定集中在本文件；
#   3. 降脂药仍按基线自报和基线前 365 天 GP 处方共同判定；
#   4. eligibility / guideline-strategy 所需信息必须真实观测；普通 outcome-model
#      协变量不再在 04 做 complete-case 排除，而是保留缺失供 06 多重插补。
#   5. 04 同时输出随后插补变量在 required-information gate 前后的人数和缺失率，
#      但这些缺失不构成纳排流程步骤。
#   6. v1.2.1 修复 required-information missingness 中标量 test 与向量 yes
#      误用 fifelse() 导致的长度不一致错误。
# ============================================================================== #

library(data.table)

# ============================================================================== #
# 0. 通用检查、计时和日志函数 ####
# ============================================================================== #
fmt_time <- function(sec) {
  if (is.null(sec) || is.na(sec)) {
    return(NA_character_)
  }
  if (sec < 60) {
    return(sprintf("%.2f sec", sec))
  }
  if (sec < 3600) {
    return(sprintf("%.2f min", sec / 60))
  }
  sprintf("%.2f hr", sec / 3600)
}

step_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

run_with_timer <- function(expr, label = NULL) {
  t0 <- Sys.time()
  if (!is.null(label)) step_message("开始：", label)
  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(label)) step_message("完成：", label, " | 耗时 = ", fmt_time(elapsed))
  if (exists("m01_record_timing", mode = "function", inherits = TRUE))
    m01_record_timing(if (is.null(label)) "unnamed_run_with_timer" else label, elapsed)
  list(value = value, elapsed_sec = elapsed)
}

# ============================================================================== #
# 1. 工具函数 ####
# ============================================================================== #
# as_binary_self_report()
#
# 说明：
#   把自报服药变量统一整理成 0 / 1 / NA，
#   以便后续直接进行向量化判定。
#
# 示例：
#   - 输入 1、"1"、1L -> 返回 1L
#   - 输入 0、"0"、0L -> 返回 0L
#   - 输入 ""、"Prefer not to answer"、NA -> 返回 NA
# ------------------------------------------------------------------------------ #
as_binary_self_report <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN", "Prefer not to answer", "Do not know")] <- NA_character_
  x_num <- suppressWarnings(as.numeric(x_chr))
  fifelse(is.na(x_num), NA_integer_, fifelse(x_num == 1, 1L, fifelse(x_num == 0, 0L, NA_integer_)))
}


# ------------------------------------------------------------------------------ #
# assert_unique_id() / merge_one_to_one_audited()
#
# 说明：
#   主表和辅助表必须预先满足一人一行；若出现重复 ID，直接停止，
#   不允许通过“保留第一条”静默掩盖数据冲突。连接后同时返回合并审计。
# ------------------------------------------------------------------------------ #
assert_unique_id <- function(dt, label, by = "Participant ID") {
  dt <- as.data.table(dt)
  require_columns(dt, by, label)
  if (anyNA(dt[[by]])) {
    stop(label, " 存在缺失 ID。", call. = FALSE)
  }
  duplicated_ids <- dt[, .N, by = by][N > 1L]
  if (nrow(duplicated_ids) > 0L) {
    stop(
      label, " 存在重复 ID（", nrow(duplicated_ids), " 个），不得静默取第一条。",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

merge_one_to_one_audited <- function(main_dt, add_dt, label, by = "Participant ID") {
  main_dt <- copy(main_dt)
  add_dt <- copy(add_dt)
  setDT(main_dt)
  setDT(add_dt)
  # CSV inference can differ between upstream outputs; IDs use character keys.
  if (by %in% names(main_dt)) set(main_dt,j=by,value=as.character(main_dt[[by]]))
  if (by %in% names(add_dt)) set(add_dt,j=by,value=as.character(add_dt[[by]]))
  assert_unique_id(main_dt, paste0(label, "：主表"), by)
  assert_unique_id(add_dt, paste0(label, "：辅助表"), by)

  matched_index <- match(main_dt[[by]], add_dt[[by]])
  overlapping_cols <- setdiff(intersect(names(main_dt), names(add_dt)), by)
  for (v in overlapping_cols) {
    main_value <- as.character(main_dt[[v]])
    add_value <- as.character(add_dt[[v]][matched_index])
    comparable <- !is.na(matched_index)
    same <- (is.na(main_value) & is.na(add_value)) | (!is.na(main_value) & !is.na(add_value) & main_value == add_value)
    if (any(comparable & !same)) {
      stop(label, " 的重名字段 ", v, " 在两个来源中不一致。", call. = FALSE)
    }
  }
  if (length(overlapping_cols) > 0L) {
    add_dt[, (overlapping_cols) := NULL]
  }

  main_dt[, `__row_order__` := .I]
  out <- merge(main_dt, add_dt, by = by, all.x = TRUE, sort = FALSE)
  setorder(out, `__row_order__`)
  out[, `__row_order__` := NULL]
  assert_unique_id(out, paste0(label, "：连接结果"), by)
  if (nrow(out) != nrow(main_dt)) {
    stop(label, " 连接后行数发生变化。", call. = FALSE)
  }

  audit <- data.table(
    merge_step = label,
    n_main = nrow(main_dt),
    n_add = nrow(add_dt),
    n_matched = sum(!is.na(matched_index)),
    n_unmatched_main = sum(is.na(matched_index)),
    n_add_not_in_main = sum(!add_dt[[by]] %in% main_dt[[by]]),
    n_output = nrow(out),
    overlapping_columns_verified = paste(overlapping_cols, collapse = ";")
  )
  list(dt = out, audit = audit)
}


# ------------------------------------------------------------------------------ #
# init_flow_log()
#
# 说明：
#   初始化纳排流程记录表，用于系统保存每一步的名称、排除人数、剩余人数。
#
# 示例：
#   返回一个空 data.table，后续每一步都可以往里追加：
#   step_no / step_name / n_before / n_excluded / n_remaining。
# ------------------------------------------------------------------------------ #
init_flow_log <- function() {
  data.table(
    step_no = integer(),
    step_name = character(),
    n_before = integer(),
    n_excluded = integer(),
    n_remaining = integer(),
    exclude_pct = numeric(),
    remaining_pct_from_start = numeric()
  )
}


# ------------------------------------------------------------------------------ #
# add_flow_step()
#
# 说明：
#   在完成某一步纳排后，把该步的名称、排除人数、剩余人数追加到 flow_log。
#   这是本脚本中“最重要的留痕表”。
#
# 示例：
#   若某一步前有 50000 人，筛后剩 42000 人，
#   则自动记录：n_before = 50000, n_excluded = 8000, n_remaining = 42000。
# ------------------------------------------------------------------------------ #
add_flow_step <- function(flow_log, step_name, n_before, n_after, n_start) {
  setDT(flow_log)
  new_row <- data.table(
    step_no = nrow(flow_log) + 1L,
    step_name = step_name,
    n_before = as.integer(n_before),
    n_excluded = as.integer(n_before - n_after),
    n_remaining = as.integer(n_after),
    exclude_pct = ifelse(n_before > 0, round(100 * (n_before - n_after) / n_before, 2), NA_real_),
    remaining_pct_from_start = ifelse(n_start > 0, round(100 * n_after / n_start, 2), NA_real_)
  )
  rbind(flow_log, new_row, fill = TRUE)
}


# ------------------------------------------------------------------------------ #
# apply_filter_with_log()
#
# 说明：
#   对一个数据集执行一次筛选，并同步记录这一步的纳排名称、排除人数和剩余人数。
#   这样就不会漏记任何一步。
#
# 示例：
#   调用一次后，会同时得到：
#   1. 筛选后的 dt
#   2. 更新后的 flow_log
#   3. 控制台打印“该步排除了多少，还剩多少”
# ------------------------------------------------------------------------------ #
apply_filter_with_log <- function(dt, i, flow_log, step_name, n_start) {
  n_before <- nrow(dt)

  idx <- eval(substitute(i), envir = dt, enclos = parent.frame())
  dt_after <- dt[idx]

  n_after <- nrow(dt_after)

  flow_log <- add_flow_step(
    flow_log = flow_log,
    step_name = step_name,
    n_before = n_before,
    n_after = n_after,
    n_start = n_start
  )

  step_message(
    "纳排步骤：", step_name,
    " | 排除 = ", format(n_before - n_after, big.mark = ","),
    " | 剩余 = ", format(n_after, big.mark = ",")
  )

  list(dt = dt_after, flow_log = flow_log)
}


# ------------------------------------------------------------------------------ #
# require_columns()
#
# 说明：
#   对关键纳排字段做显式检查。关键字段缺失时直接报错，避免静默不排除。
# ------------------------------------------------------------------------------ #
require_columns <- function(dt, cols, context = "数据表") {
  missing_cols <- setdiff(cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      context, " 缺少必要字段：",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ------------------------------------------------------------------------------ #
# flag_is_one()
#
# 说明：
#   把 1/0、TRUE/FALSE、字符型 1/0 统一转成是否为 1 的逻辑向量。
# ------------------------------------------------------------------------------ #
flag_is_one <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  x_num <- suppressWarnings(as.integer(as.character(x)))
  !is.na(x_num) & x_num == 1L
}


# ------------------------------------------------------------------------------ #
# clean_character_missing_values()
#
# 说明：
#   把字符型协变量中的非信息性回答统一转成 NA。
#   这样 complete.cases() 才能把 "Prefer not to answer" 等回答视为缺失。
# ------------------------------------------------------------------------------ #
clean_character_missing_values <- function(x,
                                           missing_strings = c(
                                             "",
                                             "NA",
                                             "N/A",
                                             "NaN",
                                             "NULL",
                                             "Do not know",
                                             "Prefer not to answer"
                                           )) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% missing_strings] <- NA_character_
  x
}


# ------------------------------------------------------------------------------ #
# read_baseline_diabetes_self_report()
#
# 说明：
#   从 UKB Field 20002 的 Instance 0 数组识别基线自报糖尿病。
#   仅使用 1220（diabetes）、1222（type 1 diabetes）和 1223（type 2 diabetes），
#   后续 Instance 不用于回填基线状态。
# ------------------------------------------------------------------------------ #
read_baseline_diabetes_self_report <- function(
    file,
    column_prefix = "Non-cancer illness code, self-reported | Instance 0") {
  if (!file.exists(file)) stop("Field 20002 自报文件不存在：", file, call. = FALSE)
  header <- fread(file, nrows = 0L)
  value_cols <- names(header)[startsWith(names(header), column_prefix)]
  require_columns(header, "Participant ID", "Field 20002 自报文件")
  if (length(value_cols) == 0L) {
    stop("Field 20002 文件中未找到 Instance 0 数组。", call. = FALSE)
  }
  raw <- fread(file, select = c("Participant ID", value_cols))
  assert_unique_id(raw, "Field 20002 基线糖尿病自报表")
  accepted <- c("1220", "1222", "1223", "diabetes", "type 1 diabetes", "type 2 diabetes")
  normalize <- function(x) tolower(trimws(as.character(x)))
  flag <- Reduce(`|`, lapply(raw[, ..value_cols], function(x) normalize(x) %in% accepted))
  out <- raw[, .(`Participant ID`)]
  out[, diabetes_self_report_baseline := as.integer(flag)]
  assert_unique_id(out, "Field 20002 基线糖尿病自报派生表")
  out[]
}


# ------------------------------------------------------------------------------ #
# prevalent_from_time_or_flag()
#
# 说明：
#   根据“基线前/基线当日事件时间”或已经整理好的 prevalent 标记判断基线既往史。
# ------------------------------------------------------------------------------ #
prevalent_from_time_or_flag <- function(dt,
                                        time_col = NULL,
                                        flag_col = NULL,
                                        baseline_col = "InCohort_time") {
  setDT(dt)
  require_columns(dt, baseline_col, "prevalence 判定")

  baseline_date <- as.IDate(dt[[baseline_col]])
  out <- rep(FALSE, nrow(dt))

  if (!is.null(flag_col) && flag_col %in% names(dt)) {
    out <- out | flag_is_one(dt[[flag_col]])
  }

  if (!is.null(time_col) && time_col %in% names(dt)) {
    event_time <- as.IDate(dt[[time_col]])
    out <- out | (!is.na(event_time) & !is.na(baseline_date) & event_time <= baseline_date)
  }

  out[is.na(out)] <- FALSE
  out
}


# ------------------------------------------------------------------------------ #
# calculate_egfr_2021()
#
# 说明：
#   使用 CKD-EPI 2021 race-free 肌酐公式计算基线 eGFR。
#   UKB 血肌酐字段单位为 umol/L，因此先除以 88.4 转成 mg/dL。
# ------------------------------------------------------------------------------ #
calculate_egfr_2021 <- function(creatinine_umol_l, age, sex) {
  creatinine_umol_l <- suppressWarnings(as.numeric(creatinine_umol_l))
  age <- suppressWarnings(as.numeric(age))
  sex_clean <- tolower(trimws(as.character(sex)))
  sex_clean[sex_clean %in% c("m", "male")] <- "male"
  sex_clean[sex_clean %in% c("f", "female")] <- "female"

  scr <- creatinine_umol_l / 88.4
  valid <- !is.na(scr) & scr > 0 & !is.na(age) & age >= 0 & sex_clean %in% c("male", "female")
  out <- rep(NA_real_, length(scr))

  female <- valid & sex_clean == "female"
  if (any(female)) {
    k <- 0.7
    alpha <- -0.241
    out[female] <- 142 *
      (pmin(scr[female] / k, 1)^alpha) *
      (pmax(scr[female] / k, 1)^(-1.200)) *
      (0.9938^age[female]) *
      1.012
  }

  male <- valid & sex_clean == "male"
  if (any(male)) {
    k <- 0.9
    alpha <- -0.302
    out[male] <- 142 *
      (pmin(scr[male] / k, 1)^alpha) *
      (pmax(scr[male] / k, 1)^(-1.200)) *
      (0.9938^age[male])
  }

  out
}


# ------------------------------------------------------------------------------ #
# derive_documented_ascvd_baseline_status()
#
# 说明：
#   使用 02 脚本中专门为主分析一级预防排除生成的 inclusion_ascvd。
#   该变量由 IHD、stroke、TIA、PAD/AAA、冠脉血运重建等 documented ASCVD
#   组分合并得到；死亡结局不用于反向定义基线 ASCVD。
# ------------------------------------------------------------------------------ #
derive_documented_ascvd_baseline_status <- function(dt,
                                                    baseline_col = "InCohort_time") {
  dt <- copy(dt)
  setDT(dt)
  require_columns(dt, baseline_col, "documented ASCVD 判定")
  if (!any(c("prevalent_inclusion_ascvd", "inclusion_ascvd_time") %in% names(dt))) {
    stop("缺少 documented ASCVD 纳排字段：prevalent_inclusion_ascvd / inclusion_ascvd_time", call. = FALSE)
  }

  dt[, documented_ascvd_baseline := as.integer(prevalent_from_time_or_flag(
    dt,
    time_col = "inclusion_ascvd_time",
    flag_col = "prevalent_inclusion_ascvd",
    baseline_col = baseline_col
  ))]
  dt[, documented_ascvd_time := if ("inclusion_ascvd_time" %in% names(dt)) as.IDate(inclusion_ascvd_time) else as.IDate(NA)]
  dt[]
}


# ------------------------------------------------------------------------------ #
# derive_diabetes_baseline_status()
#
# 说明：
#   合并四类基线糖尿病证据，并同时标记 HbA1c 是否具有真实观测值：
#   1) 02 脚本生成的 E10-E14 ICD first occurrence；
#   2) Field 20002 Instance 0 基线自报糖尿病；
#   3) 基线 HbA1c 按 UKB 原始单位 mmol/mol 固定使用 >=48 mmol/mol；
#   4) 03 脚本整理的基线降糖药/胰岛素证据。
#
# 主分析纳排规则：
#   - ICD、HbA1c 或降糖药任一提示糖尿病时，排除；
#   - HbA1c 缺失不在本函数中排除，统一交给最后 required-information gate。
# ------------------------------------------------------------------------------ #
derive_diabetes_baseline_status <- function(dt,
                                            baseline_col = "InCohort_time",
                                            hba1c_col = "Glycated_haemoglobin_(HbA1c)_0",
                                            glucose_med_col = "glucose_lowering_med_baseline",
                                            self_report_col = "diabetes_self_report_baseline") {
  dt <- copy(dt)
  setDT(dt)
  require_columns(dt, c(baseline_col, hba1c_col, glucose_med_col, self_report_col), "糖尿病判定")

  diabetes_time_cols <- intersect(
    c("Diabetes_ICD_time", "T1D_time", "T2D_time"),
    names(dt)
  )
  if (length(diabetes_time_cols) == 0L) {
    stop("缺少糖尿病 ICD 时间字段：Diabetes_ICD_time / T1D_time / T2D_time", call. = FALSE)
  }

  dt[, diabetes_icd_baseline := 0L]
  for (v in diabetes_time_cols) {
    event_time <- as.IDate(dt[[v]])
    baseline_date <- as.IDate(dt[[baseline_col]])
    dt[, diabetes_icd_baseline := as.integer(
      diabetes_icd_baseline == 1L |
        (!is.na(event_time) & !is.na(baseline_date) & event_time <= baseline_date)
    )]
  }

  hba1c_value <- suppressWarnings(as.numeric(dt[[hba1c_col]]))
  dt[, hba1c_baseline_available := as.integer(!is.na(hba1c_value))]
  dt[, diabetes_hba1c_threshold := 48]
  dt[, diabetes_hba1c_baseline := as.integer(!is.na(hba1c_value) & hba1c_value >= 48)]

  dt[, glucose_lowering_med_baseline := as.integer(flag_is_one(get(glucose_med_col)))]
  dt[, diabetes_self_report_baseline := as.integer(flag_is_one(get(self_report_col)))]
  dt[, diabetes_baseline := as.integer(
    diabetes_icd_baseline == 1L |
      diabetes_self_report_baseline == 1L |
      diabetes_hba1c_baseline == 1L |
      glucose_lowering_med_baseline == 1L
  )]

  dt[, diabetes_baseline_source := mapply(
    function(icd, self_report, hba1c, medicine) {
      sources <- c("ICD_E10_E14", "self_report_20002", "HbA1c", "glucose_lowering_med")[
        c(icd, self_report, hba1c, medicine) == 1L
      ]
      if (length(sources) == 0L) "no_evidence" else paste(sources, collapse = ";")
    },
    diabetes_icd_baseline,
    diabetes_self_report_baseline,
    diabetes_hba1c_baseline,
    glucose_lowering_med_baseline,
    USE.NAMES = FALSE
  )]
  dt[, diabetes_main_exclusion_reason := fcase(
    diabetes_icd_baseline == 1L, "diabetes_icd",
    diabetes_self_report_baseline == 1L, "diabetes_self_report",
    diabetes_hba1c_baseline == 1L, "diabetes_hba1c",
    glucose_lowering_med_baseline == 1L, "diabetes_glucose_lowering_med",
    default = "eligible"
  )]
  dt[]
}


# ------------------------------------------------------------------------------ #
# derive_ckd_baseline_status()
#
# 说明：
#   合并 ICD N18 既往史和基线 eGFR <60 来定义基线 CKD。
# ------------------------------------------------------------------------------ #
derive_ckd_baseline_status <- function(dt,
                                       baseline_col = "InCohort_time",
                                       creatinine_col = "Creatinine_0",
                                       age_col = "Age_at_recruitment",
                                       sex_col = "Sex") {
  dt <- copy(dt)
  setDT(dt)
  require_columns(dt, c(baseline_col, creatinine_col, age_col, sex_col), "CKD 判定")

  dt[, ckd_icd_baseline := as.integer(prevalent_from_time_or_flag(
    dt,
    time_col = "CKD_time",
    flag_col = "prevalent_CKD_ICD10_N18",
    baseline_col = baseline_col
  ))]

  dt[, egfr_2021_baseline := calculate_egfr_2021(
    creatinine_umol_l = get(creatinine_col),
    age = get(age_col),
    sex = get(sex_col)
  )]
  dt[, ckd_egfr_baseline := as.integer(!is.na(egfr_2021_baseline) & egfr_2021_baseline < 60)]
  dt[, CKD_baseline := as.integer(ckd_icd_baseline == 1L | ckd_egfr_baseline == 1L)]
  dt[, CKD_baseline_source := fcase(
    ckd_icd_baseline == 1L, "ICD_N18",
    ckd_egfr_baseline == 1L, "eGFR_lt60",
    default = "no_evidence"
  )]
  dt[]
}


# ------------------------------------------------------------------------------ #
# derive_markedly_elevated_risk_status()
#
# 说明：
#   标记指南中的 markedly elevated single risk factors：
#   TC >8 mmol/L、LDL-C >4.9 mmol/L、SBP >=180 或 DBP >=110。
# ------------------------------------------------------------------------------ #
derive_markedly_elevated_risk_status <- function(dt,
                                                 tc_col = "Cholesterol_0",
                                                 ldl_col = "LDL_direct_0",
                                                 sbp_col = "SystolicBP0",
                                                 dbp_col = "DiastolicBP0") {
  dt <- copy(dt)
  setDT(dt)
  require_columns(dt, c(tc_col, ldl_col, sbp_col, dbp_col), "markedly elevated single risk factors 判定")

  tc_value <- suppressWarnings(as.numeric(dt[[tc_col]]))
  ldl_value <- suppressWarnings(as.numeric(dt[[ldl_col]]))
  sbp_value <- suppressWarnings(as.numeric(dt[[sbp_col]]))
  dbp_value <- suppressWarnings(as.numeric(dt[[dbp_col]]))

  dt[, markedly_high_tc := as.integer(!is.na(tc_value) & tc_value > 8)]
  dt[, markedly_high_ldl := as.integer(!is.na(ldl_value) & ldl_value > 4.9)]
  dt[, markedly_high_bp := as.integer(
    (!is.na(sbp_value) & sbp_value >= 180) |
      (!is.na(dbp_value) & dbp_value >= 110)
  )]
  dt[, markedly_elevated_single_risk_factor := as.integer(
    markedly_high_tc == 1L |
      markedly_high_ldl == 1L |
      markedly_high_bp == 1L
  )]
  dt[]
}


# ------------------------------------------------------------------------------ #
# classify_lipid_baseline_status()
#
# 说明：
#   按“基线自报为基准，365 天窗口 GP 处方为辅助”的规则，
#   判定基线时降脂药是否正在使用。
#
# 主分析规则：
#   1. 基线自报 = 1                       -> treated
#   2. 基线自报 = 0 且 365 天内 GP = 1      -> treated（证据冲突时从严视作 treated）
#   3. 基线自报 = 0 且 365 天内 GP = 0/NA   -> untreated
#   4. 基线自报缺失 且 365 天内 GP = 1      -> treated
#   5. 基线自报缺失 且 365 天内 GP = 0/NA   -> NA（信息不足，主分析排除）
#
# 示例：
#   - 某人基线自报服药=1，哪怕 365 天内 GP 没抓到处方，也记为 treated。
#   - 某人基线自报未服药=0，但 365 天内 GP 有处方，也记为 treated。
#   - 某人基线自报缺失，且 365 天内 GP 无处方证据，则记为 NA。
#
# 这里生成的两个核心变量是：
#   - Lipid_Medicine_baseline：基线时是否服用降脂药
#   - Lipid_Medicine_time：与这个“是否服药”判定对应的时间变量
# ------------------------------------------------------------------------------ #
classify_lipid_baseline_status <- function(dt,
                                           self_col = "Lipid_medication_self_report_0",
                                           gp365_col = "gp_lipid_rx_in_365d",
                                           gp_last_col = "last_rx_time",
                                           baseline_col = "InCohort_time") {
  dt <- copy(dt)
  setDT(dt)

  dt[, lipid_self_report_0_bin := as_binary_self_report(get(self_col))]

  if (!gp365_col %in% names(dt)) dt[, (gp365_col) := NA_integer_]
  dt[, (gp365_col) := fifelse(is.na(get(gp365_col)), NA_integer_, as.integer(get(gp365_col)))]

  if (!gp_last_col %in% names(dt)) {
    dt[, (gp_last_col) := as.IDate(NA)]
  } else {
    dt[, (gp_last_col) := as.IDate(get(gp_last_col))]
  }

  dt[, baseline_date_tmp := as.IDate(get(baseline_col))]

  dt[, Lipid_Medicine_baseline := fcase(
    lipid_self_report_0_bin == 1L, 1L,
    lipid_self_report_0_bin == 0L & get(gp365_col) == 1L, 1L,
    lipid_self_report_0_bin == 0L & (is.na(get(gp365_col)) | get(gp365_col) == 0L), 0L,
    is.na(lipid_self_report_0_bin) & get(gp365_col) == 1L, 1L,
    default = NA_integer_
  )]

  dt[, lipid_status_source := fcase(
    # lipid_self_report_0_bin：把 self_col 转成二元后的基线自报结果
    lipid_self_report_0_bin == 1L & (is.na(get(gp365_col)) | get(gp365_col) == 0L), "self_report_0_only",
    lipid_self_report_0_bin == 1L & get(gp365_col) == 1L, "self_report_0+gp_365d",
    lipid_self_report_0_bin == 0L & get(gp365_col) == 1L, "conflict_self0_gp1_treated",
    lipid_self_report_0_bin == 0L & (is.na(get(gp365_col)) | get(gp365_col) == 0L), "self_report_0+no_gp365d",
    is.na(lipid_self_report_0_bin) & get(gp365_col) == 1L, "gp_365d_only",
    default = "insufficient_info"
  )]

  dt[, Lipid_Medicine := Lipid_Medicine_baseline]
  dt[, Lipid_Medicine_time := as.IDate(NA)]
  # 时间判定 情况 A：基线判定为服药，且依据是“基线自报服药”
  dt[Lipid_Medicine_baseline == 1L & lipid_self_report_0_bin == 1L, Lipid_Medicine_time := baseline_date_tmp]
  # 时间判定 情况 B：基线判定为服药，但不是靠基线自报，而是靠 GP 记录
  dt[Lipid_Medicine_baseline == 1L & lipid_self_report_0_bin != 1L & !is.na(get(gp_last_col)), Lipid_Medicine_time := get(gp_last_col)]
  # 时间判定 情况 C：基线判定为未服药，那么时间也记为基线日期
  dt[Lipid_Medicine_baseline == 0L, Lipid_Medicine_time := baseline_date_tmp]
  dt[, baseline_date_tmp := NULL]
  dt[]
}


# ------------------------------------------------------------------------------ #
# classify_hyper_baseline_status()
#
# 说明：
#   仅根据基线自报（Instance 0）判定基线时是否使用降压药。
#   本研究主分析不再使用 GP 降压药处方作为判定依据。
#
# 示例：
#   - 基线自报=1 -> hyper_Medicine_baseline = 1
#   - 基线自报=0 -> hyper_Medicine_baseline = 0
#   - 基线自报缺失 -> NA；04 不因该变量单独排除，后续作为 MI 变量处理
# ------------------------------------------------------------------------------ #
classify_hyper_baseline_status <- function(dt,
                                           self_col = "hyper_medication_self_report_0",
                                           baseline_col = "InCohort_time") {
  dt <- copy(dt)
  setDT(dt)
  dt[, hyper_self_report_0_bin := as_binary_self_report(get(self_col))]
  dt[, hyper_Medicine_baseline := hyper_self_report_0_bin]
  dt[, hyper_Medicine := hyper_Medicine_baseline]
  dt[, hyper_Medicine_time := as.IDate(fifelse(!is.na(hyper_Medicine_baseline), as.IDate(get(baseline_col)), as.IDate(NA)))]
  dt[]
}


# ------------------------------------------------------------------------------ #
# apply_required_observed_information_screen()
#
# 说明：
#   仅对“没有真实观测就无法确定 eligibility 或 guideline strategy”的变量进行
#   最终 required-information gate。该步骤不是传统 complete-case analysis：
#   ethnicity、TDI、BMI、SES 组成变量、raw MET、TG、降压药、家族史和慢性病等
#   outcome-model / subgroup 变量允许缺失，并保留给 06 multiple imputation。
#
# 返回：
#   - dt：通过 required-information gate 的最终 analytic cohort；
#   - n_before / n_after；
#   - vars_required_observed；
#   - missing_summary：gate 前各 required variable 的缺失情况。
# ------------------------------------------------------------------------------ #
apply_required_observed_information_screen <- function(
    dt,
    vars_required_observed,
    context = "eligibility / guideline-strategy required-information gate") {

  dt <- copy(as.data.table(dt))
  vars_required_observed <- unique(
    vars_required_observed[!is.na(vars_required_observed) & nzchar(vars_required_observed)]
  )
  require_columns(dt, vars_required_observed, context)

  n_before <- nrow(dt)

  char_vars <- vars_required_observed[vapply(
    vars_required_observed,
    function(v) is.character(dt[[v]]) || is.factor(dt[[v]]),
    logical(1)
  )]
  if (length(char_vars) > 0L) {
    dt[, (char_vars) := lapply(.SD, clean_character_missing_values), .SDcols = char_vars]
  }

  missing_summary <- data.table(
    variable = vars_required_observed,
    n_missing = vapply(vars_required_observed, function(v) sum(is.na(dt[[v]])), integer(1))
  )
  # n_before 是单个标量，而 n_missing 是按变量形成的向量。
  # 这里不能使用 data.table::fifelse()：fifelse 要求 test 与 yes/no 的长度
  # 相同（或 yes/no 长度为 1）。当 test 长度为 1、yes 为多变量向量时会报错。
  # 使用普通 if/else 先判断分母，再一次性生成整列缺失比例。
  missing_summary[, `:=`(
    denominator_n = n_before,
    missing_pct = if (n_before > 0L) 100 * n_missing / n_before else NA_real_
  )]
  setorder(missing_summary, -n_missing, variable)

  dt[, required_observed_information_complete :=
       as.integer(complete.cases(.SD)), .SDcols = vars_required_observed]
  dt_after <- dt[required_observed_information_complete == 1L]

  list(
    dt = dt_after[],
    n_before = n_before,
    n_after = nrow(dt_after),
    vars_required_observed = vars_required_observed,
    missing_summary = missing_summary[]
  )
}


# ------------------------------------------------------------------------------ #
# build_imputation_missingness_summary()
#
# 说明：
#   对“后续会在 06 插补”的变量进行描述性缺失报告。这里的缺失绝不用于排除。
#   同时报告：
#     1) 临床排除条件完成后、required-information gate 前；
#     2) 最终 analytic cohort（通过 required-information gate 后）。
#
#   额外添加 "Any subsequently imputed variable" 行，表示至少一个待插补变量缺失。
# ------------------------------------------------------------------------------ #
build_imputation_missingness_summary <- function(
    pre_gate_dt,
    final_dt,
    variable_methods) {

  pre <- copy(as.data.table(pre_gate_dt))
  fin <- copy(as.data.table(final_dt))
  vm <- copy(as.data.table(variable_methods))

  require_columns(vm, c("variable", "imputation_method"), "待插补变量方法表")
  vars <- unique(vm$variable)
  require_columns(pre, vars, "required-information gate 前数据")
  require_columns(fin, vars, "最终 analytic cohort")

  one_var <- function(v, method) {
    data.table(
      variable = v,
      imputation_method = method,
      n_pre_required_gate = nrow(pre),
      n_missing_pre_required_gate = sum(is.na(pre[[v]])),
      pct_missing_pre_required_gate = if (nrow(pre) > 0L) 100 * mean(is.na(pre[[v]])) else NA_real_,
      n_final_analytic_cohort = nrow(fin),
      n_missing_final_analytic_cohort = sum(is.na(fin[[v]])),
      pct_missing_final_analytic_cohort = if (nrow(fin) > 0L) 100 * mean(is.na(fin[[v]])) else NA_real_
    )
  }

  out <- rbindlist(
    lapply(seq_len(nrow(vm)), function(i) one_var(vm$variable[i], vm$imputation_method[i])),
    use.names = TRUE
  )

  pre_any <- if (length(vars)) !complete.cases(pre[, ..vars]) else rep(FALSE, nrow(pre))
  fin_any <- if (length(vars)) !complete.cases(fin[, ..vars]) else rep(FALSE, nrow(fin))

  any_row <- data.table(
    variable = "Any subsequently imputed variable",
    imputation_method = "multiple",
    n_pre_required_gate = nrow(pre),
    n_missing_pre_required_gate = sum(pre_any),
    pct_missing_pre_required_gate = if (nrow(pre) > 0L) 100 * mean(pre_any) else NA_real_,
    n_final_analytic_cohort = nrow(fin),
    n_missing_final_analytic_cohort = sum(fin_any),
    pct_missing_final_analytic_cohort = if (nrow(fin) > 0L) 100 * mean(fin_any) else NA_real_
  )

  rbindlist(list(out, any_row), use.names = TRUE)[]
}


# ------------------------------------------------------------------------------ #
# validate_flow_log() / build_eligibility_invariant_audit()
# ------------------------------------------------------------------------------ #
validate_flow_log <- function(flow_log, label = "纳排流程") {
  flow_log <- as.data.table(flow_log)
  require_columns(
    flow_log,
    c("step_no", "n_before", "n_excluded", "n_remaining"),
    label
  )
  arithmetic_ok <- flow_log$n_before - flow_log$n_excluded == flow_log$n_remaining
  sequence_ok <- c(TRUE, flow_log$n_before[-1L] == flow_log$n_remaining[-nrow(flow_log)])
  if (nrow(flow_log) > 0L && (!all(arithmetic_ok) || !all(sequence_ok))) {
    stop(label, " 的人数算术或相邻步骤衔接不一致。", call. = FALSE)
  }
  invisible(TRUE)
}

build_eligibility_invariant_audit <- function(dt) {
  dt <- as.data.table(dt)
  assert_unique_id(dt, "04 最终 eligible 数据")
  checks <- data.table(
    check = c(
      "one_row_per_id",
      "documented_ascvd_absent",
      "diabetes_absent",
      "ckd_absent",
      "lipid_treatment_absent",
      "markedly_elevated_risk_absent",
      "required_observed_information_complete"
    ),
    n_violations = c(
      nrow(dt) - uniqueN(dt[["Participant ID"]]),
      sum(is.na(dt$documented_ascvd_baseline) | dt$documented_ascvd_baseline != 0L),
      sum(is.na(dt$diabetes_baseline) | dt$diabetes_baseline != 0L),
      sum(is.na(dt$CKD_baseline) | dt$CKD_baseline != 0L),
      sum(is.na(dt$Lipid_Medicine_baseline) | dt$Lipid_Medicine_baseline != 0L),
      sum(is.na(dt$markedly_elevated_single_risk_factor) | dt$markedly_elevated_single_risk_factor != 0L),
      sum(
        is.na(dt$required_observed_information_complete) |
          dt$required_observed_information_complete != 1L
      )
    )
  )
  checks[, passed := n_violations == 0L]
  if (any(!checks$passed)) {
    stop("04 最终纳排不变量存在违规。", call. = FALSE)
  }
  checks[]
}


# ------------------------------------------------------------------------------ #
# run_eligibility_unit_tests()
# ------------------------------------------------------------------------------ #
run_eligibility_unit_tests <- function() {
  diabetes_test <- data.table(
    `Participant ID` = 1:4,
    InCohort_time = as.IDate(rep("2020-01-01", 4L)),
    Diabetes_ICD_time = as.IDate(c(NA, NA, "2020-01-02", NA)),
    T1D_time = as.IDate(NA),
    T2D_time = as.IDate(NA),
    `Glycated_haemoglobin_(HbA1c)_0` = c(47.9, 48, 30, 30),
    glucose_lowering_med_baseline = c(0L, 0L, 0L, 0L),
    diabetes_self_report_baseline = c(0L, 0L, 0L, 1L)
  )
  diabetes_result <- derive_diabetes_baseline_status(diabetes_test)
  stopifnot(
    identical(diabetes_result$diabetes_baseline, c(0L, 1L, 0L, 1L)),
    all(diabetes_result$diabetes_hba1c_threshold == 48)
  )

  gate_test <- data.table(
    `Participant ID` = 1:3,
    age = c(50, 51, 52),
    smoking = c("Never", NA, "Current"),
    optional_covariate = c(NA_real_, 1, NA_real_)
  )
  gate_result <- apply_required_observed_information_screen(
    gate_test,
    c("age", "smoking")
  )
  stopifnot(
    nrow(gate_result$dt) == 2L,
    all(gate_result$dt$`Participant ID` %in% c(1L, 3L)),
    all(is.na(gate_result$dt$optional_covariate))
  )

  missingness_test <- build_imputation_missingness_summary(
    pre_gate_dt = gate_test,
    final_dt = gate_result$dt,
    variable_methods = data.table(
      variable = "optional_covariate",
      imputation_method = "pmm"
    )
  )
  stopifnot(
    missingness_test[variable == "optional_covariate", n_missing_pre_required_gate] == 2L,
    missingness_test[variable == "optional_covariate", n_missing_final_analytic_cohort] == 2L
  )

  flow_test <- init_flow_log()
  flow_test <- add_flow_step(flow_test, "a", 10L, 8L, 10L)
  flow_test <- add_flow_step(flow_test, "b", 8L, 7L, 10L)
  validate_flow_log(flow_test, "单元测试流程")
  step_message("04 边界单元测试通过：HbA1c 48 固定阈值；required-information gate 仅删除策略/资格必要信息缺失；普通 MI 协变量缺失保留；流程算术一致。")
  invisible(TRUE)
}
