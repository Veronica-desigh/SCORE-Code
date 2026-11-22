# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_03_Clean_medicine_M01.R
# 任务：供03调用的用药识别与GP归并函数
# 输入：自报和处方原表
# 输出：一人一行用药状态
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function_03_Clean_medicine_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   为 03_Clean_medicine_M01.R 提供用药信息清理函数。
#
# 设计原则：
#   1. 03 脚本只清理“是否有用药证据”，不计算剂量、频次或处方数量；
#   2. 路径、字段名、代码集和药名规则在主脚本中声明，函数只负责读取、校验和派生；
#   3. 自报多选题按选项 token 精确识别，避免枚举所有组合导致漏项；
#   4. GP scripts 同一条处方若同时有 BNF、dm+d、Read v2 编码，按“任一编码命中”判定；
#   5. 主分析以 Instance 0 和基线前 0-365 天 GP 证据判定；另保留 90 天辅助标记；
#   6. 本应一人一行的来源表发现重复 ID 时立即停止，GP 多条处方则按规则聚合。
#   7. 降压药自报缺失必须保留为 NA；不得把“无记录”自动当成“未使用”。
#      该变量不承担 SCORE/SCORE2 或 eligibility assignment，后续由 MI 处理。
# ============================================================================== #

library(data.table)

# ============================================================================== #
# 0. 通用检查、计时和格式化函数 ####
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
  if (!is.null(label)) {
    step_message("开始：", label)
  }

  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(label)) {
    step_message("完成：", label, " | 耗时 = ", fmt_time(elapsed))
  }

  if (exists("m01_record_timing", mode = "function", inherits = TRUE))
    m01_record_timing(if (is.null(label)) "unnamed_run_with_timer" else label, elapsed)
  list(value = value, elapsed_sec = elapsed)
}

assert_file_exists <- function(path, label = "文件") {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(label, "不存在：", path, call. = FALSE)
  }
  invisible(path)
}

assert_required_columns <- function(dt, cols, label = "数据表") {
  dt <- as.data.table(dt)
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  missing_cols <- setdiff(cols, names(dt))

  if (length(missing_cols) > 0L) {
    stop(
      label, "缺少关键列：", paste(missing_cols, collapse = ", "),
      "\n当前前 20 个列名：", paste(utils::head(names(dt), 20), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_unique_id <- function(dt, label = "数据表", id_col = "Participant ID") {
  dt <- as.data.table(dt)
  assert_required_columns(dt, id_col, label)
  if (anyNA(dt[[id_col]])) stop(label, "存在缺失 ID。", call. = FALSE)
  duplicate_ids <- dt[, .N, by = id_col][N > 1L]
  if (nrow(duplicate_ids) > 0L) {
    stop(
      label, "存在 ", nrow(duplicate_ids), " 个重复 ID；不得按文件顺序静默取第一条。",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

merge_one_to_one_checked <- function(main_dt,
                                     add_dt,
                                     label,
                                     by = "Participant ID") {
  main_dt <- copy(as.data.table(main_dt))
  add_dt <- copy(as.data.table(add_dt))
  assert_unique_id(main_dt, paste0(label, "主表"), by)
  assert_unique_id(add_dt, paste0(label, "辅助表"), by)
  n_before <- nrow(main_dt)
  out <- merge(main_dt, add_dt, by = by, all.x = TRUE, sort = FALSE)
  if (nrow(out) != n_before || uniqueN(out[[by]]) != n_before) {
    stop(label, "改变了主表行数或唯一 ID 数。", call. = FALSE)
  }
  out[]
}

clean_character_missing <- function(x) {
  x <- trimws(as.character(x))
  missing_strings <- c(
    "",
    "NA",
    "N/A",
    "NaN",
    "NULL",
    "Do not know",
    "Prefer not to answer"
  )
  x[is.na(x) | x %in% missing_strings] <- NA_character_
  x
}

standardize_code_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "N/A", "NaN", "NULL")] <- NA_character_
  x
}

normalize_drug_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | x %in% c("", "na", "n/a", "nan", "null")] <- NA_character_
  x
}

standardize_output_id <- function(dt) {
  dt <- copy(as.data.table(dt))

  id_candidates <- c(
    "Participant ID",
    "Participant ID(gp_scripts - eid)",
    "Participant ID(participant - eid)",
    "Participant_ID",
    "eid"
  )
  hit <- intersect(id_candidates, names(dt))

  if (length(hit) == 0L) {
    stop(
      "当前数据集中未找到可标准化的 ID 列。现有列名为：",
      paste(names(dt), collapse = ", "),
      call. = FALSE
    )
  }

  if (hit[1] != "Participant ID") {
    setnames(dt, hit[1], "Participant ID")
  }

  dt[]
}

var_info_print <- function(data, vars) {
  data <- as.data.table(data)
  vars <- intersect(vars, names(data))

  if (length(vars) == 0L) {
    message("未找到需要检查的变量。")
    return(invisible(NULL))
  }

  out <- data.table(
    variable = vars,
    class = vapply(data[, ..vars], function(x) class(x)[1], character(1)),
    n_missing = vapply(data[, ..vars], function(x) sum(is.na(x)), integer(1)),
    n_unique = vapply(data[, ..vars], function(x) uniqueN(x, na.rm = TRUE), integer(1))
  )

  print(out)
  invisible(out)
}

# ============================================================================== #
# 1. 基线人群信息 ####
# ============================================================================== #

read_baseline_medicine_base <- function(data_exposure_path,
                                        id_col = "Participant ID",
                                        baseline_col = "InCohort_time",
                                        sex_col = "Sex") {
  assert_file_exists(data_exposure_path, "M01 基线暴露清理结果")

  keep_cols <- c(id_col, baseline_col, sex_col)
  dt <- fread(data_exposure_path, select = keep_cols)
  setDT(dt)
  assert_required_columns(dt, keep_cols, "M01 基线暴露清理结果")

  out <- dt[, .(
    `Participant ID` = get(id_col),
    baseline_time = as.IDate(get(baseline_col)),
    Sex = as.character(get(sex_col))
  )]

  if (anyNA(out[["Participant ID"]]) || anyNA(out$baseline_time)) {
    stop("M01 基线暴露清理结果存在缺失 ID 或基线日期。", call. = FALSE)
  }
  assert_unique_id(out, "M01 基线用药主表")
  setkey(out, `Participant ID`)
  out[]
}

# ============================================================================== #
# 2. Touchscreen 自报用药 ####
# ============================================================================== #

token_flag_from_multiselect <- function(x, token) {
  x <- clean_character_missing(x)
  out <- rep(NA_integer_, length(x))
  idx <- !is.na(x)

  if (any(idx)) {
    values <- unique(x[idx])
    flags <- as.integer(vapply(
      strsplit(values, "|", fixed = TRUE),
      function(tokens) token %in% trimws(tokens),
      logical(1)
    ))
    out[idx] <- flags[match(x[idx], values)]
  }

  out
}

build_touchscreen_med_self_report <- function(touchscreen_file,
                                              baseline_dt,
                                              male_columns,
                                              female_columns) {
  assert_file_exists(touchscreen_file, "Touchscreen 自报用药文件")

  baseline_dt <- copy(as.data.table(baseline_dt))
  assert_required_columns(baseline_dt, c("Participant ID", "Sex"), "基线人群表")

  select_cols <- unique(c("Participant ID", unname(male_columns), unname(female_columns)))
  raw <- fread(touchscreen_file, select = select_cols)
  setDT(raw)
  assert_required_columns(raw, select_cols, "Touchscreen 自报用药文件")
  assert_unique_id(baseline_dt, "Touchscreen 合并基线表")
  assert_unique_id(raw, "Touchscreen 自报用药文件")

  dt <- merge_one_to_one_checked(
    baseline_dt[, .(`Participant ID`, Sex)],
    raw,
    label = "Touchscreen 自报用药合并"
  )
  setDT(dt)

  for (instance in names(male_columns)) {
    male_col <- male_columns[[instance]]
    female_col <- female_columns[[instance]]
    med_value <- rep(NA_character_, nrow(dt))

    male_idx <- dt$Sex == "Male" & !is.na(dt$Sex)
    female_idx <- dt$Sex == "Female" & !is.na(dt$Sex)

    med_value[male_idx] <- clean_character_missing(dt[[male_col]][male_idx])
    med_value[female_idx] <- clean_character_missing(dt[[female_col]][female_idx])

    dt[, (paste0("touchscreen_med_raw_", instance)) := med_value]
    dt[, (paste0("Lipid_medication_self_report_", instance)) :=
      token_flag_from_multiselect(med_value, "Cholesterol lowering medication")]
    dt[, (paste0("hyper_medication_self_report_", instance)) :=
      token_flag_from_multiselect(med_value, "Blood pressure medication")]
    dt[, (paste0("insulin_self_report_", instance)) :=
      token_flag_from_multiselect(med_value, "Insulin")]
  }

  lipid_cols <- paste0("Lipid_medication_self_report_", names(male_columns))
  hyper_cols <- paste0("hyper_medication_self_report_", names(male_columns))
  insulin_cols <- paste0("insulin_self_report_", names(male_columns))

  lipid_dt <- dt[, c("Participant ID", lipid_cols), with = FALSE]
  hyper_dt <- dt[, c("Participant ID", hyper_cols), with = FALSE]
  insulin_dt <- dt[, c("Participant ID", insulin_cols), with = FALSE]

  setkey(lipid_dt, `Participant ID`)
  setkey(hyper_dt, `Participant ID`)
  setkey(insulin_dt, `Participant ID`)

  list(
    Lipid_med_self_report = lipid_dt[],
    hyper_med_self_report = hyper_dt[],
    insulin_self_report = insulin_dt[]
  )
}

# ============================================================================== #
# 3. Field 20003 verbal interview 药物名称识别 ####
# ============================================================================== #

match_drug_targets <- function(x, target_names = character(), target_patterns = character()) {
  x_clean <- normalize_drug_name(x)
  out <- rep(FALSE, length(x_clean))

  if (length(target_names) > 0L) {
    out <- out | (x_clean %chin% normalize_drug_name(target_names))
  }

  idx <- !is.na(x_clean)
  if (any(idx) && length(target_patterns) > 0L) {
    pattern_hit <- rep(FALSE, sum(idx))
    x_idx <- x_clean[idx]

    for (pat in target_patterns) {
      pattern_hit <- pattern_hit | grepl(pat, x_idx, perl = TRUE)
    }

    out[idx] <- out[idx] | pattern_hit
  }

  out[is.na(x_clean)] <- FALSE
  out
}

build_20003_medication_flag <- function(medication_file,
                                        instance = 0L,
                                        max_array = 47L,
                                        target_names = character(),
                                        target_patterns = character(),
                                        flag_prefix = "glucose_lowering_med_20003") {
  assert_file_exists(medication_file, "Field 20003 verbal interview 用药文件")

  header <- names(fread(medication_file, nrows = 0L))
  med_cols_expected <- paste0(
    "Treatment/medication code | Instance ", instance, " | Array ", 0:max_array
  )
  med_cols <- intersect(med_cols_expected, header)

  if (length(med_cols) == 0L) {
    stop(
      "Field 20003 文件中未找到 Instance ", instance,
      " 的 Treatment/medication code 数组列。",
      call. = FALSE
    )
  }

  if (length(setdiff(med_cols_expected, med_cols)) > 0L) {
    warning(
      "Field 20003 文件缺少部分数组列：",
      paste(setdiff(med_cols_expected, med_cols), collapse = ", "),
      call. = FALSE
    )
  }

  dt <- fread(medication_file, select = c("Participant ID", med_cols))
  setDT(dt)
  assert_required_columns(dt, c("Participant ID", med_cols), "Field 20003 verbal interview 用药文件")
  assert_unique_id(dt, "Field 20003 verbal interview 用药文件")

  out <- dt[, .(`Participant ID`)]
  out[, (paste0(flag_prefix, "_", instance)) := 0L]
  out[, (paste0(flag_prefix, "_n_match")) := 0L]
  out[, (paste0(flag_prefix, "_names")) := NA_character_]

  hit_long <- rbindlist(
    lapply(med_cols, function(col) {
      hit <- match_drug_targets(
        dt[[col]],
        target_names = target_names,
        target_patterns = target_patterns
      )

      if (!any(hit)) {
        return(NULL)
      }

      data.table(
        `Participant ID` = dt[hit, `Participant ID`],
        drug_name_clean = normalize_drug_name(dt[[col]][hit])
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )

  if (!is.null(hit_long) && nrow(hit_long) > 0L) {
    hit_summary <- hit_long[, .(
      flag = 1L,
      n_match = .N,
      names = paste(sort(unique(drug_name_clean)), collapse = "; ")
    ), by = `Participant ID`]

    setnames(
      hit_summary,
      c("flag", "n_match", "names"),
      c(
        paste0(flag_prefix, "_", instance),
        paste0(flag_prefix, "_n_match"),
        paste0(flag_prefix, "_names")
      )
    )

    out <- merge(out[, .(`Participant ID`)], hit_summary, by = "Participant ID", all.x = TRUE)
    out[is.na(get(paste0(flag_prefix, "_", instance))), (paste0(flag_prefix, "_", instance)) := 0L]
    out[is.na(get(paste0(flag_prefix, "_n_match"))), (paste0(flag_prefix, "_n_match")) := 0L]
  }

  setkey(out, `Participant ID`)
  out[]
}

build_glucose_med_baseline <- function(baseline_dt,
                                       glucose_20003_dt,
                                       insulin_self_report_dt,
                                       instance = 0L) {
  baseline_dt <- copy(as.data.table(baseline_dt))
  glucose_20003_dt <- copy(as.data.table(glucose_20003_dt))
  insulin_self_report_dt <- copy(as.data.table(insulin_self_report_dt))

  assert_required_columns(baseline_dt, "Participant ID", "基线人群表")
  assert_required_columns(glucose_20003_dt, c("Participant ID", paste0("glucose_lowering_med_20003_", instance)), "20003 降糖用药表")
  assert_required_columns(insulin_self_report_dt, c("Participant ID", paste0("insulin_self_report_", instance)), "Touchscreen 胰岛素自报表")

  assert_unique_id(baseline_dt, "降糖用药基线表")
  assert_unique_id(glucose_20003_dt, "20003 降糖用药表")
  assert_unique_id(insulin_self_report_dt, "Touchscreen 胰岛素自报表")

  dt <- merge_one_to_one_checked(
    baseline_dt[, .(`Participant ID`)],
    glucose_20003_dt,
    label = "20003 降糖用药合并"
  )
  dt <- merge_one_to_one_checked(
    dt,
    insulin_self_report_dt,
    label = "Touchscreen 胰岛素合并"
  )

  med20003_col <- paste0("glucose_lowering_med_20003_", instance)
  insulin_col <- paste0("insulin_self_report_", instance)

  dt[is.na(get(med20003_col)), (med20003_col) := 0L]

  med20003_yes <- dt[[med20003_col]] == 1L
  insulin_yes <- dt[[insulin_col]] == 1L
  med20003_yes[is.na(med20003_yes)] <- FALSE
  insulin_yes[is.na(insulin_yes)] <- FALSE

  dt[, glucose_lowering_med_baseline := as.integer(med20003_yes | insulin_yes)]

  dt[, glucose_lowering_med_source := fcase(
    med20003_yes & insulin_yes,
    "20003+touchscreen_insulin",
    med20003_yes & !insulin_yes,
    "20003_only",
    !med20003_yes & insulin_yes,
    "touchscreen_insulin_only",
    default = "no_evidence"
  )]

  setkey(dt, `Participant ID`)
  dt[]
}

# ============================================================================== #
# 4. GP scripts 降脂药处方 ####
# ============================================================================== #

clean_gp_scripts <- function(file,
                             lipid_prefilter_bnf = NULL,
                             lipid_prefilter_dmd = NULL,
                             lipid_prefilter_read2 = NULL) {
  assert_file_exists(file, "GP scripts 原始处方文件")

  keep_cols <- c(
    "Participant ID(gp_scripts - eid)",
    "Date prescription was issued",
    "Drug name",
    "BNF code",
    "dm+d code",
    "Read v2 code"
  )

  dt <- fread(file, select = keep_cols)
  setDT(dt)
  assert_required_columns(dt, keep_cols, "GP scripts 原始处方文件")

  dt <- dt[
    !is.na(`Participant ID(gp_scripts - eid)`) &
      `Participant ID(gp_scripts - eid)` != "" &
      !is.na(`Date prescription was issued`) &
      `Date prescription was issued` != ""
  ]

  # Optional exact-code early filter BEFORE converting millions of nonlipid
  # prescription dates. This uses exactly the same normalization/%chin%
  # membership as filter_lipid_gp; the latter still verifies all three codes.
  use_prefilter <- !is.null(lipid_prefilter_bnf) ||
    !is.null(lipid_prefilter_dmd) || !is.null(lipid_prefilter_read2)
  if (use_prefilter) {
    if (is.null(lipid_prefilter_bnf) || is.null(lipid_prefilter_dmd) ||
        is.null(lipid_prefilter_read2))
      stop("GP early-filter requires all three code dictionaries")
    keep <- standardize_code_value(dt[["BNF code"]]) %chin%
              standardize_code_value(lipid_prefilter_bnf) |
            standardize_code_value(dt[["dm+d code"]]) %chin%
              standardize_code_value(lipid_prefilter_dmd) |
            standardize_code_value(dt[["Read v2 code"]]) %chin%
              standardize_code_value(lipid_prefilter_read2)
    if (anyNA(keep)) stop("GP early-filter produced NA membership")
    if (exists("m01_checkpoint", mode = "function", inherits = TRUE))
      m01_checkpoint("03_GP_preprocessing", 1L, 3L,
        sprintf("exact code prefilter retained %s/%s prescriptions",
          format(sum(keep), big.mark = ","), format(nrow(dt), big.mark = ",")))
    dt <- dt[keep]
  }

  dt[, Prescription_time := as.IDate(`Date prescription was issued`)]
  dt <- dt[!is.na(Prescription_time)]

  for (col in c("Drug name", "BNF code", "dm+d code", "Read v2 code")) {
    dt[, (col) := standardize_code_value(get(col))]
  }

  dt[]
}

filter_lipid_gp <- function(dt,
                            bnf_codes,
                            dmd_codes,
                            read2_codes) {
  dt <- as.data.table(dt)

  bnf_codes <- standardize_code_value(bnf_codes)
  dmd_codes <- standardize_code_value(dmd_codes)
  read2_codes <- standardize_code_value(read2_codes)

  bnf_hit <- standardize_code_value(dt[["BNF code"]]) %chin% bnf_codes
  dmd_hit <- standardize_code_value(dt[["dm+d code"]]) %chin% dmd_codes
  read2_hit <- standardize_code_value(dt[["Read v2 code"]]) %chin% read2_codes
  keep <- bnf_hit | dmd_hit | read2_hit

  lipid_dt <- copy(dt[keep])
  lipid_dt[, lipid_match_bnf := bnf_hit[keep]]
  lipid_dt[, lipid_match_dmd := dmd_hit[keep]]
  lipid_dt[, lipid_match_read2 := read2_hit[keep]]
  lipid_dt[, lipid_match_source := ""]
  lipid_dt[lipid_match_bnf == TRUE, lipid_match_source := paste0(lipid_match_source, "BNF;")]
  lipid_dt[lipid_match_dmd == TRUE, lipid_match_source := paste0(lipid_match_source, "dm+d;")]
  lipid_dt[lipid_match_read2 == TRUE, lipid_match_source := paste0(lipid_match_source, "Read v2;")]
  lipid_dt[, lipid_match_source := sub(";$", "", lipid_match_source)]

  lipid_dt[]
}

build_lipid_history <- function(lipid_gp_dt) {
  dt <- copy(as.data.table(lipid_gp_dt))
  setorder(dt, `Participant ID(gp_scripts - eid)`, Prescription_time)
  dt[]
}

build_nearest_prescription_prebaseline <- function(history_dt, baseline_dt) {
  history_dt <- copy(as.data.table(history_dt))
  baseline_dt <- copy(as.data.table(baseline_dt))
  assert_unique_id(baseline_dt, "GP 处方基线表")

  dt <- merge(
    history_dt,
    baseline_dt,
    by.x = "Participant ID(gp_scripts - eid)",
    by.y = "Participant ID",
    all.x = FALSE,
    all.y = FALSE
  )
  setDT(dt)

  dt <- dt[Prescription_time <= baseline_time]

  if (nrow(dt) == 0L) {
    return(data.table())
  }

  dt[, days_before_baseline := as.integer(baseline_time - Prescription_time)]
  dt[, nearest_days__ := min(days_before_baseline), by = `Participant ID(gp_scripts - eid)`]
  nearest_dt <- dt[days_before_baseline == nearest_days__, .(
    Prescription_time = max(Prescription_time),
    days_before_baseline = nearest_days__[1L],
    gp_lipid_rx_prebaseline = 1L,
    `Drug name` = paste(sort(unique(na.omit(`Drug name`))), collapse = ";"),
    `BNF code` = paste(sort(unique(na.omit(`BNF code`))), collapse = ";"),
    `dm+d code` = paste(sort(unique(na.omit(`dm+d code`))), collapse = ";"),
    `Read v2 code` = paste(sort(unique(na.omit(`Read v2 code`))), collapse = ";"),
    lipid_match_bnf = any(lipid_match_bnf, na.rm = TRUE),
    lipid_match_dmd = any(lipid_match_dmd, na.rm = TRUE),
    lipid_match_read2 = any(lipid_match_read2, na.rm = TRUE),
    lipid_match_source = paste(sort(unique(na.omit(lipid_match_source))), collapse = ";")
  ), by = `Participant ID(gp_scripts - eid)`]
  for (v in c("Drug name", "BNF code", "dm+d code", "Read v2 code", "lipid_match_source")) {
    nearest_dt[get(v) == "", (v) := NA_character_]
  }
  assert_unique_id(nearest_dt, "最近基线前 GP 处方汇总", "Participant ID(gp_scripts - eid)")
  nearest_dt[]
}

build_window_summary <- function(history_dt, baseline_dt) {
  history_dt <- copy(as.data.table(history_dt))
  baseline_dt <- copy(as.data.table(baseline_dt))
  assert_unique_id(baseline_dt, "GP 时间窗基线表")

  dt <- merge(
    history_dt[, .(`Participant ID(gp_scripts - eid)`, Prescription_time)],
    baseline_dt,
    by.x = "Participant ID(gp_scripts - eid)",
    by.y = "Participant ID",
    all.x = FALSE,
    all.y = FALSE
  )
  setDT(dt)

  dt <- dt[Prescription_time <= baseline_time]

  if (nrow(dt) == 0L) {
    return(data.table())
  }

  dt[, diff_days := as.integer(baseline_time - Prescription_time)]

  # 边界冻结：基线当天为 0 天，第 365 天纳入，第 366 天不纳入。
  # n_lipid_rx_total 是处方记录数，仅用于质控和兼容既有 04 脚本，不代表用药剂量。
  window_dt <- dt[, .(
    gp_lipid_rx_in_90d = as.integer(any(diff_days >= 0L & diff_days <= 90L)),
    gp_lipid_rx_in_180d = as.integer(any(diff_days >= 0L & diff_days <= 180L)),
    gp_lipid_rx_in_365d = as.integer(any(diff_days >= 0L & diff_days <= 365L)),
    gp_lipid_rx_in_2y = as.integer(any(diff_days >= 0L & diff_days <= 730L)),
    gp_lipid_rx_in_5y = as.integer(any(diff_days >= 0L & diff_days <= 1825L)),
    first_rx_time = min(Prescription_time, na.rm = TRUE),
    last_rx_time = max(Prescription_time, na.rm = TRUE),
    n_lipid_rx_total = .N
  ), by = `Participant ID(gp_scripts - eid)`]

  assert_unique_id(window_dt, "GP 降脂药时间窗汇总", "Participant ID(gp_scripts - eid)")

  window_dt[]
}

# Efficient, contract-compatible joint summary for GP Rx; one join instead of
# running build_nearest_prescription_prebaseline + build_window_summary separately.
build_prebaseline_gp_summaries <- function(history_dt, baseline_dt) {
  history_dt <- data.table::as.data.table(history_dt)
  baseline_dt <- data.table::copy(data.table::as.data.table(baseline_dt))
  assert_unique_id(baseline_dt, "GP Rx reference cohort")
  id <- "Participant ID(gp_scripts - eid)"
  joined <- merge(history_dt, baseline_dt,
    by.x = id, by.y = "Participant ID", all = FALSE, sort = FALSE)
  joined <- joined[!is.na(Prescription_time) &
                   !is.na(baseline_time) & Prescription_time <= baseline_time]
  if (!nrow(joined)) return(list(nearest = data.table::data.table(),
                                 window = data.table::data.table()))
  joined[, days_before_baseline := as.integer(baseline_time - Prescription_time)]
  # Keep all ties on nearest date for correct union of BNF/dm+d/Read codes.
  joined[, nearest_days__ := min(days_before_baseline), by = id]
  nearest <- joined[days_before_baseline == nearest_days__, .(
    Prescription_time = max(Prescription_time),
    days_before_baseline = nearest_days__[1L],
    gp_lipid_rx_prebaseline = 1L,
    `Drug name` = paste(sort(unique(na.omit(`Drug name`))), collapse = ";"),
    `BNF code` = paste(sort(unique(na.omit(`BNF code`))), collapse = ";"),
    `dm+d code` = paste(sort(unique(na.omit(`dm+d code`))), collapse = ";"),
    `Read v2 code` = paste(sort(unique(na.omit(`Read v2 code`))), collapse = ";"),
    lipid_match_bnf = any(lipid_match_bnf, na.rm = TRUE),
    lipid_match_dmd = any(lipid_match_dmd, na.rm = TRUE),
    lipid_match_read2 = any(lipid_match_read2, na.rm = TRUE),
    lipid_match_source = paste(sort(unique(na.omit(lipid_match_source))), collapse = ";")
  ), by = id]
  for (v in c("Drug name", "BNF code", "dm+d code", "Read v2 code", "lipid_match_source"))
    nearest[get(v) == "", (v) := NA_character_]
  windows <- joined[, .(
    gp_lipid_rx_in_90d = as.integer(any(days_before_baseline <= 90L)),
    gp_lipid_rx_in_180d = as.integer(any(days_before_baseline <= 180L)),
    gp_lipid_rx_in_365d = as.integer(any(days_before_baseline <= 365L)),
    gp_lipid_rx_in_2y = as.integer(any(days_before_baseline <= 730L)),
    gp_lipid_rx_in_5y = as.integer(any(days_before_baseline <= 1825L)),
    first_rx_time = min(Prescription_time),
    last_rx_time = max(Prescription_time),
    n_lipid_rx_total = .N
  ), by = id]
  assert_unique_id(nearest, "Combined GP nearest", id)
  assert_unique_id(windows, "Combined GP windows", id)
  list(nearest = nearest[], window = windows[])
}

build_lipid_gp_main <- function(nearest_dt, window_dt, baseline_dt) {
  baseline_dt <- copy(as.data.table(baseline_dt))
  nearest_dt <- copy(as.data.table(nearest_dt))
  window_dt <- copy(as.data.table(window_dt))
  assert_unique_id(baseline_dt, "GP 降脂药主表基线数据")
  if (nrow(nearest_dt) > 0L) {
    assert_unique_id(nearest_dt, "GP 最近处方数据", "Participant ID(gp_scripts - eid)")
  }
  if (nrow(window_dt) > 0L) {
    assert_unique_id(window_dt, "GP 时间窗数据", "Participant ID(gp_scripts - eid)")
  }

  main_dt <- copy(baseline_dt)
  setnames(main_dt, "Participant ID", "Participant ID(participant - eid)")

  if (nrow(nearest_dt) > 0L) {
    tmp1 <- nearest_dt[, .(
      `Participant ID(participant - eid)` = `Participant ID(gp_scripts - eid)`,
      Prescription_time,
      days_before_baseline,
      gp_lipid_rx_prebaseline,
      `Drug name`,
      `BNF code`,
      `dm+d code`,
      `Read v2 code`,
      lipid_match_bnf,
      lipid_match_dmd,
      lipid_match_read2,
      lipid_match_source
    )]

    main_dt <- merge(
      main_dt,
      tmp1,
      by = "Participant ID(participant - eid)",
      all.x = TRUE
    )
  } else {
    main_dt[, `:=`(
      Prescription_time = as.IDate(NA),
      days_before_baseline = NA_integer_,
      gp_lipid_rx_prebaseline = NA_integer_,
      `Drug name` = NA_character_,
      `BNF code` = NA_character_,
      `dm+d code` = NA_character_,
      `Read v2 code` = NA_character_,
      lipid_match_bnf = NA,
      lipid_match_dmd = NA,
      lipid_match_read2 = NA,
      lipid_match_source = NA_character_
    )]
  }

  if (nrow(window_dt) > 0L) {
    tmp2 <- copy(window_dt)
    setnames(
      tmp2,
      "Participant ID(gp_scripts - eid)",
      "Participant ID(participant - eid)"
    )

    main_dt <- merge(
      main_dt,
      tmp2,
      by = "Participant ID(participant - eid)",
      all.x = TRUE
    )
  }

  window_flag_cols <- c(
    "gp_lipid_rx_in_90d",
    "gp_lipid_rx_in_180d",
    "gp_lipid_rx_in_365d",
    "gp_lipid_rx_in_2y",
    "gp_lipid_rx_in_5y"
  )

  missing_window_cols <- setdiff(window_flag_cols, names(main_dt))
  for (v in missing_window_cols) {
    main_dt[, (v) := 0L]
  }

  for (v in window_flag_cols) {
    main_dt[is.na(get(v)), (v) := 0L]
  }

  if (!"first_rx_time" %in% names(main_dt)) {
    main_dt[, first_rx_time := as.IDate(NA)]
  }
  if (!"last_rx_time" %in% names(main_dt)) {
    main_dt[, last_rx_time := as.IDate(NA)]
  }
  if (!"n_lipid_rx_total" %in% names(main_dt)) {
    main_dt[, n_lipid_rx_total := 0L]
  }

  main_dt[is.na(n_lipid_rx_total), n_lipid_rx_total := 0L]
  main_dt[is.na(gp_lipid_rx_prebaseline), gp_lipid_rx_prebaseline := 0L]

  # 兼容旧脚本接口：Prescription 表示基线前/当天是否有 GP 降脂药处方证据。
  main_dt[, Prescription := as.integer(gp_lipid_rx_prebaseline == 1L)]
  assert_unique_id(main_dt, "GP 降脂药主表", "Participant ID(participant - eid)")

  main_dt[]
}

# ============================================================================== #
# 5. 一人一行基线用药证据和来源审计 ####
# ============================================================================== #

as_medication_binary <- function(x) {
  x_chr <- clean_character_missing(x)
  x_num <- suppressWarnings(as.integer(x_chr))
  fifelse(is.na(x_num), NA_integer_, fifelse(x_num == 1L, 1L, fifelse(x_num == 0L, 0L, NA_integer_)))
}

build_medication_evidence_master <- function(baseline_dt,
                                             lipid_self_report_dt,
                                             hyper_self_report_dt,
                                             glucose_baseline_dt,
                                             lipid_gp_main_dt) {
  baseline_dt <- copy(as.data.table(baseline_dt))
  lipid_self_report_dt <- copy(as.data.table(lipid_self_report_dt))
  hyper_self_report_dt <- copy(as.data.table(hyper_self_report_dt))
  glucose_baseline_dt <- copy(as.data.table(glucose_baseline_dt))
  lipid_gp_main_dt <- standardize_output_id(lipid_gp_main_dt)

  required <- list(
    lipid_self_report = c("Participant ID", "Lipid_medication_self_report_0"),
    hyper_self_report = c("Participant ID", "hyper_medication_self_report_0"),
    glucose_baseline = c(
      "Participant ID", "glucose_lowering_med_20003_0", "insulin_self_report_0",
      "glucose_lowering_med_baseline", "glucose_lowering_med_source"
    ),
    lipid_gp = c("Participant ID", "gp_lipid_rx_in_90d", "gp_lipid_rx_in_365d")
  )
  assert_required_columns(lipid_self_report_dt, required$lipid_self_report, "降脂药自报表")
  assert_required_columns(hyper_self_report_dt, required$hyper_self_report, "降压药自报表")
  assert_required_columns(glucose_baseline_dt, required$glucose_baseline, "基线降糖用药表")
  assert_required_columns(lipid_gp_main_dt, required$lipid_gp, "GP 降脂药主表")

  dt <- baseline_dt[, .(`Participant ID`, baseline_time)]
  dt <- merge_one_to_one_checked(
    dt,
    lipid_self_report_dt[, required$lipid_self_report, with = FALSE],
    "统一用药表：降脂药自报"
  )
  dt <- merge_one_to_one_checked(
    dt,
    hyper_self_report_dt[, required$hyper_self_report, with = FALSE],
    "统一用药表：降压药自报"
  )
  dt <- merge_one_to_one_checked(
    dt,
    glucose_baseline_dt[, required$glucose_baseline, with = FALSE],
    "统一用药表：降糖用药"
  )
  dt <- merge_one_to_one_checked(
    dt,
    lipid_gp_main_dt[, required$lipid_gp, with = FALSE],
    "统一用药表：GP 降脂药"
  )

  dt[, lipid_self_report_instance0 := as_medication_binary(Lipid_medication_self_report_0)]
  dt[, lipid_gp_90d := as_medication_binary(gp_lipid_rx_in_90d)]
  dt[, lipid_gp_365d := as_medication_binary(gp_lipid_rx_in_365d)]
  dt[is.na(lipid_gp_90d), lipid_gp_90d := 0L]
  dt[is.na(lipid_gp_365d), lipid_gp_365d := 0L]
  dt[, lipid_lowering_med_baseline := fcase(
    lipid_self_report_instance0 == 1L | lipid_gp_365d == 1L, 1L,
    lipid_self_report_instance0 == 0L & lipid_gp_365d == 0L, 0L,
    default = NA_integer_
  )]
  dt[, lipid_lowering_med_source := fcase(
    lipid_self_report_instance0 == 1L & lipid_gp_365d == 1L, "touchscreen;gp_365d",
    lipid_self_report_instance0 == 1L, "touchscreen",
    lipid_gp_365d == 1L, "gp_365d",
    lipid_self_report_instance0 == 0L & lipid_gp_365d == 0L, "no_evidence",
    default = "unknown"
  )]

  # 降压药只使用 Instance 0 touchscreen 自报。这里故意保留 NA：
  #   1 = reported use; 0 = reported non-use; NA = baseline status unknown.
  # NA 不得在清洗阶段被改写成 0；最终 analytic cohort 中该变量可进入 MI，
  # 并在 MI 后与 SBP/DBP 一起被动派生 hypertension subgroup。
  dt[, antihypertensive_self_report_instance0 := as_medication_binary(hyper_medication_self_report_0)]
  dt[, antihypertensive_med_baseline := antihypertensive_self_report_instance0]
  dt[, glucose_lowering_verbal_instance0 := as_medication_binary(glucose_lowering_med_20003_0)]
  dt[, insulin_self_report_instance0 := as_medication_binary(insulin_self_report_0)]
  dt[, glucose_lowering_med_baseline := as_medication_binary(glucose_lowering_med_baseline)]

  keep <- c(
    "Participant ID", "baseline_time",
    "lipid_self_report_instance0", "lipid_gp_90d", "lipid_gp_365d",
    "lipid_lowering_med_baseline", "lipid_lowering_med_source",
    "antihypertensive_self_report_instance0", "antihypertensive_med_baseline",
    "glucose_lowering_verbal_instance0", "insulin_self_report_instance0",
    "glucose_lowering_med_baseline", "glucose_lowering_med_source"
  )
  out <- dt[, ..keep]
  assert_unique_id(out, "统一基线用药证据表")
  out[]
}

build_medication_source_audit <- function(medication_dt) {
  dt <- copy(as.data.table(medication_dt))
  assert_unique_id(dt, "用药来源审计输入")

  source_vars <- c(
    lipid_touchscreen = "lipid_self_report_instance0",
    lipid_gp_90d = "lipid_gp_90d",
    lipid_gp_365d = "lipid_gp_365d",
    antihypertensive_touchscreen = "antihypertensive_self_report_instance0",
    glucose_verbal_20003 = "glucose_lowering_verbal_instance0",
    insulin_touchscreen = "insulin_self_report_instance0"
  )
  source_rows <- rbindlist(lapply(names(source_vars), function(source_name) {
    value <- dt[[source_vars[[source_name]]]]
    medication <- sub("_(touchscreen|gp_90d|gp_365d|verbal_20003)$", "", source_name)
    data.table(
      audit_type = "source_count",
      medication = medication,
      evidence_source = source_name,
      category = c("positive", "negative", "missing"),
      n = c(sum(value == 1L, na.rm = TRUE), sum(value == 0L, na.rm = TRUE), sum(is.na(value)))
    )
  }))

  lipid_pattern <- dt[, fcase(
    lipid_self_report_instance0 == 1L & lipid_gp_365d == 1L, "touchscreen_and_gp_365d",
    lipid_self_report_instance0 == 1L, "touchscreen_only",
    lipid_gp_365d == 1L, "gp_365d_only",
    lipid_self_report_instance0 == 0L & lipid_gp_365d == 0L, "neither",
    default = "touchscreen_missing_no_gp"
  )]
  glucose_pattern <- dt[, fcase(
    glucose_lowering_verbal_instance0 == 1L & insulin_self_report_instance0 == 1L, "verbal_and_insulin",
    glucose_lowering_verbal_instance0 == 1L, "verbal_only",
    insulin_self_report_instance0 == 1L, "insulin_only",
    glucose_lowering_verbal_instance0 == 0L & insulin_self_report_instance0 == 0L, "neither",
    default = "one_or_more_sources_missing"
  )]
  intersection_rows <- rbindlist(list(
    data.table(category = lipid_pattern)[, .(n = .N), by = category][, `:=`(
      audit_type = "source_intersection", medication = "lipid", evidence_source = "touchscreen_x_gp_365d"
    )],
    data.table(category = glucose_pattern)[, .(n = .N), by = category][, `:=`(
      audit_type = "source_intersection", medication = "glucose", evidence_source = "verbal_x_insulin"
    )]
  ), use.names = TRUE, fill = TRUE)

  setcolorder(intersection_rows, names(source_rows))
  out <- rbindlist(list(source_rows, intersection_rows), use.names = TRUE, fill = TRUE)
  setorder(out, audit_type, medication, evidence_source, category)
  out[]
}

run_medication_unit_tests <- function() {
  baseline <- data.table(
    `Participant ID` = 1:2,
    baseline_time = as.IDate(c("2020-01-01", "2020-01-01"))
  )
  history <- data.table(
    `Participant ID(gp_scripts - eid)` = c(1L, 1L, 2L, 2L),
    Prescription_time = as.IDate("2020-01-01") - c(365L, 366L, 366L, 366L)
  )
  window <- build_window_summary(history, baseline)
  stopifnot(
    window[`Participant ID(gp_scripts - eid)` == 1L, gp_lipid_rx_in_365d] == 1L,
    window[`Participant ID(gp_scripts - eid)` == 2L, gp_lipid_rx_in_365d] == 0L,
    uniqueN(window[["Participant ID(gp_scripts - eid)"]]) == nrow(window)
  )

  boundary_history <- data.table(
    `Participant ID(gp_scripts - eid)` = c(1L,2L),
    Prescription_time = as.IDate("2020-01-01") - c(90L,91L))
  boundary_window <- build_window_summary(boundary_history,baseline)
  stopifnot(boundary_window[`Participant ID(gp_scripts - eid)`==1L,gp_lipid_rx_in_90d]==1L,
            boundary_window[`Participant ID(gp_scripts - eid)`==2L,gp_lipid_rx_in_90d]==0L)
  # Day 91 is outside 90d, but must remain IN the 365d main-analysis window.
  stopifnot(boundary_window[`Participant ID(gp_scripts - eid)`==2L,gp_lipid_rx_in_365d]==1L)
  # Day 365 is included, day 366 excluded (checked above); same-day (day 0) included.
  same_day <- data.table(`Participant ID(gp_scripts - eid)`=1L,
                         Prescription_time=as.IDate("2020-01-01"))
  stopifnot(build_window_summary(same_day, baseline)[
    `Participant ID(gp_scripts - eid)`==1L,gp_lipid_rx_in_365d]==1L)

  lipid_self <- data.table(
    `Participant ID` = 1:2,
    Lipid_medication_self_report_0 = c(0L, NA_integer_),
    Lipid_medication_self_report_1 = c(1L, 1L)
  )
  hyper_self <- data.table(`Participant ID` = 1:2, hyper_medication_self_report_0 = c(0L, 0L))
  glucose <- data.table(
    `Participant ID` = 1:2,
    glucose_lowering_med_20003_0 = c(0L, 0L),
    insulin_self_report_0 = c(0L, 0L),
    glucose_lowering_med_baseline = c(0L, 0L),
    glucose_lowering_med_source = c("no_evidence", "no_evidence")
  )
  gp <- data.table(`Participant ID` = 1:2, gp_lipid_rx_in_90d = c(0L, 0L), gp_lipid_rx_in_365d = c(1L, 0L))
  master <- build_medication_evidence_master(baseline, lipid_self, hyper_self, glucose, gp)
  stopifnot(
    master[`Participant ID` == 1L, lipid_lowering_med_baseline] == 1L,
    is.na(master[`Participant ID` == 2L, lipid_lowering_med_baseline])
  )
  # An individual with a GP prescription between days 91 and 365 is treated.
  stopifnot(master[`Participant ID` == 1L, lipid_lowering_med_source] == "gp_365d")
  audit <- build_medication_source_audit(master)
  stopifnot(any(audit$evidence_source == "touchscreen_x_gp_365d"))
  invisible(TRUE)
}

# High-specificity ART: established combination products or a complete regimen.
# Individual antivirals (including tenofovir/lamivudine) never qualify alone.
build_hiv_art <- function(path) {
 h<-names(fread(path,nrows=0));cols<-grep("^Treatment/medication code \\| Instance 0 \\| Array",h,value=TRUE)
 if(!length(cols))stop("No baseline medication arrays for ART audit")
 d<-fread(path,select=c("Participant ID",cols),colClasses="character")
 text<-tolower(do.call(paste,c(d[,..cols],sep=" | ")))
 combo<-grepl("\\b(atripla|trizivir|combivir|truvada|eviplera|complera|stribild|genvoya|triumeq|odefsey|biktarvy|symtuza|juluca|delstrigo|dovato)\\b",text,perl=TRUE)
 # Truvada/Descovy alone can be PrEP; require a third HIV treatment agent.
 third<-grepl("\\b(efavirenz|nevirapine|etravirine|rilpivirine|raltegravir|dolutegravir|elvitegravir|darunavir|atazanavir|lopinavir|indinavir|saquinavir|nelfinavir|fosamprenavir|kaletra)\\b",text,perl=TRUE)
 dual<-grepl("\\b(tenofovir|truvada|descovy)\\b",text,perl=TRUE)&grepl("\\b(emtricitabine|lamivudine|truvada|descovy)\\b",text,perl=TRUE)
 combo<-combo & (!grepl("truvada",text) | third | grepl("atripla|trizivir|triumeq|stribild|genvoya|eviplera|complera|odefsey|biktarvy|symtuza|juluca|delstrigo|dovato",text))
 # Combivir/Trizivir are HIV-specific products, not non-specific antivirals.
 known<-rowSums(!is.na(as.matrix(d[,..cols])))>0L
 data.table(`Participant ID`=d[[1]],hiv_specific_art=as.integer(combo | (dual & third)),
   hiv_art_ascertainment_known=as.integer(known))
}

