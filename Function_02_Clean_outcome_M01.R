# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_02_Clean_outcome_M01.R
# 任务：供02调用的日期/ICD/结局构造函数
# 输入：住院/死亡/自报事件
# 输出：基线疾病史和ASCVD/竞争死亡组件
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function_02_Clean_outcome_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   为 02_Clean_outcome_M01_v1.0.R 提供结局清理函数。
#
# 设计原则：
#   1. 文件路径、关键列名和代码集在主脚本各板块开头声明；
#   2. 函数只负责读取、校验、解析和派生；
#   3. 保留标准化完整 ICD-10，并按“范围纳入后精确排除”实现定义；
#   4. first occurrence 与详细住院 ICD-10 共同作为候选来源；
#   5. nonfatal 事件统一要求事件日期严格早于死亡日期；
#   6. 基线慢性共病仅使用 Instance 0 自报或基线前住院记录。
#   7. 本文件只负责结局/疾病史的原始定义；05 可基于这里冻结的 follow-up 和
#      competing-risk status 计算 MI auxiliary cumulative hazards，但不得重定义事件。
# ============================================================================== #

library(data.table)
library(stringr)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================================== #
# 0. 通用检查和格式化函数 ####
# ============================================================================== #

var_info_print <- function(data, vars) {
  data <- as.data.table(data)
  vars <- intersect(vars, names(data))

  if (length(vars) == 0) {
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

assert_file_exists <- function(path, label = "文件") {
  if (length(path) != 1 || is.na(path) || !file.exists(path)) {
    stop(label, "不存在：", path, call. = FALSE)
  }
  invisible(path)
}

assert_required_columns <- function(dt, cols, label = "数据表") {
  dt <- as.data.table(dt)
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  missing_cols <- setdiff(cols, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      label, "缺少关键列：", paste(missing_cols, collapse = ", "),
      "\n当前前 20 个列名：", paste(utils::head(names(dt), 20), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_no_name_conflicts <- function(dt) {
  nm <- names(dt)

  exact_dup <- unique(nm[duplicated(nm)])
  if (length(exact_dup) > 0) {
    stop("输出列名存在完全重复：", paste(exact_dup, collapse = ", "), call. = FALSE)
  }

  lower_nm <- tolower(nm)
  case_dup_key <- unique(lower_nm[duplicated(lower_nm)])
  if (length(case_dup_key) > 0) {
    detail <- vapply(case_dup_key, function(k) {
      paste(unique(nm[lower_nm == k]), collapse = " / ")
    }, character(1))
    stop("输出列名存在大小写冲突：", paste(detail, collapse = "; "), call. = FALSE)
  }

  invisible(TRUE)
}

assert_unique_id <- function(dt, label = "数据表") {
  dt <- as.data.table(dt)
  assert_required_columns(dt, "Participant ID", label)
  if (anyNA(dt[["Participant ID"]])) {
    stop(label, "存在缺失 Participant ID。", call. = FALSE)
  }
  duplicated_n <- dt[, .N, by = `Participant ID`][N > 1L, sum(N - 1L)]
  if (length(duplicated_n) == 0L || is.na(duplicated_n)) duplicated_n <- 0L
  if (duplicated_n > 0L) {
    stop(label, "存在 ", duplicated_n, " 条重复 Participant ID 记录。", call. = FALSE)
  }
  invisible(TRUE)
}

standardize_id_column <- function(dt) {
  dt <- as.data.table(dt)

  if ("Participant ID" %in% names(dt)) {
    set(dt,j="Participant ID",value=as.character(dt[["Participant ID"]]))
    return(dt)
  }

  id_candidates <- names(dt)[
    str_detect(names(dt), regex("participant.*id|participant.*eid|\\beid\\b", ignore_case = TRUE))
  ]

  if (length(id_candidates) == 0) {
    stop(
      "无法识别 Participant ID 列。当前前 20 个列名：",
      paste(utils::head(names(dt), 20), collapse = ", "),
      call. = FALSE
    )
  }

  setnames(dt, id_candidates[1], "Participant ID")
  set(dt,j="Participant ID",value=as.character(dt[["Participant ID"]]))
  dt
}

as_date_safely <- function(x) {
  # IDate 同时继承 Date，因此必须先判断 IDate。
  if (inherits(x, "IDate")) {
    return(as.Date(
      as.integer(x),
      origin = "1970-01-01"
    ))
  }
  
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  
  if (is.numeric(x)) {
    return(as.Date(
      x,
      origin = "1970-01-01"
    ))
  }
  
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NaN", "NULL")] <- NA_character_
  
  # UKB arrays are sparse with repeated dates: parse unique strings once.
  values <- unique(x)
  suppressWarnings(as.Date(values))[match(x, values)]
}

min_date_across_cols <- function(dt, cols) {
  dt <- as.data.table(dt)
  cols <- intersect(cols, names(dt))

  if (length(cols) == 0) {
    return(rep(as.Date(NA), nrow(dt)))
  }

  tmp <- dt[, ..cols]
  tmp[, (cols) := lapply(.SD, as_date_safely), .SDcols = cols]

  out <- do.call(pmin, c(tmp, list(na.rm = TRUE)))
  out_num <- suppressWarnings(as.numeric(out))
  out_num[is.infinite(out_num)] <- NA_real_
  as.Date(out_num, origin = "1970-01-01")
}

standardize_missing_values <- function(dt) {
  stopifnot(data.table::is.data.table(dt))

  missing_strings <- c("", "NA", "N/A", "NaN", "NULL", "Do not know", "Prefer not to answer")

  for (col in names(dt)) {
    x <- dt[[col]]
    if (is.character(x)) {
      x_trim <- trimws(x)
      x_trim[is.na(x_trim) | x_trim %in% missing_strings] <- NA_character_
      set(dt, j = col, value = x_trim)
    } else if (is.numeric(x)) {
      idx_nan <- which(is.nan(x))
      if (length(idx_nan) > 0) {
        set(dt, i = idx_nan, j = col, value = NA_real_)
      }
    }
  }

  dt[]
}

# ============================================================================== #
# 1. 编码提取函数 ####
# ============================================================================== #

extract_icd10_code <- function(x) {
  x <- toupper(trimws(as.character(x)))
  raw <- str_extract(x, "[A-Z][0-9]{2}(?:\\.?[0-9A-Z]{1,4})?")
  key <- gsub("\\.", "", raw)
  ifelse(
    !is.na(key) & nchar(key) > 3L,
    paste0(substr(key, 1L, 3L), ".", substr(key, 4L, nchar(key))),
    key
  )
}

extract_icd3_code <- function(x) {
  code <- extract_icd10_code(x)
  ifelse(!is.na(code), substr(code, 1, 3), NA_character_)
}

icd10_key <- function(x) {
  gsub("\\.", "", extract_icd10_code(x))
}

match_icd10_definition <- function(x,
                                   include_families = character(),
                                   include_exact = character(),
                                   exclude_families = character(),
                                   exclude_exact = character()) {
  key <- icd10_key(x)
  family <- ifelse(!is.na(key), substr(key, 1L, 3L), NA_character_)
  include_families <- unique(substr(icd10_key(include_families), 1L, 3L))
  include_exact <- unique(icd10_key(include_exact))
  exclude_families <- unique(substr(icd10_key(exclude_families), 1L, 3L))
  exclude_exact <- unique(icd10_key(exclude_exact))

  included <- (!is.na(family) & family %in% include_families) |
    (!is.na(key) & key %in% include_exact)
  excluded <- (!is.na(family) & family %in% exclude_families) |
    (!is.na(key) & key %in% exclude_exact)
  included & !excluded
}

extract_opcs4_code <- function(x) {
  x <- toupper(trimws(as.character(x)))
  str_extract(x, "[A-Z][0-9]{2}(?:\\.[0-9A-Z]+)?")
}

extract_opcs3_code <- function(x) {
  code <- extract_opcs4_code(x)
  ifelse(!is.na(code), substr(code, 1, 3), NA_character_)
}

icd_range <- function(prefix, numbers) {
  paste0(prefix, sprintf("%02d", numbers))
}

opcs_range <- function(prefix, numbers) {
  paste0(prefix, sprintf("%02d", numbers))
}

parse_array_index <- function(x) {
  out <- str_match(as.character(x), regex("Array\\s*([0-9]+)", ignore_case = TRUE))[, 2]
  suppressWarnings(as.integer(out))
}

# ============================================================================== #
# 2. 基线日期、失访和死亡数据读取 ####
# ============================================================================== #

read_assessment_and_censoring <- function(assessment_file,
                                          assessment_columns,
                                          population_file,
                                          population_columns,
                                          admin_censor_date = as.Date("2022-11-30")) {
  assert_file_exists(assessment_file, "基线评估文件")
  assert_file_exists(population_file, "失访信息文件")

  assessment_raw <- fread(assessment_file)
  assessment_raw <- standardize_id_column(assessment_raw)
  assert_required_columns(
    assessment_raw,
    c(assessment_columns$assessment_date, assessment_columns$assessment_centre),
    "基线评估文件"
  )

  assessment <- assessment_raw[, .(
    `Participant ID`,
    InCohort_time = as_date_safely(get(assessment_columns$assessment_date)),
    Region = fifelse(
      get(assessment_columns$assessment_centre) %in% c(
        "Barts", "Birmingham", "Bristol", "Bury", "Croydon", "Hounslow",
        "Leeds", "Liverpool", "Manchester", "Middlesborough", "Newcastle",
        "Nottingham", "Oxford", "Reading", "Stockport (pilot)", "Stoke",
        "Sheffield"
      ),
      "England",
      fifelse(
        get(assessment_columns$assessment_centre) %in% c("Cardiff", "Swansea", "Wrexham"),
        "Wales",
        "Scotland"
      )
    )
  )]

  pop <- fread(population_file)
  pop <- standardize_id_column(pop)
  assert_required_columns(
    pop,
    c(population_columns$loss_date, population_columns$loss_reason),
    "失访信息文件"
  )

  censoring <- pop[, .(
    `Participant ID`,
    Loss_to_followup_time = as_date_safely(get(population_columns$loss_date)),
    Censored_reason = as.character(get(population_columns$loss_reason))
  )]

  out <- merge(assessment, censoring, by = "Participant ID", all.x = TRUE)
  out[, Admin_censor_time := admin_censor_date]
  out[, Censored_time := pmin(Loss_to_followup_time, Admin_censor_time, na.rm = TRUE)]
  out[is.infinite(as.numeric(Censored_time)), Censored_time := Admin_censor_time]
  out[, Censored := fifelse(!is.na(Loss_to_followup_time) & Loss_to_followup_time <= Admin_censor_time, 1L, 0L)]

  out[]
}

read_death_records <- function(death_file, death_columns) {
  assert_file_exists(death_file, "死亡日期文件")

  death_raw <- fread(death_file)
  death_raw <- standardize_id_column(death_raw)

  date_cols <- names(death_raw)[
    str_detect(names(death_raw), regex(death_columns$date_pattern, ignore_case = TRUE))
  ]

  if (length(date_cols) == 0) {
    stop("死亡日期文件中未找到日期列：", death_columns$date_pattern, call. = FALSE)
  }

  death_raw[, (date_cols) := lapply(.SD, as_date_safely), .SDcols = date_cols]

  death <- death_raw[, .(`Participant ID`)]
  death[, Death_time := min_date_across_cols(death_raw, date_cols)]
  death[, Death := fifelse(!is.na(Death_time), 1L, 0L)]
  death[]
}

read_death_causes <- function(death_cause_file,
                              death_cause_columns,
                              fatal_cvd_include_families,
                              fatal_cvd_exclude_families = character(),
                              fatal_cvd_exclude_exact = character()) {
  assert_file_exists(death_cause_file, "死因文件")

  raw <- fread(death_cause_file)
  raw <- standardize_id_column(raw)
  assert_required_columns(
    raw,
    c(death_cause_columns$classification, death_cause_columns$cause),
    "死因文件"
  )

  raw <- raw[get(death_cause_columns$classification) == death_cause_columns$primary_label]

  if (nrow(raw) == 0) {
    warning("死因文件中未找到 primary cause 记录。", call. = FALSE)
    return(data.table(
      `Participant ID` = character(),
      primary_death_cause_available = integer(),
      fatal_cvd_death_cause = integer(),
      fatal_cvd_death_code = character()
    ))
  }

  raw[, deathcause_full := extract_icd10_code(get(death_cause_columns$cause))]
  raw[, fatal_cvd_death_cause := as.integer(match_icd10_definition(
    deathcause_full,
    include_families = fatal_cvd_include_families,
    exclude_families = fatal_cvd_exclude_families,
    exclude_exact = fatal_cvd_exclude_exact
  ))]

  out <- raw[, .(
    primary_death_cause_available = 1L,
    fatal_cvd_death_cause = as.integer(any(fatal_cvd_death_cause == 1L, na.rm = TRUE)),
    fatal_cvd_death_code = paste(unique(na.omit(deathcause_full[fatal_cvd_death_cause == 1L])), collapse = ";")
  ), by = `Participant ID`]

  out[fatal_cvd_death_code == "", fatal_cvd_death_code := NA_character_]
  out[]
}

derive_death_classification <- function(dt) {
  dt <- as.data.table(dt)

  if (!"fatal_cvd_death_cause" %in% names(dt)) {
    dt[, fatal_cvd_death_cause := NA_integer_]
  }
  if (!"primary_death_cause_available" %in% names(dt)) {
    dt[, primary_death_cause_available := NA_integer_]
  }
  if (!"fatal_cvd_death_code" %in% names(dt)) {
    dt[, fatal_cvd_death_code := NA_character_]
  }

  dt[is.na(Death), Death := 0L]
  dt[, fatal_cvd_death := fifelse(Death == 1L & fatal_cvd_death_cause == 1L, 1L, 0L)]
  dt[, fatal_cvd_death_time := fifelse(fatal_cvd_death == 1L, Death_time, as.Date(NA))]
  dt[, fatal_cvd_death_source := fifelse(fatal_cvd_death == 1L, "death_registry", NA_character_)]

  # competing_death 仅标记已完成死因判定且不属于 fatal CVD 的死亡；
  # 主要死因缺失者保持 NA，并在主脚本中单独 QA，不擅自改变冻结分类。
  dt[, competing_death := fifelse(Death == 1L & fatal_cvd_death != 1L, 1L, 0L)]
  dt[, competing_death_time := fifelse(competing_death == 1L, Death_time, as.Date(NA))]

  dt[]
}

# ============================================================================== #
# 3. first occurrence 读取和 ICD 组分派生 ####
# ============================================================================== #

read_first_occurrence_file <- function(file_path, component_label = basename(file_path), relevant_codes = NULL) {
  assert_file_exists(file_path, paste0(component_label, " first occurrence 文件"))

  header <- names(fread(file_path,nrows=0L))
  id_cols <- header[grepl("^Participant ID",header)]
  date_cols <- header[
    str_detect(header, regex("^Date\\s+[A-Z][0-9]{2}\\s+first reported", ignore_case = TRUE))
  ]

  if (!is.null(relevant_codes)) {
    codes <- toupper(str_match(date_cols,regex("^Date\\s+([A-Z][0-9]{2})\\s+first reported",ignore_case=TRUE))[,2])
    date_cols <- date_cols[codes %in% toupper(relevant_codes)]
  }
  if (length(id_cols)!=1L) stop("First occurrence参与者ID列不唯一：",file_path)
  message("读取first occurrence：",component_label,"；日期列：",length(date_cols))
  raw <- fread(file_path,select=c(id_cols,date_cols))
  raw <- standardize_id_column(raw)
  if (length(date_cols) == 0) {
    if (is.null(relevant_codes)) warning("未在 first occurrence 文件中识别到日期列：", component_label, call. = FALSE)
    return(raw[, .(`Participant ID`)])
  }

  code_from_col <- str_match(
    date_cols,
    regex("^Date\\s+([A-Z][0-9]{2})\\s+first reported", ignore_case = TRUE)
  )[, 2]
  code_from_col <- toupper(code_from_col)

  out <- raw[, .(`Participant ID`)]

  for (code in unique(code_from_col)) {
    cols <- date_cols[code_from_col == code]
    raw[, (cols) := lapply(.SD, as_date_safely), .SDcols = cols]
    # UKB first occurrence 特殊占位日期不能被误作真实的基线前事件日期。
    raw[, (cols) := lapply(.SD, function(x) {
      x[x %in% as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))] <- as.Date(NA)
      x
    }), .SDcols = cols]
    out[, paste0(code, "_time") := min_date_across_cols(raw, cols)]
    out[, (code) := fifelse(!is.na(get(paste0(code, "_time"))), 1L, 0L)]
  }

  setcolorder(out, c("Participant ID", sort(setdiff(names(out), "Participant ID"))))
  out[]
}

read_first_occurrence_files <- function(first_occurrence_files, relevant_codes = NULL) {
  if (length(first_occurrence_files) == 0) {
    stop("first_occurrence_files 为空。", call. = FALSE)
  }

  first_occurrence_list <- vector("list", length(first_occurrence_files))
  source_labels <- names(first_occurrence_files)
  for (k in seq_along(first_occurrence_files)) {
    nm <- if (!is.null(source_labels) && nzchar(source_labels[[k]]))
      source_labels[[k]] else basename(first_occurrence_files[[k]])
    first_occurrence_list[[k]] <- read_first_occurrence_file(
      first_occurrence_files[[k]], component_label = nm,
      relevant_codes = relevant_codes)
    if (exists("m01_checkpoint", mode = "function", inherits = TRUE))
      m01_checkpoint("02_first_occurrence_files", k,
                     length(first_occurrence_files), nm)
  }

  Reduce(function(x, y) merge(x, y, by = "Participant ID", all = TRUE), first_occurrence_list)
}

derive_icd_component <- function(dt, codes, component_name, source_label = "first_occurrence") {
  dt <- as.data.table(dt)
  time_cols <- intersect(paste0(codes, "_time"), names(dt))

  event_name <- component_name
  time_name <- paste0(component_name, "_time")

  if (length(time_cols) == 0) {
    warning("未找到该 ICD 组分对应的 first occurrence 列：", component_name, call. = FALSE)
    dt[, (time_name) := as.Date(NA)]
    dt[, (event_name) := 0L]
    dt[, paste0(component_name, "_code") := NA_character_]
    dt[, paste0(component_name, "_source") := NA_character_]
    return(dt[])
  }

  dt[, (time_name) := min_date_across_cols(.SD, time_cols), .SDcols = time_cols]
  dt[, (event_name) := fifelse(!is.na(get(time_name)), 1L, 0L)]
  candidate_codes <- sub("_time$", "", time_cols)
  date_matrix <- do.call(cbind, lapply(dt[, ..time_cols], function(x) as.numeric(as_date_safely(x))))
  if (!is.matrix(date_matrix)) date_matrix <- matrix(date_matrix, ncol = length(time_cols))
  date_matrix[is.na(date_matrix)] <- Inf
  earliest_index <- max.col(-date_matrix, ties.method = "first")
  earliest_code <- candidate_codes[earliest_index]
  earliest_code[dt[[event_name]] != 1L] <- NA_character_
  dt[, paste0(component_name, "_code") := earliest_code]
  dt[, paste0(component_name, "_source") := fifelse(get(event_name) == 1L, source_label, NA_character_)]
  dt[]
}

derive_nonfatal_icd_component <- function(dt,
                                          codes,
                                          component_name,
                                          death_time_var = "Death_time") {
  dt <- as.data.table(dt)
  time_cols <- intersect(paste0(codes, "_time"), names(dt))

  event_name <- component_name
  time_name <- paste0(component_name, "_time")

  if (length(time_cols) == 0) {
    warning("未找到该 nonfatal ICD 组分对应的 first occurrence 列：", component_name, call. = FALSE)
    dt[, (time_name) := as.Date(NA)]
    dt[, (event_name) := 0L]
    return(dt[])
  }

  candidate_time <- min_date_across_cols(dt, time_cols)
  death_time <- if (death_time_var %in% names(dt)) dt[[death_time_var]] else rep(as.Date(NA), nrow(dt))
  death_time <- as_date_safely(death_time)

  # first occurrence 与死亡同日时，不直接判为 nonfatal；
  # 该日期应优先由 fatal death 或 competing death 处理。
  nonfatal_flag <- !is.na(candidate_time) & (is.na(death_time) | candidate_time < death_time)
  event_time <- candidate_time
  event_time[!nonfatal_flag] <- as.Date(NA)

  dt[, (time_name) := event_time]
  dt[, (event_name) := fifelse(!is.na(get(time_name)), 1L, 0L)]

  same_day_n <- sum(!is.na(candidate_time) & !is.na(death_time) & candidate_time == death_time)
  after_death_n <- sum(!is.na(candidate_time) & !is.na(death_time) & candidate_time > death_time)
  if (same_day_n > 0 || after_death_n > 0) {
    message(
      component_name,
      "：first occurrence 同日死亡 ", same_day_n,
      " 例；晚于死亡 ", after_death_n,
      " 例；均未计入 nonfatal 组分。"
    )
  }

  dt[]
}

read_hospital_icd10_events <- function(hospital_file,
                                       code_column = "Diagnoses - ICD10",
                                       date_prefix = "Date of first in-patient diagnosis - ICD10 | Array ",
                                       relevant_families = character(), raw_input = NULL, date_parts = NULL) {
  assert_file_exists(hospital_file, "住院 ICD-10 文件")
  if (is.null(raw_input)) {
    header <- fread(hospital_file, nrows = 0L)
    date_cols <- names(header)[startsWith(names(header), date_prefix)]
    raw <- fread(hospital_file, select = c("Participant ID", code_column, date_cols))
  } else raw <- copy(raw_input)
  # Use both 41280 files; overlapping arrays must agree. IDs are matched, not rows.
  if (!is.null(date_parts)) for (part in date_parts) {
    ix <- match(as.character(raw[["Participant ID"]]),as.character(part[["Participant ID"]]))
    for (col in names(part)[startsWith(names(part),date_prefix)]) {
      value <- part[[col]][ix]
      if (col %in% names(raw)) {
        old <- as_date_safely(raw[[col]]); new <- as_date_safely(value)
        if (any(!is.na(old) & !is.na(new) & old!=new)) stop("Conflicting 41280 array: ",col)
        fill <- which(is.na(old)); set(raw,i=fill,j=col,value=value[fill])
      } else set(raw,j=col,value=value)
    }
  }
  date_cols <- names(raw)[startsWith(names(raw),date_prefix)]
  if (!length(date_cols)) stop("住院 ICD-10 文件未找到41280日期数组")
  raw <- standardize_id_column(raw)
  assert_required_columns(raw, c("Participant ID", code_column), "住院 ICD-10 文件")
  if (length(relevant_families) > 0L) {
    relevant_families <- unique(substr(icd10_key(relevant_families), 1L, 3L))
    family_pattern <- paste0(
      "(^|\\|)\\s*(?:", paste(relevant_families, collapse = "|"), ")",
      "(?:\\.|[0-9A-Z]|\\s|$)"
    )
    raw <- raw[grepl(family_pattern, get(code_column), ignore.case = TRUE)]
  }
  raw[, row_index__ := .I]

  split_codes <- strsplit(as.character(raw[[code_column]]), "\\|", fixed = FALSE)
  lens <- lengths(split_codes)
  keep_rows <- lens > 0L
  events <- data.table(
    row_index__ = rep(raw$row_index__[keep_rows], lens[keep_rows]),
    `Participant ID` = rep(raw[["Participant ID"]][keep_rows], lens[keep_rows]),
    array_index = sequence(lens[keep_rows]) - 1L,
    icd10_raw = unlist(split_codes[keep_rows], use.names = FALSE)
  )
  events[, icd10_code := extract_icd10_code(icd10_raw)]
  events[, icd10_family := extract_icd3_code(icd10_code)]
  if (length(relevant_families) > 0L) {
    events <- events[icd10_family %in% relevant_families]
  }

  events[, event_date := as.Date(NA)]
  date_index <- parse_array_index(date_cols)
  if (length(setdiff(events$array_index,date_index))) stop("住院诊断缺少对应41280数组")
  # Index events ONCE by the actual array index: repeated which() scans cost
  # O(n_events * n_arrays) on large UKB hospital exports.
  # All date mappings remain exactly the same as the original loop.
  event_rows_by_array <- split(seq_len(nrow(events)), events$array_index)
  n_array_groups <- length(event_rows_by_array)
  progress_every <- max(1L, n_array_groups %/% 10L)
  array_done <- 0L
  for (key in names(event_rows_by_array)) {
    idx <- as.integer(key)
    date_col <- date_cols[match(idx, date_index)]
    target <- event_rows_by_array[[key]]
    data.table::set(events, i = target, j = "event_date",
      value = as_date_safely(raw[[date_col]][events$row_index__[target]]))
    array_done <- array_done + 1L
    if ((array_done %% progress_every == 0L || array_done == n_array_groups) &&
        exists("m01_checkpoint", mode = "function", inherits = TRUE))
      m01_checkpoint("02_hospital_ICD_date_arrays", array_done, n_array_groups)
  }

  events <- events[!is.na(icd10_code), .(
    `Participant ID`, array_index, icd10_code, icd10_family, event_date
  )]
  setkey(events, `Participant ID`)
  events[]
}

summarize_hospital_icd_component <- function(ids,
                                             hospital_events,
                                             component_name,
                                             include_families = character(),
                                             include_exact = character(),
                                             exclude_families = character(),
                                             exclude_exact = character()) {
  ids <- unique(standardize_id_column(copy(ids))[, .(`Participant ID`)])
  events <- standardize_id_column(copy(hospital_events))
  matched <- events[match_icd10_definition(
    icd10_code,
    include_families = include_families,
    include_exact = include_exact,
    exclude_families = exclude_families,
    exclude_exact = exclude_exact
  ) & !is.na(event_date)]
  if (nrow(matched) > 0L) {
    setorder(matched, `Participant ID`, event_date, icd10_code)
    summary <- matched[, .(
      event_time__ = event_date[1L],
      event_code__ = paste(unique(icd10_code[event_date == event_date[1L]]), collapse = ";")
    ), by = `Participant ID`]
  } else {
    summary <- data.table(`Participant ID` = ids[["Participant ID"]][0], event_time__ = as.Date(character()), event_code__ = character())
  }
  out <- merge(ids, summary, by = "Participant ID", all.x = TRUE)
  setnames(out, c("event_time__", "event_code__"), c(paste0(component_name, "_time"), paste0(component_name, "_code")))
  out[, (component_name) := as.integer(!is.na(get(paste0(component_name, "_time"))))]
  out[, paste0(component_name, "_source") := fifelse(get(component_name) == 1L, "hospital", NA_character_)]
  out[]
}

combine_source_components <- function(dt, source_components, output_name) {
  dt <- as.data.table(dt)
  source_components <- source_components[source_components %in% names(dt)]
  if (length(source_components) == 0L) stop("没有可合并的来源组分：", output_name, call. = FALSE)
  time_cols <- paste0(source_components, "_time")
  dt[, paste0(output_name, "_time") := min_date_across_cols(.SD, time_cols), .SDcols = time_cols]
  dt[, (output_name) := as.integer(!is.na(get(paste0(output_name, "_time"))))]

  collect_values <- function(suffix) {
    result <- rep(NA_character_, nrow(dt))
    for (comp in source_components) {
      value_col <- paste0(comp, suffix)
      if (!value_col %in% names(dt)) next
      value <- as.character(dt[[value_col]])
      use <- !is.na(dt[[comp]]) & dt[[comp]] == 1L & !is.na(value) & nzchar(value)
      first <- use & is.na(result)
      result[first] <- value[first]
      append <- use & !is.na(result) & result != value
      result[append] <- paste0(result[append], ";", value[append])
    }
    result
  }
  dt[, paste0(output_name, "_source") := collect_values("_source")]
  dt[, paste0(output_name, "_code") := collect_values("_code")]
  dt[]
}

derive_multisource_icd_component <- function(dt,
                                             hospital_events,
                                             first_occurrence_codes,
                                             component_name,
                                             include_families,
                                             include_exact = character(),
                                             exclude_families = character(),
                                             exclude_exact = character()) {
  dt <- standardize_id_column(dt)
  fo_name <- paste0(component_name, "__first_occurrence")
  hosp_name <- paste0(component_name, "__hospital")
  dt <- derive_icd_component(dt, first_occurrence_codes, fo_name, "first_occurrence")
  hospital <- summarize_hospital_icd_component(
    dt[, .(`Participant ID`)], hospital_events, hosp_name,
    include_families, include_exact, exclude_families, exclude_exact
  )
  dt <- merge(dt, hospital, by = "Participant ID", all.x = TRUE)
  dt <- combine_source_components(dt, c(fo_name, hosp_name), component_name)
  dt[]
}

apply_nonfatal_rule <- function(dt, component_name, death_time_var = "Death_time") {
  dt <- as.data.table(dt)
  time_name <- paste0(component_name, "_time")
  code_name <- paste0(component_name, "_code")
  source_name <- paste0(component_name, "_source")
  assert_required_columns(dt, c(component_name, time_name), component_name)
  death_time <- if (death_time_var %in% names(dt)) as_date_safely(dt[[death_time_var]]) else rep(as.Date(NA), nrow(dt))
  event_time <- as_date_safely(dt[[time_name]])
  invalid <- !is.na(event_time) & !is.na(death_time) & event_time >= death_time
  same_day_n <- sum(!is.na(event_time) & !is.na(death_time) & event_time == death_time)
  after_death_n <- sum(!is.na(event_time) & !is.na(death_time) & event_time > death_time)
  if (any(invalid)) {
    dt[invalid, (component_name) := 0L]
    dt[invalid, (time_name) := as.Date(NA)]
    if (code_name %in% names(dt)) dt[invalid, (code_name) := NA_character_]
    if (source_name %in% names(dt)) dt[invalid, (source_name) := NA_character_]
  }
  message(component_name, "：同日死亡排除 ", same_day_n, " 例；死亡后排除 ", after_death_n, " 例。")
  dt[]
}

read_baseline_self_report_flags <- function(self_report_file,
                                            definitions,
                                            column_prefix = "Non-cancer illness code, self-reported | Instance 0",
                                            raw_input = NULL) {
  assert_file_exists(self_report_file, "Field 20002 自报文件")
  if (is.null(raw_input)) {
    header <- fread(self_report_file,nrows=0L)
    value_cols <- names(header)[startsWith(names(header),column_prefix)]
    raw <- fread(self_report_file,select=c("Participant ID",value_cols))
  } else {
    raw <- copy(raw_input)
    value_cols <- names(raw)[startsWith(names(raw),column_prefix)]
  }
  if (!length(value_cols)) stop("Field 20002文件中未找到Instance 0数组")
  raw <- standardize_id_column(raw)
  out <- raw[, .(`Participant ID`)]
  normalize_value <- function(x) tolower(trimws(as.character(x)))
  for (component_name in names(definitions)) {
    accepted <- unique(normalize_value(definitions[[component_name]]))
    flags <- Reduce(`|`, lapply(raw[, ..value_cols], function(x) normalize_value(x) %in% accepted))
    out[, paste0(component_name, "_self_report") := as.integer(flags)]
  }
  assert_unique_id(out, "Field 20002 基线自报派生表")
  out[]
}

derive_baseline_hospital_comorbidity <- function(ids_with_baseline,
                                                 hospital_events,
                                                 component_name,
                                                 include_families = character(),
                                                 include_exact = character(),
                                                 exclude_families = character(),
                                                 exclude_exact = character(),
                                                 baseline_time_var = "InCohort_time") {
  ids <- standardize_id_column(copy(ids_with_baseline))[, .(`Participant ID`, baseline_time__ = as_date_safely(get(baseline_time_var)))]
  events <- standardize_id_column(copy(hospital_events))[match_icd10_definition(
    icd10_code, include_families, include_exact, exclude_families, exclude_exact
  ) & !is.na(event_date)]
  events <- merge(events, ids, by = "Participant ID", all = FALSE)
  events <- events[!is.na(baseline_time__) & event_date <= baseline_time__]
  if (nrow(events) > 0L) {
    setorder(events, `Participant ID`, event_date, icd10_code)
    summary <- events[, .(
      hospital_time__ = event_date[1L],
      hospital_code__ = paste(unique(icd10_code[event_date == event_date[1L]]), collapse = ";")
    ), by = `Participant ID`]
  } else {
    summary <- data.table(`Participant ID` = ids[["Participant ID"]][0], hospital_time__ = as.Date(character()), hospital_code__ = character())
  }
  out <- merge(ids[, .(`Participant ID`)], summary, by = "Participant ID", all.x = TRUE)
  setnames(out, c("hospital_time__", "hospital_code__"), c(paste0(component_name, "_hospital_time"), paste0(component_name, "_hospital_code")))
  out[, paste0(component_name, "_hospital") := as.integer(!is.na(get(paste0(component_name, "_hospital_time"))))]
  out[]
}

combine_baseline_comorbidity_sources <- function(dt, component_name) {
  dt <- as.data.table(dt)
  self_var <- paste0(component_name, "_self_report")
  hospital_var <- paste0(component_name, "_hospital")
  hospital_time_var <- paste0(component_name, "_hospital_time")
  dt[is.na(get(self_var)), (self_var) := 0L]
  dt[is.na(get(hospital_var)), (hospital_var) := 0L]
  dt[, (component_name) := as.integer(get(self_var) == 1L | get(hospital_var) == 1L)]
  dt[, paste0(component_name, "_time") := as_date_safely(get(hospital_time_var))]
  dt[, paste0(component_name, "_source") := fcase(
    get(self_var) == 1L & get(hospital_var) == 1L, "self_report;hospital",
    get(self_var) == 1L, "self_report",
    get(hospital_var) == 1L, "hospital",
    default = NA_character_
  )]
  dt[]
}

derive_ckd_icd10_n18_component <- function(dt,
                                           ckd_code = "N18",
                                           baseline_time_var = "InCohort_time") {
  dt <- as.data.table(dt)

  # 当前 CKD 仅使用 ICD-10 N18 first occurrence 判定，不包含 eGFR 标准。
  dt <- derive_icd_component(dt, ckd_code, "CKD_ICD10_N18")
  dt <- derive_prevalent_component(
    dt,
    "CKD_ICD10_N18",
    baseline_time_var = baseline_time_var
  )

  # 保留旧版 CKD 变量名作为兼容别名；含义同 CKD_ICD10_N18。
  dt[, `:=`(
    CKD = CKD_ICD10_N18,
    CKD_time = CKD_ICD10_N18_time
  )]

  dt[]
}

# ============================================================================== #
# 4. OPCS-4 手术记录解析 ####
# ============================================================================== #

build_operation_date_map <- function(operation_tables, date_prefix) {
  parts <- list()

  for (i in seq_along(operation_tables)) {
    raw <- operation_tables[[i]]
    date_cols <- names(raw)[startsWith(names(raw), date_prefix)]

    if (length(date_cols) == 0) {
      next
    }

    tmp <- melt(
      raw[, c("Participant ID", date_cols), with = FALSE],
      id.vars = "Participant ID",
      variable.name = "date_var",
      value.name = "op_date",
      variable.factor = FALSE
    )

    tmp[, array_index := parse_array_index(date_var)]
    tmp[, op_date := as_date_safely(op_date)]
    tmp <- tmp[!is.na(array_index) & !is.na(op_date), .(`Participant ID`, array_index, op_date)]

    parts[[length(parts) + 1]] <- tmp
  }

  if (length(parts) == 0) {
    return(data.table(`Participant ID` = character(), array_index = integer(), op_date = as.Date(character())))
  }

  date_map <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  date_map[, .(op_date = min(op_date, na.rm = TRUE)), by = .(`Participant ID`, array_index)]
}

derive_coronary_revascularization <- function(operation_files,
                                              operation_columns,
                                              opcs_codes,
                                              ids,
                                              component_name = "coronary_revascularization") {
  ids <- as.data.table(ids)[, .(`Participant ID`)]
  operation_files <- operation_files[!is.na(operation_files) & file.exists(operation_files)]

  if (length(operation_files) == 0) {
    out <- copy(ids)
    out[, paste0(component_name, "_time") := as.Date(NA)]
    out[, (component_name) := 0L]
    out[, paste0(component_name, "_code") := NA_character_]
    out[, paste0(component_name, "_source") := NA_character_]
    warning("未找到手术文件；冠脉血运重建组分设为 0。", call. = FALSE)
    return(out[])
  }

  operation_tables <- lapply(operation_files, function(file_path) {
    raw <- fread(file_path)
    standardize_id_column(raw)
  })

  code_columns <- operation_columns$code_columns
  date_prefix_by_code <- operation_columns$date_prefix_by_code
  code_pattern <- paste0("(^|\\|)\\s*(", paste(opcs_codes, collapse = "|"), ")(\\.|\\s|$)")

  long_parts <- list()

  for (code_col in code_columns) {
    source_index <- which(vapply(operation_tables, function(x) code_col %in% names(x), logical(1)))[1]
    if (is.na(source_index)) {
      next
    }

    date_prefix <- date_prefix_by_code[[code_col]]
    if (is.null(date_prefix) || is.na(date_prefix)) {
      next
    }

    raw_codes <- operation_tables[[source_index]][, .(
      `Participant ID`,
      op_codes = as.character(get(code_col))
    )]

    raw_codes <- raw_codes[
      !is.na(op_codes) &
        nzchar(trimws(op_codes)) &
        grepl(code_pattern, op_codes, ignore.case = TRUE)
    ]

    if (nrow(raw_codes) == 0) {
      next
    }

    split_codes <- strsplit(raw_codes$op_codes, "\\|")
    lens <- lengths(split_codes)

    if (sum(lens) == 0) {
      next
    }

    codes_long <- data.table(
      `Participant ID` = rep(raw_codes[["Participant ID"]], lens),
      array_index = sequence(lens) - 1L,
      opcs_raw = unlist(split_codes, use.names = FALSE)
    )

    codes_long[, opcs4 := extract_opcs4_code(opcs_raw)]
    codes_long[, opcs3 := extract_opcs3_code(opcs_raw)]
    codes_long <- codes_long[opcs3 %in% opcs_codes]

    if (nrow(codes_long) == 0) {
      next
    }

    date_map <- build_operation_date_map(operation_tables, date_prefix)
    if (nrow(date_map) == 0) {
      next
    }

    long <- merge(
      codes_long[, .(`Participant ID`, array_index, opcs4, opcs3)],
      date_map,
      by = c("Participant ID", "array_index"),
      all = FALSE
    )

    long_parts[[length(long_parts) + 1]] <- long[!is.na(op_date)]
  }

  out <- copy(ids)
  time_name <- paste0(component_name, "_time")

  if (length(long_parts) == 0) {
    out[, (time_name) := as.Date(NA)]
    out[, (component_name) := 0L]
    out[, paste0(component_name, "_code") := NA_character_]
    out[, paste0(component_name, "_source") := NA_character_]
    return(out[])
  }

  events <- unique(rbindlist(long_parts, use.names = TRUE, fill = TRUE))
  setorder(events, `Participant ID`, op_date, opcs4)
  tmp <- events[, .(
    event_time = op_date[1L],
    event_code = paste(unique(opcs4[op_date == op_date[1L]]), collapse = ";")
  ), by = `Participant ID`]
  tmp[is.infinite(as.numeric(event_time)), event_time := as.Date(NA)]

  out <- merge(out, tmp, by = "Participant ID", all.x = TRUE)
  setnames(out, c("event_time", "event_code"), c(time_name, paste0(component_name, "_code")))
  out[, (component_name) := fifelse(!is.na(get(time_name)), 1L, 0L)]
  out[, paste0(component_name, "_source") := fifelse(get(component_name) == 1L, "opcs", NA_character_)]

  message(component_name, "：识别到 ", sum(out[[component_name]] == 1L, na.rm = TRUE), " 人。")
  out[]
}

# ============================================================================== #
# 5. 复合终点和随访变量 ####
# ============================================================================== #

combine_components <- function(dt, component_time_vars, output_name) {
  dt <- as.data.table(dt)
  component_time_vars <- intersect(component_time_vars, names(dt))
  time_name <- paste0(output_name, "_time")

  dt[, (time_name) := min_date_across_cols(.SD, component_time_vars), .SDcols = component_time_vars]
  dt[, (output_name) := fifelse(!is.na(get(time_name)), 1L, 0L)]
  dt[]
}

derive_prevalent_component <- function(dt, component_name, baseline_time_var = "InCohort_time") {
  dt <- as.data.table(dt)
  time_var <- paste0(component_name, "_time")
  prevalent_var <- paste0("prevalent_", component_name)

  if (!time_var %in% names(dt)) {
    stop("缺少组分时间变量：", time_var, call. = FALSE)
  }
  if (!baseline_time_var %in% names(dt)) {
    stop("缺少基线时间变量：", baseline_time_var, call. = FALSE)
  }

  dt[, (prevalent_var) := fifelse(
    !is.na(get(time_var)) &
      !is.na(get(baseline_time_var)) &
      get(time_var) <= get(baseline_time_var),
    1L,
    0L
  )]
  dt[]
}

make_analysis_endpoint <- function(dt,
                                   endpoint_name,
                                   event_time_var,
                                   competing_time_var = "competing_death_time",
                                   baseline_time_var = "InCohort_time",
                                   censor_time_var = "Censored_time") {
  dt <- as.data.table(dt)

  required <- c(event_time_var, competing_time_var, baseline_time_var, censor_time_var)
  missing_required <- setdiff(required, names(dt))
  if (length(missing_required) > 0) {
    stop("生成分析结局 ", endpoint_name, " 时缺少变量：", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  event_date_var <- paste0(endpoint_name, "_analysis_date")
  status_var <- paste0(endpoint_name, "_status")
  event_var <- paste0(endpoint_name, "_event")
  comp_var <- paste0(endpoint_name, "_competing_event")
  fup_days_var <- paste0(endpoint_name, "_followup_days")
  fup_years_var <- paste0(endpoint_name, "_followup_years")

  event_after_baseline <- paste0(endpoint_name, "__target_event_date")
  comp_after_baseline <- paste0(endpoint_name, "__competing_event_date")

  dt[, (event_after_baseline) := fifelse(
    !is.na(get(event_time_var)) &
      !is.na(get(baseline_time_var)) &
      !is.na(get(censor_time_var)) &
      get(event_time_var) > get(baseline_time_var) &
      get(event_time_var) <= get(censor_time_var),
    get(event_time_var),
    as.Date(NA)
  )]

  dt[, (comp_after_baseline) := fifelse(
    !is.na(get(competing_time_var)) &
      !is.na(get(baseline_time_var)) &
      !is.na(get(censor_time_var)) &
      get(competing_time_var) > get(baseline_time_var) &
      get(competing_time_var) <= get(censor_time_var),
    get(competing_time_var),
    as.Date(NA)
  )]

  dt[, (event_date_var) := pmin(get(event_after_baseline), get(comp_after_baseline), get(censor_time_var), na.rm = TRUE)]
  dt[is.infinite(as.numeric(get(event_date_var))), (event_date_var) := as.Date(NA)]

  # 同一天同时满足目标事件和竞争死亡时，目标事件优先。
  dt[, (status_var) := fcase(
    !is.na(get(event_after_baseline)) & get(event_after_baseline) == get(event_date_var), 1L,
    !is.na(get(comp_after_baseline)) & get(comp_after_baseline) == get(event_date_var), 2L,
    !is.na(get(event_date_var)), 0L,
    default = NA_integer_
  )]

  dt[, (event_var) := fifelse(get(status_var) == 1L, 1L, 0L)]
  dt[, (comp_var) := fifelse(get(status_var) == 2L, 1L, 0L)]
  dt[, (fup_days_var) := as.numeric(get(event_date_var) - get(baseline_time_var))]
  dt[, (fup_years_var) := get(fup_days_var) / 365.25]

  dt[, c(event_after_baseline, comp_after_baseline) := NULL]
  dt[]
}

print_endpoint_qa <- function(dt, vars) {
  dt <- as.data.table(dt)
  vars <- intersect(vars, names(dt))

  for (v in vars) {
    cat("\n", v, "\n", sep = "")
    print(dt[, .N, by = v][order(get(v))])
  }

  invisible(NULL)
}



# ============================================================================== #
# 文件模块说明 | Function_RM_sources_M01.R
# 任务：RM来源路径、定义、已观测分类与总体RM派生
# 输入：01/02 RM原始来源，05已标准化字段
# 输出：RM源信息及0/1/NA分类
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# Shared RM source cleaning and classification for 01, 02, 05 and 06.
# Explicit three-valued OR: a known positive wins, two negatives yield zero, otherwise NA.
rm01_tri_or <- function(a,b) {
  if (length(a) != length(b)) stop("RM binary OR length mismatch")
  ifelse((!is.na(a) & a == 1L) | (!is.na(b) & b == 1L), 1L,
         ifelse(!is.na(a) & !is.na(b) & a == 0L & b == 0L,0L,NA_integer_))
}
rm01_version <- "M01_RM_v1.8_20261101"
rm01_source_registry <- NULL
rm01_configure_paths <- function(workbook, root = NULL) {
  if(!file.exists(workbook)) {
    mode <- Sys.getenv("M01_RM_SOURCE_MODE", "workbook")
    if (!identical(mode,"scan")) stop(
      "RM workbook missing: ", workbook,
      ". Set M01_RM_PATH_WORKBOOK or explicitly set M01_RM_SOURCE_MODE=scan ",
      "to resolve uniquely named CSVs below M01_UKB_DIR; duplicates will stop.")
    if (is.null(root) || !dir.exists(root)) stop("M01_UKB_DIR not found: ",root)
    candidates <- list.files(root, pattern="[.]csv$", recursive=TRUE, full.names=TRUE)
    candidates <- candidates[!grepl("协变量数据集", candidates)]
    if(!length(candidates)) stop("No RM CSV source files under: ",root)
    registry<-data.table::data.table(file_name=basename(candidates),
                                      path=normalizePath(candidates,winslash="/",mustWork=TRUE))
    rm01_source_registry <<- registry
    return(registry)
  }
  if(!requireNamespace("readxl",quietly=TRUE))stop("RM workbook paths require readxl")
  dictionary<-as.data.frame(readxl::read_excel(workbook,sheet="字段路径字典",col_types="text"))
  column<-"完整本地路径"
  if(!column %in% names(dictionary))stop("RM path workbook missing path column")
  paths<-unique(trimws(dictionary[[column]]))
  paths<-paths[!is.na(paths)&nzchar(paths)]
  paths<-normalizePath(paths,winslash="/",mustWork=FALSE)
  if(any(!file.exists(paths)))stop("Workbook source files missing: ",paste(paths[!file.exists(paths)],collapse="; "))
  registry<-data.table::data.table(file_name=basename(paths),path=paths)
  if(anyDuplicated(registry$file_name))stop("Ambiguous file names in RM path workbook")
  rm01_source_registry <<- registry
  registry
}
rm01_resolve <- function(root,name) {
  if(!is.null(rm01_source_registry)) {
    path<-rm01_source_registry$path[rm01_source_registry$file_name==name]
    if(length(path)!=1L)stop("Source not uniquely specified in RM workbook: ",name)
    return(path)
  }
  paths<-list.files(root,pattern="\\.csv$",full.names=TRUE,recursive=TRUE)
  paths<-paths[basename(paths)==name & !grepl("协变量数据集",paths)]
  if(length(paths)!=1L)stop("Expected one RM source: ",name)
  paths
}
rm01_num <- function(x) {
  z <- suppressWarnings(as.numeric(x)); z[!is.finite(z) | z < 0] <- NA_real_; z
}
rm01_binary <- function(x) {
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c("yes", "1"), 1L, ifelse(z %in% c("no", "0"), 0L, NA_integer_))
}
rm01_unique <- function(d, id) {
  if (!id %in% names(d) || anyNA(d[[id]]) || any(!nzchar(trimws(as.character(d[[id]])))) || anyDuplicated(d[[id]])) stop("RM source ID invalid")
}
rm01_read <- function(root, name, columns = NULL) {
  paths <- rm01_resolve(root,name)
  h <- names(data.table::fread(paths, nrows = 0L))
  if (!is.null(columns)) {
    absent <- setdiff(columns, h)
    if (length(absent)) stop(name, " missing columns: ", paste(absent, collapse = ", "))
  }
  d <- data.table::fread(paths, select = columns, colClasses = "character", na.strings = c("", "NA"))
  rm01_unique(d, "Participant ID"); data.table::setattr(d, "source_path", paths); d
}
rm01_align <- function(d, ids) d[match(as.character(ids), as.character(d[["Participant ID"]]))]
rm01_codes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  toupper(gsub("\\.", "", sub(" .*", "", strsplit(x, "|", fixed = TRUE)[[1L]])))
}
rm01_record_definitions <- function() list(
  inflammatory = list(icd="^(M05|M06|M07|M32|L40|K50|K51)", self=c("rheumatoid arthritis","systemic lupus erythematosis/sle","psoriasis","psoriatic arthropathy","inflammatory bowel disease","crohns disease","ulcerative colitis")),
  psychiatric = list(icd="^(F20|F31)|^(F322|F323|F332|F333)$", self=c("schizophrenia","mania/bipolar disorder/manic depression")),
  af = list(icd="^I48(0|1|2|9|11|19|20|21)$", self="atrial fibrillation"),
  hdp = list(icd="^(O11|O13|O14|O15|O16)", self="gestational hypertension/pre-eclampsia"),
  hiv = list(icd="^(B20|B21|B22|B23|B24|Z21)", self="hiv/aids"),
  sleep_apnoea = list(icd="^G473$", self="sleep apnoea"),
  fatty_liver = list(icd="^K760$", self=character())
)
rm01_diagnoses <- function(ids, baseline, self, hes, dates) {
  # Order in 41270's pipe-separated list is the UKB array index.
  n <- length(ids); out <- data.table::data.table(`Participant ID`=as.character(ids))
  sc <- grep("^Non-cancer illness code, self-reported \\| Instance 0 \\| Array", names(self), value=TRUE)
  if (!length(sc)) stop("Missing baseline 20002 arrays")
  self <- rm01_align(self, ids); hes <- rm01_align(hes, ids)
  dates <- lapply(dates,rm01_align,ids=ids)
  self_known <- rowSums(vapply(sc,function(v) {
    z<-tolower(trimws(self[[v]]));!is.na(z)&nzchar(z)&!grepl("do not know|prefer not|-1$|-3$",z)
  },logical(n)))>0L
  htxt <- hes[["Diagnoses - ICD10"]]
  hknown <- !is.na(htxt) & nzchar(htxt)
  lists <- lapply(htxt, rm01_codes)
  bd <- as.Date(baseline)
  for (key in names(rm01_record_definitions())) {
    def <- rm01_record_definitions()[[key]]
    spos <- if (length(def$self)) rowSums(vapply(sc, function(c) {
      tolower(self[[c]]) %in% def$self
    }, logical(n))) > 0L else rep(FALSE,n)
    hits <- lapply(lists, function(codes) which(grepl(def$icd,codes))-1L)
    pos <- unknown <- rep(FALSE,n)
    # Build sparse array-to-row mapping once instead of scanning 502k IDs per array.
    hit_rows <- rep(seq_len(n),lengths(hits))
    hit_arrays <- unlist(hits,use.names=FALSE)
    rows_by_array <- split(hit_rows,hit_arrays)
    for (array_name in names(rows_by_array)) {
      a <- as.integer(array_name)
      rows <- rows_by_array[[array_name]]
      col <- paste0("Date of first in-patient diagnosis - ICD10 | Array ",a)
      dt <- as.Date(rep(NA_character_,length(rows)))
      found <- FALSE
      for (part in dates) if (col %in% names(part)) {
        v <- as.Date(part[[col]][rows])
        v[v %in% as.Date(c("1900-01-01","1901-01-01","1902-02-02","1903-03-03","1909-09-09","2037-07-07"))] <- as.Date(NA)
        conflict <- !is.na(dt) & !is.na(v) & dt != v
        if (any(conflict)) stop("Conflicting 41280 date arrays: ", col)
        dt[is.na(dt)] <- v[is.na(dt)]; found <- TRUE
      }
      if (!found) stop("Missing 41280 array for target ICD: ",a)
      unknown[rows] <- unknown[rows] | is.na(dt) | is.na(bd[rows]) | dt < as.Date("1900-01-01")
      good <- !is.na(dt) & !is.na(bd[rows]) & dt >= as.Date("1900-01-01") & dt <= bd[rows]
      pos[rows[good]] <- TRUE
    }
    covered <- if (key=="fatty_liver") hknown else self_known | hknown
    value <- ifelse(spos | pos,1L,ifelse(covered & !unknown,0L,NA_integer_))
    out[, (paste0("rm_source_",key)) := as.integer(value)]
  }
  out
}
rm02_read_records <- function(root) {
  self_name <- "3_Assessment centre_4_Verbal interyiew__3_Medical conditions5.csv"
  h <- names(data.table::fread(rm01_resolve(root, self_name), nrows=0L))
  self <- rm01_read(root, self_name, c("Participant ID", grep(
    "^Non-cancer illness code, self-reported \\| Instance 0 \\| Array", h, value=TRUE)))
  dates <- lapply(1:2, function(k) {
    name <- paste0("8_Health-related outcomes_3_Hospital inpatient_2_Summary Diagnoses",k,".csv")
    h <- names(data.table::fread(rm01_resolve(root,name),nrows=0L))
    cols <- c("Participant ID", intersect("Diagnoses - ICD10",h),
      grep("^Date of first in-patient diagnosis - ICD10",h,value=TRUE))
    rm01_read(root,name,cols)
  })
  list(self=self, hes=dates[[1L]], dates=dates)
}
rm01_attach <- function(base, sources, id) {
  rm01_unique(base,id); rm01_unique(sources,"Participant ID")
  ix <- match(as.character(base[[id]]),as.character(sources[["Participant ID"]]))
  if(anyNA(ix)) stop("RM sources do not cover analytic cohort")
  overlap <- intersect(names(base),setdiff(names(sources),"Participant ID"))
  if(length(overlap)) stop("RM source columns already present: ",paste(overlap,collapse=", "))
  for(v in setdiff(names(sources),"Participant ID")) data.table::set(base,j=v,value=sources[[v]][ix])
  base
}

# Descriptive record audit. Observed maxima are never treated as verified coverage.
audit_record_cutoffs <- function(data_date,data_outcome,records,admin_date) {
 ids<-as.character(data_date[["Participant ID"]]);regions<-data_date$Region
 safe_max<-function(x) {x<-as.Date(unique(x));if(all(is.na(x)))as.Date(NA) else max(x,na.rm=TRUE)}
 result<-list();k<-0L
 for(g in sort(unique(regions))) {
  use<-which(regions==g)
  for(src in c("hospital_diagnosis","death")) {
   dates<-as.Date(character())
   if(src=="death")dates<-data_outcome$Death_time[match(ids[use],as.character(data_outcome[["Participant ID"]]))]
   else for(part in records$dates) {
    ix<-match(ids[use],as.character(part[["Participant ID"]]))
    cs<-grep("^Date of first in-patient diagnosis - ICD10",names(part),value=TRUE)
    for(v in cs) dates<-c(dates,safe_max(part[[v]][ix]))
   }
   k<-k+1L
   result[[k]]<-data.table(region=g,source=src,observed_latest_record=as.character(safe_max(dates)),
    adopted_admin_date=as.character(admin_date),verified_export_coverage_date=NA_character_,
    coverage_status="Unverified for local export; observed maxima do not establish completeness",
    region_basis="Baseline assessment centre; not record-provider geography",
    reference_url="https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates.html")
  }
 }
 audit<-rbindlist(result)
 audit[,admin_minus_latest_days:=as.integer(as.Date(adopted_admin_date)-as.Date(observed_latest_record))]
 audit[,admin_beyond_observed_max:=admin_minus_latest_days>0L]
 audit
}

