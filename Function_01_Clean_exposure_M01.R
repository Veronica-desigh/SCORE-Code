# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_01_Clean_exposure_M01.R
# 任务：供01调用的UKB暴露清理函数
# 输入：UKB字段向量/表
# 输出：标准化基线暴露及检查结果
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function__01_exposure_M01.R
# ------------------------------------------------------------------------------ #
# Purpose:
#   Helper functions used by 01_Clean_exposure_M01.R for the M01 ASCVD
#   baseline exposure/covariate cleaning pipeline.
#
# How this file maps to 01_Clean_exposure_M01.R:
#   0. General QA helpers
#      - var_info_print()
#   1.2 Touchscreen sociodemographics and SES components
#      - classify_education()
#      - classify_household_income()
#      - classify_employment_status()
#      - classify_ethnicity()
#      - classify_physical_activity()
#   3.3 Alcohol score
#      - encode_alcohol()
#      - alcohol_fun_one_instance()
#   3.4 Diet score
#      - process_*_data()
#      - calc_diet_vars()
#   5.1 Blood biomarkers
#      - standardize_variable_names()
#   6.1 Family history of CVD
#      - classify_family_history_cvd()
#
# Design notes:
#   - These functions do not set local paths and are intended to be called by the
#     01 main script after data tables have been read from disk.
#   - SES components are baseline variables only. The main script passes
#     Qualifications | Instance 0, Current employment status | Instance 0, and
#     Average total household income before tax | Instance 0.
#   - Final SES is intentionally NOT fitted in Function_01/01. The formal
#     2-class LCA is fitted after MI in 06 using completed SES components.
#   - Missing/non-informative UKB responses are generally returned as typed NA.
# ============================================================================== #

library(tidyverse)
library(dplyr)
library(data.table)


# ============================================================================== #
# 0. General quality-control helper used across the 01 main script
# ============================================================================== #

# ---- 0.1 Print variable type and missingness ----
# Called after reading or deriving key variables in 01_Clean_exposure_M01.R.
# It prints the first class, number of missing values, missing proportion, and
# total length for each variable. This is only a diagnostic helper; it does not
# modify the input data.
var_info_print <- function(data, var_names) {
  for (var_name in var_names) {
    x <- data[[var_name]]
    cat(
      "变量:", var_name,
      "| 类型:", class(x)[1],
      "| 缺失值:", sum(is.na(x)),
      "| 缺失比例:", round(mean(is.na(x)) * 100, 1), "%",
      "| 总数:", length(x), "\n"
    )
  }
}

# ---- 0.2 Assert one row per participant ----
# This is deliberately strict: duplicate or missing IDs must be resolved by the
# component-specific cleaning logic, never silently dropped before a merge.
assert_one_row_per_id <- function(data,
                                  data_name,
                                  id_col = "Participant ID") {
  dt <- data.table::as.data.table(data)

  if (!id_col %in% names(dt)) {
    stop(data_name, " is missing ID column: ", id_col, call. = FALSE)
  }

  id <- dt[[id_col]]
  missing_id <- is.na(id)
  if (is.character(id)) {
    missing_id <- missing_id | trimws(id) == ""
  }
  if (any(missing_id)) {
    stop(
      data_name, " contains ", sum(missing_id), " missing/blank ", id_col,
      " value(s).",
      call. = FALSE
    )
  }

  duplicate_rows <- duplicated(id) | duplicated(id, fromLast = TRUE)
  if (any(duplicate_rows)) {
    duplicate_ids <- unique(id[duplicate_rows])
    preview <- paste(utils::head(duplicate_ids, 10L), collapse = ", ")
    stop(
      data_name, " contains ", length(duplicate_ids), " duplicated ", id_col,
      " value(s). Example(s): ", preview,
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---- 0.3 Audited left join ----
# The first argument is always the explicit master table. The helper preserves
# its row order and ID set, rejects overlapping non-ID columns, and appends one
# row to the merge audit for every component join.
audited_left_join <- function(master,
                              component,
                              component_name,
                              audit_log = NULL,
                              id_col = "Participant ID") {
  master_dt <- data.table::as.data.table(master)
  component_dt <- data.table::as.data.table(component)

  assert_one_row_per_id(master_dt, "master table", id_col)

  if (!id_col %in% names(component_dt)) {
    stop(component_name, " is missing ID column: ", id_col, call. = FALSE)
  }

  component_id <- component_dt[[id_col]]
  component_missing_id <- is.na(component_id)
  if (is.character(component_id)) {
    component_missing_id <- component_missing_id | trimws(component_id) == ""
  }
  component_duplicate_ids <- unique(component_id[
    !component_missing_id &
      (duplicated(component_id) | duplicated(component_id, fromLast = TRUE))
  ])

  if (any(component_missing_id)) {
    stop(
      component_name, " contains ", sum(component_missing_id),
      " missing/blank ", id_col, " value(s).",
      call. = FALSE
    )
  }
  if (length(component_duplicate_ids) > 0L) {
    preview <- paste(utils::head(component_duplicate_ids, 10L), collapse = ", ")
    stop(
      component_name, " contains ", length(component_duplicate_ids),
      " duplicated ", id_col, " value(s). Example(s): ", preview,
      call. = FALSE
    )
  }

  overlapping_cols <- intersect(
    setdiff(names(master_dt), id_col),
    setdiff(names(component_dt), id_col)
  )
  if (length(overlapping_cols) > 0L) {
    stop(
      component_name, " has non-ID column(s) already present in the master: ",
      paste(overlapping_cols, collapse = ", "),
      call. = FALSE
    )
  }

  master_ids_before <- master_dt[[id_col]]
  component_ids <- component_dt[[id_col]]
  matched_ids <- sum(master_ids_before %in% component_ids)

  row_order_col <- ".m01_master_row_order"
  if (row_order_col %in% c(names(master_dt), names(component_dt))) {
    stop("Reserved merge-order column already exists: ", row_order_col, call. = FALSE)
  }

  result <- data.table::copy(master_dt)
  result[, (row_order_col) := .I]
  result <- merge(
    result,
    data.table::copy(component_dt),
    by = id_col,
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setorderv(result, row_order_col)
  result[, (row_order_col) := NULL]
  data.table::setcolorder(
    result,
    c(names(master_dt), setdiff(names(component_dt), id_col))
  )

  if (nrow(result) != nrow(master_dt) ||
      !identical(result[[id_col]], master_ids_before)) {
    stop(
      "Left join with ", component_name,
      " changed the master row count or participant order.",
      call. = FALSE
    )
  }

  audit_row <- data.table::data.table(
    join_order = if (is.null(audit_log)) 1L else nrow(audit_log) + 1L,
    component = as.character(component_name),
    master_n_before = nrow(master_dt),
    master_unique_ids_before = data.table::uniqueN(master_ids_before),
    component_n = nrow(component_dt),
    component_unique_ids = data.table::uniqueN(component_ids),
    component_duplicate_ids = length(component_duplicate_ids),
    matched_ids = matched_ids,
    unmatched_ids = nrow(master_dt) - matched_ids,
    master_n_after = nrow(result),
    master_unique_ids_after = data.table::uniqueN(result[[id_col]])
  )

  audit_out <- data.table::rbindlist(
    list(audit_log, audit_row),
    use.names = TRUE,
    fill = TRUE
  )

  list(data = result, audit = audit_out)
}

# ---- 0.4 Final critical-column type checks ----
assert_column_types <- function(data, type_spec) {
  dt <- data.table::as.data.table(data)
  missing_cols <- setdiff(names(type_spec), names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Final exposure data are missing critical column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  for (column_name in names(type_spec)) {
    expected <- type_spec[[column_name]]
    value <- dt[[column_name]]
    valid <- switch(
      expected,
      numeric = is.numeric(value),
      character = is.character(value),
      Date = inherits(value, "Date"),
      integer = is.integer(value),
      logical = is.logical(value),
      stop("Unsupported expected type: ", expected, call. = FALSE)
    )
    if (!isTRUE(valid)) {
      stop(
        "Column ", column_name, " must be ", expected,
        "; actual class is ", paste(class(value), collapse = "/"), ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# ============================================================================== #
# 1.2 Touchscreen sociodemographics and individual-level SES helpers
# ============================================================================== #

# ---- 1.2.1 Baseline educational attainment ----
# Main-script source field:
#   Qualifications | Instance 0
# Main-script output variable:
#   Education
#
# Classification rule:
#   High education   = College or University degree
#   Medium education = A levels/AS levels, O levels/GCSEs, CSEs,
#                      NVQ or HND or HNC, Other professional qualifications
#   Low education    = None of the above
#
# Important revision:
#   O levels/GCSEs and CSEs are classified as Medium education, not as Less than
#   high school, to keep the UKB SES classification closer to high-impact UKB
#   SES/LCA applications.
classify_education <- function(qualifications) {
  # Field 6138: Qualifications.
  # Baseline individual-level SES component. Use the highest available
  # qualification within Instance 0. In line with UKB SES/LCA applications,
  # O levels/GCSEs and CSEs are classified as medium education rather than
  # less than high school.
  qualifications <- trimws(as.character(qualifications))
  q_lookup <- data.table::fifelse(is.na(qualifications), "", qualifications)
  
  fcase(
    is.na(qualifications) |
      qualifications == "" |
      qualifications %in% c("Prefer not to answer", "Do not know"),
    NA_character_,
    
    grepl("College or University degree", q_lookup, ignore.case = TRUE),
    "High education",
    
    grepl(
      paste(
        c(
          "A levels/AS levels",
          "O levels/GCSEs",
          "CSEs",
          "NVQ or HND or HNC",
          "Other professional qualifications"
        ),
        collapse = "|"
      ),
      q_lookup,
      ignore.case = TRUE
    ),
    "Medium education",
    
    grepl("None of the above", q_lookup, ignore.case = TRUE),
    "Low education",
    
    default = NA_character_
  )
}

# ---- 1.2.2 Baseline household income ----
# Main-script source field:
#   Average total household income before tax | Instance 0
# Main-script output variable:
#   Household_income
#
# Field 738 contains ordered household-income bands. For LCA, the original UKB
# income bands are retained rather than collapsed to low/medium/high, because
# LCA can use the full ordinal pattern across income, education, and employment.
# Non-informative responses are set to NA_character_.
classify_household_income <- function(income) {
  income <- trimws(as.character(income))
  income_lookup <- data.table::fifelse(is.na(income), "", income)
  
  fcase(
    is.na(income) |
      income == "" |
      income %in% c("Do not know", "Prefer not to answer"),
    NA_character_,
    
    grepl("Less than.*18", income_lookup, ignore.case = TRUE),
    "Less than 18,000",
    
    grepl("18[, ]*000.*30[, ]*999", income_lookup, ignore.case = TRUE),
    "18,000 to 30,999",
    
    grepl("31[, ]*000.*51[, ]*999", income_lookup, ignore.case = TRUE),
    "31,000 to 51,999",
    
    grepl("52[, ]*000.*100[, ]*000", income_lookup, ignore.case = TRUE),
    "52,000 to 100,000",
    
    grepl("Greater than.*100", income_lookup, ignore.case = TRUE),
    "Greater than 100,000",
    
    default = NA_character_
  )
}


# ---- 1.2.3 Baseline current employment status ----
# Main-script source field:
#   Current employment status | Instance 0
# Main-script output variable:
#   employment_status
#
# Naming rationale:
#   UKB Field 6142 is current employment status, not a true occupation code.
#   Therefore the derived variable is named employment_status rather than
#   Occupation.
#
# Three-category classification:
#   Employed/self-employed = paid employment or self-employment
#   Retired                = retired
#   Other not employed     = looking after home/family, unable to work,
#                            unemployed, unpaid/voluntary work, or student
#
# This three-category version is intended for the LCA-derived individual-level
# SES variable requested in the main analysis update.
classify_employment_status <- function(employment) {
  employment <- trimws(as.character(employment))
  emp_lookup <- data.table::fifelse(is.na(employment), "", employment)
  
  fcase(
    is.na(employment) |
      employment == "" |
      grepl("Prefer not to answer|None of the above", emp_lookup, ignore.case = TRUE),
    NA_character_,
    
    grepl("In paid employment|self-employed", emp_lookup, ignore.case = TRUE),
    "Employed/self-employed",
    
    grepl("Retired", emp_lookup, ignore.case = TRUE),
    "Retired",
    
    grepl(
      paste(
        c(
          "Looking after home",
          "Unable to work",
          "Unemployed",
          "unpaid",
          "voluntary",
          "Full or part-time student"
        ),
        collapse = "|"
      ),
      emp_lookup,
      ignore.case = TRUE
    ),
    "Other not employed",
    
    default = NA_character_
  )
}

# ---- 1.2.4 Ethnicity ----
# Main-script source fields:
#   Ethnic background | Instance 0-3
# Main-script output variables:
#   Ethnic0, Ethnic1, Ethnic2, Ethnic3, Ethnic
#
# The function collapses ethnicity to white/non-white. The main script derives
# the baseline Ethnic variable from Instance 0 only; Instances 1-3 are retained
# as separate follow-up variables and never backfill the baseline value.
classify_ethnicity <- function(ethnic_var) {
  fcase(
    ethnic_var %in% c("British", "Irish", "Any other white background"), "white",
    ethnic_var %in% c("Prefer not to answer", "Do not know", ""), NA_character_,
    !is.na(ethnic_var), "non-white", # 其他非缺失值
    default = NA_character_
  )
}

# ---- 1.2.5 Physical activity ----
# Main-script source fields:
#   Summed MET minutes per week for all activity | Instance 0-3
# Main-script output variables:
#   MET_minutes_week0-3 (raw numeric values are created in 01)
#   PhysicalActivity0-3 (deterministic categories)
#
# Missing raw MET must remain missing. 06 will impute the raw baseline MET value
# and re-apply classify_physical_activity() after imputation.
#
# MET-minutes are grouped into four ordered activity levels. Negative values and
# missing values are treated as missing.
classify_physical_activity <- function(met_value) {
  fcase(
    is.na(met_value) | met_value < 0, NA_character_, # 处理缺失值和负值
    met_value <= 600, "level1", # 低活动
    met_value <= 3000, "level2", # 中活动
    met_value <= 6000, "level3", # 高活动
    met_value > 6000, "level4", # 极高活动
    default = NA_character_
  )
}



# ============================================================================== #
# 6.1 Family history of CVD helper
# ============================================================================== #
# Main-script source fields:
#   Illnesses of father, Illnesses of mother, Illnesses of siblings
# Main-script output variable:
#   Family_history_CVD
#
# Definition:
#   Positive if father, mother, or siblings reported heart disease or stroke.
#   Non-informative responses are treated as missing. If all available responses
#   are missing/non-informative, the derived variable is set to NA_integer_.
classify_family_history_cvd <- function(father, mother, siblings) {
  df <- data.table(
    father = as.character(father),
    mother = as.character(mother),
    siblings = as.character(siblings)
  )
  
  non_informative <- c("", "Do not know", "Prefer not to answer")
  
  for (v in names(df)) {
    df[[v]] <- trimws(df[[v]])
    df[[v]][df[[v]] %in% non_informative] <- NA_character_
  }
  
  positive_mat <- sapply(df, function(x) {
    grepl("heart disease|stroke", x, ignore.case = TRUE)
  })
  informative_mat <- sapply(df, function(x) {
    !is.na(x)
  })
  
  data.table::fifelse(
    rowSums(positive_mat, na.rm = TRUE) > 0,
    1L,
    data.table::fifelse(
      rowSums(informative_mat, na.rm = TRUE) > 0,
      0L,
      NA_integer_
    )
  )
}

# ============================================================================== #
# 3.3 Alcohol-score helpers
# ============================================================================== #
# Main-script source fields:
#   Alcohol intake frequency | Instance 0-3
#   Average monthly/weekly intake fields for red wine, white wine/champagne,
#   beer/cider, spirits, and fortified wine
# Main-script output variables:
#   alcohol_score_0, alcohol_score_1, alcohol_score_2, alcohol_score_3
#
# Overall logic:
#   For each instance, the main script calls alcohol_fun_one_instance(). The
#   helper converts beverage-specific monthly/weekly servings to pure alcohol
#   g/day and then derives a binary score: 1 if <=10 g/day, 0 if >10 g/day.
encode_alcohol <- function(status) {
  fcase(
    status == "Never", "0",
    status == "Previous", "1",
    status == "Current", "2",
    default = NA_character_
  )
}

# ---- 3.3.1 Derive alcohol score for one UKB instance ----
# Arguments:
#   data = the wide raw alcohol dataset read in the main script
#   k    = UKB instance number: 0, 1, 2, or 3
# Returns:
#   Participant ID, pure_alcohol_<k>, and alcohol_score_<k>
alcohol_fun_one_instance <- function(data, k) {
  k_chr <- as.character(k)
  
  ## 2.1 这一轮用到的列名（严格照 UKB 原名）
  freq_name <- paste0("Alcohol intake frequency. | Instance ", k_chr)
  
  monthly_vars <- c(
    "Average monthly red wine intake | Instance ",
    "Average monthly champagne plus white wine intake | Instance ",
    "Average monthly beer plus cider intake | Instance ",
    "Average monthly spirits intake | Instance ",
    "Average monthly fortified wine intake | Instance "
  ) %>% paste0(k_chr)
  
  weekly_vars <- c(
    "Average weekly red wine intake | Instance ",
    "Average weekly champagne plus white wine intake | Instance ",
    "Average weekly beer plus cider intake | Instance ",
    "Average weekly spirits intake | Instance ",
    "Average weekly fortified wine intake | Instance "
  ) %>% paste0(k_chr)
  
  ## 2.2 只抽取这一轮相关变量，重编码 alcohol_frequency，并转为 numeric
  df_k <- data %>%
    select(
      `Participant ID`,
      all_of(freq_name),
      all_of(monthly_vars),
      all_of(weekly_vars)
    ) %>%
    # 重编码 alcohol_frequency：0=从不，1=偶尔/月饮，2=周饮及以上
    mutate(
      alcohol_frequency = case_when(
        .data[[freq_name]] == "Never" ~ 0L,
        .data[[freq_name]] %in% c("Special occasions only", "One to three times a month") ~ 1L,
        .data[[freq_name]] %in% c("Once or twice a week", "Three or four times a week", "Daily or almost daily") ~ 2L,
        TRUE ~ NA_integer_
      )
    ) %>%
    # 每月 / 每周饮酒量转成 numeric
    mutate(
      across(
        .cols = c(all_of(monthly_vars), all_of(weekly_vars)),
        .fns  = ~ suppressWarnings(as.numeric(.x))
      )
    )
  
  ## 2.3 按饮酒频率分组，先算每周 g/week
  # 频率缺失 → 全 NA
  data_na <- df_k[is.na(df_k$alcohol_frequency), ] %>%
    mutate(
      red_wine2             = NA_real_,
      white_wine_champagne2 = NA_real_,
      beer_cider2           = NA_real_,
      spirits2              = NA_real_,
      fortified_wine2       = NA_real_
    )
  
  # 从不饮酒 → 全 0
  data_never <- df_k %>%
    filter(alcohol_frequency == 0) %>%
    mutate(
      red_wine2             = 0,
      white_wine_champagne2 = 0,
      beer_cider2           = 0,
      spirits2              = 0,
      fortified_wine2       = 0
    )
  
  # 偶尔 / 每月饮酒：用月杯数，/4 → 周，再按系数换算 g/week
  data_less <- df_k %>%
    filter(alcohol_frequency == 1) %>%
    mutate(
      red_wine2 = (.[[monthly_vars[1]]] / 4) * 0.85 * 8,
      white_wine_champagne2 = (.[[monthly_vars[2]]] / 4) * 0.85 * 8,
      beer_cider2 = (.[[monthly_vars[3]]] / 4) * 1.28 * 8,
      spirits2 = (.[[monthly_vars[4]]] / 4) * 0.57 * 8,
      fortified_wine2 = (.[[monthly_vars[5]]] / 4) * 0.56 * 8
    )
  
  # 常规 / 每周饮酒：直接用周杯数换算 g/week
  data_regular <- df_k %>%
    filter(alcohol_frequency == 2) %>%
    mutate(
      red_wine2 = (.[[weekly_vars[1]]]) * 0.85 * 8,
      white_wine_champagne2 = (.[[weekly_vars[2]]]) * 0.85 * 8,
      beer_cider2 = (.[[weekly_vars[3]]]) * 1.28 * 8,
      spirits2 = (.[[weekly_vars[4]]]) * 0.57 * 8,
      fortified_wine2 = (.[[weekly_vars[5]]]) * 0.56 * 8
    )
  
  ## 2.4 合并四类人群 → 算每周总量，再转每天 g/day
  out_k <- bind_rows(data_na, data_never, data_less, data_regular) %>%
    mutate(
      # 判断是否 5 种酒全 NA
      all_na = if_else(
        is.na(red_wine2) &
          is.na(white_wine_champagne2) &
          is.na(beer_cider2) &
          is.na(spirits2) &
          is.na(fortified_wine2),
        TRUE, FALSE
      ),
      # 总纯酒精 g/week：部分 NA 视为 0，全 NA 则 NA
      pure_alcohol_week = if_else(
        all_na,
        NA_real_,
        rowSums(
          cbind(
            red_wine2,
            white_wine_champagne2,
            beer_cider2,
            spirits2,
            fortified_wine2
          ),
          na.rm = TRUE
        )
      ),
      # 每天 g/day
      pure_alcohol_day = pure_alcohol_week / 7
    ) %>%
    transmute(
      `Participant ID`,
      !!paste0("pure_alcohol_", k_chr) := pure_alcohol_day,
      !!paste0("alcohol_score_", k_chr) := case_when(
        !is.na(pure_alcohol_day) & pure_alcohol_day <= 10 ~ 1L,
        !is.na(pure_alcohol_day) & pure_alcohol_day > 10 ~ 0L, # 10 - 按照指南里的标准来
        TRUE ~ NA_integer_
      )
    )
  
  return(out_k)
}


# ============================================================================== #
# 3.4 Diet-score helpers
# ============================================================================== #
# Main-script source fields:
#   Baseline and repeat dietary touchscreen fields across Instances 0-3
# Main-script output variables:
#   diet_score_w0, diet_score_w1, diet_score_w2, diet_score_w3
#
# Overall logic:
#   Each process_*_data() function cleans one dietary component for one instance.
#   calc_diet_vars() combines component-level indicators into the final healthy
#   diet score for each wave/instance.

# ---- 3.4.1 Fruit intake ----
# Cleans fresh fruit and dried fruit intake for one UKB instance.
# Returns component-level fruit variables used by calc_diet_vars().
process_fruit_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 处理新鲜水果
      !!paste0("fre_fruit_w", instance) := as.numeric(
        ifelse(
          get(paste0("Fresh.fruit.intake", suffix)) == "Less than one", 0.5,
          ifelse(
            get(paste0("Fresh.fruit.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
            NA,
            get(paste0("Fresh.fruit.intake", suffix))
          )
        )
      ),
      # 处理干果
      !!paste0("dry_fruit_w", instance) := as.numeric(
        ifelse(
          get(paste0("Dried.fruit.intake", suffix)) == "Less than one", 0.5,
          ifelse(
            get(paste0("Dried.fruit.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
            NA,
            get(paste0("Dried.fruit.intake", suffix))
          )
        )
      ),
      
      # 计算总水果摄入量 (5份干果为1 serving，一份新鲜水果为1 serving)
      !!paste0("fruit_w", instance) := ifelse(
        is.na(get(paste0("fre_fruit_w", instance))) & is.na(get(paste0("dry_fruit_w", instance))),
        NA_real_,
        replace_na(get(paste0("fre_fruit_w", instance)), 0) +
          (replace_na(get(paste0("dry_fruit_w", instance)), 0) / 5)
      ),
      
      # 创建水果评分
      !!paste0("fruit_score_w", instance) := case_when(
        replace_na(get(paste0("fre_fruit_w", instance)), 0) >= 3 ~ 1,
        replace_na(get(paste0("dry_fruit_w", instance)) / 2, 0) >= 3 ~ 1,
        replace_na(get(paste0("fruit_w", instance)), 0) >= 3 ~ 1,
        is.na(get(paste0("fruit_w", instance))) &
          is.na(get(paste0("fre_fruit_w", instance))) &
          is.na(get(paste0("dry_fruit_w", instance))) ~ NA_real_,
        TRUE ~ 0
      )
    )
}


# vegetable
# ---- 3.4.2 Vegetable intake ----
# Cleans cooked vegetable and salad/raw vegetable intake for one UKB instance.
process_vegetable_data <- function(data, instance) {
  # 构造变量名后缀
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 处理熟蔬菜
      !!paste0("cook_vegetable_w", instance) := as.numeric(
        ifelse(
          get(paste0("Cooked.vegetable.intake", suffix)) == "Less than one", 0.5,
          ifelse(
            get(paste0("Cooked.vegetable.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
            NA,
            get(paste0("Cooked.vegetable.intake", suffix))
          )
        )
      ),
      
      # 处理生蔬菜
      !!paste0("raw_vegetable_w", instance) := as.numeric(
        ifelse(
          get(paste0("Salad...raw.vegetable.intake", suffix)) == "Less than one", 0.5,
          ifelse(
            get(paste0("Salad...raw.vegetable.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
            NA,
            get(paste0("Salad...raw.vegetable.intake", suffix))
          )
        )
      ),
      
      # 计算总蔬菜摄入量(serving)，3份生或者熟为1serving
      !!paste0("vegetable_w", instance) := ifelse(
        is.na(get(paste0("cook_vegetable_w", instance))) &
          is.na(get(paste0("raw_vegetable_w", instance))),
        NA_real_,
        round(
          (replace_na(get(paste0("cook_vegetable_w", instance)), 0) +
             replace_na(get(paste0("raw_vegetable_w", instance)), 0)) / 3,
          2
        )
      ),
      
      # 创建蔬菜评分
      !!paste0("vegetable_score_w", instance) := case_when(
        get(paste0("cook_vegetable_w", instance)) / 3 >= 3 ~ 1,
        get(paste0("raw_vegetable_w", instance)) / 3 >= 3 ~ 1,
        get(paste0("vegetable_w", instance)) >= 3 ~ 1,
        is.na(get(paste0("vegetable_w", instance))) &
          is.na(get(paste0("cook_vegetable_w", instance))) &
          is.na(get(paste0("raw_vegetable_w", instance))) ~ NA_real_,
        TRUE ~ 0
      )
    )
}


# Whole grains、Refined grains
# ---- 3.4.3 Whole/refined grains ----
# Derives bread/cereal-related variables and classifies grain quality for one
# UKB instance.
process_whole_Refined_grains_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    # 处理面包摄入量
    mutate(
      !!paste0("bread_w", instance) := ifelse(
        get(paste0("Bread.intake", suffix)) == "Less than one", 0.5,
        ifelse(
          get(paste0("Bread.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
          NA,
          get(paste0("Bread.intake", suffix))
        )
      ),
      !!paste0("bread_w", instance, "t") := as.numeric(get(paste0("bread_w", instance))) / 7
    ) %>%
    # 处理面包类型
    mutate(
      !!paste0("bread_type_w", instance) := case_when(
        get(paste0("Bread.type", suffix)) == "White" ~ "White",
        get(paste0("Bread.type", suffix)) %in% c("Brown", "Other type of bread", "Wholemeal or wholegrain") ~ "other"
      )
    ) %>%
    # 计算非白面包摄入量
    mutate(
      !!paste0("withoutwhite_bread_w", instance) :=
        ifelse(get(paste0("bread_type_w", instance)) == "White", 0, get(paste0("bread_w", instance, "t")))
    ) %>%
    # 处理谷物摄入量
    mutate(
      !!paste0("cereal_w", instance) := as.numeric(
        ifelse(
          get(paste0("Cereal.intake", suffix)) == "Less than one", 0.5,
          ifelse(
            get(paste0("Cereal.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
            NA,
            get(paste0("Cereal.intake", suffix))
          )
        )
      ),
      !!paste0("cereal_w", instance, "t") :=
        get(paste0("cereal_w", instance)) / 7 * 2
    ) %>%
    # 计算全谷物评分
    mutate(
      !!paste0("Wholegrains_score_w", instance) := case_when(
        is.na(get(paste0("withoutwhite_bread_w", instance))) & is.na(get(paste0("cereal_w", instance, "t"))) & is.na(get(paste0("bread_w", instance, "t"))) ~ NA_real_,
        replace_na(get(paste0("withoutwhite_bread_w", instance)), 0) + replace_na(get(paste0("cereal_w", instance, "t")), 0) < 3 ~ 0,
        replace_na(get(paste0("withoutwhite_bread_w", instance)), 0) + replace_na(get(paste0("cereal_w", instance, "t")), 0) >= 3 ~ 1,
        get(paste0("withoutwhite_bread_w", instance)) >= 3 ~ 1,
        get(paste0("cereal_w", instance, "t")) >= 3 ~ 1,
        is.na(get(paste0("bread_type_w", instance))) &
          replace_na(get(paste0("bread_w", instance, "t")), 0) + replace_na(get(paste0("cereal_w", instance, "t")), 0) < 3 ~ 0
      )
    ) %>%
    # 计算白面包评分
    mutate(
      !!paste0("whitebread_score_w", instance) := case_when(
        get(paste0("bread_w", instance, "t")) <= 2 & get(paste0("bread_type_w", instance)) == "White" ~ 1,
        get(paste0("bread_w", instance, "t")) > 2 & get(paste0("bread_type_w", instance)) == "White" ~ 0,
        # 只要不吃白面包（选择全麦/棕色/其他），无论多少都直接得1分
        get(paste0("bread_type_w", instance)) == "other" ~ 1,
        # 当面包类型未知时，如果谷物摄入≥3份/天且面包≤2份，得1分
        get(paste0("cereal_w", instance, "t")) >= 3 &
          is.na(get(paste0("bread_type_w", instance))) &
          get(paste0("bread_w", instance, "t")) <= 2 ~ 1,
        # 当面包类型未知时，如果总谷物摄入<3份/天且面包≤2份，得1分
        is.na(get(paste0("bread_type_w", instance))) &
          get(paste0("bread_w", instance, "t")) <= 2 &
          replace_na(get(paste0("bread_w", instance, "t")), 0) + replace_na(get(paste0("cereal_w", instance, "t")), 0) < 3 ~ 1
      )
    )
}

# fish
# ---- 3.4.4 Fish intake ----
# Cleans oily and non-oily fish intake for one UKB instance.
process_fish_data <- function(data, instance) {
  # 构造变量名后缀
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 处理油性鱼
      !!paste0("oily_fish_w", instance) := case_when(
        get(paste0("Oily.fish.intake", suffix)) == "Never" ~ 0,
        get(paste0("Oily.fish.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Oily.fish.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Oily.fish.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Oily.fish.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Oily.fish.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 处理非油性鱼
      !!paste0("non_oily_fish_w", instance) := case_when(
        get(paste0("Non.oily.fish.intake", suffix)) == "Never" ~ 0,
        get(paste0("Non.oily.fish.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Non.oily.fish.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Non.oily.fish.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Non.oily.fish.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Non.oily.fish.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 计算总鱼类摄入量
      !!paste0("fish_w", instance) := ifelse(
        is.na(get(paste0("oily_fish_w", instance))) & is.na(get(paste0("non_oily_fish_w", instance))),
        NA_real_,
        replace_na(get(paste0("oily_fish_w", instance)), 0) +
          replace_na(get(paste0("non_oily_fish_w", instance)), 0)
      ),
      
      # 计算鱼类评分
      !!paste0("fish_score_w", instance) := case_when(
        get(paste0("oily_fish_w", instance)) >= 2 ~ 1,
        get(paste0("non_oily_fish_w", instance)) >= 2 ~ 1,
        get(paste0("fish_w", instance)) >= 2 ~ 1,
        is.na(get(paste0("fish_w", instance))) &
          is.na(get(paste0("oily_fish_w", instance))) &
          is.na(get(paste0("non_oily_fish_w", instance))) ~ NA_real_,
        TRUE ~ 0
      )
    )
}

# processed_meat
# ---- 3.4.5 Processed meat intake ----
# Cleans processed-meat frequency for one UKB instance.
process_processed_meat_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 处理加工肉摄入量
      !!paste0("processed_meat_w", instance) := case_when(
        get(paste0("Processed.meat.intake", suffix)) == "Never" ~ 0,
        get(paste0("Processed.meat.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Processed.meat.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Processed.meat.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Processed.meat.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Processed.meat.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 计算加工肉评分
      !!paste0("processed_meat_score_w", instance) := case_when(
        get(paste0("processed_meat_w", instance)) <= 1 ~ 1,
        get(paste0("processed_meat_w", instance)) > 1 ~ 0,
        TRUE ~ NA_real_
      )
    )
}

# unprocessed_meat
# ---- 3.4.6 Red meat intake ----
# Cleans beef, lamb/mutton, and pork intake for one UKB instance.
process_red_meat_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 处理牛肉摄入
      !!paste0("beef_w", instance) := case_when(
        get(paste0("Beef.intake", suffix)) == "Never" ~ 0,
        get(paste0("Beef.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Beef.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Beef.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Beef.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Beef.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 处理猪肉摄入
      !!paste0("pork_w", instance) := case_when(
        get(paste0("Pork.intake", suffix)) == "Never" ~ 0,
        get(paste0("Pork.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Pork.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Pork.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Pork.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Pork.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 处理羊肉摄入
      !!paste0("lamb_mutton_w", instance) := case_when(
        get(paste0("Lamb.mutton.intake", suffix)) == "Never" ~ 0,
        get(paste0("Lamb.mutton.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Lamb.mutton.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Lamb.mutton.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Lamb.mutton.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Lamb.mutton.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # 计算红肉总摄入量
      !!paste0("red_meat_w", instance) := ifelse(
        is.na(get(paste0("beef_w", instance))) &
          is.na(get(paste0("pork_w", instance))) &
          is.na(get(paste0("lamb_mutton_w", instance))),
        NA_real_,
        round(
          as.numeric(
            replace_na(get(paste0("beef_w", instance)), 0) +
              replace_na(get(paste0("pork_w", instance)), 0) +
              replace_na(get(paste0("lamb_mutton_w", instance)), 0)
          ), 2
        )
      ),
      
      # 计算红肉评分
      !!paste0("red_meat_score_w", instance) := case_when(
        get(paste0("red_meat_w", instance)) <= 2 ~ 1,
        is.na(get(paste0("beef_w", instance))) & is.na(get(paste0("pork_w", instance))) & is.na(get(paste0("lamb_mutton_w", instance))) & is.na(get(paste0("red_meat_w", instance))) ~ NA_real_,
        TRUE ~ 0
      )
    )
}

# Dairy
# ---- 3.4.7 Dairy intake ----
# Cleans milk type and cheese intake for one UKB instance.
process_Dairy_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # cheese 摄入量处理
      !!paste0("cheese_w", instance) := case_when(
        get(paste0("Cheese.intake", suffix)) == "Never" ~ 0,
        get(paste0("Cheese.intake", suffix)) == "Less than once a week" ~ 0.5,
        get(paste0("Cheese.intake", suffix)) == "Once a week" ~ 1,
        get(paste0("Cheese.intake", suffix)) == "2-4 times a week" ~ 3,
        get(paste0("Cheese.intake", suffix)) == "5-6 times a week" ~ 5.5,
        get(paste0("Cheese.intake", suffix)) == "Once or more daily" ~ 7,
        TRUE ~ NA_real_
      ),
      
      # cheese 每日摄入量标准化
      !!paste0("cheese_w", instance, "t") := get(paste0("cheese_w", instance)) / 7,
      
      # milk 摄入量处理
      !!paste0("milk_w", instance) := as.numeric(
        case_when(
          get(paste0("Milk.type.used", suffix)) == "Full cream" ~ 1,
          get(paste0("Milk.type.used", suffix)) == "Semi-skimmed" ~ 1,
          get(paste0("Milk.type.used", suffix)) == "Skimmed" ~ 1,
          get(paste0("Milk.type.used", suffix)) == "Soya" ~ 1,
          get(paste0("Milk.type.used", suffix)) == "Other type of milk" ~ 1,
          get(paste0("Milk.type.used", suffix)) == "Never/rarely have milk" ~ 0,
          get(paste0("Milk.type.used", suffix)) %in% c("Do not know", "Prefer not to answer", "") ~ NA_real_,
          TRUE ~ NA_real_
        )
      ),
      
      # Dairy 总摄入量计算
      !!paste0("Dairy_w", instance) := ifelse(
        is.na(get(paste0("cheese_w", instance, "t"))) &
          is.na(get(paste0("milk_w", instance))),
        NA_real_,
        round(
          as.numeric(
            replace_na(get(paste0("cheese_w", instance, "t")), 0) +
              replace_na(get(paste0("milk_w", instance)), 0)
          ),
          2
        )
      ),
      
      # Dairy 评分
      !!paste0("Dairy_score_w", instance) := case_when(
        is.na(get(paste0("Dairy_w", instance))) ~ NA_real_,
        get(paste0("Dairy_w", instance)) >= 2 ~ 1,
        TRUE ~ 0
      )
    )
}

# vegetable oils
# ---- 3.4.8 Vegetable-oil indicator ----
# Cleans spread/oil-related dietary fields for one UKB instance.
process_veg_oil_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  data %>%
    mutate(
      # 计算每天吃的面包片数
      !!paste0("bread_w", instance) := ifelse(
        get(paste0("Bread.intake", suffix)) == "Less than one", 0.5,
        ifelse(
          get(paste0("Bread.intake", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
          NA,
          get(paste0("Bread.intake", suffix))
        )
      ),
      !!paste0("bread_slices_perday_w", instance) := as.numeric(get(paste0("bread_w", instance))) / 7,
      
      
      # 判断是否为植物油类涂抹酱
      !!paste0("is_veg_oil_spread_w", instance) :=
        ifelse(
          get(paste0("Spread.type", suffix)) %in% c("Do not know", "Prefer not to answer", ""),
          NA_real_,
          ifelse(
            get(paste0("Spread.type", suffix)) == "Butter/spreadable butter",
            0,
            1 # 除黄油之外都是植物油
          )
        ),
      
      # 计算植物油评分
      !!paste0("veg_oil_score_w", instance) := case_when(
        is.na(get(paste0("is_veg_oil_spread_w", instance))) &
          is.na(get(paste0("bread_slices_perday_w", instance))) ~ NA_real_,
        get(paste0("is_veg_oil_spread_w", instance)) == 1 &
          get(paste0("bread_slices_perday_w", instance)) >= 4 ~ 1, # 2 servings
        TRUE ~ 0
      )
    )
}

# sugar
# ---- 3.4.9 Sugar-sweetened beverage / added sugar indicator ----
# Cleans sugar-related dietary variables for one UKB instance.
process_sugar_data <- function(data, instance) {
  suffix <- paste0("...Instance.", instance)
  
  # 需要判断的取值
  sugar_values <- c(
    "Wheat products|Sugar or foods/drinks containing sugar",
    "Sugar or foods/drinks containing sugar",
    "Eggs or foods containing eggs|Wheat products|Sugar or foods/drinks containing sugar",
    "Eggs or foods containing eggs|Sugar or foods/drinks containing sugar",
    "Eggs or foods containing eggs|Dairy products|Wheat products|Sugar or foods/drinks containing sugar",
    "Eggs or foods containing eggs|Dairy products|Sugar or foods/drinks containing sugar",
    "Dairy products|Wheat products|Sugar or foods/drinks containing sugar",
    "Dairy products|Sugar or foods/drinks containing sugar"
  )
  
  data %>%
    mutate(
      !!paste0("sugar_score_w", instance) := case_when(
        get(paste0("Never.eat.eggs..dairy..wheat..sugar", suffix)) %in% sugar_values ~ 1,
        get(paste0("Never.eat.eggs..dairy..wheat..sugar", suffix)) %in% c("Prefer not to answer", "") ~ NA_real_,
        TRUE ~ 0
      )
    )
}

# 定义函数：给定 wave 生成 sumdiet、dietna、diet_score
# ---- 3.4.10 Final diet score for one wave/instance ----
# Combines dietary component indicators into diet_score_w<wave>. The main script
# calls this function for waves 0-3 and then keeps diet_score_w0-diet_score_w3.
calc_diet_vars <- function(data, wave) {
  diet_var <- c(
    "vegetable_score",
    "fruit_score",
    "fish_score",
    "red_meat_score",
    "processed_meat_score",
    "whitebread_score",
    "Wholegrains_score",
    "Dairy_score",
    "sugar_score",
    "veg_oil_score"
  )
  
  diet_var_wave <- paste0(diet_var, "_w", wave)
  
  data %>%
    mutate(
      !!paste0("sumdiet_w", wave) := rowSums(across(all_of(diet_var_wave)), na.rm = TRUE),
      !!paste0("dietna_w", wave) := rowSums(is.na(across(all_of(diet_var_wave)))),
      !!paste0("diet_score_w", wave) := case_when(
        rowSums(!is.na(across(all_of(diet_var_wave)))) == 0 ~ NA_real_,
        .data[[paste0("sumdiet_w", wave)]] >= 5 ~ 1,
        .data[[paste0("sumdiet_w", wave)]] < 5 ~ 0,
        .data[[paste0("sumdiet_w", wave)]] + .data[[paste0("dietna_w", wave)]] <= 3 ~ 0
      )
    )
}


# 4.1 转成长数据集的函数 -----------------------------------------------------------

# ============================================================================== #
# 5.1 Biomarker-name helper
# ============================================================================== #
# Main-script use:
#   Standardize UKB biomarker column names after replacing " | Instance " with
#   an underscore and converting remaining spaces to underscores. This helper is
#   kept here for compatibility with the biomarker-cleaning section.
standardize_variable_names <- function(df) {
  col_names <- colnames(df)
  
  # 遍历所有列名
  for (i in 1:length(col_names)) {
    original_name <- col_names[i]
    
    # 检查是否是时间点变量（以数字结尾但没有下划线分隔）
    if (grepl("^(.*[^0-9_])(0|1|2|3)$", original_name)) {
      # 提取变量名和时间点
      matches <- regexec("^(.*[^0-9_])(0|1|2|3)$", original_name)
      match_result <- regmatches(original_name, matches)[[1]]
      
      if (length(match_result) == 3) {
        base_name <- match_result[2]
        time_point <- match_result[3]
        
        # 构建新的变量名
        new_name <- paste0(base_name, "_", time_point)
        
        # 检查新名称是否已存在（避免重复）
        if (!new_name %in% col_names) {
          col_names[i] <- new_name
          cat("标准化:", original_name, "->", new_name, "\n")
        }
      }
    }
  }
  
  colnames(df) <- col_names
  return(df)
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
rm01_family <- function(d) {
  cs <- paste0("Illnesses of ", c("father", "mother", "siblings"), " | Instance 0")
  positive <- known <- matrix(FALSE, nrow(d), 3L)
  for (j in seq_along(cs)) {
    s <- tolower(as.character(d[[cs[j]]]))
    positive[,j] <- !is.na(s) & grepl("heart disease|stroke", s)
    # Group 2 answers do not establish a negative Group 1 (CVD) answer.
    known[,j] <- !is.na(s) & !grepl("do not know.*group 1|prefer not.*group 1", s) &
      grepl("heart disease|stroke|high blood pressure|bronchitis|alzheimer|diabetes|none of the above.*group 1|^none of the above$", s)
  }
  ifelse(rowSums(positive) > 0L, 1L, ifelse(rowSums(known) == 3L, 0L, NA_integer_))
}
rm01_build_exposure_sources <- function(root) {
  pop <- rm01_read(root,"2_Population characteristics.csv",c("Participant ID","Sex","Age at recruitment"))
  ids <- pop[["Participant ID"]]
  socname <- "3_Assessment centre_2_Touchscreen_2_Sociodemographics.csv"
  soc <- rm01_align(rm01_read(root,socname,c("Participant ID","Date of attending assessment centre | Instance 0","Ethnic background | Instance 0","At or above moderate/vigorous recommendation | Instance 0")),ids)
  out <- data.table::data.table(`Participant ID`=ids)
  ethnic <- soc[["Ethnic background | Instance 0"]]
  out[, rm_source_south_asian := ifelse(ethnic %in% c("Indian","Pakistani","Bangladeshi"),1L,
    ifelse(is.na(ethnic) | ethnic %in% c("Do not know","Prefer not to answer"),NA_integer_,0L))]
  out[, rm_source_inactive := 1L-rm01_binary(soc[["At or above moderate/vigorous recommendation | Instance 0"]])]
  family <- rm01_align(rm01_read(root,"3_Assessment centre_2_Touchscreen_4_Family history.csv",c("Participant ID",paste0("Illnesses of ",c("father","mother","siblings")," | Instance 0"))),ids)
  out[, rm_source_family_proxy := rm01_family(family)]
  stress <- rm01_align(rm01_read(root,"3_Assessment centre_2_Touchscreen_5_Psychosocial factors1.csv",c("Participant ID","Illness, injury, bereavement, stress in last 2 years | Instance 0")),ids)[[2L]]
  valid_stress <- c("Serious illness, injury or assault to yourself","Serious illness, injury or assault of a close relative","Death of a spouse or partner","Death of your spouse or partner","Death of a close relative","Marital separation/divorce","Financial difficulties")
  out[, rm_source_stressor_proxy := vapply(stress,function(s) {
    if (is.na(s)) return(NA_integer_)
    tokens <- strsplit(s,"|",fixed=TRUE)[[1L]]
    if (any(tokens %in% valid_stress)) 1L else if (identical(tokens,"None of the above")) 0L else NA_integer_
  },integer(1))]
  bi <- list()
  for (k in 1:3) {
    name <- paste0("4_Biological samples_1_Blood assays_2_Blood biochemistry",k,".csv")
    path <- rm01_resolve(root,name)
    h <- names(data.table::fread(path,nrows=0L))
    wanted <- intersect(h,paste0(c("C-reactive protein","Lipoprotein A reportability","C-reactive protein reportability","Apolipoprotein B reportability","Triglycerides reportability")," | Instance 0"))
    if (length(wanted)) bi[[k]] <- rm01_align(rm01_read(root,name,c("Participant ID",wanted)),ids)
  }
  flagfields <- c(crp="C-reactive protein",lpa_flag="Lipoprotein A reportability",crp_flag="C-reactive protein reportability",apob_flag="Apolipoprotein B reportability",tg_flag="Triglycerides reportability")
  for (v in names(flagfields)) {
    col <- paste0(flagfields[[v]]," | Instance 0"); d <- Filter(function(d) col %in% names(d),bi)
    if (length(d)!=1L) stop("Missing/ambiguous RM biochemistry: ",col)
    out[, (paste0("rm_source_",v)) := if(v=="crp") rm01_num(d[[1L]][[col]]) else d[[1L]][[col]]]
  }
  sx <- c(menopause="Had menopause",menopause_age="Age at menopause (last menstrual period)",hysterectomy="Ever had hysterectomy (womb removed)",hysterectomy_age="Age at hysterectomy",oophorectomy="Bilateral oophorectomy (both ovaries removed)",oophorectomy_age="Age at bilateral oophorectomy (both ovaries removed)")
  sexd <- rm01_align(rm01_read(root,"3_Assessment centre_2_Touchscreen_7_Sex-specific factors.csv",c("Participant ID",paste0(sx," | Instance 0"))),ids)
  for(v in names(sx)) out[, (paste0("rm_source_",v)) := if(grepl("age$",v)) rm01_num(sexd[[paste0(sx[[v]]," | Instance 0")]]) else rm01_binary(sexd[[paste0(sx[[v]]," | Instance 0")]])]
  natural_invalid <- out$rm_source_hysterectomy==1L &
    (is.na(out$rm_source_hysterectomy_age) | is.na(out$rm_source_menopause_age) |
       out$rm_source_hysterectomy_age <= out$rm_source_menopause_age)
  out[which(natural_invalid), `:=`(rm_source_menopause=NA_integer_,rm_source_menopause_age=NA_real_)]
  out[, rm_definition_version := rm01_version]
  out
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
rm01_manifest <- function() {
  path <- file.path(code_dir,"config","RM_domains_v1.7.csv")
  if (!file.exists(path)) stop("Missing RM domain manifest: ",path)
  x <- data.table::fread(path, encoding="UTF-8")
  required <- c("rm_id","information_type","include_2019","include_2025","rule")
  if (!all(required %in% names(x)) || nrow(x)!=16L || anyDuplicated(x$rm_id))
    stop("Invalid or incomplete RM domains manifest; expected 16 measurable domains")
  for (v in c("include_2019","include_2025")) {
    x[, (v) := as.logical(get(v))]
    if (anyNA(x[[v]])) stop("Invalid RM flag ",v)
  }
  x[, version := rm01_version]
  x
}
rm01_derive <- function(d, cutoff, preimputed_premature_menopause=FALSE) {
  if(length(cutoff)!=1L || !is.finite(cutoff)) stop("Invalid frozen Townsend Q80")
  out <- data.table::copy(d); female <- out$sex=="female"
  manifest <- rm01_manifest()
  for(i in seq_len(nrow(manifest))) {
    value <- eval(parse(text=manifest$rule[i]),envir=out,enclos=environment())
    out[, (paste0("rm_",manifest$rm_id[i])) := as.integer(value)]
  }
  # Reportability can establish Lp(a) status without inventing a concentration.
  flag <- as.character(out$rm_source_lpa_flag)
  low <- !is.na(flag)&(flag %in% c("2","4")|grepl("too low",flag,ignore.case=TRUE))
  high <- !is.na(flag)&(flag %in% c("3","5")|grepl("too high",flag,ignore.case=TRUE))
  out[is.na(rm_lpa)&low,rm_lpa:=0L];out[is.na(rm_lpa)&high,rm_lpa:=1L]
  if("rm_lpa_status" %in% names(out)) out[!is.na(rm_lpa_status),rm_lpa:=as.integer(as.character(rm_lpa_status))]
  early <- (out$rm_source_menopause==1L & out$rm_source_menopause_age<40) |
    (out$rm_source_oophorectomy==1L & out$rm_source_oophorectomy_age<40)
  out[, rm_menopause_uncertain := as.integer(female & is.na(early))]
  out[, rm_premature_menopause := as.integer(female & !is.na(early) & early)]
  out[!female,rm_hdp:=0L]
  for(year in c("2019","2025")) for(domain in c("H","T","any")) {
    use <- manifest[[paste0("include_",year)]] & (domain=="any" | manifest$information_type==domain)
    columns <- paste0("rm_",manifest$rm_id[use]);a<-as.matrix(out[,..columns])
    value<-ifelse(rowSums(a==1L,na.rm=TRUE)>0L,1L,ifelse(rowSums(is.na(a))==0L,0L,NA_integer_))
    out[, (paste0("rm_",tolower(domain),"_",year)) := as.integer(value)]
  }
  out
}
