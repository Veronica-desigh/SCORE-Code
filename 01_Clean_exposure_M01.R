# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 01_Clean_exposure_M01.R
# ------------------------------------------------------------------------------ #
# Purpose:
#   Clean baseline exposure and covariate variables for the M01 ASCVD analysis.
#
# Scope of this revision:
#   1. Preserve the validated baseline cleaning rules and one-row-per-participant contract.
#   2. Keep Education, household income, and employment status as separate SES components;
#      do NOT derive the final SES class in 01. SES is re-derived after MI in 06.
#   3. Preserve raw summed MET-min/week (Instance 0-3) in addition to the legacy
#      physical-activity categories. The baseline raw MET value is the variable to
#      be imputed later; the category is deterministically re-derived after MI.
#   4. Standardize missing values before final output without converting genuine
#      missingness to a reference/negative category.
#   5. Existing validated input paths and core output filenames remain unchanged.
# 模块任务：清理血液生化指标
# 模块输入：基线血液生化数据
# 模块输出：脂质、糖代谢、肌酐等
# ============================================================================== #

# ------------------------------------------------------------------------------ #
# 文件接口契约（M01 2026-10-15 reporting版；不改变临床与统计方法）
# 输入：项目路径及前序步骤输出；请按编号01→06依次运行。
# 输出数据：M01_ASCVD 下标准化CSV/RDS，06完整数据位于MI_M01/MI_completed。
# 输出报告：A=纳排，B=分析基础数据，C=RM专属报告，MI_M01/quality_control=MI通用QC。
# 特别约束：报告不可替代数据集；不将M套MI记录堆叠计算比例。
# ------------------------------------------------------------------------------ #

rm(list = ls())
gc()

# ============================================================================== #
# 0. Environment, project paths, packages, and helper functions
# 模块任务：配置读写路径、加载清理函数并建立基线主表
# 模块输入：UKB路径、原始assessment master、01函数
# 模块输出：受检者ID主表、清理环境
# ============================================================================== #

# ---- 0.1 Project paths ----
project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
data_dir <- Sys.getenv("M01_UKB_DIR", file.path(project_root, "DATA", "UKB"))
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
rm_report_dir <- file.path(data_clean_dir,"风险修饰因素定义与处理记录_M01")
quality_control_dir <- file.path(data_clean_dir, "quality_control")
data_intervention_dir <- file.path(project_root, "LJW_outcome", "data_intervention", "M01_ASCVD")

# 在任何读取、计算或输出前创建本模块的根目录，避免依赖手工建 M01_ASCVD。
for (path in c(data_clean_dir, quality_control_dir, data_intervention_dir, rm_report_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ---- 0.2 Packages ----
library(tidyverse)
library(dplyr)
library(purrr)
library(data.table)


setwd(project_root)

# ---- 0.3 Helper functions ----
source(file.path(code_dir, "Function_01_Clean_exposure_M01.R"))

source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("01", data_clean_dir)

# 单独运行默认12线程；并行入口默认每进程4线程，可由M01_DT_THREADS调整。
m01_configure_threads()
message("data.table 线程数：", data.table::getDTthreads())

# ---- 0.4 Missing-value standardization helpers ----
# Purpose:
#   Standardize common missing-value encodings before final output.
#
# Rules:
#   1) Numeric columns: NaN is converted to NA_real_.
#   2) Integer columns: missing values are kept as typed integer NA.
#   3) Character columns: common UKB/non-standard missing strings are converted
#      to NA_character_ after trimming leading/trailing whitespace.
#   4) Logical and Date columns: missing values are kept as their corresponding
#      typed NA values.
#
# Note:
#   R stores missing values with column-specific types. Therefore, NA_real_ is
#   the correct standardized missing value for numeric columns; character and
#   Date columns must use typed NA values to preserve their classes.
standardize_missing_values <- function(dt) {
  stopifnot(data.table::is.data.table(dt))
  
  missing_strings <- c(
    "",
    "NA",
    "N/A",
    "NaN",
    "NULL",
    "Do not know",
    "Prefer not to answer"
  )
  
  for (col in names(dt)) {
    x <- dt[[col]]
    
    if (is.character(x)) {
      x_trim <- trimws(x)
      idx_missing <- is.na(x_trim) | x_trim %in% missing_strings
      x_trim[idx_missing] <- NA_character_
      data.table::set(dt, j = col, value = x_trim)
    } else if (is.numeric(x)) {
      idx_nan <- which(is.nan(x))
      if (length(idx_nan) > 0) {
        data.table::set(dt, i = idx_nan, j = col, value = NA_real_)
      }
    }
  }
  
  invisible(dt)
}

# Row-wise mean for repeated readings. If all readings are missing, return NA_real_
# rather than NaN.
row_mean_standard_na <- function(...) {
  out <- rowMeans(cbind(...), na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

# Select columns from a data.table/data.frame and return a data.table.
# This avoids the common data.table pitfall that DT[, c("a", "b")] returns
# a character vector rather than selecting columns.
select_cols <- function(dt, cols, object_name = deparse(substitute(dt))) {
  dt <- as.data.table(dt)
  missing_cols <- setdiff(cols, names(dt))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing expected columns in ", object_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dt[, ..cols]
}

# ---- 0.5 Explicit baseline assessment master ----
# Every downstream component is joined to this one-row-per-participant table.
# This prevents the participant set from depending on component availability or
# on the order in which component tables are joined.
assessment_source_file <- file.path(
  data_dir,
  "3_Assessment center/3_Assessment centre_2_Touchscreen_2_Lifestyle and environment1.csv"
)
assessment_raw <- fread(
  assessment_source_file,
  select = c(
    "Participant ID",
    "Date of attending assessment centre | Instance 0",
    "UK Biobank assessment centre | Instance 0"
  )
)

data_assessment <- assessment_raw[, .(
  `Participant ID`,
  InCohort_time = as.Date(
    `Date of attending assessment centre | Instance 0`,
    format = "%Y-%m-%d"
  ),
  Region = fifelse(
    `UK Biobank assessment centre | Instance 0` %in% c(
      "Barts", "Birmingham", "Bristol", "Bury", "Croydon", "Hounslow",
      "Leeds", "Liverpool", "Manchester", "Middlesborough", "Newcastle",
      "Nottingham", "Oxford", "Reading", "Stockport (pilot)", "Stoke",
      "Sheffield"
    ),
    "England",
    fifelse(
      `UK Biobank assessment centre | Instance 0` %in% c(
        "Cardiff", "Swansea", "Wrexham"
      ),
      "Wales",
      fifelse(`UK Biobank assessment centre | Instance 0` %in% c(
        "Edinburgh", "Glasgow"
      ), "Scotland", NA_character_)
    )
  )
)]

assert_one_row_per_id(data_assessment, "baseline assessment master")
baseline_master_ids <- data_assessment[["Participant ID"]]
baseline_master_n <- nrow(data_assessment)
merge_audit <- data.table::data.table()

m01_checkpoint("01_pipeline", 1L, 8L, "baseline assessment master")
# ============================================================================ #
# 1. Demographic variables
# 模块任务：整理年龄、性别、种族及社会经济组成信息
# 模块输入：UKB baseline人口学与问卷
# 模块输出：data_demographic与SES原始组成变量
# ============================================================================ #
# Variables cleaned in this section: Age_at_recruitment, Sex, Ethnic, and baseline assessment information.


# ---- 1.1 Age at recruitment and sex ----
# Source: Population characteristics.
# Output variables: Age_at_recruitment, Sex.

data_demographic <- m01_timed("01_population_csv", fread(file.path(data_dir, "协变量数据集/2_Population characteristics.csv")))
colnames(data_demographic)

data_demographic[, Age_at_recruitment := as.numeric(`Age at recruitment`)]
summary(data_demographic$Age_at_recruitment)

data_demographic <- select_cols(data_demographic, c(
  "Participant ID",
  "Age_at_recruitment",
  "Sex"
))

var_info_print(data_demographic, c(
  "Participant ID",
  "Age_at_recruitment",
  "Sex"
))

data_demographic[, Sex := fcase(
  Sex == "Female", "Female",
  Sex == "Male", "Male",
  default = NA
)]
table(data_demographic$Sex)

# ---- 1.2 Ethnicity, SES components, and physical-activity source table ----
# Source: Touchscreen sociodemographics.
# This section uses baseline touchscreen sociodemographic fields for ethnicity,
# education, household income, current employment status, and physical activity.
# Only the three SES components are cleaned here; the final 2-class SES is NOT
# fitted in 01 because missing SES components will be imputed in 06 before LCA.

data_sociodemographic <- fread(file.path(data_dir, "协变量数据集/3_Assessment centre_2_Touchscreen_2_Sociodemographics.csv"))
colnames(data_sociodemographic)

var_info_print(data_sociodemographic, c(
  "Participant ID",
  "Qualifications | Instance 0",
  "Current employment status | Instance 0",
  "Ethnic background | Instance 0",
  "Ethnic background | Instance 1",
  "Ethnic background | Instance 2",
  "Ethnic background | Instance 3",
  "Summed MET minutes per week for all activity | Instance 0",
  "Summed MET minutes per week for all activity | Instance 1",
  "Summed MET minutes per week for all activity | Instance 2",
  "Summed MET minutes per week for all activity | Instance 3"
))

# Household income is stored in the assessment-centre touchscreen
# sociodemographics file. The filename in the exported UKB folder may or may not
# contain a blank before ".csv", so two equivalent path variants are checked.
income_file_candidates <- c(
  file.path(data_dir, "3_Assessment center", "3_Assessment centre_2_Touchscreen_2_Sociodemographics.csv")
)
income_file <- income_file_candidates[file.exists(income_file_candidates)][1]
if (is.na(income_file)) {
  stop(
    "Cannot find household income source file. Checked: ",
    paste(income_file_candidates, collapse = " ; ")
  )
}

data_income <- fread(
  income_file,
  select = c(
    "Participant ID",
    "Average total household income before tax | Instance 0"
  )
)

var_info_print(data_income, c(
  "Participant ID",
  "Average total household income before tax | Instance 0"
))

data_income[, Household_income := classify_household_income(
  get("Average total household income before tax | Instance 0")
)]
table(data_income$Household_income)

data_income <- select_cols(data_income, c(
  "Participant ID",
  "Household_income"
))

# Individual-level SES component 1: baseline educational attainment.
# Only Instance 0 is used; O levels/GCSEs and CSEs are classified as medium
# education in classify_education().
data_sociodemographic[, Education := classify_education(get("Qualifications | Instance 0"))]
table(data_sociodemographic$Education)

data_sociodemographic[, `:=`(
  Ethnic0 = classify_ethnicity(get("Ethnic background | Instance 0")),
  Ethnic1 = classify_ethnicity(get("Ethnic background | Instance 1")),
  Ethnic2 = classify_ethnicity(get("Ethnic background | Instance 2")),
  Ethnic3 = classify_ethnicity(get("Ethnic background | Instance 3"))
)]

# Ethnic is a baseline covariate and must come from Instance 0 only. Follow-up
# instance variables are retained under their existing names but never backfill
# the baseline Ethnic value.
data_sociodemographic[, Ethnic := Ethnic0]
table(data_sociodemographic$Ethnic)

table(data_sociodemographic$`Current employment status | Instance 0`)

# Individual-level SES component 2: baseline current employment status.
# Only Instance 0 is used. This is not an occupation code, so the derived
# variable is named employment_status rather than Occupation.
data_sociodemographic[, employment_status := classify_employment_status(
  get("Current employment status | Instance 0")
)]
table(data_sociodemographic$employment_status)

# Add household income to the touchscreen sociodemographics component using an
# audited left join anchored on the sociodemographics participant set.
join_result <- audited_left_join(
  master = data_sociodemographic,
  component = data_income,
  component_name = "household_income_to_sociodemographic",
  audit_log = merge_audit
)
data_sociodemographic <- join_result$data
merge_audit <- join_result$audit

# Individual-level SES component 3 is baseline household income, cleaned above.
#
# IMPORTANT:
#   01 deliberately stops at the three observed SES components:
#     - Education
#     - Household_income
#     - employment_status
#   The final Low/High SES class is fitted as a 2-class LCA only after multiple
#   imputation in 06. This prevents participants with incomplete SES components
#   from being excluded or assigned SES before missing-data handling.

m01_checkpoint("01_pipeline", 2L, 8L, "demographic source ready")
# ============================================================================ #
# 2. Socioeconomic variables
# 模块任务：清理Townsend deprivation index
# 模块输入：UKB recruitment deprivation字段
# 模块输出：连续TDI
# ============================================================================ #
# Variables cleaned in this section: Townsend deprivation index.
# Education, household income, and employment status are cleaned above because
# they share touchscreen sociodemographics source tables. Final SES is deferred to 06.


# ---- 2.1 Townsend deprivation index ----
# Source: Population characteristics.
# Output variable: TDI.
# Note: TDI remains the output variable name to preserve downstream compatibility.

# 更新说明：本轮删除原先按 Region 选择的 region-specific Index of Multiple
# Deprivation，改用 UKB 的 Townsend deprivation index at recruitment。
data_tdi <- fread(
  file.path(data_dir, "2_Population characteristics.csv"),
  select = c(
    "Participant ID",
    "Townsend deprivation index at recruitment"
  )
)

data_tdi[, TDI := suppressWarnings(as.numeric(`Townsend deprivation index at recruitment`))]

data_tdi <- select_cols(data_tdi, c(
  "Participant ID",
  "TDI"
))

m01_checkpoint("01_pipeline", 3L, 8L, "socioeconomic components ready")
# ============================================================================ #
# 3. Lifestyle variables
# 模块任务：清理运动、吸烟、饮酒及饮食
# 模块输入：基线问卷Instance 0
# 模块输出：生活方式协变量，保留缺失
# ============================================================================ #
# Variables cleaned in this section: raw MET-min/week, PhysicalActivity0-3, Smoke0-3,
# alcohol_score_0-3, and diet_score_w0-w3.


# ---- 3.1 Physical activity ----
# Source: Touchscreen sociodemographics.
#
# Two representations are intentionally retained:
#   1. MET_minutes_week0-3: raw summed MET-min/week. Baseline MET_minutes_week0 is
#      the quantity used by the future MI module (06) when missing.
#   2. PhysicalActivity0-3: deterministic category derived from the raw MET value.
#      These categories are retained for QC/backward compatibility only. The
#      formal post-MI category must be re-derived from imputed raw MET in 06.
#
# Missing raw MET remains NA; it must never be converted to a physical-activity
# reference category before MI.

data_sociodemographic[, `:=`(
  MET_minutes_week0 = suppressWarnings(as.numeric(get("Summed MET minutes per week for all activity | Instance 0"))),
  MET_minutes_week1 = suppressWarnings(as.numeric(get("Summed MET minutes per week for all activity | Instance 1"))),
  MET_minutes_week2 = suppressWarnings(as.numeric(get("Summed MET minutes per week for all activity | Instance 2"))),
  MET_minutes_week3 = suppressWarnings(as.numeric(get("Summed MET minutes per week for all activity | Instance 3")))
)]

data_sociodemographic[, `:=`(
  PhysicalActivity0 = classify_physical_activity(MET_minutes_week0),
  PhysicalActivity1 = classify_physical_activity(MET_minutes_week1),
  PhysicalActivity2 = classify_physical_activity(MET_minutes_week2),
  PhysicalActivity3 = classify_physical_activity(MET_minutes_week3)
)]

var_info_print(data_sociodemographic, c(
  "MET_minutes_week0", "MET_minutes_week1", "MET_minutes_week2", "MET_minutes_week3",
  "PhysicalActivity0", "PhysicalActivity1", "PhysicalActivity2", "PhysicalActivity3"
))

data_sociodemographic <- select_cols(data_sociodemographic, c(
  "Participant ID",
  "Education",
  "Household_income",
  "employment_status",
  "Ethnic",
  "Ethnic0",
  "Ethnic1",
  "Ethnic2",
  "Ethnic3",
  "MET_minutes_week0",
  "MET_minutes_week1",
  "MET_minutes_week2",
  "MET_minutes_week3",
  "PhysicalActivity0",
  "PhysicalActivity1",
  "PhysicalActivity2",
  "PhysicalActivity3"
))

# ---- 3.2 Smoking status ----
# Source: Touchscreen lifestyle and environment.
# Output variables: Smoke0, Smoke1, Smoke2, Smoke3.

data_smoking <- fread(file.path(data_dir, "3_Assessment center/3_Assessment centre_2_Touchscreen_2_Lifestyle and environment2.csv"))
colnames(data_smoking)

var_info_print(data_smoking, c(
  "Participant ID",
  "Smoking status | Instance 0",
  "Smoking status | Instance 1",
  "Smoking status | Instance 2",
  "Smoking status | Instance 3"
))

data_smoking[, `:=`(
  Smoke0 = get("Smoking status | Instance 0"),
  Smoke1 = get("Smoking status | Instance 1"),
  Smoke2 = get("Smoking status | Instance 2"),
  Smoke3 = get("Smoking status | Instance 3")
)]

data_smoking <- select_cols(data_smoking, c(
  "Participant ID",
  "Smoke0",
  "Smoke1",
  "Smoke2",
  "Smoke3"
))

# ---- 3.3 Alcohol score ----
# Source: Alcohol frequency and alcohol intake fields.
# Output variables: alcohol_score_0, alcohol_score_1, alcohol_score_2, alcohol_score_3.

drink_month <- fread(
  file.path(data_dir, "3_Assessment center/3_Assessment centre_2_Cognitive function_12_Pairs matching.csv"),
  select = c(
    "Participant ID",
    # 4 次随访的饮酒频率
    "Alcohol intake frequency. | Instance 0",
    "Alcohol intake frequency. | Instance 1",
    "Alcohol intake frequency. | Instance 2",
    "Alcohol intake frequency. | Instance 3",
    # 4 次随访的每月饮酒量
    "Average monthly red wine intake | Instance 0",
    "Average monthly red wine intake | Instance 1",
    "Average monthly red wine intake | Instance 2",
    "Average monthly red wine intake | Instance 3",
    "Average monthly champagne plus white wine intake | Instance 0",
    "Average monthly champagne plus white wine intake | Instance 1",
    "Average monthly champagne plus white wine intake | Instance 2",
    "Average monthly champagne plus white wine intake | Instance 3",
    "Average monthly beer plus cider intake | Instance 0",
    "Average monthly beer plus cider intake | Instance 1",
    "Average monthly beer plus cider intake | Instance 2",
    "Average monthly beer plus cider intake | Instance 3",
    "Average monthly spirits intake | Instance 0",
    "Average monthly spirits intake | Instance 1",
    "Average monthly spirits intake | Instance 2",
    "Average monthly spirits intake | Instance 3",
    "Average monthly fortified wine intake | Instance 0",
    "Average monthly fortified wine intake | Instance 1",
    "Average monthly fortified wine intake | Instance 2",
    "Average monthly fortified wine intake | Instance 3"
  )
)

drink_week <- fread(
  file.path(data_dir, "3_Assessment center/3_Assessment centre_2_Touchscreen_4_Family history.csv"),
  select = c(
    "Participant ID",
    # 0-3 次随访的各类酒每周摄入量
    "Average weekly red wine intake | Instance 0",
    "Average weekly red wine intake | Instance 1",
    "Average weekly red wine intake | Instance 2",
    "Average weekly red wine intake | Instance 3",
    "Average weekly champagne plus white wine intake | Instance 0",
    "Average weekly champagne plus white wine intake | Instance 1",
    "Average weekly champagne plus white wine intake | Instance 2",
    "Average weekly champagne plus white wine intake | Instance 3",
    "Average weekly beer plus cider intake | Instance 0",
    "Average weekly beer plus cider intake | Instance 1",
    "Average weekly beer plus cider intake | Instance 2",
    "Average weekly beer plus cider intake | Instance 3",
    "Average weekly spirits intake | Instance 0",
    "Average weekly spirits intake | Instance 1",
    "Average weekly spirits intake | Instance 2",
    "Average weekly spirits intake | Instance 3",
    "Average weekly fortified wine intake | Instance 0",
    "Average weekly fortified wine intake | Instance 1",
    "Average weekly fortified wine intake | Instance 2",
    "Average weekly fortified wine intake | Instance 3"
  )
)

# 合并 0–3 次所有饮酒原始信息。两张来源表均左连接到固定基线主表，
# 缺少任一来源只形成缺失值，不得删除参与者。
join_result <- audited_left_join(
  master = data_assessment[, .(`Participant ID`)],
  component = drink_month,
  component_name = "alcohol_monthly_intake",
  audit_log = merge_audit
)
drinking_wide_raw <- join_result$data
merge_audit <- join_result$audit

join_result <- audited_left_join(
  master = drinking_wide_raw,
  component = drink_week,
  component_name = "alcohol_weekly_intake",
  audit_log = merge_audit
)
merge_audit <- join_result$audit
drinking_wide_raw <- join_result$data %>%
  mutate(
    across(
      .cols = -`Participant ID`,
      .fns  = ~ ifelse(.x %in% c("", "Do not know", "Prefer not to answer"), NA, .x)
    )
  )

# 基础 ID 列沿用固定基线主表顺序。
alcohol_score_wide <- data.table::as.data.table(drinking_wide_raw)[
  , .(`Participant ID`)
]

# 依次左连接 0、1、2、3 次的得分，并记录每次连接审计。
for (k in 0:3) {
  join_result <- audited_left_join(
    master = alcohol_score_wide,
    component = alcohol_fun_one_instance(drinking_wide_raw, k),
    component_name = paste0("alcohol_score_instance_", k),
    audit_log = merge_audit
  )
  alcohol_score_wide <- join_result$data
  merge_audit <- join_result$audit
}

# 只保留最后二分类的变量
data_alcohol <- select_cols(alcohol_score_wide, c(
  "Participant ID",
  "alcohol_score_0",
  "alcohol_score_1",
  "alcohol_score_2",
  "alcohol_score_3"
))

# ---- 3.4 Diet score ----
# Source: Touchscreen dietary fields.
# Output variables: diet_score_w0, diet_score_w1, diet_score_w2, diet_score_w3.

diet <- m01_timed("01_diet_read_and_40_component_transformations",
  data.table::fread(file.path(data_dir, "3_Assessment center/3_Assessment centre_2_Cognitive function_12_Pairs matching.csv"),
                    data.table = FALSE, check.names = TRUE) %>%
  process_fruit_data(0) %>%
  process_fruit_data(1) %>%
  process_fruit_data(2) %>%
  process_fruit_data(3) %>%
  process_vegetable_data(0) %>%
  process_vegetable_data(1) %>%
  process_vegetable_data(2) %>%
  process_vegetable_data(3) %>%
  process_whole_Refined_grains_data(0) %>%
  process_whole_Refined_grains_data(1) %>%
  process_whole_Refined_grains_data(2) %>%
  process_whole_Refined_grains_data(3) %>%
  process_fish_data(0) %>%
  process_fish_data(1) %>%
  process_fish_data(2) %>%
  process_fish_data(3) %>%
  process_processed_meat_data(0) %>%
  process_processed_meat_data(1) %>%
  process_processed_meat_data(2) %>%
  process_processed_meat_data(3) %>%
  process_red_meat_data(0) %>%
  process_red_meat_data(1) %>%
  process_red_meat_data(2) %>%
  process_red_meat_data(3) %>%
  process_Dairy_data(0) %>%
  process_Dairy_data(1) %>%
  process_Dairy_data(2) %>%
  process_Dairy_data(3) %>%
  process_veg_oil_data(0) %>%
  process_veg_oil_data(1) %>%
  process_veg_oil_data(2) %>%
  process_veg_oil_data(3) %>%
  process_sugar_data(0) %>%
  process_sugar_data(1) %>%
  process_sugar_data(2) %>%
  process_sugar_data(3))
m01_checkpoint("01_pipeline", 4L, 8L, "diet variables ready")

# 膳食得分汇总
# 对 w0,w1,w2,w3 依次生成变量
diet_all <- diet %>%
  select(
    Participant.ID,
    paste0(c(
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
    ), "_w0"),
    paste0(c(
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
    ), "_w1"),
    paste0(c(
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
    ), "_w2"),
    paste0(c(
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
    ), "_w3")
  ) %>%
  calc_diet_vars(0) %>%
  calc_diet_vars(1) %>%
  calc_diet_vars(2) %>%
  calc_diet_vars(3)

data_diet <- select_cols(diet_all, c(
  "Participant.ID",
  "diet_score_w0",
  "diet_score_w1",
  "diet_score_w2",
  "diet_score_w3"
))
setnames(data_diet, "Participant.ID", "Participant ID")

# ============================================================================ #
# 4. Clinical measurements
# 模块任务：清理体格及血压测量
# 模块输入：基线体测、重复测量
# 模块输出：BMI、腰围、血压
# ============================================================================ #
# Variables cleaned in this section: anthropometric measures and blood pressure.


# ---- 4.1 Anthropometric measures ----
# Output variables: WaistCircumference0-3, Weight0-3, BMI0-3, HipCircumference0-3, Height0-3.

data_anthropometry <- fread(file.path(data_dir, "3_Assessment center/3_Assessment centre_5_Physical measures_6_ Anthropometry_5_Body size measures1.csv"))
colnames(data_anthropometry)

var_info_print(data_anthropometry, c(
  "Participant ID",
  "Waist circumference | Instance 0",
  "Waist circumference | Instance 1",
  "Waist circumference | Instance 2",
  "Waist circumference | Instance 3",
  "Weight | Instance 0",
  "Weight | Instance 1",
  "Weight | Instance 2",
  "Weight | Instance 3",
  "Body mass index (BMI) | Instance 0",
  "Body mass index (BMI) | Instance 1",
  "Body mass index (BMI) | Instance 2",
  "Body mass index (BMI) | Instance 3",
  "Hip circumference | Instance 0",
  "Hip circumference | Instance 1",
  "Hip circumference | Instance 2",
  "Hip circumference | Instance 3",
  "Standing height | Instance 0",
  "Standing height | Instance 1",
  "Standing height | Instance 2",
  "Standing height | Instance 3"
))

data_anthropometry[, `:=`(
  WaistCircumference0 = get("Waist circumference | Instance 0"),
  WaistCircumference1 = get("Waist circumference | Instance 1"),
  WaistCircumference2 = get("Waist circumference | Instance 2"),
  WaistCircumference3 = get("Waist circumference | Instance 3")
)]

data_anthropometry[, `:=`(
  Weight0 = get("Weight | Instance 0"),
  Weight1 = get("Weight | Instance 1"),
  Weight2 = get("Weight | Instance 2"),
  Weight3 = get("Weight | Instance 3")
)]

data_anthropometry[, `:=`(
  BMI0 = get("Body mass index (BMI) | Instance 0"),
  BMI1 = get("Body mass index (BMI) | Instance 1"),
  BMI2 = get("Body mass index (BMI) | Instance 2"),
  BMI3 = get("Body mass index (BMI) | Instance 3")
)]

data_anthropometry[, `:=`(
  HipCircumference0 = get("Hip circumference | Instance 0"),
  HipCircumference1 = get("Hip circumference | Instance 1"),
  HipCircumference2 = get("Hip circumference | Instance 2"),
  HipCircumference3 = get("Hip circumference | Instance 3")
)]

data_anthropometry[, `:=`(
  Height0 = get("Standing height | Instance 0"),
  Height1 = get("Standing height | Instance 1"),
  Height2 = get("Standing height | Instance 2"),
  Height3 = get("Standing height | Instance 3")
)]

data_anthropometry <- select_cols(data_anthropometry, c(
  "Participant ID",
  "WaistCircumference0",
  "WaistCircumference1",
  "WaistCircumference2",
  "WaistCircumference3",
  "Weight0",
  "Weight1",
  "Weight2",
  "Weight3",
  "BMI0",
  "BMI1",
  "BMI2",
  "BMI3",
  "HipCircumference0",
  "HipCircumference1",
  "HipCircumference2",
  "HipCircumference3",
  "Height0",
  "Height1",
  "Height2",
  "Height3"
))

# ---- 4.2 Blood pressure ----
# Output variables: DiastolicBP0-3 and SystolicBP0-3.
# Rule: automated array 0 and array 1 readings are averaged within each instance.

data_blood_pressure <- fread(file.path(data_dir, "3_Assessment center/3_Assessment centre_5_Physical measures_4_Hearing test1.csv"))
colnames(data_blood_pressure)

var_info_print(data_blood_pressure, c(
  "Participant ID",
  "Diastolic blood pressure, automated reading | Instance 0 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 0 | Array 1",
  "Diastolic blood pressure, automated reading | Instance 1 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 1 | Array 1",
  "Diastolic blood pressure, automated reading | Instance 2 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 2 | Array 1",
  "Diastolic blood pressure, automated reading | Instance 3 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 3 | Array 1",
  "Systolic blood pressure, automated reading | Instance 0 | Array 0",
  "Systolic blood pressure, automated reading | Instance 1 | Array 0",
  "Systolic blood pressure, automated reading | Instance 0 | Array 1",
  "Systolic blood pressure, automated reading | Instance 1 | Array 1",
  "Systolic blood pressure, automated reading | Instance 2 | Array 0",
  "Systolic blood pressure, automated reading | Instance 2 | Array 1",
  "Systolic blood pressure, automated reading | Instance 3 | Array 0",
  "Systolic blood pressure, automated reading | Instance 3 | Array 1"
))

data_blood_pressure[, `:=`(
  # 舒张压
  DiastolicBP0 = row_mean_standard_na(
    get("Diastolic blood pressure, automated reading | Instance 0 | Array 0"),
    get("Diastolic blood pressure, automated reading | Instance 0 | Array 1")
  ),
  DiastolicBP1 = row_mean_standard_na(
    get("Diastolic blood pressure, automated reading | Instance 1 | Array 0"),
    get("Diastolic blood pressure, automated reading | Instance 1 | Array 1")
  ),
  DiastolicBP2 = row_mean_standard_na(
    get("Diastolic blood pressure, automated reading | Instance 2 | Array 0"),
    get("Diastolic blood pressure, automated reading | Instance 2 | Array 1")
  ),
  DiastolicBP3 = row_mean_standard_na(
    get("Diastolic blood pressure, automated reading | Instance 3 | Array 0"),
    get("Diastolic blood pressure, automated reading | Instance 3 | Array 1")
  ),
  
  # 收缩压
  SystolicBP0 = row_mean_standard_na(
    get("Systolic blood pressure, automated reading | Instance 0 | Array 0"),
    get("Systolic blood pressure, automated reading | Instance 0 | Array 1")
  ),
  SystolicBP1 = row_mean_standard_na(
    get("Systolic blood pressure, automated reading | Instance 1 | Array 0"),
    get("Systolic blood pressure, automated reading | Instance 1 | Array 1")
  ),
  SystolicBP2 = row_mean_standard_na(
    get("Systolic blood pressure, automated reading | Instance 2 | Array 0"),
    get("Systolic blood pressure, automated reading | Instance 2 | Array 1")
  ),
  SystolicBP3 = row_mean_standard_na(
    get("Systolic blood pressure, automated reading | Instance 3 | Array 0"),
    get("Systolic blood pressure, automated reading | Instance 3 | Array 1")
  )
)]

data_blood_pressure <- select_cols(data_blood_pressure, c(
  "Participant ID",
  "DiastolicBP0",
  "DiastolicBP1",
  "DiastolicBP2",
  "DiastolicBP3",
  "SystolicBP0",
  "SystolicBP1",
  "SystolicBP2",
  "SystolicBP3"
))

m01_checkpoint("01_pipeline", 5L, 8L, "clinical measurements ready")
# ============================================================================ #
# 5. Blood biomarkers
# ============================================================================ #
# Variables cleaned in this section: lipids, apolipoproteins, glucose/HbA1c, and creatinine.


# ---- 5.1 Lipid, glucose, HbA1c, and creatinine variables ----
# Source: Blood biochemistry files.
# Output variables are generated by standardizing UKB instance field names.

blood_biochemistry_1 <- m01_timed("01_blood_lipids_part2", fread(file.path(data_dir, "4_Biological sample/4_1_2Blood biochemistry/4_Biological samples_1_Blood assays_2_Blood biochemistry2.csv")))
blood_biochemistry_2 <- m01_timed("01_blood_lipids_part3", fread(file.path(data_dir, "4_Biological sample/4_1_2Blood biochemistry/4_Biological samples_1_Blood assays_2_Blood biochemistry3.csv")))

join_result <- audited_left_join(
  master = data_assessment[, .(`Participant ID`)],
  component = blood_biochemistry_1,
  component_name = "blood_biochemistry_part_1",
  audit_log = merge_audit
)
data_biomarkers <- join_result$data
merge_audit <- join_result$audit

join_result <- audited_left_join(
  master = data_biomarkers,
  component = blood_biochemistry_2,
  component_name = "blood_biochemistry_part_2",
  audit_log = merge_audit
)
data_biomarkers <- join_result$data
merge_audit <- join_result$audit
colnames(data_biomarkers)

Blood_vars <- c(
  "Participant ID",
  
  # 1. 血脂相关
  "Cholesterol | Instance 0",
  "Cholesterol | Instance 1",
  "LDL direct | Instance 0",
  "LDL direct | Instance 1",
  "HDL cholesterol | Instance 0",
  "HDL cholesterol | Instance 1",
  "Triglycerides | Instance 0",
  "Triglycerides | Instance 1",
  
  # 2. 载脂蛋白相关
  "Lipoprotein A | Instance 0",
  "Lipoprotein A | Instance 1",
  "Apolipoprotein A | Instance 0",
  "Apolipoprotein A | Instance 1",
  "Apolipoprotein B | Instance 0",
  "Apolipoprotein B | Instance 1",
  
  # 3. 血糖相关
  "Glucose | Instance 0",
  "Glucose | Instance 1",
  "Glycated haemoglobin (HbA1c) | Instance 0",
  "Glycated haemoglobin (HbA1c) | Instance 1",
  
  # 4. CKD相关
  "Creatinine | Instance 0",
  "Creatinine | Instance 1"
)

var_info_print(data_biomarkers, Blood_vars) # 查看原始值缺失情况

data_biomarkers <- data_biomarkers[, ..Blood_vars]

# 批量重命名：先替换" | Instance "为"_"，再替换剩余的空格为"_"
new_names <- gsub(" ", "_", gsub(" | Instance ", "_", Blood_vars, fixed = TRUE), fixed = TRUE)
setnames(data_biomarkers, old = Blood_vars, new = new_names)

setnames(data_biomarkers, "Participant_ID", "Participant ID")

m01_checkpoint("01_pipeline", 6L, 8L, "blood biochemistry ready")
# ============================================================================ #
# 6. Disease history, family history, and baseline assessment information
# 模块任务：整理家族史和评估日期
# 模块输入：UKB基线家族史及评估
# 模块输出：家族史与基线时间
# ============================================================================ #
# Variables cleaned in this section: Family_history_CVD, InCohort_time, and Region.


# ---- 6.1 Family history of CVD ----
# Source: Family history fields.
# Output variable: Family_history_CVD.

family_history_raw <- fread(
  file.path(data_dir, "3_Assessment center/3_Assessment centre_2_Touchscreen_4_Family history.csv"),
  select = c(
    "Participant ID",
    "Illnesses of father | Instance 0",
    "Illnesses of mother | Instance 0",
    "Illnesses of siblings | Instance 0"
  )
)

family_history_raw[, Family_history_CVD := classify_family_history_cvd(
  father   = `Illnesses of father | Instance 0`,
  mother   = `Illnesses of mother | Instance 0`,
  siblings = `Illnesses of siblings | Instance 0`
)]

data_family_history <- select_cols(family_history_raw, c(
  "Participant ID",
  "Family_history_CVD"
))

# ---- 6.2 Assessment date and assessment-centre region ----
# Source: Baseline assessment date and assessment centre.
# Output variables: InCohort_time, Region.
# data_assessment was created and validated in section 0.5 because it is the
# explicit master table for every component merge.

m01_checkpoint("01_pipeline", 7L, 8L, "family history ready")
# ============================================================================ #
# 7. Final merge, variable check, and output
# 模块任务：按ID汇总并验证输出
# 模块输入：以上全部基线组件
# 模块输出：M01_ASCVD/data_exposure_M01.csv及QA
# ============================================================================ #
# All cleaned components are left-joined to the explicit baseline assessment
# master. Component order may change column order, but can never change the
# participant set or row order.
# Output file: data_exposure_M01.csv.

exposure_components <- list(
  population_demographic = data_demographic,
  touchscreen_sociodemographic = data_sociodemographic,
  anthropometry = data_anthropometry,
  blood_pressure = data_blood_pressure,
  smoking = data_smoking,
  alcohol = data_alcohol,
  townsend_deprivation = data_tdi,
  family_history = data_family_history,
  diet = data_diet,
  blood_biomarkers = data_biomarkers
)

data_exposure <- data.table::copy(data_assessment)
for (component_name in names(exposure_components)) {
  join_result <- audited_left_join(
    master = data_exposure,
    component = exposure_components[[component_name]],
    component_name = component_name,
    audit_log = merge_audit
  )
  data_exposure <- join_result$data
  merge_audit <- join_result$audit
}

# ============================================================================== #
# RM基线源变量：族裔、活动、家族史、生活事件、生殖史、CRP和尿液指标
# 输入：路径工作簿“字段路径字典”的完整本地CSV路径；疾病史由02处理。
# 输出：并入data_exposure_M01.csv；连续量及检测限标记保留至06。
# ============================================================================== #
rm_workbook <- Sys.getenv("M01_RM_PATH_WORKBOOK", file.path(project_root,
  "LJW_code", "Risk modifiers修改方案", "RM的指南定义和数据定义", "RM_2019_2025_缺失核查_20261010.xlsx"))
rm_path_registry <- rm01_configure_paths(rm_workbook, root = data_dir)
print(rm_path_registry)
rm_exposure_sources <- rm01_build_exposure_sources(data_dir)
data_exposure <- rm01_attach(data_exposure, rm_exposure_sources, "Participant ID")
fwrite(rm_path_registry, file.path(rm_report_dir,"RM_source_path_audit_M01.csv"),bom=TRUE)
fwrite(data.table(workbook=if(file.exists(rm_workbook)) normalizePath(rm_workbook,winslash="/") else NA_character_,
  md5=if(file.exists(rm_workbook)) unname(tools::md5sum(rm_workbook)) else NA_character_),
  file.path(rm_report_dir,"RM_path_workbook_provenance_M01.csv"),bom=TRUE)
fwrite(rm01_manifest(),file.path(rm_report_dir,"RM_definition_manifest_M01.csv"),bom=TRUE)
fwrite(data.table(variable=names(rm_exposure_sources),n_total=nrow(rm_exposure_sources),
  missing=vapply(rm_exposure_sources,function(x)sum(is.na(x)),integer(1)),
  pct_missing=vapply(rm_exposure_sources,function(x)100*mean(is.na(x)),numeric(1))),
  file.path(rm_report_dir,"RM_exposure_source_missingness_M01.csv"),bom=TRUE)

# Standardize missing values across the final exposure dataset.
# This safeguards against both numeric NaN values and non-standard character
# missing strings before writing the cleaned exposure file.
standardize_missing_values(data_exposure)

assert_one_row_per_id(data_exposure, "final exposure data")
if (nrow(data_exposure) != baseline_master_n ||
    !identical(data_exposure[["Participant ID"]], baseline_master_ids)) {
  stop(
    "Final exposure data do not preserve the baseline assessment master cohort.",
    call. = FALSE
  )
}

assert_column_types(data_exposure, c(
  Age_at_recruitment = "numeric",
  Sex = "character",
  Ethnic = "character",
  InCohort_time = "Date",
  Region = "character",
  TDI = "numeric",
  BMI0 = "numeric",
  DiastolicBP0 = "numeric",
  SystolicBP0 = "numeric"
))

if (any(merge_audit$master_n_before != merge_audit$master_n_after)) {
  stop("At least one audited merge changed the master row count.", call. = FALSE)
}

data.table::fwrite(
  merge_audit,
  file.path(quality_control_dir, "merge_audit_M01.csv")
)

m01_timed("01_write_exposure_master_csv",
  data.table::fwrite(data_exposure, file.path(data_clean_dir, "data_exposure_M01.csv")))
m01_checkpoint("01_pipeline", 8L, 8L, "exposure master saved")

m01_runtime_finish()
