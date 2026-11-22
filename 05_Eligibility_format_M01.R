# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 05_Eligibility_format_M01.R
# ------------------------------------------------------------------------------ #
# 用途：
#   基于 04 输出的一人一行最终 analytic cohort，建立“允许普通分析协变量缺失”
#   的 M01 ASCVD analysis base，为两条后续支路提供共同起点：
#     A. 07 observed SCORE/SCORE2 + guideline-strategy branch；
#     B. 06 multiple-imputation branch。
#
# 本次重构的核心原则：
#   1. 05 不再执行 complete-case deletion；样本量和 participant_id 顺序必须与 04
#      完全一致。
#   2. 只有 eligibility / guideline-strategy required variables 必须完整；
#      outcome-model / subgroup covariates 可缺失，并由变量字典标记是否在 06 MI。
#   3. 05 不生成最终 SES、不填补 physical activity、不生成最终 hypertension。
#      SES 在 06 completed datasets 中由 Education + income + employment 直接拟合
#      2-class LCA；physical activity 由 imputed raw MET 被动派生；hypertension 由
#      observed SBP/DBP + imputed antihypertensive medication 被动派生。
#   4. 结局直接复用 02 的冻结 competing-risk 定义；05 仅标准化接口，并构造
#      cause-specific Nelson-Aalen cumulative hazards 作为 MICE auxiliary predictors。
#
# 主要输出：
#   - data_analysis_base_M01.csv                  （兼容既有路径；允许 MI 变量缺失）
#   - data_analysis_base_incomplete_M01.rds       （推荐 06 读取，保留 R 类型）
#   - 分析基础数据与随访审计_M01/variable_dictionary_M01.csv
#   - B. .../analysis_base_missingness_M01.csv
#   - B. .../mi_variable_manifest_M01.csv
#   - B. .../post_mi_derivation_manifest_M01.csv
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
# 0. 环境、路径、R 包和函数 ####
# 模块任务：确定分析基础数据输入输出路径
# 模块输入：04合格样本CSV，05 helper
# 模块输出：B目录字典和M01数据根目录
# ============================================================================== #

project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
description_dir <- file.path(
  data_clean_dir,
  "分析基础数据与随访审计_M01"
)

input_path <- file.path(data_clean_dir, "data_eligible_M01.csv")
output_csv_path <- file.path(data_clean_dir, "data_analysis_base_M01.csv")
output_rds_path <- file.path(data_clean_dir, "data_analysis_base_incomplete_M01.rds")
rm_report_dir <- file.path(data_clean_dir,"风险修饰因素定义与处理记录_M01")


# 仅供开发阶段端到端验证。正式运行默认 FALSE。
validation_mode <- identical(Sys.getenv("M01_VALIDATION_MODE", unset = "0"), "1")
if (validation_mode) {
  validation_root <- file.path(tempdir(), "M01_05_validation")
  description_dir <- file.path(validation_root, "分析基础数据与随访审计_M01")
  output_csv_path <- file.path(validation_root, "data_analysis_base_M01.csv")
  output_rds_path <- file.path(validation_root, "data_analysis_base_incomplete_M01.rds")
  rm_report_dir <- file.path(validation_root,"风险修饰因素定义与处理记录_M01")
}

dir.create(data_clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(description_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rm_report_dir, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(parallel)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("05 需要 R 包 survival 来构造 MI Nelson-Aalen auxiliary predictors。", call. = FALSE)
}

setwd(project_root)
source(file.path(code_dir, "Function_05_Eligibility_format_M01.R"))
source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("05", data_clean_dir)

run_analysis_base_unit_tests()

m01_configure_threads()

step_message("data.table 线程数：", data.table::getDTthreads())
step_message("运行模式：", if (validation_mode) "validation" else "formal")


# ============================================================================== #
# 1. 纳排审计规则和变量角色字典 ####
# 模块任务：定义变量角色与字典
# 模块输入：原始UKB字段映射
# 模块输出：required/impute/derive变量字典
# ============================================================================== #

# 04 已完成临床纳排和 required-information gate；这些变量在最终队列中应为常数。
# 注意：不再存在 baseline_covariates_complete。
expected_eligibility_values <- list(
  documented_ascvd_baseline = 0L,
  diabetes_baseline = 0L,
  CKD_baseline = 0L,
  Lipid_Medicine_baseline = 0L,
  markedly_elevated_single_risk_factor = 0L,
  required_observed_information_complete = 1L,
  lipid_main_analysis_included = 1L
)

# ------------------------------------------------------------------------------
# 变量字典说明
#
# required_observed:
#   TRUE = 04 final gate 后必须无缺失；通常属于 eligibility / strategy / outcome
#          interface 的硬条件。
#
# impute:
#   TRUE = 06 对该“原始分析变量”进行 MI。
#
# derive_after_mi:
#   TRUE = 不直接插补；由已插补/已观测组成变量在每套 completed dataset 中重新派生。
#   本文件中 physical_activity 仅保留 observed baseline category 作为 QC/兼容字段；
#   06 必须根据 completed raw MET 覆盖重建。
# ------------------------------------------------------------------------------

variable_dictionary <- rbindlist(list(
  # ID
  new_variable_spec(
    "Participant_ID", "participant_id", "Participant identifier", NA,
    "identifier", "character",
    required_column = TRUE, required_observed = TRUE
  ),

  # Demographics
  new_variable_spec(
    "Age_at_recruitment", "age", "Age at recruitment", "years",
    "demographic", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_strategy = TRUE,
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),
  new_variable_spec(
    "Sex", "sex", "Sex", NA,
    "demographic", "character",
    required_observed = TRUE,
    used_for_strategy = TRUE, used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),
  new_variable_spec(
    "Region", "region", "Assessment-centre region", NA,
    "demographic", "character",
    required_observed = FALSE, impute = TRUE, imputation_method = "polyreg",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "Ethnic", "ethnicity", "High-risk South Asian ethnicity (0/1)", NA,
    "demographic", "character",
    required_observed = FALSE, impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),

  # SES components: final SES is NOT generated in 05.
  new_variable_spec(
    "TDI", "townsend_deprivation_index", "Townsend deprivation index", "index",
    "socioeconomic", "numeric",
    impute = TRUE, imputation_method = "pmm",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "Education", "education", "Educational attainment", NA,
    "socioeconomic_component", "character",
    impute = TRUE, imputation_method = "polr",
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),
  new_variable_spec(
    "Household_income", "household_income", "Annual household income before tax", "GBP_band",
    "socioeconomic_component", "character",
    impute = TRUE, imputation_method = "polr",
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),
  new_variable_spec(
    "employment_status", "employment_status", "Current employment status", NA,
    "socioeconomic_component", "character",
    impute = TRUE, imputation_method = "polyreg",
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),

  # Lifestyle
  new_variable_spec(
    "BMI0", "bmi", "Body mass index", "kg/m2",
    "lifestyle", "numeric",
    impute = TRUE, imputation_method = "pmm",
    used_for_outcome_model = FALSE
  ),
  new_variable_spec(
    "Smoke0", "smoking_status", "Smoking status", NA,
    "lifestyle", "character",
    required_observed = TRUE,
    used_for_strategy = TRUE, used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "MET_minutes_week0", "physical_activity_met_minutes_week",
    "Summed MET minutes per week for all activity", "MET-min/week",
    "lifestyle_component", "numeric",
    impute = FALSE,
    used_for_outcome_model = FALSE
  ),
  new_variable_spec(
    "PhysicalActivity0", "physical_activity",
    "Physical inactivity indicator (Field 22035, 1=inactive)", NA,
    "post_mi_derived", "integer",
    impute = FALSE, derive_after_mi = TRUE,
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "alcohol_score_0", "alcohol_score", "Alcohol score", "0/1",
    "lifestyle", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "diet_score_w0", "diet_score", "Weighted diet score", "0/1",
    "lifestyle", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),

  # Anthropometry (retained for description/QC; not in frozen outcome model)
  new_variable_spec(
    "WaistCircumference0", "waist_circumference", "Waist circumference", "cm",
    "anthropometry", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "HipCircumference0", "hip_circumference", "Hip circumference", "cm",
    "anthropometry", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Height0", "height", "Standing height", "cm",
    "anthropometry", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Weight0", "weight", "Body weight", "kg",
    "anthropometry", "numeric", required_column = TRUE
  ),

  # Blood pressure: strategy required
  new_variable_spec(
    "SystolicBP0", "systolic_blood_pressure", "Systolic blood pressure", "mmHg",
    "clinical", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_strategy = TRUE,
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),
  new_variable_spec(
    "DiastolicBP0", "diastolic_blood_pressure", "Diastolic blood pressure", "mmHg",
    "clinical", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_strategy = TRUE,
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),

  # Biomarkers
  new_variable_spec(
    "Cholesterol_0", "total_cholesterol", "Total cholesterol", "mmol/L",
    "biomarker", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_strategy = TRUE
  ),
  new_variable_spec(
    "LDL_direct_0", "ldl_cholesterol", "Direct LDL cholesterol", "mmol/L",
    "biomarker", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_strategy = TRUE, used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "HDL_cholesterol_0", "hdl_cholesterol", "HDL cholesterol", "mmol/L",
    "biomarker", "numeric",
    required_observed = TRUE,
    used_for_strategy = TRUE, used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "Triglycerides_0", "triglycerides", "Triglycerides", "mmol/L",
    "biomarker", "numeric",
    impute = TRUE, imputation_method = "pmm",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "Lipoprotein_A_0", "lipoprotein_a", "Lipoprotein(a)", "nmol/L",
    "biomarker", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Apolipoprotein_A_0", "apolipoprotein_a", "Apolipoprotein A", "g/L",
    "biomarker", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Apolipoprotein_B_0", "apolipoprotein_b", "Apolipoprotein B", "g/L",
    "biomarker", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Creatinine_0", "creatinine", "Serum creatinine", "umol/L",
    "biomarker", "numeric",
    required_observed = TRUE, used_for_eligibility = TRUE
  ),
  new_variable_spec(
    "Glucose_0", "glucose", "Serum glucose", "mmol/L",
    "biomarker", "numeric", required_column = TRUE
  ),
  new_variable_spec(
    "Glycated_haemoglobin_(HbA1c)_0", "hba1c", "Glycated haemoglobin", "mmol/mol",
    "biomarker", "numeric",
    required_observed = TRUE, used_for_eligibility = TRUE
  ),
  new_variable_spec(
    "egfr_2021_baseline", "egfr_2021", "eGFR using 2021 CKD-EPI equation", "mL/min/1.73m2",
    "clinical", "numeric",
    required_observed = TRUE,
    used_for_eligibility = TRUE, used_for_outcome_model = TRUE
  ),

  # Comorbidity / history / medication: retained and imputed where missing.
  new_variable_spec(
    "chronic_respiratory_disease", "chronic_respiratory_disease",
    "Chronic respiratory disease at baseline", "0/1",
    "comorbidity", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "chronic_liver_disease", "chronic_liver_disease",
    "Chronic liver disease at baseline", "0/1",
    "comorbidity", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "Family_history_CVD", "family_history_cvd",
    "Family history of cardiovascular disease", "0/1",
    "history", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE
  ),
  new_variable_spec(
    "hyper_Medicine_baseline", "antihypertensive_medication_baseline",
    "Antihypertensive medication at baseline", "0/1",
    "treatment", "integer",
    impute = TRUE, imputation_method = "logreg",
    used_for_outcome_model = TRUE, used_for_subgroup = TRUE
  ),

  # Follow-up and frozen outcome interface from 02.
  new_variable_spec(
    "InCohort_time", "baseline_date", "Baseline assessment date", "date",
    "followup", "date", required_observed = TRUE
  ),
  new_variable_spec(
    "Loss_to_followup_time", "loss_to_followup_date", "Loss-to-follow-up date", "date",
    "followup", "date", required_column = TRUE
  ),
  new_variable_spec(
    "Admin_censor_time", "administrative_censor_date", "Administrative censoring date", "date",
    "followup", "date", required_observed = TRUE
  ),
  new_variable_spec(
    "Censored", "loss_to_followup_event", "Loss-to-follow-up indicator", "0/1",
    "followup", "integer", required_observed = TRUE
  ),
  new_variable_spec(
    "Censored_time", "censor_date", "Individual censoring date", "date",
    "followup", "date", required_observed = TRUE
  ),
  new_variable_spec(
    "Censored_reason", "censor_reason", "Censoring reason", NA,
    "followup", "character", required_column = TRUE
  ),
  new_variable_spec(
    "primary_ascvd_analysis_date", "followup_end_date", "Primary ASCVD analysis end date", "date",
    "outcome", "date", required_observed = TRUE
  ),
  new_variable_spec(
    "primary_ascvd_status", "event_status", "Primary ASCVD competing-risk status", "0/1/2",
    "outcome", "integer", required_observed = TRUE
  ),
  new_variable_spec(
    "primary_ascvd_event", "ascvd_event", "Incident primary ASCVD event", "0/1",
    "outcome", "integer", required_observed = TRUE
  ),
  new_variable_spec(
    "primary_ascvd_time", "ascvd_source_date", "Source date of primary ASCVD", "date",
    "outcome_source", "date", required_column = TRUE
  ),
  new_variable_spec(
    "competing_death_time", "competing_death_source_date", "Source date of competing death", "date",
    "outcome_source", "date", required_column = TRUE
  ),
  new_variable_spec(
    "primary_ascvd_followup_years", "source_followup_years", "Primary ASCVD follow-up from 02", "years",
    "outcome_source", "numeric", required_observed = TRUE
  )
), use.names = TRUE, fill = TRUE)

validate_variable_dictionary(variable_dictionary)


# ============================================================================== #
# 2. 读取 04 最终人群并审计纳排常数 ####
# 模块任务：验证04纳排不变量
# 模块输入：data_eligible_M01.csv
# 模块输出：一致性审计
# ============================================================================== #

res_read <- run_with_timer(
  read_eligible_analysis_base(input_path, id_col = "Participant_ID"),
  "读取 04 最终 analytic cohort"
)

data_eligible <- res_read$value
input_n <- nrow(data_eligible)
input_ids <- as.character(data_eligible$Participant_ID)

step_message("04 输入样本量：", format(input_n, big.mark = ","))

res_eligibility_audit <- run_with_timer(
  audit_expected_values(
    data_eligible,
    expected_values = expected_eligibility_values,
    stop_on_failure = TRUE
  ),
  "核验 04 clinical eligibility + required-information gate"
)

eligibility_audit <- res_eligibility_audit$value
print(eligibility_audit)


# ============================================================================== #
# 3. 按变量字典选择变量并标准化取值 ####
# 模块任务：标准化数据并冻结观测RM
# 模块输入：04合格队列及RM来源
# 模块输出：0/1/NA观测RM、C目录RM定义表
# ============================================================================== #

# RM source columns pass through the same role dictionary; eligibility is unchanged.
rm_source_names <- grep("^rm_source_", names(data_eligible), value=TRUE)
if (!length(rm_source_names)) stop("04缺少RM源变量，请先运行01、02和04")
extra <- lapply(rm_source_names, function(v) new_variable_spec(
  v,v,v,NA_character_,"risk_modifier_source",
  if(grepl("flag$",v)) "character" else if(v %in% rm06_raw_binary_sources()) "integer" else "numeric",
  impute=v %in% rm06_imputation_targets(),
  imputation_method=if(v %in% rm06_binary_sources()) "logreg" else "pmm"
))
variable_dictionary <- rbindlist(c(list(variable_dictionary),extra),fill=TRUE)
# 2025 ethnicity is already imputed once as the high-risk ethnicity variable.
# rm_source_south_asian is a deterministic mirror, never a second MICE target.
variable_dictionary[standard_name == "rm_source_south_asian",
                    `:=`(impute=FALSE,imputation_method=NA_character_)]
for(v in intersect(c("waist_circumference"),variable_dictionary$standard_name))
  variable_dictionary[standard_name==v, `:=`(impute=TRUE,imputation_method="pmm")]
res_standardize <- run_with_timer(
  {
    dt <- apply_variable_dictionary(data_eligible, variable_dictionary)
    dt <- standardize_analysis_values(dt)
    assert_unique_id(dt, "participant_id", "标准化后的 analysis base")
    dt
  },
  "按变量角色字典建立 incomplete analysis base"
)

data_analysis_base <- res_standardize$value
# Freeze the definition of the single ethnicity predictor used by the expanded model.
# The old Ethnic field was broader (white/non-white) and is not the new covariate.
if (!"rm_source_south_asian" %in% names(data_analysis_base))
  stop("Missing baseline South-Asian RM input from 01")
data_analysis_base[, ethnicity := as.character(rm_source_south_asian)]
# Observed physical inactivity is Field 22035, NOT the legacy MET 4-level category.
data_analysis_base[, physical_activity := as.integer(rm_source_inactive)]
if (any(!is.na(data_analysis_base$ethnicity) &
        !data_analysis_base$ethnicity %in% c("0","1")))
  stop("High-risk ethnicity must be binary 0/1/NA")
rm_q80 <- as.numeric(quantile(data_analysis_base$townsend_deprivation_index,.8,na.rm=TRUE,type=7))
if (!is.finite(rm_q80)) stop("No observed Townsend to freeze Q80")
data_analysis_base[, rm_townsend_q80 := rm_q80]
observed_rm <- rm01_derive(data_analysis_base,rm_q80)
data_analysis_base[, rm_lpa_status := observed_rm$rm_lpa]
data_analysis_base[, rm_menopause_uncertain := observed_rm$rm_menopause_uncertain]
data_analysis_base[, family_history_cvd := as.integer(rm_source_family_proxy)]
variable_dictionary[standard_name %in% c("lipoprotein_a","apolipoprotein_b","family_history_cvd"), `:=`(impute=FALSE,imputation_method=NA_character_)]
variable_dictionary <- rbindlist(list(variable_dictionary,new_variable_spec(
 "derived_05_rm_lpa_status","rm_lpa_status","Lp(a) >105 nmol/L status","0/1",
 "risk_modifier_binary","integer",required_column=FALSE,impute=TRUE,imputation_method="logreg")),fill=TRUE)
# Frozen study classification: confirmed early onset=1; otherwise=0.
# Unknown onset is a study-rule zero; reproductive sources are never imputed.
data_analysis_base[, rm_premature_menopause := as.integer(observed_rm$rm_premature_menopause)]
variable_dictionary <- rbindlist(list(variable_dictionary, new_variable_spec(
  "derived_05_rm_premature_menopause", "rm_premature_menopause",
  "Premature menopause (<40 years), fixed confirmed-onset study rule",
  "0/1", "risk_modifier_binary", "integer", required_column=FALSE,
  impute=FALSE, imputation_method=NA_character_)), fill=TRUE)
# Do not allow any reproductive component or urinary ACR source in the MI target list.
forbidden_rm_targets <- c(rm06_reproductive_components(),
                           "rm_source_urine_albumin", "rm_source_urine_creatinine")
if (any(variable_dictionary[standard_name %in% forbidden_rm_targets & impute==TRUE, .N] > 0L))
  stop("Forbidden component imputation target remains in 05 dictionary")
observed_columns <- m01qa_rm_vars()
for(v in observed_columns) data_analysis_base[, (paste0(v,"_observed")) := observed_rm[[v]]]
variable_dictionary <- rbindlist(c(list(variable_dictionary),list(new_variable_spec(
 "rm_townsend_q80","rm_townsend_q80","Frozen cohort Townsend Q80",NA_character_,"RM_parameter","numeric",required_observed=TRUE)),
 lapply(observed_columns,function(v)new_variable_spec(paste0("derived_05_",v,"_observed"),paste0(v,"_observed"),v,NA_character_,"RM_observed","integer",required_column=FALSE))),fill=TRUE)
fwrite(rm01_manifest(),file.path(rm_report_dir,"RM_definition_manifest_M01.csv"),bom=TRUE)
fwrite(data.table(variable=observed_columns,positive=vapply(observed_columns,function(v)sum(observed_rm[[v]]==1,na.rm=TRUE),integer(1)),unknown=vapply(observed_columns,function(v)sum(is.na(observed_rm[[v]])),integer(1))),file.path(rm_report_dir,"RM_observed_classification_M01.csv"),bom=TRUE)
# 模块 3A｜插补前RM来源和0/1/NA描述
# 输入：05已标准化的data_analysis_base，已观测RM结果observed_rm。
# 输出：C目录中的源变量缺失表、二分类RM分布；未知绝不自动归零。
rm_source_pre_vars <- intersect(unique(c(
  grep("^rm_source_",names(data_analysis_base),value=TRUE),
  rm06_continuous_sources(),
  "triglycerides","bmi","townsend_deprivation_index")),names(data_analysis_base))
rm_source_missingness <- m01qa_missingness(data_analysis_base,rm_source_pre_vars,"preMI_05")
rm_source_missingness[,imputation_target := variable %in% rm06_imputation_targets()]
fwrite(rm_source_missingness,
  file.path(rm_report_dir,"RM_preMI_source_missingness_M01.csv"),bom=TRUE)
fwrite(m01qa_distribution(data_analysis_base,rm_source_pre_vars,"preMI_05"),
  file.path(rm_report_dir,"RM_preMI_source_distribution_M01.csv"),bom=TRUE)
fwrite(m01qa_rm_binary(observed_rm,observed_columns,"preMI_observed_05"),
  file.path(rm_report_dir,"RM_preMI_binary_distribution_M01.csv"),bom=TRUE)


# 05 明确不拟合 SES：
#   education / household_income / employment_status 原样保留；
#   ses_class 将由 06 在每套 completed dataset 中重新生成。
if ("ses_class" %in% names(data_analysis_base)) {
  stop("05 不应包含 ses_class；SES 必须延迟到 06 post-MI derivation。", call. = FALSE)
}


# ============================================================================== #
# 4. 整理 competing-risk 结局和随访接口 ####
# 模块任务：整理主要ASCVD随访接口
# 模块输入：02冻结结局
# 模块输出：事件和时间
# ============================================================================== #

res_followup <- run_with_timer(
  standardize_primary_ascvd_followup(data_analysis_base),
  "整理 primary ASCVD / competing death 随访接口"
)

data_analysis_base <- res_followup$value


# ============================================================================== #
# 5. 构造 MI auxiliary Nelson-Aalen cumulative hazards ####
# 模块任务：计算MI辅助风险信息
# 模块输入：ASCVD/竞争死亡及随访
# 模块输出：Nelson–Aalen辅助预测变量
# ============================================================================== #

res_na <- run_with_timer(
  add_mi_nelson_aalen_predictors(data_analysis_base),
  "构造 cause-specific Nelson-Aalen MI auxiliary predictors"
)

data_analysis_base <- res_na$value


# ============================================================================== #
# 6. 最终变量顺序和派生变量字典 ####
# 模块任务：整理完整字典与顺序
# 模块输入：变量角色映射及派生字段
# 模块输出：最终字段顺序
# ============================================================================== #

final_variable_order <- c(
  "participant_id",
  "age", "sex", "region", "ethnicity",
  "townsend_deprivation_index",
  "education", "household_income", "employment_status",
  "bmi", "smoking_status",
  "physical_activity_met_minutes_week", "physical_activity",
  "alcohol_score", "diet_score",
  "waist_circumference", "hip_circumference", "height", "weight",
  "systolic_blood_pressure", "diastolic_blood_pressure",
  "total_cholesterol", "ldl_cholesterol", "hdl_cholesterol", "triglycerides",
  "lipoprotein_a", "apolipoprotein_a", "apolipoprotein_b",
  "creatinine", "glucose", "hba1c", "egfr_2021",
  "chronic_respiratory_disease", "chronic_liver_disease",
  "family_history_cvd", "antihypertensive_medication_baseline",
  "baseline_date", "loss_to_followup_date", "administrative_censor_date",
  "loss_to_followup_event", "censor_date", "censor_reason",
  "followup_end_date", "followup_days", "followup_years",
  "event_status", "ascvd_event", "ascvd_event_date",
  "competing_death_event", "competing_death_date",
  "na_hazard_ascvd", "na_hazard_death"
)

assert_required_columns(data_analysis_base, final_variable_order, "05 最终 incomplete analysis base")
setcolorder(data_analysis_base, final_variable_order)

# 05 内部新派生变量也纳入变量字典。
# source_name 必须保持唯一：这里使用 derived_05_<standard_name> 作为审计来源标签，
# 而不是让所有派生变量共用同一个 source_name = "derived_05"。
# 这些 source_name 只是 provenance/字典标识，不表示 04 中真实存在同名原始列。
derived_dictionary <- rbindlist(list(
  new_variable_spec(
    "derived_05_followup_days", "followup_days", "Observed analysis follow-up", "days",
    "followup", "numeric", required_observed = TRUE
  ),
  new_variable_spec(
    "derived_05_followup_years", "followup_years", "Observed analysis follow-up", "years",
    "followup", "numeric", required_observed = TRUE
  ),
  new_variable_spec(
    "derived_05_ascvd_event_date", "ascvd_event_date", "Observed primary ASCVD event date", "date",
    "outcome", "date"
  ),
  new_variable_spec(
    "derived_05_competing_death_event", "competing_death_event", "Observed competing death event", "0/1",
    "outcome", "integer", required_observed = TRUE
  ),
  new_variable_spec(
    "derived_05_competing_death_date", "competing_death_date", "Observed competing death date", "date",
    "outcome", "date"
  ),
  new_variable_spec(
    "derived_05_na_hazard_ascvd", "na_hazard_ascvd",
    "Marginal cause-specific Nelson-Aalen cumulative hazard for ASCVD", "cumulative_hazard",
    "mi_auxiliary", "numeric", required_observed = TRUE
  ),
  new_variable_spec(
    "derived_05_na_hazard_death", "na_hazard_death",
    "Marginal cause-specific Nelson-Aalen cumulative hazard for non-ASCVD competing death",
    "cumulative_hazard", "mi_auxiliary", "numeric", required_observed = TRUE
  )
), use.names = TRUE, fill = TRUE)

internal_vars <- c("ascvd_source_date", "competing_death_source_date", "source_followup_years")
variable_dictionary_output <- rbindlist(list(
  variable_dictionary[!standard_name %in% internal_vars],
  derived_dictionary
), use.names = TRUE, fill = TRUE)

variable_dictionary_output <- variable_dictionary_output[
  standard_name %in% names(data_analysis_base)
]
variable_dictionary_output[, output_order := match(standard_name, names(data_analysis_base))]
setorder(variable_dictionary_output, output_order)
variable_dictionary_output[, output_order := NULL]

# 合并原始变量与 05 派生变量后再次做唯一性 QA。
# 这里若失败，说明 source_name 或 standard_name 的 provenance 设计发生冲突，
# 不应通过放宽 validate_variable_dictionary() 来掩盖。
validate_variable_dictionary(variable_dictionary_output)


# ============================================================================== #
# 7. 缺失、随访和数据契约 QA ####
# 模块任务：插补前描述及数据质量审计
# 模块输入：analysis base、字典、观测RM
# 模块输出：B目录缺失/分布与C目录RM分布
# ============================================================================== #

res_description <- run_with_timer(
  list(
    followup_outcome_summary = build_followup_outcome_summary(data_analysis_base),
    missingness_summary = build_missingness_summary(
      data_analysis_base,
      dictionary = variable_dictionary_output
    )
  ),
  "生成随访、结局和变量缺失描述"
)

followup_outcome_summary <- res_description$value$followup_outcome_summary
missingness_summary <- res_description$value$missingness_summary

# 模块 7A｜05插补前全变量分布快照；仅聚合统计，无个人数据。
analysis_preMI_distribution <- m01qa_distribution(data_analysis_base,
  vars=setdiff(names(data_analysis_base),"participant_id"),
  stage="preMI_05",dataset_id=0L)


validate_final_analysis_base(
  data_analysis_base,
  input_n = input_n,
  input_ids = input_ids
)

row_identity_audit <- build_row_identity_audit(input_ids, data_analysis_base$participant_id)
type_contract_audit <- build_type_contract_audit(
  data_analysis_base,
  variable_dictionary_output
)
required_observed_audit <- build_required_observed_audit(
  data_analysis_base,
  variable_dictionary_output
)
value_domain_audit <- build_value_domain_audit(data_analysis_base)
followup_event_consistency_audit <- build_followup_event_consistency_audit(data_analysis_base)

# MI manifest 只列出实际需要 06 处理的原始变量，以及 post-MI 派生占位变量。
mi_variable_manifest <- variable_dictionary_output[
  impute == TRUE | derive_after_mi == TRUE,
  .(
    variable = standard_name,
    role,
    type,
    impute,
    imputation_method,
    derive_after_mi,
    used_for_outcome_model,
    used_for_subgroup
  )
]

# 06 必须在每套 completed dataset 中按相同规则重新派生这些变量。
# 这里单独保存 manifest，避免把“直接插补”和“插补后确定性/模型派生”混在一起。
post_mi_derivation_manifest <- data.table(
  derived_variable = c("physical_activity", "ses_class", "hypertension"),
  source_variables = c(
    "physical_activity_met_minutes_week",
    "education + household_income + employment_status",
    "systolic_blood_pressure + diastolic_blood_pressure + antihypertensive_medication_baseline"
  ),
  derivation_rule = c(
    "Re-apply classify_physical_activity() using completed raw MET-min/week",
    "Fit 2-class LCA within each completed MI dataset; align classes by socioeconomic profile",
    "1 if SBP>=140 or DBP>=90 or antihypertensive medication==1; otherwise 0"
  ),
  direct_imputation = FALSE,
  used_for_outcome_model = c(TRUE, TRUE, FALSE),
  used_for_subgroup = c(FALSE, TRUE, TRUE)
)

# 明确检查：普通 MI 变量存在缺失时不得停止；这正是 05 输出的设计。
n_any_imputation_missing <- missingness_summary[
  variable %in% mi_variable_manifest[impute == TRUE, variable],
  sum(n_missing > 0L)
]
step_message("存在缺失的待插补变量数：", n_any_imputation_missing)
step_message("随访与结局摘要：")
print(followup_outcome_summary)


# ============================================================================== #
# 8. 保存 incomplete analysis base 和辅助结果 ####
# 模块任务：保存数据与报告
# 模块输入：所有核验通过的05输出
# 模块输出：M01根目录CSV/RDS，B/C目录报告
# ============================================================================== #

res_save <- run_with_timer(
  {
    # CSV 保留既有路径，方便人工检查和后续兼容；RDS 是 06 推荐读取来源。
    fwrite(
      data_analysis_base,
      file = output_csv_path,
      row.names = FALSE,
      bom = TRUE
    )
    saveRDS(data_analysis_base, file = output_rds_path, compress = FALSE)

    fwrite(
      variable_dictionary_output,
      file = file.path(description_dir, "variable_dictionary_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      mi_variable_manifest,
      file = file.path(description_dir, "mi_variable_manifest_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      post_mi_derivation_manifest,
      file = file.path(description_dir, "post_mi_derivation_manifest_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      eligibility_audit,
      file = file.path(description_dir, "eligibility_invariant_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      missingness_summary,
      file = file.path(description_dir, "analysis_base_missingness_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(analysis_preMI_distribution,
      file.path(description_dir,"analysis_base_preMI_distribution_M01.csv"),bom=TRUE)

    fwrite(
      row_identity_audit,
      file = file.path(description_dir, "analysis_base_row_identity_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      type_contract_audit,
      file = file.path(description_dir, "analysis_base_type_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      required_observed_audit,
      file = file.path(description_dir, "analysis_base_required_observed_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      value_domain_audit,
      file = file.path(description_dir, "analysis_base_value_domain_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      followup_event_consistency_audit,
      file = file.path(description_dir, "followup_event_consistency_check_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )
    fwrite(
      followup_outcome_summary,
      file = file.path(description_dir, "followup_outcome_summary_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    TRUE
  },
  "保存 incomplete analysis base 与 B 类辅助结果"
)

step_message("05 运行完成：样本量未改变，普通分析协变量缺失已保留供 06 MI。")

followup_support <- build_followup_support_audit(data_analysis_base)
for(nm in names(followup_support)) fwrite(followup_support[[nm]],file.path(description_dir,paste0(nm,"_M01.csv")),bom=TRUE)
fwrite(data.table(parameter=c("RM_version","lipid_window_days","townsend_q80","q80_n_observed","q80_algorithm","risk_set_convention","annual_interval","year_length_days"),
 value=c("1.8","365",as.character(rm_q80),sum(!is.na(data_analysis_base$townsend_deprivation_index)),"R quantile type 7; eligible cohort observed values","Y(t-): followup_years >= t","[start,end)","365.25")),
 file.path(description_dir,"analysis_parameters_M01.csv"),bom=TRUE)

m01_runtime_finish()
