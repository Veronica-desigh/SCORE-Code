# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 02_Clean_outcome_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   清理 M01 ASCVD 分析所需的疾病史、死亡、fatal / nonfatal 结局和随访变量。
#
# v1.1 修改要点：
#   1. fatal CVD 保留完整 ICD-10，并落实家族与精确排除；
#   2. cardiovascular 组分整合 first occurrence 与详细住院 ICD-10；
#   3. documented ASCVD/CVD 不再包含 death-cause 组分；
#   4. 慢性呼吸病/肝病按 Field 20002 Instance 0 或基线前 41270/41280 定义；
#   5. 所有 nonfatal 组分统一要求 event_date < death_date；
#   6. 输出组件计数、来源贡献和日期一致性 QA 表。
#   7. 保留 competing-risk 分析所需的 primary ASCVD status、随访时间和竞争死亡
#      信息；05 将直接复用这些结局构造 MI 的 cause-specific Nelson-Aalen
#      cumulative-hazard auxiliary predictors，不在 05 重新定义结局。
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
# 模块任务：配置、加载ICD和RM共病解析工具
# 模块输入：UKB来源路径、RM来源路径字典
# 模块输出：运行配置
# ============================================================================== #

project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
data_dir <- Sys.getenv("M01_UKB_DIR", file.path(project_root, "DATA", "UKB"))
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
rm_report_dir <- file.path(data_clean_dir,"风险修饰因素定义与处理记录_M01")
quality_control_dir <- file.path(data_clean_dir, "quality_control")

dir.create(data_clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(quality_control_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rm_report_dir, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(stringr)

setwd(project_root)

source(file.path(code_dir, "Function_02_Clean_outcome_M01.R"))
source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("02", data_clean_dir)

rm_workbook <- Sys.getenv("M01_RM_PATH_WORKBOOK",file.path(project_root,
  "LJW_code","Risk modifiers修改方案","RM的指南定义和数据定义","RM_2019_2025_缺失核查_20261010.xlsx"))
rm_path_registry <- rm01_configure_paths(rm_workbook, root = data_dir)


# 单独运行默认12线程；并行入口默认每进程4线程，可由M01_DT_THREADS调整。
m01_configure_threads()
message("data.table 线程数：", data.table::getDTthreads())

# 行政删失日期：与研究方案中随访截止日期保持一致。
admin_censor_date <- as.Date("2022-11-30")

# ============================================================================== #
# 1. 基线评估日期、地区和失访信息 ####
# 模块任务：建立基线时点和失访信息
# 模块输入：assessment和loss-to-follow-up
# 模块输出：基线与删失信息
# ============================================================================== #

assessment_file <- file.path(
  data_dir,
  "3_Assessment center",
  "3_Assessment centre_2_Touchscreen_2_Lifestyle and environment1.csv"
)

assessment_columns <- list(
  assessment_date = "Date of attending assessment centre | Instance 0",
  assessment_centre = "UK Biobank assessment centre | Instance 0"
)

population_file <- file.path(data_dir, "2_Population characteristics.csv")

population_columns <- list(
  loss_date = "Date lost to follow-up",
  loss_reason = "Reason lost to follow-up"
)

data_date <- read_assessment_and_censoring(
  assessment_file = assessment_file,
  assessment_columns = assessment_columns,
  population_file = population_file,
  population_columns = population_columns,
  admin_censor_date = admin_censor_date
)

var_info_print(data_date, c(
  "Participant ID", 
  "InCohort_time", 
  "Region", 
  "Loss_to_followup_time",
  "Admin_censor_time", 
  "Censored_time", 
  "Censored", 
  "Censored_reason"
))

# ============================================================================== #
# 2. 死亡日期、主要死因和竞争死亡 ####
# 模块任务：确认死亡、死因及竞争死亡
# 模块输入：death registry
# 模块输出：死亡状态/日期/竞争事件
# ============================================================================== #

death_file <- file.path(
  data_dir,
  "8_Health-related outcomes",
  "8_4Death Record",
  "8_Health-related outcomes_4_ Death Record.csv"
)

death_columns <- list(
  date_pattern = "date of death"
)

death_cause_file <- file.path(
  data_dir,
  "8_Health-related outcomes",
  "8_4Death Record",
  "8_Health-related outcomes_4_ Death Record_data_death_cause.csv"
)

death_cause_columns <- list(
  classification = "Classification of cause of death",
  cause = "Cause of death - ICD-10",
  primary_label = "Primary cause of death"
)

# 方案定义：先按 ICD-10 家族范围纳入，再排除 I60/I62 全家族和四个精确代码。
fatal_cvd_include_families <- unique(c(
  icd_range("I", 10:16),
  icd_range("I", 20:25),
  icd_range("I", 46:52),
  icd_range("I", 60:69),
  icd_range("I", 70:73),
  "R96"
))
fatal_cvd_exclude_families <- c("I60", "I62")
fatal_cvd_exclude_exact <- c("I51.4", "I67.1", "I67.5", "I68.2")

data_death <- read_death_records(
  death_file = death_file,
  death_columns = death_columns
)

data_death_cause <- read_death_causes(
  death_cause_file = death_cause_file,
  death_cause_columns = death_cause_columns,
  fatal_cvd_include_families = fatal_cvd_include_families,
  fatal_cvd_exclude_families = fatal_cvd_exclude_families,
  fatal_cvd_exclude_exact = fatal_cvd_exclude_exact
)

data_base <- merge(data_date, data_death, by = "Participant ID", all.x = TRUE)
data_base <- merge(data_base, data_death_cause, by = "Participant ID", all.x = TRUE)
data_base <- derive_death_classification(data_base)

# Review-only endpoint flag: I48 (AF/flutter) deaths remain in the prespecified
# fatal cardiovascular outcome, but their contribution requires clinical review.
if ("fatal_cvd_death_code" %in% names(data_base)) {
  i48_review <- data_base[, .(
    fatal_cvd_code = as.character(fatal_cvd_death_code),
    status = as.integer(fatal_cvd_death)
  )][!is.na(fatal_cvd_code) & grepl("^I48",fatal_cvd_code)]
  fwrite(i48_review[, .N,by=.(fatal_cvd_code,status)],
    file.path(quality_control_dir,"fatal_cvd_I48_code_review_M01.csv"),bom=TRUE)
}


cat("\n死亡和死因 QA：\n")
cat("有死亡日期人数：", data_base[Death == 1L, .N], "\n", sep = "")
cat(
  "死亡但缺主要死因人数：",
  data_base[Death == 1L & is.na(primary_death_cause_available), .N],
  "\n",
  sep = ""
)
print(data_base[, .N, by = .(Death, fatal_cvd_death, competing_death)][order(Death, fatal_cvd_death, competing_death)])

var_info_print(data_base, c(
  "Participant ID", 
  "Death", 
  "Death_time", 
  "fatal_cvd_death",
  "fatal_cvd_death_time",
  "fatal_cvd_death_code",
  "fatal_cvd_death_source",
  "competing_death", 
  "competing_death_time"
))

# ============================================================================== #
# 3. first occurrence 文件读取 ####
# 模块任务：读取first-occurrence疾病记录
# 模块输入：UKB first occurrences
# 模块输出：疾病候选记录
# ============================================================================== #

first_occurrence_files <- c(
  endocrine = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_3_Endocrine, nutritional and metabolic diseases.csv"
  ),
  nervous = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_5_Nervous system disorders.csv"
  ),
  circulatory = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_8_Circulatory system disorders.csv"
  ),
  respiratory = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_9_Respiratory system disorders.csv"
  ),
  digestive = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_10_Digestive system disorders.csv"
  ),
  genitourinary = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_7First occurrences",
    "8_Health-related outcomes_7_First occurrences_13_Genitourinary system disorders.csv"
  )
)

diabetes_codes <- icd_range("E", 10:14)
ckd_codes <- "N18"
ihd_codes <- icd_range("I", 20:25)
mi_codes <- icd_range("I", 21:23)
stroke_include_families <- icd_range("I", 60:69)
stroke_exclude_families <- c("I60", "I62")
stroke_exclude_exact <- c("I67.1", "I67.5", "I68.2")
# First occurrence 只有三位码；I67/I68 无法执行精确排除，因此仅使用无歧义家族。
stroke_first_occurrence_codes <- c("I61", icd_range("I", 63:66), "I69")
tia_codes <- "G45"
vascular_codes <- icd_range("I", 70:73)
heart_failure_codes <- "I50"
cardiac_arrest_codes <- "I46"
atrial_fibrillation_codes <- "I48"

# 只读取随后使用的first-occurrence家族；RM疾病由20002和住院来源定义。
first_occurrence_codes <- unique(c(diabetes_codes,ckd_codes,ihd_codes,
  stroke_first_occurrence_codes,tia_codes,vascular_codes,heart_failure_codes,
  cardiac_arrest_codes,atrial_fibrillation_codes))
data_first_occurrence <- m01_timed("02_first_occurrence_relevant_codes",
  read_first_occurrence_files(first_occurrence_files,
    relevant_codes=first_occurrence_codes))

# 仅保留有基线评估信息的参与者。
data_outcome <- merge(data_base, data_first_occurrence, by = "Participant ID", all.x = TRUE)

cat("\nfirst occurrence 合并后人数：", nrow(data_outcome), "\n", sep = "")

# ============================================================================== #
# 4. 详细住院 ICD-10、疾病史和 nonfatal 结局组分 ####
# 模块任务：整合住院和自报疾病史及ASCVD
# 模块输入：41270/41280和20002
# 模块输出：基线疾病史和结局组分
# ============================================================================== #

hospital_diagnosis_file <- file.path(
  data_dir,
  "8_Health-related outcomes",
  "8_3Hospital inpatient",
  "8_Health-related outcomes_3_Hospital inpatient_2_Summary Diagnoses1.csv"
)
hospital_diagnosis_columns <- list(
  code = "Diagnoses - ICD10",
  date_prefix = "Date of first in-patient diagnosis - ICD10 | Array "
)
relevant_hospital_families <- unique(c(
  diabetes_codes, ckd_codes, ihd_codes, stroke_include_families, tia_codes,
  vascular_codes, heart_failure_codes, cardiac_arrest_codes, atrial_fibrillation_codes,
  icd_range("J", 41:47), "J84", icd_range("J", 60:67), "J96",
  "B18", "K70", "K71", "K72", "K73", "K74", "K75", "K76", "I85", "E83"
))

# 共享输入：20002基线数组、41270代码和两个41280日期文件；路径来自上传工作簿。
# 一次读取，供原结局/共病定义及RM疾病定义共同使用。
rm_records <- m01_timed("02_RM_hospital_self_report_sources", rm02_read_records(data_dir))
hospital_diagnosis_file <- rm01_resolve(data_dir,basename(hospital_diagnosis_file))
print(rm_path_registry[file_name %in% c(basename(hospital_diagnosis_file),
  "8_Health-related outcomes_3_Hospital inpatient_2_Summary Diagnoses2.csv",
  "3_Assessment centre_4_Verbal interyiew__3_Medical conditions5.csv")])
hospital_icd10_events <- m01_timed("02_hospital_ICD_code_date_alignment", read_hospital_icd10_events(
  hospital_file = hospital_diagnosis_file,
  code_column = hospital_diagnosis_columns$code,
  date_prefix = hospital_diagnosis_columns$date_prefix,
  relevant_families = relevant_hospital_families,
  raw_input = rm_records$hes, date_parts = rm_records$dates
))
cat("\n住院 ICD-10 候选记录数：", nrow(hospital_icd10_events), "\n", sep = "")

# 常规疾病史变量；First Occurrence 均有真实日期时才形成记录。
data_outcome <- derive_icd_component(data_outcome, "E10", "T1D")
data_outcome <- derive_icd_component(data_outcome, "E11", "T2D")
data_outcome <- derive_icd_component(data_outcome, diabetes_codes, "Diabetes_ICD")
data_outcome <- derive_ckd_icd10_n18_component(data_outcome, ckd_code = ckd_codes, baseline_time_var = "InCohort_time")
data_outcome <- derive_icd_component(data_outcome, ihd_codes, "Coronary_Heart_Disease")
data_outcome <- derive_icd_component(data_outcome, stroke_first_occurrence_codes, "Stroke")
data_outcome <- derive_icd_component(data_outcome, atrial_fibrillation_codes, "Atrial_Fibrillation")
data_outcome <- derive_icd_component(data_outcome, heart_failure_codes, "Heart_Failure_history")

# 随访组件整合 first occurrence 与详细住院来源；stroke 的精确排除仅由完整住院码判定。
component_definitions <- list(
  nonfatal_ihd = list(fo = ihd_codes, include = ihd_codes),
  nonfatal_mi = list(fo = mi_codes, include = mi_codes),
  nonfatal_stroke = list(
    fo = stroke_first_occurrence_codes, include = stroke_include_families,
    exclude_families = stroke_exclude_families, exclude_exact = stroke_exclude_exact
  ),
  tia = list(fo = tia_codes, include = tia_codes),
  atherosclerosis_aortic_aneurysm_pvd = list(fo = vascular_codes, include = vascular_codes),
  heart_failure = list(fo = heart_failure_codes, include = heart_failure_codes),
  cardiac_arrest = list(fo = cardiac_arrest_codes, include = cardiac_arrest_codes)
)

for (component_name in names(component_definitions)) {
  definition <- component_definitions[[component_name]]
  data_outcome <- derive_multisource_icd_component(
    data_outcome,
    hospital_events = hospital_icd10_events,
    first_occurrence_codes = definition$fo,
    component_name = component_name,
    include_families = definition$include,
    exclude_families = definition$exclude_families %||% character(),
    exclude_exact = definition$exclude_exact %||% character()
  )
  data_outcome <- apply_nonfatal_rule(data_outcome, component_name)
}

# ============================================================================== #
# 5. 基线慢性呼吸病和慢性肝病 ####
# 模块任务：确定慢性呼吸/肝病
# 模块输入：基线自报、住院记录
# 模块输出：基线慢性共病
# ============================================================================== #

self_report_file <- file.path(
  data_dir,
  "3_Assessment center",
  "3_Assessment centre_4_Verbal interyiew__3_Medical conditions5.csv"
)
self_report_definitions <- list(
  chronic_respiratory_disease = c(
    "1111", "asthma", "1112", "chronic obstructive airways disease/copd",
    "1113", "emphysema/chronic bronchitis", "1472", "emphysema",
    "1114", "bronchiectasis", "1115", "interstitial lung disease",
    "1121", "pulmonary fibrosis", "1122", "fibrosing alveolitis/unspecified alveolitis",
    "1120", "asbestosis"
  ),
  chronic_liver_disease = c(
    "1579", "hepatitis b", "1580", "hepatitis c", "1581", "hepatitis d",
    "1604", "alcoholic liver disease / alcoholic cirrhosis",
    "1157", "non-infective hepatitis", "1158", "liver failure/cirrhosis",
    "1475", "sclerosing cholangitis", "1506", "primary biliary cirrhosis",
    "1141", "oesophageal varicies", "1507", "haemochromatosis"
  )
)
self_report_file <- rm01_resolve(data_dir,basename(self_report_file))
self_report_flags <- read_baseline_self_report_flags(self_report_file, self_report_definitions,
  raw_input=rm_records$self)

# 基线RM疾病：炎症性疾病、精神疾病、严格AF、HDP、HIV、睡眠呼吸暂停、脂肪肝。
# 严格AF不等同于原I48房颤/房扑病史；K76.0代理不等同于原慢性肝病。
# ICD10+日期按数组对齐；日期未知保留NA，自报阳性可确定为1。
rm_disease_sources <- m01_timed("02_RM_diagnostic_source_derivation",
  rm01_diagnoses(data_outcome[["Participant ID"]],
    data_outcome$InCohort_time,rm_records$self,rm_records$hes,rm_records$dates))
data_outcome <- rm01_attach(data_outcome,rm_disease_sources,"Participant ID")
fwrite(data.table(variable=names(rm_disease_sources),n_total=nrow(rm_disease_sources),
  missing=vapply(rm_disease_sources,function(x)sum(is.na(x)),integer(1)),
  pct_missing=vapply(rm_disease_sources,function(x)100*mean(is.na(x)),numeric(1))),
  file.path(rm_report_dir,"RM_disease_source_missingness_M01.csv"),bom=TRUE)
record_cutoff_audit <- audit_record_cutoffs(data_date,data_outcome,rm_records,admin_censor_date)
fwrite(record_cutoff_audit,file.path(quality_control_dir,"record_coverage_by_country_source_M01.csv"),bom=TRUE)
rm(rm_records)
gc(FALSE)
data_outcome <- merge(data_outcome, self_report_flags, by = "Participant ID", all.x = TRUE)

chronic_definitions <- list(
  chronic_respiratory_disease = list(
    include_families = c(icd_range("J", 41:46), "J47", "J84", icd_range("J", 60:67)),
    include_exact = "J96.1"
  ),
  chronic_liver_disease = list(
    include_families = c("B18", "K70", "K73", "K74", "I85"),
    include_exact = c("K75.4", "K71.5", "K71.6", "K71.7", "K72.1", "K76.0", "K76.6", "K76.7", "E83.0", "E83.1"),
    exclude_families = c(icd_range("B", 15:17), "B19", icd_range("K", 80:82)),
    exclude_exact = c("K72.0", "K75.0", "K83.0")
  )
)

for (component_name in names(chronic_definitions)) {
  definition <- chronic_definitions[[component_name]]
  hospital_component <- derive_baseline_hospital_comorbidity(
    data_outcome[, .(`Participant ID`, InCohort_time)],
    hospital_icd10_events,
    component_name = component_name,
    include_families = definition$include_families,
    include_exact = definition$include_exact %||% character(),
    exclude_families = definition$exclude_families %||% character(),
    exclude_exact = definition$exclude_exact %||% character()
  )
  data_outcome <- merge(data_outcome, hospital_component, by = "Participant ID", all.x = TRUE)
  data_outcome <- combine_baseline_comorbidity_sources(data_outcome, component_name)
}

# 兼容旧脚本使用的简称，但含义已更新为文档冻结定义。
data_outcome[, `:=`(
  CPD = chronic_respiratory_disease,
  CPD_time = chronic_respiratory_disease_time,
  CLD = chronic_liver_disease,
  CLD_time = chronic_liver_disease_time
)]

# ============================================================================== #
# 6. OPCS-4 冠脉血运重建组分 ####
# 模块任务：读取冠脉血运重建操作
# 模块输入：OPCS4数据
# 模块输出：冠脉操作结局组分
# ============================================================================== #

operation_files <- c(
  operations1 = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_3Hospital inpatient",
    "8_Health-related outcomes_3_Hospital inpatient_3_Summary Operations1.csv"
  ),
  operations2 = file.path(
    data_dir,
    "8_Health-related outcomes",
    "8_3Hospital inpatient",
    "8_Health-related outcomes_3_Hospital inpatient_3_Summary Operations2.csv"
  )
)

operation_columns <- list(
  code_columns = c(
    "Operative procedures - OPCS4",
    "Operative procedures - main OPCS4",
    "Operative procedures - secondary OPCS4",
    "Operative procedures - OPCS3",
    "Operative procedures - main OPCS3",
    "Operative procedures - secondary OPCS3"
  ),
  date_prefix_by_code = c(
    "Operative procedures - OPCS4" = "Date of first operative procedure - OPCS4 | Array ",
    "Operative procedures - main OPCS4" = "Date of first operative procedure - main OPCS4 | Array ",
    "Operative procedures - secondary OPCS4" = "Date of first operative procedure - secondary OPCS4 | Array ",
    "Operative procedures - OPCS3" = "Date of first operative procedure - OPCS3 | Array ",
    "Operative procedures - main OPCS3" = "Date of first operative procedure - main OPCS3 | Array ",
    "Operative procedures - secondary OPCS3" = "Date of first operative procedure - secondary OPCS3 | Array "
  )
)

coronary_revascularization_opcs_codes <- unique(c(
  opcs_range("K", 40:46),
  "K49", "K50", "K75"
))

procedure_ids <- data_outcome[, .(`Participant ID`)]

data_revascularization <- derive_coronary_revascularization(
  operation_files = operation_files,
  operation_columns = operation_columns,
  opcs_codes = coronary_revascularization_opcs_codes,
  ids = procedure_ids,
  component_name = "coronary_revascularization"
)

data_outcome <- merge(data_outcome, data_revascularization, by = "Participant ID", all.x = TRUE)
data_outcome <- apply_nonfatal_rule(data_outcome, "coronary_revascularization")

# ============================================================================== #
# 7. 纳排变量和复合结局变量 ####
# 模块任务：组合结局和纳排状态
# 模块输入：前述临床事件
# 模块输出：ASCVD/CVD派生状态
# ============================================================================== #

# inclusion_ascvd：用于主分析一级预防人群 documented ASCVD 排除。
inclusion_ascvd_components <- c(
  "nonfatal_ihd_time",
  "nonfatal_stroke_time",
  "tia_time",
  "atherosclerosis_aortic_aneurysm_pvd_time",
  "coronary_revascularization_time"
)

data_outcome <- combine_components(data_outcome, inclusion_ascvd_components, "inclusion_ascvd")
data_outcome <- derive_prevalent_component(data_outcome, "inclusion_ascvd", baseline_time_var = "InCohort_time")

# inclusion_cvd：用于 secondary CVD 分析的 documented CVD 排除/参考。
inclusion_cvd_components <- c(
  "inclusion_ascvd_time",
  "heart_failure_time",
  "cardiac_arrest_time"
)

data_outcome <- combine_components(data_outcome, inclusion_cvd_components, "inclusion_cvd")
data_outcome <- derive_prevalent_component(data_outcome, "inclusion_cvd", baseline_time_var = "InCohort_time")

# primary_ascvd：fatal CVD + nonfatal MI + nonfatal stroke。
primary_ascvd_components <- c(
  "fatal_cvd_death_time",
  "nonfatal_mi_time",
  "nonfatal_stroke_time"
)

data_outcome <- combine_components(data_outcome, primary_ascvd_components, "primary_ascvd")
data_outcome <- derive_prevalent_component(data_outcome, "primary_ascvd", baseline_time_var = "InCohort_time")

# secondary_cvd：fatal CVD + 更宽泛的 nonfatal CVD 组分。
secondary_cvd_components <- c(
  "fatal_cvd_death_time",
  "nonfatal_ihd_time",
  "nonfatal_stroke_time",
  "heart_failure_time",
  "cardiac_arrest_time",
  "atherosclerosis_aortic_aneurysm_pvd_time",
  "tia_time",
  "coronary_revascularization_time"
)

data_outcome <- combine_components(data_outcome, secondary_cvd_components, "secondary_cvd")
data_outcome <- derive_prevalent_component(data_outcome, "secondary_cvd", baseline_time_var = "InCohort_time")

# ============================================================================== #
# 8. 随访分析状态和随访时间 ####
# 模块任务：定义随访和竞争风险接口
# 模块输入：事件日期/删失日期
# 模块输出：冻结随访与事件状态
# ============================================================================== #

# 状态编码：0 = 删失；1 = 目标结局；2 = 竞争死亡。
data_outcome <- make_analysis_endpoint(
  dt = data_outcome,
  endpoint_name = "primary_ascvd",
  event_time_var = "primary_ascvd_time",
  competing_time_var = "competing_death_time",
  baseline_time_var = "InCohort_time",
  censor_time_var = "Censored_time"
)

data_outcome <- make_analysis_endpoint(
  dt = data_outcome,
  endpoint_name = "secondary_cvd",
  event_time_var = "secondary_cvd_time",
  competing_time_var = "competing_death_time",
  baseline_time_var = "InCohort_time",
  censor_time_var = "Censored_time"
)

# ============================================================================== #
# 9. QA 检查与可审计输出 ####
# 模块任务：核查事件构成和日期
# 模块输入：上游组分
# 模块输出：结局QA表
# ============================================================================== #

endpoint_check_vars <- c(
  "fatal_cvd_death", "competing_death",
  "inclusion_ascvd", "prevalent_inclusion_ascvd",
  "inclusion_cvd", "prevalent_inclusion_cvd",
  "primary_ascvd", "prevalent_primary_ascvd", "primary_ascvd_event", "primary_ascvd_status",
  "secondary_cvd", "prevalent_secondary_cvd", "secondary_cvd_event", "secondary_cvd_status"
)

print_endpoint_qa(data_outcome, endpoint_check_vars)

cat("\n随访时间 summary：primary_ascvd_followup_years\n")
print(summary(data_outcome$primary_ascvd_followup_years))
cat("\n随访时间 summary：secondary_cvd_followup_years\n")
print(summary(data_outcome$secondary_cvd_followup_years))

endpoint_component_map <- list(
  inclusion_ascvd = sub("_time$", "", inclusion_ascvd_components),
  inclusion_cvd = sub("_time$", "", inclusion_cvd_components),
  primary_ascvd = sub("_time$", "", primary_ascvd_components),
  secondary_cvd = sub("_time$", "", secondary_cvd_components)
)
endpoint_component_qa <- rbindlist(lapply(names(endpoint_component_map), function(endpoint) {
  rbindlist(lapply(endpoint_component_map[[endpoint]], function(component) {
    component_time <- paste0(component, "_time")
    endpoint_time <- paste0(endpoint, "_time")
    data.table(
      endpoint = endpoint,
      component = component,
      n_component = sum(data_outcome[[component]] == 1L, na.rm = TRUE),
      n_at_composite_earliest_date = sum(
        !is.na(data_outcome[[component_time]]) &
          data_outcome[[component_time]] == data_outcome[[endpoint_time]],
        na.rm = TRUE
      )
    )
  }))
}))

source_components <- c(
  "fatal_cvd_death", names(component_definitions), "coronary_revascularization",
  "chronic_respiratory_disease", "chronic_liver_disease"
)
source_labels <- c("death_registry", "first_occurrence", "hospital", "opcs", "self_report")
source_contribution_qa <- rbindlist(lapply(source_components, function(component) {
  source_var <- paste0(component, "_source")
  if (!source_var %in% names(data_outcome)) return(NULL)
  rbindlist(lapply(source_labels, function(source) {
    data.table(
      component = component,
      source = source,
      n_participants = sum(
        grepl(paste0("(^|;)", source, "(;|$)"), data_outcome[[source_var]]),
        na.rm = TRUE
      )
    )
  }))
}))

nonfatal_components <- c(names(component_definitions), "coronary_revascularization")
nonfatal_on_or_after_death <- sum(vapply(nonfatal_components, function(component) {
  event_time <- data_outcome[[paste0(component, "_time")]]
  sum(!is.na(event_time) & !is.na(data_outcome$Death_time) & event_time >= data_outcome$Death_time)
}, integer(1)))
date_consistency_qa <- data.table(
  check = c(
    "nonfatal_event_on_or_after_death",
    "chronic_respiratory_hospital_after_baseline",
    "chronic_liver_hospital_after_baseline",
    "primary_negative_followup",
    "secondary_negative_followup",
    "inclusion_ascvd_contains_death_component_definition"
  ),
  n_violations = c(
    nonfatal_on_or_after_death,
    data_outcome[!is.na(chronic_respiratory_disease_hospital_time) & chronic_respiratory_disease_hospital_time > InCohort_time, .N],
    data_outcome[!is.na(chronic_liver_disease_hospital_time) & chronic_liver_disease_hospital_time > InCohort_time, .N],
    data_outcome[primary_ascvd_followup_years < 0, .N],
    data_outcome[secondary_cvd_followup_years < 0, .N],
    as.integer("fatal_cvd_death_time" %in% inclusion_ascvd_components)
  )
)

fwrite(endpoint_component_qa, file.path(quality_control_dir, "endpoint_component_count_M01.csv"))
fwrite(source_contribution_qa, file.path(quality_control_dir, "source_contribution_M01.csv"))
fwrite(date_consistency_qa, file.path(quality_control_dir, "date_consistency_M01.csv"))
if (any(date_consistency_qa$n_violations > 0L)) {
  stop(
    "02 日期一致性 QA 存在违规，请检查：",
    file.path(quality_control_dir, "date_consistency_M01.csv"),
    call. = FALSE
  )
}

# ============================================================================== #
# 10. 选择最终变量并输出 ####
# 模块任务：保存最终结局数据
# 模块输入：全部结局与RM疾病来源
# 模块输出：M01_ASCVD/data_outcome_M01.csv
# ============================================================================== #

base_vars <- c(
  "Participant ID",
  "InCohort_time",
  "Region",
  "Loss_to_followup_time",
  "Admin_censor_time",
  "Censored",
  "Censored_time",
  "Censored_reason",
  "Death",
  "Death_time",
  "fatal_cvd_death",
  "fatal_cvd_death_time",
  "fatal_cvd_death_code",
  "fatal_cvd_death_source",
  "competing_death",
  "competing_death_time"
)

history_vars <- c(
  "T1D", "T1D_time",
  "T2D", "T2D_time",
  "Diabetes_ICD", "Diabetes_ICD_time",
  "CKD_ICD10_N18", "CKD_ICD10_N18_time", "prevalent_CKD_ICD10_N18",
  "CKD", "CKD_time",
  "chronic_liver_disease", "chronic_liver_disease_time", "chronic_liver_disease_source",
  "chronic_liver_disease_self_report", "chronic_liver_disease_hospital",
  "chronic_liver_disease_hospital_time", "chronic_liver_disease_hospital_code",
  "chronic_respiratory_disease", "chronic_respiratory_disease_time", "chronic_respiratory_disease_source",
  "chronic_respiratory_disease_self_report", "chronic_respiratory_disease_hospital",
  "chronic_respiratory_disease_hospital_time", "chronic_respiratory_disease_hospital_code",
  "CLD", "CLD_time",
  "CPD", "CPD_time",
  "Coronary_Heart_Disease", "Coronary_Heart_Disease_time",
  "Stroke", "Stroke_time",
  "Atrial_Fibrillation", "Atrial_Fibrillation_time",
  "Heart_Failure_history", "Heart_Failure_history_time"
)

component_vars <- c(
  "nonfatal_ihd", "nonfatal_ihd_time", "nonfatal_ihd_code", "nonfatal_ihd_source",
  "nonfatal_mi", "nonfatal_mi_time", "nonfatal_mi_code", "nonfatal_mi_source",
  "nonfatal_stroke", "nonfatal_stroke_time", "nonfatal_stroke_code", "nonfatal_stroke_source",
  "tia", "tia_time", "tia_code", "tia_source",
  "atherosclerosis_aortic_aneurysm_pvd", "atherosclerosis_aortic_aneurysm_pvd_time",
  "atherosclerosis_aortic_aneurysm_pvd_code", "atherosclerosis_aortic_aneurysm_pvd_source",
  "heart_failure", "heart_failure_time", "heart_failure_code", "heart_failure_source",
  "cardiac_arrest", "cardiac_arrest_time", "cardiac_arrest_code", "cardiac_arrest_source",
  "coronary_revascularization", "coronary_revascularization_time",
  "coronary_revascularization_code", "coronary_revascularization_source"
)

endpoint_vars <- c(
  "inclusion_ascvd", "inclusion_ascvd_time", "prevalent_inclusion_ascvd",
  "inclusion_cvd", "inclusion_cvd_time", "prevalent_inclusion_cvd",
  "primary_ascvd", "primary_ascvd_time", "prevalent_primary_ascvd",
  "primary_ascvd_analysis_date", "primary_ascvd_status", "primary_ascvd_event",
  "primary_ascvd_followup_years",
  "secondary_cvd", "secondary_cvd_time", "prevalent_secondary_cvd",
  "secondary_cvd_analysis_date", "secondary_cvd_status", "secondary_cvd_event",
  "secondary_cvd_followup_years"
)

all_needed_vars <- unique(c(
  base_vars,
  history_vars,
  component_vars,
  endpoint_vars,
  grep("^rm_source_",names(data_outcome),value=TRUE)
))
all_needed_vars <- intersect(all_needed_vars, names(data_outcome))

data_outcome_final <- data_outcome[, ..all_needed_vars]
data_outcome_final <- standardize_missing_values(data_outcome_final)
assert_no_name_conflicts(data_outcome_final)

fwrite(
  data_outcome_final,
  file = file.path(data_clean_dir, "data_outcome_M01.csv")
)

cat("\n02_Clean_outcome_M01_v1.1.R 已完成\n")
cat("主输出文件：", file.path(data_clean_dir, "data_outcome_M01.csv"), "\n", sep = "")
cat("质量控制目录：", quality_control_dir, "\n", sep = "")
cat("参与者数量：", nrow(data_outcome_final), "\n", sep = "")
cat("变量数量：", ncol(data_outcome_final), "\n", sep = "")

m01_runtime_finish()
