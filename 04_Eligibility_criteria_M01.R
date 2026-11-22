# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 04_Eligibility_criteria_M01_v1.2.R
# ------------------------------------------------------------------------------ #
# 用途：
#   基于 01 暴露/协变量、02 结局/疾病史和 03 用药清理结果，生成 M01 ASCVD
#   主分析人群，并输出主文 Figure 1 和补充材料所需的纳入排除流程记录。
#
# 本次重构要点：
#   1. 保留既有 clinical eligibility 定义：年龄、既往 ASCVD、糖尿病、CKD、
#      近期降脂药和 markedly elevated single risk factors。
#   2. 删除原“Complete baseline covariates”完整病例筛选。
#   3. 新增唯一的 required-information gate：只要求 eligibility 和
#      guideline-strategy assignment 所需变量真实观测。
#   4. Ethnicity、TDI、SES 组成变量、BMI、raw MET、TG、降压药、家族史和
#      慢性病等普通分析协变量允许缺失，并保留给 06 multiple imputation。
#   5. 纳排流程最终节点固定为：
#      Complete information required for eligibility and guideline-strategy assignment。
#   6. 随后插补变量的缺失人数/比例单独输出为 descriptive missingness table，
#      但不再构成 participant exclusion。
#   7. 所有一对一连接、流程人数和最终纳排不变量继续输出 QA 审计。
# 模块任务：输出纳排质量审计
# 模块输入：ID及规则检测
# 模块输出：QA文件
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
# 模块任务：配置纳排函数及所需路径
# 模块输入：01-03输出路径与04函数
# 模块输出：运行环境
# ============================================================================== #

project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
data_dir <- Sys.getenv("M01_UKB_DIR", file.path(project_root, "DATA", "UKB"))
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
quality_control_dir <- file.path(data_clean_dir, "quality_control")
data_intervention_dir <- file.path(project_root, "LJW_outcome", "data_intervention", "M01_ASCVD")

# 显式创建本模块输出根目录，后续敏感性分析无需预先手工建同类文件夹。
for (path in c(data_clean_dir, quality_control_dir, data_intervention_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

library(data.table)
library(parallel)

setwd(project_root)

source(file.path(code_dir, "Function_04_Eligibility_criteria_M01.R"))
source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("04", data_clean_dir)

run_eligibility_unit_tests()

# 高并行配置：默认 12 个 data.table 线程，可通过 M01_DT_THREADS 在 8–14 间调整。
m01_configure_threads()
step_message("data.table 线程数：", data.table::getDTthreads())


# ============================================================================== #
# 1. 读取并严格合并 exposure/outcome/medicine/糖尿病自报 ####
# 模块任务：按ID合并主要来源
# 模块输入：data_exposure/outcome/medicine
# 模块输出：一人一行初始数据
# ============================================================================== #
# 输入文件：
#   1. data_exposure_M01.csv：01 输出的基线暴露和协变量；
#   2. data_outcome_M01.csv：02 输出的疾病史、死亡和结局变量；
#   3. data_medicine_M01.csv：03 输出的一人一行基线用药证据；
#   4. Field 20002 Instance 0：基线自报糖尿病。
# 模块任务：排除显著升高单项危险因素
# 模块输入：血压、脂质等已观察指标
# 模块输出：符合主要临床标准人群
# ------------------------------------------------------------------------------ #

# 读取前核查全部上游输出，直接标明需要完成的清理步骤。
upstream_inputs <- data.table(
  producer = c("01_Clean_exposure_M01.R", "02_Clean_outcome_M01.R", "03_Clean_medicine_M01.R"),
  path = file.path(data_clean_dir,c("data_exposure_M01.csv", "data_outcome_M01.csv", "data_medicine_M01.csv"))
)
upstream_inputs[, readable := vapply(path,function(p) {
  if(!file.exists(p))return(FALSE)
  tryCatch({fread(p,nrows=0L);TRUE},error=function(e)FALSE)
},logical(1))]
if (any(!upstream_inputs$readable)) {
  missing_inputs <- upstream_inputs[readable==FALSE]
  stop("04缺少或无法读取上游输出：\n",
    paste(paste0(missing_inputs$producer," -> ",missing_inputs$path),collapse="\n"),
    "\n请先完成对应脚本；如使用M01_CLEAN_DIR，01–04必须设置相同输出目录。",call.=FALSE)
}

res_read <- run_with_timer(
  {
    data_exposure <- fread(file.path(data_clean_dir, "data_exposure_M01.csv"))
    data_outcome <- fread(file.path(data_clean_dir, "data_outcome_M01.csv"))
    data_medicine <- fread(file.path(data_clean_dir, "data_medicine_M01.csv"))
    diabetes_self_report <- read_baseline_diabetes_self_report(file.path(
      data_dir,
      "3_Assessment center",
      "3_Assessment centre_4_Verbal interyiew__3_Medical conditions5.csv"
    ))
    setDT(data_exposure)
    setDT(data_outcome)
    setDT(data_medicine)

    merge_audit <- list()
    joined <- merge_one_to_one_audited(data_exposure, data_outcome, "01 exposure + 02 outcome")
    merge_audit[[1L]] <- joined$audit
    joined <- merge_one_to_one_audited(joined$dt, data_medicine, "01/02 + 03 medicine")
    merge_audit[[2L]] <- joined$audit
    joined <- merge_one_to_one_audited(joined$dt, diabetes_self_report, "01/02/03 + baseline diabetes self-report")
    merge_audit[[3L]] <- joined$audit

    list(dt = joined$dt, merge_audit = rbindlist(merge_audit, fill = TRUE))
  },
  "读取并严格合并 01/02/03 与基线糖尿病自报"
)

data_01 <- res_read$value$dt
required_rm_sources <- c(paste0("rm_source_",c("inflammatory","psychiatric","af","hdp","hiv","sleep_apnoea","fatty_liver")),
  "rm_source_crp", "rm_source_south_asian")
if (!all(required_rm_sources %in% names(data_01))) stop("请重新运行01和02，生成完整RM源字段")
# Positive infection OR verified HIV-specific ART, preserving unresolved sources.
data_01[, rm_source_hiv := fifelse(rm_source_hiv==1L | hiv_specific_art==1L,1L,
  fifelse(!is.na(rm_source_hiv) & !is.na(hiv_specific_art),0L,NA_integer_))]
merge_audit_04 <- res_read$value$merge_audit
fwrite(
  merge_audit_04,
  file.path(quality_control_dir, "merge_audit_M01.csv"),
  bom = TRUE
)
if (any(merge_audit_04$n_unmatched_main > 0L)) {
  stop("04 输入表 ID 集不完整，请检查 04_merge_audit_M01.csv 的上游来源。", call. = FALSE)
}
n_start <- nrow(data_01)
flow_log_main <- init_flow_log()
flow_log_detail <- init_flow_log()

step_message("初始样本量：", format(n_start, big.mark = ","))


# ============================================================================== #
# 2. 主分析基础人群：基线、年龄和疾病史 ####
# 模块任务：应用年龄及基础疾病排除
# 模块输入：基线日期、ASCVD/糖尿病/CKD
# 模块输出：基础一级预防人群
# ============================================================================== #
# 该部分定义 SCORE/SCORE2 主分析的基础一级预防人群。
# 基线前/基线时 documented ASCVD、糖尿病或 CKD 者不进入主分析。
# ------------------------------------------------------------------------------ #

res_base_elig <- run_with_timer(
  {
    dt <- copy(data_01)
    setDT(dt)

    # ------------------------------------------------------------------------ #
    # Step 1. 可识别且具有明确基线时点 ####
    # ------------------------------------------------------------------------ #
    # 主文合并展示；补充材料拆分 ID 和基线日期两个技术性条件。
    n_before_step1 <- nrow(dt)

    n_before_s1a <- nrow(dt)
    dt <- dt[!is.na(`Participant ID`)]
    flow_log_detail <- add_flow_step(
      flow_log_detail,
      "S1a. Participant ID 未缺失",
      n_before_s1a,
      nrow(dt),
      n_start
    )

    n_before_s1b <- nrow(dt)
    dt <- dt[!is.na(InCohort_time)]
    flow_log_detail <- add_flow_step(
      flow_log_detail,
      "S1b. 基线调查日期 InCohort_time 未缺失",
      n_before_s1b,
      nrow(dt),
      n_start
    )

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 1. 可识别且具有明确基线时点（Identifiable participants with defined baseline）",
      n_before_step1,
      nrow(dt),
      n_start
    )
    step_message(
      "纳排步骤：Step 1. 可识别且具有明确基线时点",
      " | 排除 = ", format(n_before_step1 - nrow(dt), big.mark = ","),
      " | 剩余 = ", format(nrow(dt), big.mark = ",")
    )

    # ------------------------------------------------------------------------ #
    # Step 2. 年龄 40-69 岁 ####
    # ------------------------------------------------------------------------ #
    n_before <- nrow(dt)
    tmp <- apply_filter_with_log(
      dt = dt,
      i = Age_at_recruitment >= 40 & Age_at_recruitment <= 69,
      flow_log = flow_log_main,
      step_name = "Step 2. 入组年龄 40-69 岁（Age-eligible at baseline）",
      n_start = n_start
    )
    dt <- tmp$dt
    flow_log_main <- tmp$flow_log
    flow_log_detail <- add_flow_step(flow_log_detail, "S2. 入组年龄 40-69 岁", n_before, nrow(dt), n_start)

    # ------------------------------------------------------------------------ #
    # Step 3. 基线前/基线时无 documented ASCVD ####
    # ------------------------------------------------------------------------ #
    # 使用 02 脚本中为主分析生成的 inclusion_ascvd。
    n_before <- nrow(dt)
    dt <- derive_documented_ascvd_baseline_status(dt, baseline_col = "InCohort_time")
    dt <- dt[documented_ascvd_baseline == 0L]

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 3. 基线前/基线时无 documented ASCVD（Free of documented ASCVD at baseline）",
      n_before,
      nrow(dt),
      n_start
    )
    flow_log_detail <- add_flow_step(flow_log_detail, "S3. 基线前/基线时无 documented ASCVD", n_before, nrow(dt), n_start)
    step_message(
      "纳排步骤：Step 3. 基线前/基线时无 documented ASCVD",
      " | 排除 = ", format(n_before - nrow(dt), big.mark = ","),
      " | 剩余 = ", format(nrow(dt), big.mark = ",")
    )

    # ------------------------------------------------------------------------ #
    # Step 4. 基线前/基线时无糖尿病 ####
    # ------------------------------------------------------------------------ #
    # 糖尿病证据包括 ICD、Field 20002 Instance 0 自报、HbA1c 和 03 输出的降糖药/胰岛素。
    # HbA1c 缺失不在糖尿病证据函数内提前删除；最终由 Step 8 required-information gate 处理。
    n_before <- nrow(dt)
    dt <- derive_diabetes_baseline_status(dt, baseline_col = "InCohort_time")
    dt <- dt[diabetes_baseline == 0L]

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 4. 基线前/基线时无糖尿病（Free of diabetes at baseline）",
      n_before,
      nrow(dt),
      n_start
    )
    flow_log_detail <- add_flow_step(flow_log_detail, "S4. 基线前/基线时无糖尿病（ICD/自报/HbA1c/降糖药）", n_before, nrow(dt), n_start)
    step_message(
      "纳排步骤：Step 4. 基线前/基线时无糖尿病",
      " | 排除 = ", format(n_before - nrow(dt), big.mark = ","),
      " | 剩余 = ", format(nrow(dt), big.mark = ",")
    )

    # ------------------------------------------------------------------------ #
    # Step 5. 基线前/基线时无 CKD ####
    # ------------------------------------------------------------------------ #
    # CKD 证据包括 ICD N18 和基线 eGFR <60。
    n_before <- nrow(dt)
    dt <- derive_ckd_baseline_status(dt, baseline_col = "InCohort_time")
    dt <- dt[CKD_baseline == 0L]

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 5. 基线前/基线时无 CKD（Free of CKD at baseline）",
      n_before,
      nrow(dt),
      n_start
    )
    flow_log_detail <- add_flow_step(flow_log_detail, "S5. 基线前/基线时无 CKD（ICD N18 或 eGFR <60）", n_before, nrow(dt), n_start)
    step_message(
      "纳排步骤：Step 5. 基线前/基线时无 CKD",
      " | 排除 = ", format(n_before - nrow(dt), big.mark = ","),
      " | 剩余 = ", format(nrow(dt), big.mark = ",")
    )

    list(
      dt = dt,
      flow_log_main = flow_log_main,
      flow_log_detail = flow_log_detail
    )
  },
  "完成主分析基础人群筛选"
)

data_02 <- res_base_elig$value$dt
flow_log_main <- res_base_elig$value$flow_log_main
flow_log_detail <- res_base_elig$value$flow_log_detail
step_message("基础人群筛选后样本量：", format(nrow(data_02), big.mark = ","))


# ============================================================================== #
# 3. 匹配并排除近期降脂药使用者 ####
# 模块任务：排除基线前近期降脂药
# 模块输入：自报与GP处方
# 模块输出：无近期治疗人群
# ============================================================================== #
# 数据来源：data_medicine_M01.csv 中冻结的基线自报和基线前 365 天 GP 证据。
#
# 主分析规则保持原先约定：
#   1. 基线自报降脂药 = 1 -> treated，排除；
#   2. 基线自报 = 0，但基线前 365 天 GP 有降脂药处方 -> treated，排除；
#   3. 基线自报缺失，但基线前 365 天 GP 有降脂药处方 -> treated，排除；
#   4. 基线自报缺失且 365 天内无 GP 处方证据 -> 状态不可判定；
#      Step 6 暂留，最终在 Step 8 required-information gate 排除；
#   5. 最终 analytic cohort 仅保留 Lipid_Medicine_baseline == 0 的 untreated 人群。
# 模块任务：要求策略分配关键字段完整
# 模块输入：用于纳排及SCORE/SCORE2的已观察字段
# 模块输出：required-information gate + 缺失描述
# ------------------------------------------------------------------------------ #

res_lipid_screen <- run_with_timer(
  {
    dt <- copy(data_02)
    setDT(dt)

    require_columns(
      dt,
      c("lipid_lowering_med_baseline", "lipid_lowering_med_source", "lipid_gp_365d"),
      "03 统一降脂药证据"
    )
    dt[, Lipid_Medicine_baseline := as.integer(lipid_lowering_med_baseline)]
    dt[, Lipid_Medicine := Lipid_Medicine_baseline]
    dt[, Lipid_Medicine_time := as.IDate(fifelse(
      !is.na(Lipid_Medicine_baseline), as.IDate(InCohort_time), as.IDate(NA)
    ))]
    dt[, lipid_status_source := lipid_lowering_med_source]

    dt[, lipid_status_known := as.integer(!is.na(Lipid_Medicine_baseline))]

    # Step 6 只排除“已知正在使用近期降脂药”的参与者。
    # 对于 Lipid_Medicine_baseline == NA 的参与者，此处暂时保留；因为无法确认
    # guideline strategy，最终会在 Step 8 required-information gate 中排除。
    # 这样可把“临床排除”和“无法确定策略所需信息”在流程图中清楚分开。
    n_before <- nrow(dt)
    dt_after_lipid_treated_exclusion <- dt[
      is.na(Lipid_Medicine_baseline) | Lipid_Medicine_baseline == 0L
    ]

    flow_log_detail <- add_flow_step(
      flow_log_detail,
      "S6. 排除已知近期降脂药使用者；状态未知者暂留至 required-information gate",
      n_before,
      nrow(dt_after_lipid_treated_exclusion),
      n_start
    )

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 6. 排除近期降脂药使用者（Exclude recent lipid-lowering treatment）",
      n_before,
      nrow(dt_after_lipid_treated_exclusion),
      n_start
    )

    step_message(
      "纳排步骤：Step 6. 排除已知近期降脂药使用者",
      " | 排除 = ", format(n_before - nrow(dt_after_lipid_treated_exclusion), big.mark = ","),
      " | 状态未知暂留 = ",
      format(sum(is.na(dt_after_lipid_treated_exclusion$Lipid_Medicine_baseline)), big.mark = ","),
      " | 剩余 = ", format(nrow(dt_after_lipid_treated_exclusion), big.mark = ",")
    )

    list(
      dt = dt_after_lipid_treated_exclusion,
      flow_log_main = flow_log_main,
      flow_log_detail = flow_log_detail,
      lipid_source_distribution = dt[, .N, by = lipid_status_source][order(-N)]
    )
  },
  "匹配降脂药并执行用药纳排"
)

data_03 <- res_lipid_screen$value$dt
flow_log_main <- res_lipid_screen$value$flow_log_main
flow_log_detail <- res_lipid_screen$value$flow_log_detail

step_message("降脂药状态来源分布：")
print(res_lipid_screen$value$lipid_source_distribution)
step_message("降脂药纳排后样本量：", format(nrow(data_03), big.mark = ","))


# ============================================================================== #
# 4. 排除显著升高的单一危险因素 ####
# ============================================================================== #
# 排除条件：
#   - TC > 8 mmol/L；
#   - LDL-C > 4.9 mmol/L；
#   - SBP >= 180 mmHg 或 DBP >= 110 mmHg。
#
# 这部分人群后续可在敏感性分析中重新纳入。
# ------------------------------------------------------------------------------ #

res_markedly_high <- run_with_timer(
  {
    dt <- copy(data_03)
    setDT(dt)

    n_before <- nrow(dt)
    dt <- derive_markedly_elevated_risk_status(dt)
    dt <- dt[markedly_elevated_single_risk_factor == 0L]

    flow_log_main <- add_flow_step(
      flow_log_main,
      "Step 7. 排除显著升高的单一危险因素（TC>8, LDL-C>4.9, or BP>=180/110）",
      n_before,
      nrow(dt),
      n_start
    )
    flow_log_detail <- add_flow_step(flow_log_detail, "S7. 排除 TC>8、LDL-C>4.9 或 BP>=180/110", n_before, nrow(dt), n_start)

    step_message(
      "纳排步骤：Step 7. 排除 TC>8、LDL-C>4.9 或 BP>=180/110",
      " | 排除 = ", format(n_before - nrow(dt), big.mark = ","),
      " | 剩余 = ", format(nrow(dt), big.mark = ",")
    )

    list(
      dt = dt,
      flow_log_main = flow_log_main,
      flow_log_detail = flow_log_detail
    )
  },
  "排除显著升高的单一危险因素"
)

data_04 <- res_markedly_high$value$dt
flow_log_main <- res_markedly_high$value$flow_log_main
flow_log_detail <- res_markedly_high$value$flow_log_detail


# ============================================================================== #
# 5. Required-information gate + subsequently-imputed-variable missingness ####
# ============================================================================== #
# 这一节是 04 的最后一个纳排步骤。
#
# 5.1 required observed information
#   只有在缺失时无法确定 eligibility 或 SCORE/SCORE2 / LDL-C strategy 的变量，
#   才允许在这里造成排除。普通 outcome-model / subgroup 协变量不参与 gate。
#
# 5.2 variables subsequently imputed
#   Region、Ethnicity、TDI、SES 组成变量、BMI、raw MET、TG、降压药、家族史、
#   慢性呼吸病和慢性肝病等变量的缺失只做描述，保留进入 06 MI。
# ------------------------------------------------------------------------------ #

res_required_gate <- run_with_timer(
  {
    dt_pre_gate <- copy(data_04)
    setDT(dt_pre_gate)

    # 保留降压药真实观测状态，但不因其缺失而排除。
    require_columns(dt_pre_gate, "antihypertensive_med_baseline", "03 统一降压药证据")
    dt_pre_gate[, hyper_Medicine_baseline := as.integer(antihypertensive_med_baseline)]
    dt_pre_gate[, hyper_Medicine := hyper_Medicine_baseline]
    dt_pre_gate[, hyper_Medicine_time := as.IDate(fifelse(
      !is.na(hyper_Medicine_baseline), as.IDate(InCohort_time), as.IDate(NA)
    ))]
    dt_pre_gate[, hyper_status_known := as.integer(!is.na(hyper_Medicine_baseline))]

    # 这些变量决定 eligibility 或 guideline-strategy assignment。
    # 注意：Region（assessment-centre region）不等于 SCORE risk_region；当前 SCORE
    # risk_region 在 07 固定，因此 Region 仅作为 outcome-model covariate，可进入 MI。
    vars_required_observed <- c(
      "Participant ID",
      "InCohort_time",
      "Age_at_recruitment",
      "Sex",
      "Smoke0",
      "Cholesterol_0",
      "LDL_direct_0",
      "HDL_cholesterol_0",
      "SystolicBP0",
      "DiastolicBP0",
      "Glycated_haemoglobin_(HbA1c)_0",
      "Creatinine_0",
      "egfr_2021_baseline",
      "documented_ascvd_baseline",
      "diabetes_baseline",
      "CKD_baseline",
      "Lipid_Medicine_baseline",
      "markedly_elevated_single_risk_factor"
    )

    gate_res <- apply_required_observed_information_screen(
      dt = dt_pre_gate,
      vars_required_observed = vars_required_observed,
      context = "Complete information required for eligibility and guideline-strategy assignment"
    )

    flow_log_main2 <- add_flow_step(
      flow_log_main,
      "Step 8. Complete information required for eligibility and guideline-strategy assignment",
      gate_res$n_before,
      gate_res$n_after,
      n_start
    )
    flow_log_detail2 <- add_flow_step(
      flow_log_detail,
      paste0(
        "S8. eligibility / guideline-strategy 必要信息完整：",
        "Age, Sex, smoking, TC, LDL-C, HDL-C, SBP, DBP, HbA1c, creatinine/eGFR, ",
        "baseline ASCVD/DM/CKD status and recent lipid-lowering status"
      ),
      gate_res$n_before,
      gate_res$n_after,
      n_start
    )

    final_dt <- gate_res$dt

    # 06 正式插补变量。这里仅报告缺失，不作排除。
    imputation_variable_methods <- data.table(
      variable = c(
        "Region",
        "Ethnic",
        "TDI",
        "Education",
        "Household_income",
        "employment_status",
        "BMI0",
        "MET_minutes_week0",
        "alcohol_score_0",
        "diet_score_w0",
        "Triglycerides_0",
        "hyper_Medicine_baseline",
        "Family_history_CVD",
        "chronic_respiratory_disease",
        "chronic_liver_disease"
      ),
      imputation_method = c(
        "polyreg",
        "logreg",
        "pmm",
        "polr",
        "polr",
        "polyreg",
        "pmm",
        "pmm",
        "logreg",
        "logreg",
        "pmm",
        "logreg",
        "logreg",
        "logreg",
        "logreg"
      )
    )

    imputation_missingness <- build_imputation_missingness_summary(
      pre_gate_dt = dt_pre_gate,
      final_dt = final_dt,
      variable_methods = imputation_variable_methods
    )

    list(
      dt = final_dt,
      pre_gate_dt = dt_pre_gate,
      flow_log_main = flow_log_main2,
      flow_log_detail = flow_log_detail2,
      vars_required_observed = gate_res$vars_required_observed,
      required_missing_summary = gate_res$missing_summary,
      imputation_missingness = imputation_missingness,
      imputation_variable_methods = imputation_variable_methods
    )
  },
  "执行 eligibility / guideline-strategy required-information gate"
)

data_eligible <- res_required_gate$value$dt
flow_log_main <- res_required_gate$value$flow_log_main
flow_log_detail <- res_required_gate$value$flow_log_detail
required_information_missing_summary <- res_required_gate$value$required_missing_summary
covariate_missing_summary <- res_required_gate$value$imputation_missingness

# 兼容既有下游标记：最终样本中降脂药状态已知且均为未使用。
data_eligible[, lipid_main_analysis_included := 1L]
# hyper_main_analysis_included 仅表示“该人仍在主分析队列”，不再表示降压药状态完整。
data_eligible[, hyper_main_analysis_included := 1L]
rm_missing <- data.table(variable=grep("^rm_source_",names(data_eligible),value=TRUE))
rm_missing[, `:=`(n=nrow(data_eligible),
  n_missing=vapply(variable,function(v)sum(is.na(data_eligible[[v]])),integer(1)))]
rm_missing[, missing_fraction:=n_missing/n]
fwrite(rm_missing,file.path(quality_control_dir,"RM_source_missingness_eligible_M01.csv"),bom=TRUE)

validate_flow_log(flow_log_main, "04 主流程")
validate_flow_log(flow_log_detail, "04 详细流程")
eligibility_invariant_audit <- build_eligibility_invariant_audit(data_eligible)

required_gate_n_excluded <- nrow(res_required_gate$value$pre_gate_dt) - nrow(data_eligible)
step_message(
  "Required-information gate：排除 = ",
  format(required_gate_n_excluded, big.mark = ","),
  "；最终 analytic cohort = ", format(nrow(data_eligible), big.mark = ",")
)
step_message("随后插补变量缺失概况：")
print(covariate_missing_summary)



# ============================================================================== #
# 6. 纳排流程表输出 ####
# 模块任务：写出纳排流程
# 模块输入：每步计数
# 模块输出：人群纳排记录_M01
# ============================================================================== #

step_message("\n========== 主文 Figure 1 纳排流程表 ==========")
print(flow_log_main)

step_message("\n========== 补充材料详细纳排表 ==========")
print(flow_log_detail)


# ============================================================================== #
# 7. 简单质控输出 ####
# ============================================================================== #

step_message("\n==================== 质控输出 ====================")
step_message("最终 data_eligible 样本量：", format(nrow(data_eligible), big.mark = ","))
step_message("唯一 ID 数：", format(uniqueN(data_eligible$`Participant ID`), big.mark = ","))
step_message("最终降脂药 treated：", data_eligible[, format(sum(Lipid_Medicine_baseline == 1L, na.rm = TRUE), big.mark = ",")])
step_message("最终降脂药 untreated：", data_eligible[, format(sum(Lipid_Medicine_baseline == 0L, na.rm = TRUE), big.mark = ",")])
step_message("最终降压药 treated：", data_eligible[, format(sum(hyper_Medicine_baseline == 1L, na.rm = TRUE), big.mark = ",")])
step_message("最终降压药 untreated：", data_eligible[, format(sum(hyper_Medicine_baseline == 0L, na.rm = TRUE), big.mark = ",")])
step_message("最终降压药 status missing（保留供 MI）：", data_eligible[, format(sum(is.na(hyper_Medicine_baseline)), big.mark = ",")])
step_message("最终样本中 365 天内 GP 降脂药处方人数：", data_eligible[, format(sum(lipid_gp_365d == 1L, na.rm = TRUE), big.mark = ",")])

step_message("最终样本 Education 分布：")
print(data_eligible[, .N, by = Education][order(-N)])

step_message("最终样本 Household_income 分布：")
print(data_eligible[, .N, by = Household_income][order(-N)])

step_message("最终样本 employment_status 分布：")
print(data_eligible[, .N, by = employment_status][order(-N)])

step_message("最终样本 Family_history_CVD 分布：")
print(data_eligible[, .N, by = Family_history_CVD][order(Family_history_CVD)])


# ============================================================================== #
# 8. 保存重要数据 ####
# 模块任务：保存最终合格队列
# 模块输入：通过全部排除检查的样本
# 模块输出：data_eligible_M01.csv
# ============================================================================== #

res_save <- run_with_timer(
  {
    dt <- copy(data_eligible)
    setDT(dt)

    base_clean_dir_M01 <- data_clean_dir
    if (!dir.exists(base_clean_dir_M01)) {
      dir.create(base_clean_dir_M01, recursive = TRUE)
    }
    if (!dir.exists(quality_control_dir)) {
      dir.create(quality_control_dir, recursive = TRUE)
    }

    flow_dir <- file.path(base_clean_dir_M01, "人群纳排记录_M01")
    if (!dir.exists(flow_dir)) {
      dir.create(flow_dir, recursive = TRUE)
    }

    fwrite(
      flow_log_main,
      file = file.path(flow_dir, "flow_log_main_基线干预_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      flow_log_detail,
      file = file.path(flow_dir, "flow_log_detail_基线干预_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      flow_log_main,
      file = file.path(flow_dir, "flow_log_基线干预_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      covariate_missing_summary,
      file = file.path(flow_dir, "covariate_missing_summary_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      required_information_missing_summary,
      file = file.path(flow_dir, "required_information_missing_summary_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      merge_audit_04,
      file = file.path(quality_control_dir, "merge_audit_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    fwrite(
      eligibility_invariant_audit,
      file = file.path(quality_control_dir, "eligibility_invariant_audit_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    # 与后续脚本保持兼容：输出前统一把列名中的空格替换为下划线。
    setnames(dt, names(dt), gsub(" ", "_", names(dt), fixed = TRUE))

    variables_to_delete <- intersect(
      c(
        "Participant_ID(gp_scripts_-_eid)",
        "baseline_time",
        "lipid_status_known",
        "hyper_status_known"
      ),
      names(dt)
    )
    if (length(variables_to_delete) > 0L) {
      dt[, (variables_to_delete) := NULL]
    }

    fwrite(
      dt,
      file = file.path(base_clean_dir_M01, "data_eligible_M01.csv"),
      row.names = FALSE,
      bom = TRUE
    )

    TRUE
  },
  "保存 data_eligible 与 flow_log"
)

step_message("脚本运行完成。")

m01_runtime_finish()
