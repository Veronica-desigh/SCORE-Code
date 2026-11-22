# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 03_Clean_medicine_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   清理 M01 ASCVD 主分析所需的降脂、降压和降糖用药信息。
#
# v1.0 修改要点：
#   1. 将 03 专用函数拆分到 Function_03_Clean_medicine_M01.R；
#   2. 自报用药按 touchscreen 问卷多选项 token 识别，输出降脂、降压和胰岛素自报；
#   3. Field 20003 verbal interview 用药名称用于识别降糖药，和 touchscreen 胰岛素自报合并；
#   4. GP scripts 仅用于降脂药处方证据，按 BNF / dm+d / Read v2 任一编码命中判定；
#   5. 主分析基线状态使用 Instance 0 和基线前 0-365 天 GP 处方；90 天标记仅作辅助审计；
#   6. 输出统一用药证据表和来源审计，所有参与者级输出保持一人一行。
#   7. 基线降压药缺失保持为 NA，不在 03 中改写为 0，也不在 03 中决定排除；
#      该缺失状态在 04 保留到最终 required-information gate，并在符合主分析
#      纳入条件者中作为后续 MI 变量处理。
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
# 模块任务：配置与用药单元测试
# 模块输入：药物函数与工作目录
# 模块输出：可复现运行环境
# ============================================================================== #

project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
data_dir <- Sys.getenv("M01_UKB_DIR", file.path(project_root, "DATA", "UKB"))
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
quality_control_dir <- file.path(data_clean_dir, "quality_control")

dir.create(data_clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(quality_control_dir, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(parallel)

setwd(project_root)

source(file.path(code_dir, "Function_03_Clean_medicine_M01.R"))
source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("03", data_clean_dir)

run_medication_unit_tests()
step_message("03 用药边界单元测试通过：Instance 1 不回填；365天边界纳入；第366天排除；90天单独记录供审计；GP 聚合后一人一行。")

m01_configure_threads()
step_message("data.table 线程数：", data.table::getDTthreads())

# ============================================================================== #
# 1. 输入文件路径和字段声明 ####
# 模块任务：配置药名与代码映射
# 模块输入：UKB GP及自报字段
# 模块输出：识别规则
# ============================================================================== #

data_exposure_path <- file.path(data_clean_dir, "data_exposure_M01.csv")

touchscreen_med_file <- file.path(
  data_dir,
  "3_Assessment center",
  "3_Assessment centre_2_Touchscreen_5_Psychosocial factors2.csv"
)

verbal_med_file <- file.path(
  data_dir,
  "3_Assessment center",
  "3_Assessment centre_4_Verbal interyiew__4_Medications2.csv"
)

gp_scripts_file <- file.path(
  data_dir,
  "8_Health-related outcomes",
  "8_2Primary care",
  "8_Health-related outcomes_2_Primary care_data_gp_scripts.csv"
)

# Touchscreen 自报用药字段：
#   男性字段 6177：降脂药、降压药、胰岛素；
#   女性字段 6153：同上，另含 HRT 和口服避孕药选项。
touchscreen_male_columns <- setNames(
  paste0("Medication for cholesterol, blood pressure or diabetes | Instance ", 0:3),
  as.character(0:3)
)

touchscreen_female_columns <- setNames(
  paste0(
    "Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones | Instance ",
    0:3
  ),
  as.character(0:3)
)

# Field 20003 为 verbal interview 中的药物名称文本。
# 这里使用降糖药通用名/药名片段识别是否存在降糖用药证据，不区分剂量、剂型和频次。
glucose_lowering_med_names <- c(
  "metformin",
  "insulin product",
  "gliclazide",
  "glimepiride",
  "glibenclamide",
  "glipizide",
  "tolbutamide",
  "chlorpropamide",
  "acarbose",
  "pioglitazone",
  "rosiglitazone",
  "repaglinide",
  "nateglinide",
  "sitagliptin",
  "vildagliptin",
  "saxagliptin",
  "linagliptin",
  "alogliptin",
  "dapagliflozin",
  "empagliflozin",
  "canagliflozin",
  "liraglutide",
  "exenatide",
  "semaglutide",
  "dulaglutide",
  "lixisenatide",
  "rosiglitazone 1mg / metformin 500mg tablet"
)

glucose_lowering_med_patterns <- c(
  "\\binsulin\\b",
  "\\bmetformin\\b",
  "\\bgliclazide\\b",
  "\\bglimepiride\\b",
  "\\bglibenclamide\\b",
  "\\bglipizide\\b",
  "\\btolbutamide\\b",
  "\\bchlorpropamide\\b",
  "\\bacarbose\\b",
  "\\bpioglitazone\\b",
  "\\brosiglitazone\\b",
  "\\brepaglinide\\b",
  "\\bnateglinide\\b",
  "\\bsitagliptin\\b",
  "\\bvildagliptin\\b",
  "\\bsaxagliptin\\b",
  "\\blinagliptin\\b",
  "\\balogliptin\\b",
  "\\bdapagliflozin\\b",
  "\\bempagliflozin\\b",
  "\\bcanagliflozin\\b",
  "\\bliraglutide\\b",
  "\\bexenatide\\b",
  "\\bsemaglutide\\b",
  "\\bdulaglutide\\b",
  "\\blixisenatide\\b"
)

# GP scripts 降脂药结构化编码。
# 注意：同一条处方可能同时带有 BNF、dm+d 和 Read v2 编码，筛选时按任一编码命中判定。
lipid_bnf_codes <- c(
  "212",
  "02120400",
  "0212000Y0AAAAAA", "0212000Y0AAADAD", "0212000Y0AAABAB",
  "02.12.04.00.00",
  "0212000B0AAAAAA", "0212000B0AAABAB", "0212000B0AAACAC",
  "0212000X0AAAAAA", "0212000X0AAABAB", "0212000X0AAADAD",
  "1090200",
  "01.09.02.00.00",
  "02.12.01.00.00",
  "2120000",
  "02.12.02.00.00",
  "02.12.03.00.00"
)

lipid_dmd_codes <- c(
  "408060004",
  "320030001", "320031002", "320029006",
  "319996000", "319997009", "320000009", "134489001",
  "320012008", "320013003", "320014009",
  "320023007", "320025000",
  "408024009", "408036003", "408037007", "409108001",
  "4090511000001100"
)

lipid_read2_codes <- c(
  "bxd1.", "bxd1.00", "bxd2.", "bxd2.00", "bxd5.", "bxd5.00",
  "bxi1.00", "bxi3.00", "bxiz.00", "bxi2.00",
  "bxe5.00", "bxe6.00", "bxe7.00",
  "bxf1.00", "bxf2.00",
  "bxg1.00", "bxg2.00",
  "bx43.00", "bx44.00", "bx23.00"
)

# ============================================================================== #
# 2. 读取 M01 基线人群、基线日期和性别 ####
# 模块任务：建立基线用药主表
# 模块输入：01产生的基线暴露
# 模块输出：ID/性别/基线日期
# ============================================================================== #

res_base <- run_with_timer(
  {
    read_baseline_medicine_base(data_exposure_path)
  },
  "读取 M01 基线人群、基线时间和性别"
)

baseline_dt <- res_base$value
var_info_print(baseline_dt, c("Participant ID", "baseline_time", "Sex"))
step_message("M01 基线人群数：", format(nrow(baseline_dt), big.mark = ","))

# ============================================================================== #
# 3. Touchscreen 自报：降脂药、降压药和胰岛素 ####
# 模块任务：解析基线用药自报
# 模块输入：触屏问卷
# 模块输出：脂质/降压/胰岛素标记
# ============================================================================== #

res_touchscreen <- run_with_timer(
  {
    build_touchscreen_med_self_report(
      touchscreen_file = touchscreen_med_file,
      baseline_dt = baseline_dt,
      male_columns = touchscreen_male_columns,
      female_columns = touchscreen_female_columns
    )
  },
  "整理 touchscreen 自报降脂药、降压药和胰岛素"
)

Lipid_med_self_report <- res_touchscreen$value$Lipid_med_self_report
hyper_med_self_report <- res_touchscreen$value$hyper_med_self_report
insulin_self_report <- res_touchscreen$value$insulin_self_report

step_message(
  "降脂药自报 instance0 可判定人数：",
  format(Lipid_med_self_report[, sum(!is.na(Lipid_medication_self_report_0))], big.mark = ",")
)
step_message(
  "降脂药自报 instance0 用药人数：",
  format(Lipid_med_self_report[, sum(Lipid_medication_self_report_0 == 1L, na.rm = TRUE)], big.mark = ",")
)
step_message(
  "降压药自报 instance0 可判定人数：",
  format(hyper_med_self_report[, sum(!is.na(hyper_medication_self_report_0))], big.mark = ",")
)
step_message(
  "降压药自报 instance0 用药人数：",
  format(hyper_med_self_report[, sum(hyper_medication_self_report_0 == 1L, na.rm = TRUE)], big.mark = ",")
)
step_message(
  "胰岛素自报 instance0 用药人数：",
  format(insulin_self_report[, sum(insulin_self_report_0 == 1L, na.rm = TRUE)], big.mark = ",")
)

fwrite(Lipid_med_self_report, file.path(data_clean_dir, "Lipid_med_self_report_M01.csv"))
fwrite(hyper_med_self_report, file.path(data_clean_dir, "hyper_med_self_report_M01.csv"))
fwrite(insulin_self_report, file.path(data_clean_dir, "insulin_self_report_M01.csv"))

# ============================================================================== #
# 4. Field 20003 verbal interview：降糖用药 ####
# 模块任务：解析口述降糖药
# 模块输入：field20003
# 模块输出：降糖药证据
# ============================================================================== #

res_glucose_20003 <- run_with_timer(
  {
    build_20003_medication_flag(
      medication_file = verbal_med_file,
      instance = 0L,
      max_array = 47L,
      target_names = glucose_lowering_med_names,
      target_patterns = glucose_lowering_med_patterns,
      flag_prefix = "glucose_lowering_med_20003"
    )
  },
  "整理 Field 20003 verbal interview 降糖用药"
)

glucose_med_20003 <- res_glucose_20003$value

res_glucose_baseline <- run_with_timer(
  {
    build_glucose_med_baseline(
      baseline_dt = baseline_dt,
      glucose_20003_dt = glucose_med_20003,
      insulin_self_report_dt = insulin_self_report,
      instance = 0L
    )
  },
  "合并 Field 20003 降糖药和 touchscreen 胰岛素自报"
)

Glucose_med_self_report <- res_glucose_baseline$value

step_message(
  "Field 20003 降糖用药人数：",
  format(glucose_med_20003[, sum(glucose_lowering_med_20003_0 == 1L, na.rm = TRUE)], big.mark = ",")
)
step_message(
  "合并后基线降糖用药证据人数：",
  format(Glucose_med_self_report[, sum(glucose_lowering_med_baseline == 1L, na.rm = TRUE)], big.mark = ",")
)
step_message("降糖用药证据来源分布：")
print(Glucose_med_self_report[, .N, by = glucose_lowering_med_source][order(-N)])

fwrite(glucose_med_20003, file.path(data_clean_dir, "Glucose_med_20003_M01.csv"))
fwrite(Glucose_med_self_report, file.path(data_clean_dir, "Glucose_med_self_report_M01.csv"))

# ============================================================================== #
# 5. GP scripts：降脂药处方证据 ####
# 模块任务：整合GP处方
# 模块输入：GP scripts、BNF/dm+d/Read
# 模块输出：降脂药处方标记
# ============================================================================== #

res_gp_pipeline <- run_with_timer(
  {
    gp_raw <- m01_timed("03_GP_read_and_early_exact_code_prefilter",
      clean_gp_scripts(gp_scripts_file,
        lipid_prefilter_bnf = lipid_bnf_codes,
        lipid_prefilter_dmd = lipid_dmd_codes,
        lipid_prefilter_read2 = lipid_read2_codes))
    m01_checkpoint("03_GP_preprocessing", 2L, 3L, "dates parsed for qualifying Rx")

    lipid_gp_raw <- filter_lipid_gp(
      gp_raw,
      bnf_codes = lipid_bnf_codes,
      dmd_codes = lipid_dmd_codes,
      read2_codes = lipid_read2_codes
    )

    lipid_gp_history <- build_lipid_history(lipid_gp_raw)
    m01_checkpoint("03_GP_preprocessing", 3L, 3L,
      sprintf("%s lipid prescriptions identified", format(nrow(lipid_gp_history), big.mark = ",")))

    gp_summary <- m01_timed("03_GP_one_join_nearest_and_windows",
      build_prebaseline_gp_summaries(
        history_dt = lipid_gp_history,
        baseline_dt = baseline_dt
      ))
    lipid_gp_nearest_prebaseline <- gp_summary$nearest
    lipid_gp_window <- gp_summary$window

    lipid_gp_main <- build_lipid_gp_main(
      nearest_dt = lipid_gp_nearest_prebaseline,
      window_dt = lipid_gp_window,
      baseline_dt = baseline_dt
    )

    list(
      lipid_gp_history = lipid_gp_history,
      lipid_gp_nearest_prebaseline = lipid_gp_nearest_prebaseline,
      lipid_gp_window = lipid_gp_window,
      lipid_gp_main = lipid_gp_main
    )
  },
  "构建 GP scripts 降脂药处方证据"
)

lipid_gp_history <- standardize_output_id(res_gp_pipeline$value$lipid_gp_history)
lipid_gp_nearest_prebaseline <- standardize_output_id(res_gp_pipeline$value$lipid_gp_nearest_prebaseline)
lipid_gp_window <- standardize_output_id(res_gp_pipeline$value$lipid_gp_window)
lipid_gp_main <- standardize_output_id(res_gp_pipeline$value$lipid_gp_main)

fwrite(lipid_gp_history, file.path(data_clean_dir, "Lipid_med_GP_history_M01.csv"))
fwrite(lipid_gp_nearest_prebaseline, file.path(data_clean_dir, "Lipid_med_GP_nearest_prebaseline_M01.csv"))
fwrite(lipid_gp_window, file.path(data_clean_dir, "Lipid_med_GP_window_summary_M01.csv"))
fwrite(lipid_gp_main, file.path(data_clean_dir, "Lipid_med_GP_main_M01.csv"))

step_message("完整 GP 降脂药处方记录数：", format(nrow(lipid_gp_history), big.mark = ","))
step_message("有基线前/当天 GP 降脂药处方人数：", format(nrow(lipid_gp_nearest_prebaseline), big.mark = ","))
step_message("GP 降脂药主表人数：", format(nrow(lipid_gp_main), big.mark = ","))
step_message(
  "365 天窗口内有 GP 降脂药处方人数（主分析）：",
  format(lipid_gp_main[, sum(gp_lipid_rx_in_365d == 1L, na.rm = TRUE)], big.mark = ",")
)
step_message("GP 降脂药编码来源分布：")
print(lipid_gp_history[, .N, by = lipid_match_source][order(-N)])

# ============================================================================== #
# 6. 统一基线用药证据和来源审计 ####
# 模块任务：汇总来源并写出
# 模块输入：上述用药来源
# 模块输出：data_medicine_M01.csv及来源审计
# ============================================================================== #

data_medicine <- build_medication_evidence_master(
  baseline_dt = baseline_dt,
  lipid_self_report_dt = Lipid_med_self_report,
  hyper_self_report_dt = hyper_med_self_report,
  glucose_baseline_dt = Glucose_med_self_report,
  lipid_gp_main_dt = lipid_gp_main
)
art_evidence <- build_hiv_art(verbal_med_file)
art_index <- match(as.character(data_medicine[["Participant ID"]]),as.character(art_evidence[["Participant ID"]]))
if(anyNA(art_index))stop("ART source does not cover medication cohort")
for(v in setdiff(names(art_evidence),"Participant ID"))set(data_medicine,j=v,value=art_evidence[[v]][art_index])
fwrite(art_evidence,file.path(quality_control_dir,"HIV_specific_ART_evidence_M01.csv"),bom=TRUE)
medication_source_audit <- build_medication_source_audit(data_medicine)

assert_unique_id(data_medicine, "03 统一基线用药证据表")
if (nrow(data_medicine) != nrow(baseline_dt)) {
  stop("03 统一基线用药证据表未保持完整基线人群。", call. = FALSE)
}

fwrite(data_medicine, file.path(data_clean_dir, "data_medicine_M01.csv"))
fwrite(
  medication_source_audit,
  file.path(quality_control_dir, "medication_source_audit_M01.csv")
)

step_message("统一基线用药证据表人数：", format(nrow(data_medicine), big.mark = ","))
step_message("统一基线用药证据表变量数：", ncol(data_medicine))
step_message("用药来源审计目录：", quality_control_dir)

step_message("03 用药信息清理完成。")

m01_runtime_finish()
