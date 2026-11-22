# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 06_Multiple_imputation_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   Branch B：对 05 输出的 incomplete analysis base 执行 multiple imputation。
#   本轮按用户要求固定 M=5。
#
# 输出：
#   1. mids_M01.rds
#   2. MI_completed/MI_01_analysis_M01.rds ... MI_05_analysis_M01.rds
#   3. MI method / predictor matrix / missingness / observed-vs-imputed / completed-data / LCA / subgroup QA
#
# 重要：本脚本绝不计算 SCORE/SCORE2，也绝不生成 guideline strategy。
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

library(data.table)

project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
code_dir <- file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")
data_clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD"))
rm_report_dir <- file.path(data_clean_dir,"风险修饰因素定义与处理记录_M01")
base_report_dir <- file.path(data_clean_dir, "分析基础数据与随访审计_M01")
mi_dir <- file.path(data_clean_dir, "MI_M01")
mi_completed_dir <- file.path(mi_dir, "MI_completed")
mi_qc_dir <- file.path(mi_dir, "quality_control")

dir.create(mi_completed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mi_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rm_report_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(code_dir, "Function_06_Multiple_imputation_M01.R"), encoding = "UTF-8")

source(file.path(code_dir, "Function_00_Runtime_M01.R"), encoding = "UTF-8")
m01_runtime_init("06", data_clean_dir)

if (identical(Sys.getenv("M01_EXPORT_ONLY"),"1")) {
  mi06_export_workbooks(data_clean_dir,code_dir)
  quit(save="no",status=0)
}

if (!requireNamespace("mice", quietly = TRUE)) {
  stop("06 需要 R 包 mice：install.packages('mice')。", call. = FALSE)
}
if (!requireNamespace("poLCA", quietly = TRUE)) {
  stop("06 需要 R 包 poLCA：install.packages('poLCA')。", call. = FALSE)
}

# ------------------------------------------------------------------------------ #
# 0. 冻结参数
# 模块任务：配置MI参数
# 模块输入：05不完整基线数据、方法字典
# 模块输出：固定M、maxit、随机种子
# ------------------------------------------------------------------------------ #
M <- as.integer(Sys.getenv("M01_MI_M","5"))
maxit <- as.integer(Sys.getenv("M01_MI_MAXIT","5"))
mi_seed <- 20260916L
lca_nrep <- as.integer(Sys.getenv("M01_LCA_NREP","30"))
lca_maxiter <- as.integer(Sys.getenv("M01_LCA_MAXITER","5000"))
lca_seed_base <- 205000L

input_rds <- file.path(data_clean_dir, "data_analysis_base_incomplete_M01.rds")
dictionary_path <- file.path(base_report_dir, "variable_dictionary_M01.csv")
mids_path <- file.path(mi_dir, "mids_M01.rds")

mi06_assert_file(input_rds, "05 incomplete analysis base")
mi06_assert_file(dictionary_path, "05 variable dictionary")

mi06_step_message("Branch B 开始：M=", M, "；maxit=", maxit, "；seed=", mi_seed)

# ------------------------------------------------------------------------------ #
# 1. 读取 05 incomplete analysis base
# 模块任务：校验05输入和禁入变量
# 模块输入：05 RDS/CSV
# 模块输出：预插补分析数据
# ------------------------------------------------------------------------------ #
data_analysis_base <- readRDS(input_rds)
data_analysis_base <- as.data.table(data_analysis_base)
variable_dictionary <- fread(dictionary_path)
mi06_assert_unique_id(data_analysis_base)
reference_ids <- as.character(data_analysis_base$participant_id)
reference_n <- nrow(data_analysis_base)

# 禁止 strategy 变量提前进入 MI branch。
forbidden_strategy_vars <- c(
  "score_2019_percent", "score2_2025_percent",
  "score_2019_category", "score2_2025_category",
  "z_common", "ldl_target_2019", "ldl_target_2025",
  "intervention_trigger_2019", "intervention_trigger_2025",
  "ldl_2019", "ldl_2025"
)
if (any(forbidden_strategy_vars %in% names(data_analysis_base))) {
  stop(
    "05 incomplete base 中出现 strategy 变量；Branch B 必须与 Branch A 解耦：",
    paste(intersect(forbidden_strategy_vars, names(data_analysis_base)), collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------ #
# 2. MICE specification
# 模块任务：构建MICE配置及插补前QC
# 模块输入：字典、rm06_prepare、predictor matrix
# 模块输出：mice输入、方法与缺失报告
# ------------------------------------------------------------------------------ #
rm_context <- m01_timed("06_RM_preparation", rm06_prepare(data_analysis_base))
fwrite(rm_context$reportability_audit,file.path(rm_report_dir,"RM_reportability_audit_M01.csv"),bom=TRUE)
spec <- m01_timed("06_MICE_predictor_matrix_and_methods", {
  result <- mi06_build_mice_spec(rm_context$data, variable_dictionary)
  rm06_extend_spec(result, rm_context)
})
# Fail if a duplicate/stale derived RM outcome became an independent MI target.
if (any(c("rm_any_2019","rm_any_2025","rm_albuminuria") %in% names(spec$method)))
  stop("Derived composite/removed RM must not enter MICE")
worker_setting <- Sys.getenv("M01_MI_WORKERS", "auto")
physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(physical_cores) || physical_cores < 1L) physical_cores <- 2L
if (tolower(worker_setting) == "auto") {
  mi_workers <- if (physical_cores >= 4L) 2L else 1L
} else {
  mi_workers <- suppressWarnings(as.integer(worker_setting))
}
mi06_step_message("MICE worker selection: ", mi_workers,
  " (M01_MI_WORKERS=", worker_setting,
  "; PSOCK workers each duplicate a substantial part of MI data in memory)")
# Progress display: on by default; safe to disable with M01_MI_PROGRESS=0.
# The parallel master prints the progress once per completed FCS iteration.
input_gb <- as.numeric(utils::object.size(spec$mice_input)) / 1024^3
mi06_step_message(sprintf("MICE input data size: %.2f GiB; %d worker(s) hold independent R copies; monitor RAM usage.",
                          input_gb, mi_workers))
mi_progress_enabled <- !identical(Sys.getenv("M01_MI_PROGRESS", "1"), "0")
mi_progress_poll <- suppressWarnings(as.numeric(
  Sys.getenv("M01_MI_PROGRESS_INTERVAL", "1")))
mi_progress_dir <- file.path(mi_qc_dir, "progress")
if (!mi_workers %in% 1:3 || M<1L || maxit<1L) stop("Invalid MI configuration")

method_table <- data.table(
  variable = names(spec$method),
  method = unname(spec$method),
  n_missing = vapply(spec$mice_input, function(x) sum(is.na(x)), integer(1))
)
predictor_matrix_table <- data.table(variable = rownames(spec$predictor_matrix))
predictor_matrix_table <- cbind(
  predictor_matrix_table,
  as.data.table(spec$predictor_matrix)
)

fwrite(method_table, file.path(mi_qc_dir, "imputation_method_M01.csv"))
fwrite(predictor_matrix_table, file.path(mi_qc_dir, "predictor_matrix_M01.csv"))
fwrite(spec$missingness, file.path(mi_qc_dir, "missingness_before_MI_M01.csv"))

# 模块 2A｜插补前实际MICE输入（与05源数据的缺失口径区分）
# 输入：spec$mice_input / spec$method；输出：MI质量控制中的输入缺失与分布。
fwrite(m01qa_missingness(spec$mice_input, names(spec$mice_input),"preMI_mice_input"),
  file.path(mi_qc_dir,"MI_input_preMI_missingness_M01.csv"),bom=TRUE)
fwrite(m01qa_distribution(spec$mice_input,names(spec$mice_input),"preMI_mice_input"),
  file.path(mi_qc_dir,"MI_input_preMI_distribution_M01.csv"),bom=TRUE)


# ------------------------------------------------------------------------------ #
# 3. 正式执行 MICE/FCS：M=5
# 模块任务：并行MICE插补
# 模块输入：完整spec与RM上下文
# 模块输出：mids_M01.rds和链诊断
# ------------------------------------------------------------------------------ #
mi06_step_message("进入正式 MICE：", M, " 套 × ", maxit,
                  " 轮；", mi_workers, " 个 worker；实时进度=",
                  if (mi_progress_enabled) "开启" else "关闭")
signature <- list(input=mi06_md5(input_rds),dictionary=mi06_md5(dictionary_path),
 function_code=mi06_model_hash(),
 definitions=mi06_md5(file.path(code_dir,"config","RM_domains_v1.7.csv")),M=M,maxit=maxit,seed=mi_seed)
signature_path <- file.path(mi_dir,"MI_input_signature_M01.rds")
model_snapshot_path <- file.path(mi_dir,"imputation_function_snapshot_M01.R")
saved_signature <- if(file.exists(signature_path))readRDS(signature_path)else NULL
if(!is.null(saved_signature) && file.exists(model_snapshot_path) && identical(saved_signature$function_code,mi06_md5(model_snapshot_path))) {
  # Upgrade an earlier v1.7 checkpoint only after comparing actual model bodies.
  mi06_check_model_snapshot(data_clean_dir,code_dir)
  saved_signature$function_code <- mi06_model_hash()
}
if(identical(Sys.getenv("M01_MI_RESUME"),"1") && file.exists(mids_path) && identical(saved_signature,signature)) {
 imp <- readRDS(mids_path)
 saveRDS(signature,signature_path)
} else {
 if(!file.copy(file.path(code_dir,"Function_06_Multiple_imputation_M01.R"),model_snapshot_path,overwrite=TRUE))stop("Cannot save model-function snapshot")
imp <- m01_timed("06_MICE_all_parallel_chains", rm06_parallel_mice(
  spec, rm_context, M, maxit, mi_seed, mi_workers, code_dir,
  progress_dir = mi_progress_dir,
  show_progress = mi_progress_enabled,
  poll_interval = mi_progress_poll
))
 saveRDS(signature,signature_path)
}

mi06_step_message("MICE 已结束；现在保存 mids 并执行 completed-data 派生与质量控制。")

if (imp$m != M) stop("mids 对象中的 M 与冻结参数不一致。", call. = FALSE)
saveRDS(imp, mids_path, compress = FALSE)

# 保存 MICE 原生 convergence objects；避免为了 CSV 强行改变其维度语义。
saveRDS(imp$chainMean, file.path(mi_qc_dir, "mice_chain_mean_M01.rds"), compress = FALSE)
saveRDS(imp$chainVar, file.path(mi_qc_dir, "mice_chain_variance_M01.rds"), compress = FALSE)
logged_events <- if (is.null(imp$loggedEvents)) {
  data.table(
    it = integer(), im = integer(), dep = character(),
    meth = character(), out = character()
  )
} else {
  as.data.table(imp$loggedEvents)
}
fwrite(logged_events, file.path(mi_qc_dir, "mice_logged_events_M01.csv"))

# ------------------------------------------------------------------------------ #
# 4. 逐套 completed dataset 回填到原 analysis base，并重新派生 PA/SES/HTN
# 模块任务：形成每套completed与RM
# 模块输入：MICE结果/05基线
# 模块输出：MI_completed数据和逐套描述
# ------------------------------------------------------------------------------ #
membership_list <- vector("list", M)
completed_imputed_only_list <- vector("list", M)
completed_audit_list <- vector("list", M)
strategy_input_audit_list <- vector("list", M)
manifest_list <- vector("list", M)
lca_model_summary_list <- vector("list", M)
lca_class_distribution_list <- vector("list", M)
lca_component_profile_list <- vector("list", M)
distribution_summary_list <- vector("list", M)
mi_post_missingness_list <- vector("list",M)
mi_post_distribution_list <- vector("list",M)
mi_pre_post_list <- vector("list",M)
rm_binary_audit_list <- vector("list",M)
# Standardized MI input is the correct before-state for MI change QA:
# 05 observed raw MET negative codes / zero assay values may be normalized to NA
# BEFORE MICE; comparing directly to 05 would incorrectly label them as changed.
mi_pre_model <- data.table::as.data.table(data.table::copy(spec$mice_input))
mi_pre_model[,participant_id := reference_ids]
rm_post_distribution_list <- vector("list",M)
rm_pre_post_list <- vector("list",M)
rm_source_post_distribution_list <- vector("list",M)
rm_pre_vars <- m01qa_rm_vars()
m01qa_assert_vars(data_analysis_base,paste0(rm_pre_vars,"_observed"),"05 observed RM baseline")
rm_before <- data_analysis_base[,c("participant_id",paste0(rm_pre_vars,"_observed")),with=FALSE]


imputed_vars <- names(spec$method)[spec$method != ""]
summary_variables <- unique(c(
  imputed_vars,
  "physical_activity", "ses_class", "hypertension"
))

for (m in seq_len(M)) {
  mi06_step_message("完成数据处理进度：MI ", sprintf("%02d", m), "/", M,
                    "（", round(100 * (m - 1L) / M, 1),
                    "% 已完成；当前进行 LCA、RM、QA 与磁盘写入）")

  completed_mice <- mice::complete(imp, action = m)
  completed_mice <- mi06_restore_completed_types(completed_mice)

  completed <- copy(rm_context$data)
  for (v in imputed_vars) {
    if (!v %in% names(completed_mice)) stop("mice complete 缺少变量：", v, call. = FALSE)
    data.table::set(completed, j = v, value = completed_mice[[v]])
  }

  derived <- m01_timed(sprintf("06_MI_%02d_LCA_and_subgroup_derivation", m),
    mi06_post_imputation_derivation(
      completed, imputation_id = m,
      lca_nrep = lca_nrep,
      lca_maxiter = lca_maxiter,
      lca_seed_base = lca_seed_base))
  completed <- rm06_finalize(derived$data,data_analysis_base,unique(data_analysis_base$rm_townsend_q80))
for(v in grep("^rm_.*_observed$",names(data_analysis_base),value=TRUE))
  completed[, (v) := data_analysis_base[[v]]]
# Continuous values must satisfy original reportability bounds.
for(v in names(rm_context$bounds)) {
 b<-rm_context$bounds[[v]]
 if(any(completed[[v]]<b$lower | completed[[v]]>b$upper,na.rm=TRUE))stop("Completed continuous RM violates bounds: ",v)
}

  # 模块 4A｜每套完成数据全变量缺失、重点变量分布、逐人配对变化
  # 输入：05基线数据和本次m completed（ID顺序必须相同）。
  # 输出：供跨M审计汇总的长表；不在循环中重复覆盖文件。
  mi_post_missingness_list[[m]] <- m01qa_missingness(completed,
    setdiff(names(completed),"participant_id"),"postMI_completed",m)
  mi_post_distribution_list[[m]] <- m01qa_distribution(completed,
    intersect(unique(c(imputed_vars,"physical_activity","ses_class","hypertension")),names(completed)),
    "postMI_completed",m)
  mi_pre_post_list[[m]] <- m01qa_pre_post(mi_pre_model,completed,
    intersect(imputed_vars,names(mi_pre_model)),m)
  rm_binary_audit_list[[m]] <- rm06_binary_rm_audit(completed,data_analysis_base,m)
  rm_post_distribution_list[[m]] <- m01qa_rm_binary(completed,rm_pre_vars,"postMI_completed",m)
  rm_source_post_distribution_list[[m]] <- m01qa_distribution(completed,
    intersect(rm06_imputation_targets(),names(completed)),"postMI_completed",m)
  rm_pre_post_list[[m]] <- m01qa_pre_post(data_analysis_base,completed,
    rm_pre_vars,m,pre_suffix="_observed")

  audit <- mi06_validate_completed_dataset(completed, reference_ids, m)
  completed_audit_list[[m]] <- audit

  # Branch A strategy-required variables必须在 MI 前后逐人完全一致。
  strategy_input_audit_list[[m]] <- mi06_validate_strategy_inputs_unchanged(
    original_dt = data_analysis_base,
    completed_dt = completed,
    imputation_id = m
  )

  membership_list[[m]] <- completed[, .(participant_id, ses_class, hypertension)]
  completed_imputed_only_list[[m]] <- completed[
    ,
    c("participant_id", imputed_vars),
    with = FALSE
  ]

  lca_model_summary <- copy(derived$lca$model_summary)
  lca_model_summary[, imputation_id := m]
  lca_model_summary_list[[m]] <- lca_model_summary

  lca_class_distribution <- copy(derived$lca$class_distribution)
  lca_class_distribution[, imputation_id := m]
  lca_class_distribution_list[[m]] <- lca_class_distribution

  lca_component_profile <- copy(derived$lca$component_profile)
  lca_component_profile[, imputation_id := m]
  lca_component_profile_list[[m]] <- lca_component_profile

  distribution_summary_list[[m]] <- mi06_summarize_completed_distribution(
    completed, m, summary_variables
  )

  out_path <- file.path(
    mi_completed_dir,
    sprintf("MI_%02d_analysis_M01.rds", m)
  )
  
  # 保存 RDS：正式分析文件
  m01_timed(sprintf("06_MI_%02d_write_RDS", m),
            saveRDS(completed, out_path, compress = FALSE))
  
  # 同步保存 CSV：便于人工查看
  csv_path <- file.path(
    mi_completed_dir,
    sprintf("MI_%02d_analysis_M01.csv", m)
  )
  
  if (!identical(Sys.getenv("M01_MI_EXPORT_CSV", "1"), "0")) {
    m01_timed(sprintf("06_MI_%02d_write_CSV", m),
      fwrite(completed, csv_path, na = "NA", bom = TRUE))
  }
  
  m01_checkpoint("06_completed_datasets", m, M,
    sprintf("imputation %02d RDS complete", m))
  manifest_list[[m]] <- data.table(
    imputation_id = m,
    n = nrow(completed),
    n_ascvd = sum(completed$ascvd_event == 1L),
    n_competing_death = sum(completed$competing_death_event == 1L),
    n_low_ses = sum(completed$ses_class == "low"),
    n_hypertension = sum(completed$hypertension == 1L),
    file_name = basename(out_path),
    path = normalizePath(out_path, winslash = "/", mustWork = TRUE),
    md5 = mi06_md5(out_path),
    mi_seed = mi_seed,
    lca_seed = lca_seed_base + m
  )
}

# ------------------------------------------------------------------------------ #
# 5. 跨 MI QA
# 模块任务：完成跨套质量控制及汇总
# 模块输入：每套MI完成数据
# 模块输出：MI质量控制文件与C目录RM前后对照
# ------------------------------------------------------------------------------ #
completed_audit <- rbindlist(completed_audit_list, use.names = TRUE, fill = TRUE)
strategy_input_audit <- rbindlist(strategy_input_audit_list, use.names = TRUE, fill = TRUE)
manifest <- rbindlist(manifest_list, use.names = TRUE, fill = TRUE)
lca_model_summary <- rbindlist(lca_model_summary_list, use.names = TRUE, fill = TRUE)
lca_class_distribution <- rbindlist(lca_class_distribution_list, use.names = TRUE, fill = TRUE)
lca_component_profile <- rbindlist(lca_component_profile_list, use.names = TRUE, fill = TRUE)
distribution_summary <- rbindlist(distribution_summary_list, use.names = TRUE, fill = TRUE)
mi_post_missingness <- rbindlist(mi_post_missingness_list,use.names=TRUE,fill=TRUE)
mi_post_distribution <- rbindlist(mi_post_distribution_list,use.names=TRUE,fill=TRUE)
mi_pre_post <- rbindlist(mi_pre_post_list,use.names=TRUE,fill=TRUE)
rm_binary_audit <- rbindlist(rm_binary_audit_list,use.names=TRUE,fill=TRUE)
rm_post_distribution <- rbindlist(rm_post_distribution_list,use.names=TRUE,fill=TRUE)
rm_source_post_distribution <- rbindlist(rm_source_post_distribution_list,use.names=TRUE,fill=TRUE)
rm_pre_post <- rbindlist(rm_pre_post_list,use.names=TRUE,fill=TRUE)
# 临床RM不得改变原始已知分类。普通MI也不得改变原始已观测目标值。
if(any(rm_pre_post$n_changed_observed>0L | rm_pre_post$n_new_missing>0L))
  stop("RM observed classifications changed or became missing in a completed dataset")
if(any(mi_pre_post$n_changed_observed>0L | mi_pre_post$n_new_missing>0L))
  stop("Original observed MI target values were changed or became missing")

stability <- mi06_build_subgroup_stability(membership_list)
observed_imputed_diagnostic <- mi06_compare_observed_imputed(
  # Use the same normalized before-state as the paired MI QA above: assay
  # zeros and negative MET codes were missing in the actual MICE input.
  original_dt = mi_pre_model,
  completed_list = completed_imputed_only_list,
  imputed_vars = imputed_vars
)

if (nrow(manifest) != M || any(manifest$n != reference_n)) {
  stop("5 套 completed datasets 样本量不一致。", call. = FALSE)
}
if (any(!completed_audit$passed)) stop("completed dataset QA 未全部通过。", call. = FALSE)
if (any(!strategy_input_audit$passed)) {
  stop("MI 改变了 Branch A strategy-required observed inputs。", call. = FALSE)
}

fwrite(manifest, file.path(mi_dir, "MI_completed_manifest_M01.csv"))
fwrite(completed_audit, file.path(mi_qc_dir, "completed_dataset_contract_audit_M01.csv"))
fwrite(strategy_input_audit, file.path(mi_qc_dir, "strategy_required_input_invariance_M01.csv"))
fwrite(
  observed_imputed_diagnostic,
  file.path(mi_qc_dir, "observed_vs_imputed_distribution_M01.csv")
)
fwrite(lca_model_summary, file.path(mi_qc_dir, "SES_LCA_model_summary_M01.csv"))
fwrite(lca_class_distribution, file.path(mi_qc_dir, "SES_LCA_class_distribution_M01.csv"))
fwrite(lca_component_profile, file.path(mi_qc_dir, "SES_LCA_component_profile_M01.csv"))
fwrite(distribution_summary, file.path(mi_qc_dir, "completed_value_distribution_M01.csv"))
# 模块 5A｜在完成五套插补和QA后统一写出新增描述表。
fwrite(mi_post_missingness,
  file.path(mi_qc_dir,"MI_completed_missingness_by_imputation_M01.csv"),bom=TRUE)
fwrite(mi_post_distribution,
  file.path(mi_qc_dir,"MI_completed_distribution_by_imputation_M01.csv"),bom=TRUE)
fwrite(mi_pre_post,
  file.path(mi_qc_dir,"MI_pre_post_missingness_by_imputation_M01.csv"),bom=TRUE)
fwrite(m01qa_cross_imputation_missingness(mi_pre_post),
  file.path(mi_qc_dir,"MI_pre_post_missingness_summary_M01.csv"),bom=TRUE)
fwrite(rm_binary_audit,
  file.path(rm_report_dir,"RM_binary_imputation_audit_M01.csv"),bom=TRUE)
fwrite(rm_source_post_distribution,
  file.path(rm_report_dir,"RM_postMI_source_distribution_by_imputation_M01.csv"),bom=TRUE)
fwrite(mi_pre_post[variable %in% rm06_imputation_targets()],
  file.path(rm_report_dir,"RM_input_pre_post_missingness_by_imputation_M01.csv"),bom=TRUE)
fwrite(rm_post_distribution,
  file.path(rm_report_dir,"RM_postMI_binary_distribution_by_imputation_M01.csv"),bom=TRUE)
# 跨插补RM阳性数变异仅用于稳定性审计，绝不将M*N视作独立个体。
rm_positive_variation <- rm_post_distribution[,.(
  n_imputations=.N,
  n_positive_min=min(n_positive),n_positive_max=max(n_positive),
  n_positive_mean=mean(n_positive),
  pct_positive_total_min=min(pct_positive_total),
  pct_positive_total_max=max(pct_positive_total),
  pct_positive_total_mean=mean(pct_positive_total)),by=variable]
fwrite(rm_positive_variation,
  file.path(rm_report_dir,"RM_postMI_positive_variation_M01.csv"),bom=TRUE)
fwrite(rm_pre_post,
  file.path(rm_report_dir,"RM_pre_post_missingness_by_imputation_M01.csv"),bom=TRUE)
fwrite(m01qa_cross_imputation_missingness(rm_pre_post),
  file.path(rm_report_dir,"RM_pre_post_missingness_summary_M01.csv"),bom=TRUE)

fwrite(stability$summary, file.path(mi_qc_dir, "subgroup_membership_stability_summary_M01.csv"))
saveRDS(
  stability$participant,
  file.path(mi_qc_dir, "subgroup_membership_stability_participant_M01.rds"),
  compress = FALSE
)

parameter_manifest <- data.table(
  parameter = c("analysis_n", "M", "maxit", "mi_workers", "mi_seed", "lca_nrep",
                "lca_maxiter", "lca_seed_base", "mi_progress_enabled",
                "mi_progress_poll_sec"),
  value = as.character(c(reference_n, M, maxit, mi_workers, mi_seed, lca_nrep,
                         lca_maxiter, lca_seed_base, mi_progress_enabled,
                         mi_progress_poll))
)
fwrite(parameter_manifest, file.path(mi_qc_dir, "MI_parameter_manifest_M01.csv"))
parameter_manifest <- rbindlist(list(parameter_manifest,data.table(parameter="validation_label",value=Sys.getenv("M01_VALIDATION_LABEL",""))))
fwrite(parameter_manifest, file.path(mi_qc_dir, "MI_parameter_manifest_M01.csv"))

provenance <- data.table(
  component = c("input_analysis_base", "variable_dictionary", "main_06", "function_06", "RM_definitions"),
  path = c(
    input_rds,
    dictionary_path,
    file.path(code_dir, "06_Multiple_imputation_M01.R"),
    model_snapshot_path,
    file.path(code_dir, "config", "RM_domains_v1.7.csv")
  )
)
provenance[, exists := file.exists(path)]
provenance[, md5 := vapply(
  path,
  function(p) if (file.exists(p)) unname(as.character(tools::md5sum(p))) else NA_character_,
  character(1)
)]
fwrite(provenance, file.path(mi_qc_dir, "MI_code_input_provenance_M01.csv"))

capture.output(sessionInfo(), file = file.path(mi_qc_dir, "sessionInfo_M01.txt"))

mi06_step_message(
  "Branch B 完成：M=", M,
  "；每套 N=", format(reference_n, big.mark = ","),
  "；输出目录：", mi_dir
)

# 全部RM报告在01/02/05/06的C目录集中生成。
# 需要保留05 RDS数据在根目录，以及原有MI_M01/MI_completed路径不变。

if (identical(Sys.getenv("M01_EXPORT_QA_XLSX", "0"),"1"))
  mi06_export_workbooks(data_clean_dir,code_dir)

m01_runtime_finish()
