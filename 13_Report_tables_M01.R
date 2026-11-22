# =============================================================================
# 13_Report_tables_M01.R | 投稿补充表：按表格独立模块，所有路径/参数在主脚本
# 一份主脚本只 source 一份 Function_13_Report_tables_M01.R
# 数据只从 05–12 已保存的 CSV 读取，绝不重新拟合、插补或合并方差。
# =============================================================================
rm(list = ls())
options(stringsAsFactors = FALSE)


# 0. 公共路径、版式、运行参数（逐表专属参数见下方每个模块）
project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code") # 项目根路径
code_dir <- Sys.getenv("M01_CODE_DIR", file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")) # 唯一函数脚本与模板所在目录
analysis_dir <- Sys.getenv("M01_ANALYSIS_OUTPUT", file.path(project_root, "LJW_outcome", "data_intervention", "M01_ASCVD", "analysis_results")) # 07–12 pooled 结果根目录（只读）
clean_dir <- Sys.getenv("M01_CLEAN_DIR", file.path(project_root, "LJW_outcome", "data_clean", "M01_ASCVD")) # 05/06 已保存 MI QC 根目录（只读）
report_root <- Sys.getenv("M01_REPORT_OUTPUT", file.path(analysis_dir, "publication_outputs_M01", "v6.2")) # 13/14 报告输出根目录，绝不写回上游分析输出
study_config_path <- file.path(analysis_dir, "study_config.rds") # 上游冻结研究配置（只读）；仅核验时间点/情景

execution <- list(
  validate_only = FALSE, # TRUE 只检查逐表文件，不生成结果
  stop_on_missing = TRUE, # TRUE 某表缺少必需结果时停止；FALSE 记录并跳过该表
  export_word = TRUE, # TRUE 在所有表模块结束后填入 Word 原模板
  export_excel = TRUE, # TRUE 同时导出多 sheet Excel（数据检查用途）
  export_full_precision = FALSE, # TRUE 每张表保存读取的未舍入 pooled 原文件副本；源文件本身绝不改动
  strict_row_matching = TRUE, # TRUE Word 模板中任何统计占位无法精确匹配即报错
  overwrite = TRUE, # TRUE 允许覆盖本次 13 输出文件，不会覆盖 05–12 源文件
  verify_stage10_manifest = TRUE, # TRUE 校验10步统一导出批次和CSV哈希
  export_main_table1 = FALSE # 可选：只在配备主文Table 1模板时启用
)

report_style <- list(
  horizons = c(10, 15), # 允许展示的随访年数；必须已由 07–12 保存结果
  primary_horizon = 15, # 四张主图及部分补充表采用的主要随访年数
  scenario_order = c("S0", "SH", "ST", "SHT"), # 上游内部代码，不更改分析
  display_scenarios = c("S0", "SN", "ST", "SNT"), # 论文显示名；SN↔SH、SNT↔SHT
  scope_order = c(
    "overall", "age_lt_60", "age_ge_60", "female", "male", "high_ses", "low_ses",
    "no_hypertension", "hypertension"
  ), # 总体/年龄/性别/SES/高血压展示顺序
  transition_order = c(
    "unchanged_untreated", "newly_intervention", "intervention_removed",
    "treatment_intensified", "unchanged_intensity", "treatment_deintensified"
  ), # 六种治疗转换固定顺序
  digits = list(
    count = 0L, # 人数显示的小数位（必须为 0）
    proportion = 2L, # 人群比例、事件捕获比例等的小数位
    risk = 2L, # ASCVD 累积发病风险（%）及其置信区间的小数位
    rd = 2L, # 风险差/ARR（百分点）及其置信区间的小数位
    rr = 2L, # 风险比 RR 的小数位（不乘 100）
    ldl = 2L, # 人均 LDL-C 降幅（mmol/L）的小数位
    total_ldl = 2L, # 人群 LDL-C 降幅总和（mmol/L）的小数位
    p = 2L # 双侧 P 值显示的小数位
  ),
  labels = list(
    scope = c(
      overall = "Overall", age_lt_60 = "<60 years", age_ge_60 = "≥60 years",
      female = "Female", male = "Male", high_ses = "High SES", low_ses = "Low SES",
      no_hypertension = "No hypertension", hypertension = "Hypertension"
    ), # 亚组英文显示名
    transition = c(
      unchanged_untreated = "Unchanged untreated", newly_intervention = "New intervention required",
      intervention_removed = "Intervention removed", treatment_intensified = "Treatment intensified",
      unchanged_intensity = "Unchanged intensity", treatment_deintensified = "Treatment de-intensified"
    ) # 转换组英文显示名
  ),
  colors = list(
    guideline_2019 = "#2166AC", # 2019 SCORE 指南颜色：蓝色
    guideline_2025 = "#E66132", # 2025 SCORE2 指南颜色：橘色
    difference = "#773DB1", # 指南效果差异点/置信区间的紫色
    rm = list(S0 = "#888888", SH = "#2166AC", ST = "#E66132", SHT = "#31853F") # RM 情景固定配色
  ),
  table = list(
    font_family = "Times New Roman", # Word 表格英文字体
    font_size = 9 # 采用v6.2 Word原生9pt字体，R不重新设置版式
  ),
  figure = list(
    width = 12, # 单张 Figure 的输出宽度（英寸）
    height = 9, # 单张 Figure 的输出高度（英寸）；12:9=4:3
    dpi = as.integer(Sys.getenv("M01_REPORT_DPI", "600")), # TIFF 图片分辨率（dpi），PDF 保持矢量
    formats = strsplit(Sys.getenv("M01_REPORT_FORMATS", "pdf,tiff"), ",", fixed = TRUE)[[1]], # 图像导出格式；可选单独 'pdf' 或 'tiff'
    font_family = "Arial", # 图形英文字体；未安装时需改为本机可用字体
    base_size = 11, # 坐标轴标签、图例及文字基准字号
    title_size = 12, # A/B/C/D Panel 标题字号
    point_size = 2.8, # 点估计图形大小
    ci_width = 0.65, # 95%CI 误差线粗细
    heatmap_text_size = 3.6, # Figure 1 重分类矩阵内文字大小
    x_angle = 25 # Figure 3/4 横坐标标签旋转角度（度）
  )
)

build_legacy_report <- function(report_root, cache_only = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("第 13 步需要 R 包 data.table 来生成统计缓存。请在当前 R 环境安装 data.table 后重跑。", call. = FALSE)
  suppressPackageStartupMessages(library(data.table))
# v6.2: 27张补充表，新增E-values；沿用原Word母版，绝不从示例数值补数据。
source(Sys.getenv("M01_REPORT_FUNCTION",
  file.path(code_dir,"Function_13_Report_tables_M01.R")),encoding="UTF-8")
options(M01_REPORT_OVERWRITE = execution$overwrite)
cfg <- rp_load_config(study_config_path, report_style)
template_dir <- Sys.getenv("M01_REPORT_LEGACY_TEMPLATE_DIR",file.path(code_dir,"config","report_v62"))
word_template_path <- Sys.getenv("M01_LEGACY_SUPPLEMENTARY_TEMPLATE", file.path(template_dir, "Supplementary_Tables_v6.2_20261110.docx"))
main_template_path <- Sys.getenv("M01_MAIN_TABLE1_TEMPLATE",
  file.path(code_dir, "config", "report_v52", "main_table1_template.docx"))
blueprint_path <- Sys.getenv("M01_SUPPLEMENTARY_BLUEPRINT",
  file.path(template_dir,"M01_supplementary_v6.2_blueprint.csv"))
P <- function(board, name, layer = "pooled") rp_file(analysis_dir, board, name, layer)
paths <- list(
  baseline = P("baseline_and_imputation", "baseline_and_imputation_descriptive_summary", "per_imputation"),
  methods = P("baseline_and_imputation", "actual_imputation_methods"),
  main = P("baseline_and_imputation", "main_table1_by_observed_ASCVD"),
  lca = file.path(clean_dir, "MI_M01", "quality_control", "SES_LCA_component_profile_M01.csv"),
  transitions = P("treatment_transitions", "transition_counts_and_required_reductions"),
  requirements = P("intervention_requirements", "classification_and_requirements"),
  s0 = P("strategy_risks_and_effects", "strategy_risks_and_risk_differences"),
  s0_het = P("strategy_risks_and_effects", "guideline_risk_difference_heterogeneity"),
  observed = P("observed_risk_and_event_capture", "observed_population_risks"),
  transition_effects = P("treatment_transitions", "transition_risks_and_effects"),
  transition_observed = P("treatment_transitions", "observed_transition_risks"),
  contribution = P("treatment_transitions", "risk_difference_contributions"),
  domains = P("risk_classification", "domains"),
  ht_counts = P("risk_classification", "origin_category_partition_counts"),
  ht_risk = P("observed_risk_and_event_capture", "HT_stratum_risks_and_normalized_event_capture"),
  upgraded = P("observed_risk_and_event_capture", "high_risk_and_upgraded_risks_capture"),
  new_counts = P("intervention_requirements", "rm_newly_eligible_counts"),
  new_risk = P("observed_risk_and_event_capture", "rm_newly_eligible_observed_risks"),
  mean_reduction = P("treatment_transitions", "transition_mean_required_reduction_contrasts"),
  all_effects = P("strategy_risks_and_effects", "all_scenario_risks_and_effects"),
  rm_effects = P("risk_modifier_incremental_effects", "incremental_benefits_and_guideline_difference"),
  rm_het = P("risk_modifier_incremental_effects", "incremental_benefit_heterogeneity"),
  evalues = P("unmeasured_confounding_sensitivity", "risk_ratio_evalues")
)
s01_dir <- Sys.getenv("M01_S01_OUTPUT", file.path(project_root, "LJW_outcome", "data_intervention", "S01_different_RM_strategies", "analysis_results"))
s01_publication <- file.path(s01_dir, "publication_outputs_S01")
paths <- c(paths, list(
  s01_dictionary = file.path(s01_publication, "Table1_scenario_dictionary.csv"),
  s01_counts = file.path(s01_publication, "S01_counts_all_directions.csv"),
  s01_transitions = file.path(s01_publication, "S01_transitions_all_directions.csv"),
  s01_effects = file.path(s01_publication, "S01_effects_all_directions.csv"),
  s01_domains = file.path(s01_dir, "pooled", "DOMAIN", "S01_DOMAIN0_positive_domains_in_moderate_risk.csv")
))

# 1. 源文件及模板校验：文件不存在即失败，不使用Word展示数值兜底。
quality <- rp_dir(file.path(report_root, "quality_control"))
tables  <- rp_dir(file.path(report_root, "tables"))
options(M01_REPORT_OVERWRITE = execution$overwrite)
rp_require(word_template_path,'v6.2 Word模板')
blueprint <- rp62_blueprint_read(blueprint_path)
# 在读取统计结果前验证模板静态文字、行列数和合并单元格。
if (!cache_only && !requireNamespace('xml2', quietly=TRUE)) stop('需要R包xml2')
if (!cache_only) { 
.tmp <- tempfile('docx_preflight_')
dir.create(.tmp, recursive = TRUE)
tryCatch({
  utils::unzip(word_template_path, files = 'word/document.xml', exdir = .tmp)
  .doc_path <- file.path(.tmp, 'word', 'document.xml')
  if (!file.exists(.doc_path)) stop('Word模板缺少 word/document.xml')
  .doc <- xml2::read_xml(.doc_path)
  rp62_template_check(.doc, blueprint)
}, finally = {
  unlink(.tmp, recursive = TRUE)
})

}
preflight <- rp_preflight(paths)
rp_save_csv(preflight,file.path(quality,"source_manifest.csv"))
if(any(!preflight$exists))stop('MISSING_SOURCE: ',
  paste(sprintf('%s=%s',preflight$source[!preflight$exists],
    preflight$path[!preflight$exists]),collapse='; '))

cfg <- rp_load_config(study_config_path,report_style)
if(execution$verify_stage10_manifest)
  rp13_verify_stage10_manifest(analysis_dir,
    P('analysis_quality_control','classification_unified_export_manifest'))
src <- lapply(paths,rp_read)
method_audit <- rp52_method_audit(src$methods)
rp_save_csv(method_audit,file.path(quality,'imputation_method_template_audit.csv'))
if(any(!method_audit$matches_template))stop(
  'MISSING_SOURCE: Table 4插补方法与v6.2 Word冲突。查看imputation_method_template_audit.csv')
for(nm in names(src)) {
  d <- src[[nm]]
  if('m' %in% names(d) && any(d$m!=cfg$m,na.rm=TRUE))stop('MI count mismatch: ',nm)
  rp_check_ci(d,nm)
}

# 2. 只读汇总结果的加总核验。main含overall/no-observed/observed三个互斥列。
m <- src$main[,.(n=unique(n_total)),by=scope]
if(m[scope=='overall',n]!=m[scope!='overall',sum(n)])
  stop('Main Table 1的互斥列人数不满足总体之和')
ct <- src$transitions[,.(n=sum(n),den=unique(n_scope)[1],
  p=sum(proportion),k=.N),by=.(scope,scenario)]
if(any(abs(ct$n-ct$den)>1e-6 | abs(ct$p-1)>1e-8 | ct$k!=6))
  stop('六类治疗转变未能完整划分人群')
z <- merge(src$contribution[scenario=='S0',.(total=sum(estimate)),by=.(scope,time)],
  src$s0[scenario=='S0' & estimand=='rd_2025_vs_2019',
    .(scope,time,estimate)],by=c('scope','time'),all=TRUE)
if(anyNA(z)||any(abs(z$total-z$estimate)>1e-6))
  stop('治疗转换贡献之和与总体2025-2019 RD不符')
if(nrow(src$rm_het[time==report_style$primary_horizon])!=16L)
  stop('必须包含16项RM亚组异质性检验')

# 3. 动态格填值：S04、S09–S21来自M01；S22–S26来自S01；S27来自09步E-values。
base_cells <- rp61_build(src,blueprint,report_style)
s01_cells <- rp61_sensitivity_cells(src,blueprint,report_style)
eval_cells <- rp62_evalues_cells(src$evalues,report_style)
cells <- data.table::rbindlist(list(base_cells,s01_cells,eval_cells),fill=TRUE)
# 静态定义S01–S03、S05–S08不做内容覆写。
cells <- rp62_validate_cell_coverage(cells,blueprint)
provenance <- rp62_provenance(cells,paths,report_style)
rp_save_csv(provenance,file.path(quality,'supplementary_cell_provenance.csv'))
rp_save_csv(data.table(table=seq_len(27L),
  source_cells=vapply(seq_len(27L),function(i)sum(cells$table==i),integer(1)),
  demo_numeric=vapply(seq_len(27L),function(i)
    sum(blueprint$table==i & blueprint$demo_numeric),integer(1))),
  file.path(quality,'template_cell_coverage.csv'))
for(k in seq_len(27L))
  rp_save_csv(provenance[table==k],
    file.path(tables,sprintf('Supplementary_Table_%02d_cells.csv',k)))

if(!cache_only && !execution$validate_only) {
  if(execution$export_word)rp62_word_write(word_template_path,
    file.path(tables,'M01_Supplementary_Tables_v6.2.docx'),
    provenance,blueprint,execution$overwrite)
  if(execution$export_excel) {
    if(!requireNamespace('openxlsx',quietly=TRUE))stop('Excel审计输出需要R包openxlsx')
    wb<-openxlsx::createWorkbook()
    for(k in seq_len(27L)) {
      sheet<-sprintf('S%02d',k)
      openxlsx::addWorksheet(wb,sheet)
      openxlsx::writeData(wb,sheet,provenance[table==k])
    }
    openxlsx::saveWorkbook(wb,
      file.path(tables,'M01_Supplementary_Tables_v6.2_audit.xlsx'),
      overwrite=execution$overwrite)
  }
  if(execution$export_main_table1) {
    rp_require(main_template_path,'主文Table1模板')
    # 主文Table1的独立输出需要其对应的Table1填值器；不与27张补充表混写。
    rp_save_csv(rp51_main_cells(src$main),
      file.path(tables,'Main_Table_1_cell_provenance.csv'))
    rp_msg('Main Table 1 provenance已输出；主文模板格式不受本次27张补充表更新影响。')
  }
  if(execution$export_full_precision) {
    folder<-rp_dir(file.path(tables,'full_precision_sources'))
    for(nm in names(paths))file.copy(paths[[nm]],
      file.path(folder,paste0(nm,'.csv')),overwrite=execution$overwrite)
  }
}
rp_save_csv(data.table(check=c('27_tables_template_schema','all_numeric_cells_filled',
  'six_transition_partition','contribution_additivity','imputation_method_match',
  'E_values_3x4'),passed=TRUE),file.path(quality,'table_quality_checks.csv'))
rp_msg('M01 v6.2 共27张补充表；validate_only=',execution$validate_only,
  '；输出目录：',tables)

}

report_version <- Sys.getenv("M01_REPORT_VERSION", "6.3")
if (identical(report_version, "6.3")) {
  source(Sys.getenv("M01_REPORT_FUNCTION", file.path(code_dir, "Function_13_Report_tables_M01.R")),
         encoding = "UTF-8")
  cache <- Sys.getenv("M01_REPORT_LEGACY_OUTPUT", file.path(
    analysis_dir, "publication_outputs_M01", "v6.2"))
  manifest <- file.path(cache, "quality_control", "source_manifest.csv")
  required_tables <- file.path(cache, "tables", sprintf("Supplementary_Table_%02d_cells.csv", 1:27))
  fresh <- FALSE
  if (all(file.exists(c(manifest, required_tables)))) {
    sources <- read.csv(manifest, fileEncoding = "UTF-8-BOM")
    fresh <- nrow(sources) > 0L &&
      isTRUE(all(unname(tools::md5sum(sources$path)) == sources$md5))
  }
  if (!fresh) {
    rp_msg("Rebuilding supplementary statistics from saved analysis outputs")
    build_legacy_report(cache, cache_only = TRUE)
  }
  rp63_generate(project_root, code_dir, analysis_dir, execution$overwrite)
} else if (identical(report_version, "6.2")) {
  build_legacy_report(report_root)
} else stop("Unknown M01_REPORT_VERSION: ", report_version)
