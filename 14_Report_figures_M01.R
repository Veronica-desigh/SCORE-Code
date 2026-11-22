# =============================================================================
# 14_Report_figures_M01.R | 图形只读上游结果；末尾独立AJ交叉验证QC；版式修订 v4
# 一个主脚本仅对应一个 Function_14 脚本；直接读取 pooled 全精度结果。
# 本文件分为：路径/来源、共用展示规则、Figure 1/2/3/4 参数、执行与 QC。
# =============================================================================
rm(list = ls())
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
# =============================================================================
# 0. 主脚本参数和全部读取路径（路径及参数均在此，不藏在函数文件）
# 本版只导出PDF与SVG矢量图，不导出其他图像格式；数据CSV仍可按开关导出。
# =============================================================================
project_root <- Sys.getenv("M01_PROJECT_ROOT", "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code") # M01 项目根目录
code_dir <- Sys.getenv("M01_CODE_DIR", file.path(project_root, "LJW_code", "Main_code", "M01_ASCVD")) # 14 主脚本及函数文件存放目录
analysis_dir <- Sys.getenv("M01_ANALYSIS_OUTPUT", file.path(project_root, "LJW_outcome", "data_intervention", "M01_ASCVD", "analysis_results")) # 07–12 正式分析结果根目录（只读）
report_root <- file.path(analysis_dir, "publication_outputs_M01", "v6.3")
study_config_path <- file.path(analysis_dir, "study_config.rds") # 07 科学参数配置（只读）

source_paths <- list(
  reclassification = file.path(analysis_dir, "risk_classification", "pooled", "reclassification.csv"), # Figure 1 A 的 SCORE→SCORE2 四乘四矩阵
  reclassification_direction = file.path(analysis_dir, "risk_classification", "pooled", "reclassification_direction.csv"), # Figure 1 A 辅助说明的升降级比例
  hvh_observed = file.path(analysis_dir, "observed_risk_and_event_capture", "pooled", "high_risk_and_upgraded_risks_capture.csv"), # Figure 1 B/C 高危风险及捕获率
  fig1_eligibility = file.path(analysis_dir, "intervention_requirements", "pooled", "classification_and_requirements.csv"), # Figure 1 D 真实 n_intervention 和 MI_mean_proportion
  transition_counts = file.path(analysis_dir, "treatment_transitions", "pooled", "transition_counts_and_required_reductions.csv"), # Figure 2 A/B 六组人数/人均及组内 LDL-C 所需降幅
  transition_observed = file.path(analysis_dir, "treatment_transitions", "pooled", "observed_transition_risks.csv"), # Figure 2 C 六组观察风险
  transition_effects = file.path(analysis_dir, "treatment_transitions", "pooled", "transition_risks_and_effects.csv"), # Figure 2 D 六组风险差
  s0_observed = file.path(analysis_dir, "observed_risk_and_event_capture", "pooled", "observed_population_risks.csv"), # 兼容保留的历史观察风险路径；新版Figure 3不再读取
  s0_effects = file.path(analysis_dir, "strategy_risks_and_effects", "pooled", "strategy_risks_and_risk_differences.csv"), # 兼容保留的历史风险差路径；新版Figure 3不再读取
  rm_incremental = file.path(analysis_dir, "risk_modifier_incremental_effects", "pooled", "incremental_benefits_and_guideline_difference.csv"), # Figure 4 A/B RM 相对 S0 的额外获益
  all_effects = file.path(analysis_dir, "strategy_risks_and_effects", "pooled", "all_scenario_risks_and_effects.csv") # Figure 4 C RM 策略下 2025−2019 残余风险差
)

source_paths$population_mean <- file.path(analysis_dir, "intervention_requirements", "pooled", "population_mean_required_reduction_contrasts.csv")
source_paths$transition_contribution <- file.path(analysis_dir, "treatment_transitions", "pooled", "risk_difference_contributions.csv") # 12 已合并的 absolute_RD_contribution；估计值为风险单位，在14仅乘100显示为pp
source_paths$fig3_hvh <- file.path(analysis_dir, "risk_classification", "pooled", "classification_and_requirements.csv") # Figure 3A：Module 10 最终H/VH人数与RM升级人数（仅用于展示，不重新估计）
source_paths$joint_counts <- file.path(analysis_dir, "risk_classification", "pooled", "origin_category_partition_counts.csv")
source_paths$joint_risks <- file.path(analysis_dir, "observed_risk_and_event_capture", "pooled", "joint_origin_category_risks_and_event_capture.csv")

execution <- list(
  figures = c(1L, 2L, 3L, 4L, 5L, 6L, 7L), # 本次只生成补充图SN/ST/SNT，不重新生成或覆盖正文四幅图
  validate_only = FALSE, # TRUE 只核验配置/结果文件/置信区间，不生成图片
  export_figure_data = TRUE, # TRUE 同时保存每幅图完整精度的绘图数据（CSV）
  strict_internal_style = TRUE, # TRUE 检查14自身的风险时间点、分组、情景、显示精度等；不读取13文件
  overwrite = TRUE, # TRUE 允许覆盖报告目录中的同名图；不写上游 CSV
  scenario_input_mode = "allow_legacy_rename", # 仅显示层：允许把旧 SH/SHT 编码映射到 SN/SNT；设为 "strict_new" 禁用
  reject_mixed_scenario_generations = TRUE # 发现一部分 CSV 为旧编码、另一部分为新编码时停止，防止混合不同分析版本
)

# 补充图的独立统计QC参数：只用Module07/10逐插补数据核对小样本AJ，绝不改动统计结果。
# 首次运行需要安装标准R软件包cmprsk；survival是现有分析流程已使用的R包。
supplementary_small_sample_qc <- list(
  enabled = TRUE, # TRUE：在绘图结束后独立进行小样本竞争风险交叉验证
  max_cell_n = 50L, # 每个插补数据集、每个联合分类单元格人数<=50，纳入审计
  horizon = 15, # 验证15年累计ASCVD风险；必须存在于study_config$times中
  point_tolerance = 1e-6, # 与两种标准软件AJ点估计的最大允许绝对差（概率单位）
  se_relative_warning = 0.15 # 原影响函数SE与cmprsk标准SE相对差异>15%时仅作QC预警
)
# =============================================================================
# 0A. PUBLICATION DESIGN SYSTEM：全部可见字体、颜色、命名的唯一设计来源
# 改动这里会同步影响四张主图；特殊需求在 figure_layout$fig1/fig2/fig3 覆盖。
# 颜色语义：2019=蓝、2025=橙；RM=S0灰/SN蓝/ST橙/SNT绿；
# 图2转换组使用低饱和分类色；观察风险/风险热图使用同一蓝系；capture 用紫系。
# =============================================================================
publication_design <- list(
  typography = list( # 全套图形共享的字号与字体设置
    family = "Arial", # 出版图形统一字体
    base_pt = 13.0, # 正文图基础字号（pt）
    panel_title_pt = 14.0, # 正文图面板标题基准字号（pt）
    axis_title_pt = 13.0, # 正文图轴标题基准字号（pt）
    axis_tick_pt = 12.7, # 正文图坐标刻度基准字号（pt）
    legend_pt = 13.0, # 正文图图例基准字号（pt）
    fig2_export_scale = 1.12, # Figure 2字号相对全局的放大倍数
    fig2_bar_value_mm = 4.0, # Figure 2 A柱末百分比文字字号（mm）
    fig3_favors_mm = 3.8, # Figure 4倾向说明文字字号（mm）
    supp_panel_title_pt = 16.0, # 补充图面板标题字号（pt）
    supp_base_pt = 11.5, # 补充图基础字号（pt）
    supp_cell_n_mm = 3.05, # 补充图单元格人数文字字号（mm）
    supp_cell_ci_mm = 2.85, # 补充图单元格置信区间文字字号（mm）
    supp_header_mm = 3.15, # 补充图表头字号（mm）
    supp_region_mm = 3.0, # 补充图区域标签字号（mm）
    supp_footer_mm = 2.60, # 补充图底部说明字号（mm）
    supp_legend_mm = 2.8 # 补充图色标字号（mm）
  ),
  colors = list( # 统一绘图配色表
    guideline = c(`2019` = "#2166AC", `2025` = "#E66132"), # 2019蓝色与2025橙色的指南配色
    rm = c(S0 = "#858C94", SN = "#2166AC", ST = "#E66132", SNT = "#35865A"), # S0/SN/ST/SNT四种RM策略颜色
    observed = "#2166AC", # 观察风险图的主蓝色
    transitions = c( # 六种治疗转换组的柱形颜色
      unchanged_untreated = "#2166AC", # 干预不变且无需治疗组标签
      newly_intervention = "#2166AC", # 新增干预组标签
      intervention_removed = "#2166AC", # 撤销干预组标签
      treatment_intensified = "#2166AC", # 治疗加强组标签
      unchanged_intensity = "#2166AC", # 治疗强度不变组标签
      treatment_deintensified = "#2166AC" # 治疗减弱组标签
    ),
    risk_gradient = c("#F0F6FC", "#C8DDF1", "#78A6D0", "#2166AC"), # 风险矩阵连续蓝色色阶
    capture_gradient = c("#FCF8FD", "#EDDCF6", "#BC8FDB", "#773DB1"), # 事件捕获连续紫色色阶
    text = "#1F2933", # 正文文字颜色
    axis = "#303840", # 坐标轴线与刻度颜色
    subtle_text = "#58636E", # 次级说明文字颜色
    panel_guide = "#E8EDF2", # 分组水平辅助线颜色
    divider = "#D5DCE3", # 分区分隔线颜色
    zero_line = "#9DA7B1", # 风险差参考零线颜色
    header_fill = "#EAF0F5", # 矩阵表头背景色
    total_fill = "#F5F7F9", # 矩阵合计格背景色
    empty_fill = "#EEEEEE", # 无数据单元格背景色
    border = "#BFC9D3", # 常规单元格描边颜色
    strong_border = "#394955", # 重点分区描边颜色
    white = "#FFFFFF", # 纯白背景色
    bar_edge = "#343E49", # 柱形黑灰边框色
    connector = "#BBC5CE" # 杠铃图两点连接线颜色
  ),
  labels = list( # 分类、坐标或图例的英文显示名称集合
    short_subgroup = c( # 亚组简写标签
      overall = "Overall", # 总体人群标签
      age_lt_60 = "<60 years", # 年龄小于60岁标签
      age_ge_60 = "≥60 years", # 年龄不小于60岁标签
      female = "Female", # 女性标签
      male = "Male", # 男性标签
      high_ses = "High SES", # 高社会经济地位标签
      low_ses = "Low SES", # 低社会经济地位标签
      no_hypertension = "No HTN", # 非高血压标签
      hypertension = "HTN" # 高血压标签
    )
  )
)
report_style <- list(
  horizons = c(10, 15), # 可选的研究随访时点（年）
  primary_horizon = 15, # 主图使用的随访年限（年）
  ci_interval_separator = " to ", # 全部正文图及补充图95%CI上下限之间的英文连接词（保留两侧空格，避免与负号混淆）
  scenario_order = c("S0", "SN", "ST", "SNT"), # RM策略展示顺序
  scope_order = c( # 总体与亚组的展示顺序
    "overall", 
    "age_lt_60", 
    "age_ge_60", 
    "female", 
    "male", 
    "high_ses", 
    "low_ses",
    "no_hypertension", 
    "hypertension"
  ), # 总体/年龄/性别/SES/高血压展示顺序
  transition_order = c( # 六种治疗转换组的展示顺序
    "unchanged_untreated", 
    "newly_intervention", 
    "intervention_removed",
    "treatment_intensified", 
    "unchanged_intensity", 
    "treatment_deintensified"
  ), # 六种治疗转换固定顺序
  digits = list( # 所有数值标签的小数位设置
    count = 0L, # 人数的小数位数
    proportion = 2L, # 比例和事件捕获率的小数位数
    risk = 2L, # 累计风险及区间的小数位数
    rd = 2L, # 风险差及区间的小数位数
    rr = 2L, # 风险比的小数位数
    ldl = 2L, # LDL-C平均降幅的小数位数
    total_ldl = 2L, # LDL-C总降幅的小数位数
    p = 2L # P值的小数位数
  ),
  labels = list( # 分类、坐标或图例的英文显示名称集合
    scope = c( # 亚组英文名称对应表
      overall = "Overall", # 总体人群标签
      age_lt_60 = "<60 years", # 年龄小于60岁标签
      age_ge_60 = "≥60 years", # 年龄不小于60岁标签
      female = "Female", # 女性标签
      male = "Male", # 男性标签
      high_ses = "High SES", # 高社会经济地位标签
      low_ses = "Low SES", # 低社会经济地位标签
      no_hypertension = "No hypertension", # 非高血压标签
      hypertension = "Hypertension" # 高血压标签
    ), # 亚组英文显示名
    transition = c( # 治疗转换组英文名称对应表
      unchanged_untreated = "Unchanged untreated", # 干预不变且无需治疗组标签
      newly_intervention = "New intervention required", # 新增干预组标签
      intervention_removed = "Intervention removed", # 撤销干预组标签
      treatment_intensified = "Treatment intensified", # 治疗加强组标签
      unchanged_intensity = "Unchanged intensity", # 治疗强度不变组标签
      treatment_deintensified = "Treatment de-intensified" # 治疗减弱组标签
    ) # 转换组英文显示名
  ),
  colors = list( # 统一绘图配色表
    guideline_2019 = unname(publication_design$colors$guideline[["2019"]]), # 2019指南比较系列颜色
    guideline_2025 = unname(publication_design$colors$guideline[["2025"]]), # 2025指南比较系列颜色
    difference = tail(publication_design$colors$capture_gradient, 1L), # 指南差异系列颜色
    rm = as.list(publication_design$colors$rm) # S0/SN/ST/SNT四种RM策略颜色
  ),
  table = list( # 结果表格字体配置
    font_family = "Arial", # 图形字体族
    font_size = 8.5 # 结果表格英文字号（pt）
  ),
  figure = list( # 所有图像的默认导出参数
    width = 12, # 图像导出宽度（英寸）
    height = 9, # 图像导出高度（英寸）
    dpi = as.integer(Sys.getenv("M01_REPORT_DPI", "600")), # 位图输出分辨率（dpi）
    formats = c("pdf", "svg"), # 仅输出PDF和SVG两种矢量图；不生成TIFF、PNG等位图（不受旧环境变量覆盖）
    font_family = publication_design$typography$family, # 图形字体族
    base_size = publication_design$typography$base_pt, # ggplot主题基础字号（pt）
    background_color = publication_design$colors$white # 导出图像的画布背景色
  )
)
figure_layout <- list(
  design = publication_design, # 统一指向出版风格设计系统
  common = list( # 各张图共享的基础版式参数
    panel_title_size = publication_design$typography$panel_title_pt, # 面板主标题字号（pt）
    panel_prefix_separator = "  ", # Panel字母与标题之间的空格
    panel_title_gap_pt = 7, # Panel标题至绘图区间距（pt）
    axis_title_size = publication_design$typography$axis_title_pt, # 坐标轴标题字号（pt）
    axis_text_size = publication_design$typography$axis_tick_pt, # 刻度数字及类别名称字号（pt）
    group_dividers = c(1.5, 3.5, 5.5, 7.5), # 亚组大类之间的分隔线位置
    group_divider_color = publication_design$colors$divider, # 亚组分隔线颜色
    group_divider_width = 0.4, # 亚组分隔线线宽
    subgroup_short_labels = publication_design$labels$short_subgroup, # 亚组坐标轴缩写名称
    figure_margin_pt = c(8, 12, 7, 10) # 通用图形外边距：上/右/下/左（pt）
  ),
  fig1 = list( # Figure 1专用的所有参数
    export_width = 14.0, # Figure 1导出宽度（英寸）
    export_height = 11.0, # Figure 1导出高度（英寸）
    layout = list( # 面板A、B、C、D的拼图排版配置
      matrix_to_bottom_height = c(1.66, 1.24), # 矩阵A与下方B/C/D的相对高度
      summary_height_ratio = c(11.5, 1.75), # Panel A矩阵与其下方汇总带高度比例
      panel_gap_ratio = 0.035, # Panel A与下方三个面板之间的垂直留白占比
      bottom_panel_widths = c(1, 1, 1) # Panel B/C/D三列宽度比例
    ),
    guideline_colors = publication_design$colors$guideline, # 2019蓝色与2025橙色的柱状图颜色
    matrix_colors = publication_design$colors$risk_gradient, # Panel A风险比例填充色阶
    matrix_header_fill = publication_design$colors$header_fill, # Panel A行列表头及跨行列标题的背景色
    matrix_total_fill = publication_design$colors$total_fill, # Panel A合计行列背景色
    matrix_total_border = publication_design$colors$border, # Panel A合计格描边颜色
    summary_divider_color = publication_design$colors$divider, # Panel A汇总带分隔线颜色
    bar_border_color = publication_design$colors$bar_edge, # 柱状图轮廓线颜色
    panel_A = list( # Panel A：风险重分类矩阵与汇总带的配置
      letter = "A", # 该Panel左上角编号字母
      title = "Reclassification from 2019 SCORE to 2025 SCORE2", # 该Panel标题原文
      title_size_pt = 17.5, # 该Panel标题字号（pt）
      title_gap_pt = 5, # 该Panel标题到图形主体的间距（pt）
      title_lineheight = 1.02, # 该Panel标题多行文本的行高系数
      panel_margin_pt = c(2, 5, 0, 4), # 该Panel上/右/下/左留白（pt）
      matrix_side_header_label = "2019 SCORE-based\nrisk category", # 矩阵左侧2019风险分组总表头文字
      matrix_column_header = "2025 SCORE2-based risk category", # 矩阵顶部2025风险分组总表头文字
      matrix_superheader_x = c(1.504, 5.496), # 2025顶部合并表头的左右坐标范围
      matrix_superheader_y = c(5.505, 6.475), # 2025顶部合并表头的上下坐标范围
      matrix_sideheader_x = c(-.28, .504), # 2019左侧合并表头的左右坐标范围
      matrix_sideheader_y = c(.51, 4.49), # 2019左侧合并表头的上下坐标范围
      matrix_superheader_angle = 0, # 2025顶部合并表头文字旋转角度（度）
      matrix_sideheader_angle = 90, # 2019左侧合并表头文字旋转角度（度）
      matrix_merged_header_border_color = publication_design$colors$white, # 合并表头边框颜色
      matrix_merged_header_border_width = .60, # 合并表头边框线宽
      matrix_column_labels = c( # 2025年L/M/H/VH及Total列标题文本
        "Low\n(<2%)", "Moderate\n(2% to <10%)",
        "High\n(10% to <20%)", "Very high\n(≥20%)", "Total"
      ),
      matrix_row_labels = c( # 2019年L/M/H/VH及Total行标题文本
        "Low\n(<1%)", "Moderate\n(1% to <5%)",
        "High\n(5% to <10%)", "Very high\n(≥10%)", "Total"
      ),
      matrix_text_size_mm = 4.55, # 重分类格和总计格内文字字号（mm）
      matrix_header_size_mm = 4.45, # 矩阵L/M/H/VH行列标题字号（mm）
      matrix_top_header_size_mm = 5.00, # 2025顶部合并表头字号（mm）
      matrix_side_header_size_mm = 5.00, # 2019左侧合并表头字号（mm）
      matrix_cell_lineheight = .98, # 矩阵内部数值标签行距系数
      matrix_header_lineheight = .96, # 普通行列标题行距系数
      matrix_side_header_lineheight = .98, # 2019侧表头行距系数
      matrix_superheader_lineheight = .98, # 2025顶表头行距系数
      matrix_x_limits = c(-.31, 6.51), # 矩阵显示的横向坐标范围
      matrix_y_limits = c(-.52, 6.51), # 矩阵显示的纵向坐标范围
      matrix_tile_width = .992, # 矩阵单元格宽度（绘图坐标单位）
      matrix_tile_height = .98, # 矩阵单元格高度（绘图坐标单位）
      matrix_header_border_width = .60, # 普通表头单元格边框线宽
      matrix_cell_border_width = .85, # 风险矩阵单元格边框线宽
      matrix_total_border_width = .45, # 合计格描边线宽
      matrix_fill_max = 100, # 矩阵蓝色色阶上界（行比例，%）
      matrix_dark_text_threshold = 60, # 达到此比例时矩阵格文字改成白色（%）
      matrix_proportion_digits = 1L, # 矩阵比例数值显示小数位数
      show_row_total = TRUE, # 是否保留右侧合计列
      show_col_total = TRUE, # 是否保留下方合计行
      show_direction_counts = TRUE, # 是否显示降级、不变、升级三项人数和比例
      summary_color_title = "Row proportion (%)", # 矩阵汇总带色标标题
      summary_direction_labels = c( # 汇总带降级/不变/升级的标签集合
        lower = "Reclassified to lower risk", # 风险下降重分类的英文标签
        unchanged = "Unchanged", # 风险分类不变的英文标签
        higher = "Reclassified to higher risk" # 风险上升重分类的英文标签
      ),
      summary_text_size_mm = 4.45, # 汇总带类别标题字号（mm）
      summary_value_size_mm = 4.45, # 汇总带人数/比例字号（mm）
      summary_legend_tick_size_mm = 4.45, # 色标刻度数值字号（mm）
      summary_proportion_digits = 1L, # 汇总带比例显示小数位数
      summary_gradient_steps = 90L, # 色标渐变分段数量
      summary_divider = TRUE, # 是否显示汇总带竖向分隔线
      summary_divider_width = .35, # 汇总带分隔线宽度
      summary_margin_pt = c(2, 5, 1, 4), # 汇总带外边距：上/右/下/左（pt）
      summary_x_limits = c(0, 4), # 汇总带横向绘图坐标范围
      summary_y_limits = c(.15, .91), # 汇总带纵向绘图坐标范围
      summary_col_centers = c(1.5, 2.5, 3.5), # 降级、不变、升级三列中心横坐标
      summary_title_y = .74, # 汇总带类别标题纵坐标
      summary_value_y = .37, # 汇总带人数与比例标签纵坐标
      summary_legend_title_x = .50, # 汇总带色标标题横坐标
      summary_legend_title_y = .78, # 汇总带色标标题纵坐标
      summary_legend_y = .22, # 汇总带色标刻度数值纵坐标
      summary_gradient_x = c(.12, .88), # 色标渐变条左右边界
      summary_gradient_y = c(.34, .49), # 色标渐变条上下边界
      summary_divider_x = c(1, 2, 3), # 汇总带三条分隔线横坐标
      summary_divider_y = c(.18, .85) # 汇总带分隔线纵向起止坐标
    ),
    panel_B = list( # Panel B：高危观察风险柱状图配置
      letter = "B", # 该Panel左上角编号字母
      title = "Observed 15-year ASCVD risk\nin high/very-high risk", # 该Panel标题原文
      y_title = "Observed ASCVD risk (%)", # 该Panel纵轴标题文本
      x_title = NULL, # 该Panel横轴标题文本
      x_labels = c(`2019` = "2019 SCORE\n(≥5%)", `2025` = "2025 SCORE2\n(≥10%)"), # 横轴2019/2025两组类别标签
      title_size_pt = 16.2, # 该Panel标题字号（pt）
      title_gap_pt = 6, # 该Panel标题到图形主体的间距（pt）
      title_lineheight = 1.03, # 该Panel标题多行文本的行高系数
      axis_title_size_pt = 13.8, # 该Panel坐标轴标题字号（pt）
      axis_tick_size_pt = 13.0, # 该Panel坐标轴刻度字号（pt）
      axis_title_gap_pt = 5, # Panel纵轴标题与刻度文字间距（pt）
      x_tick_gap_pt = 3, # Panel横轴刻度数字与轴线间距（pt）
      panel_margin_pt = c(5, 5, 5, 5), # 该Panel上/右/下/左留白（pt）
      bar_width = .54, # 柱宽（图形坐标单位）
      bar_border_width = .55, # 柱轮廓线宽
      show_ci = TRUE, # 是否绘制置信区间误差线
      errorbar_width = .022, # 误差线端帽宽度
      errorbar_linewidth = .70, # 误差线线宽
      label_mode = "ci", # 柱顶标签模式（ci=百分比与95%CI；n=人数与比例）
      label_digits = 2L, # 柱顶百分比显示小数位数
      label_size_mm = 4.50, # 柱顶数值标签字号（mm）
      label_lineheight = 1.08, # 柱顶双行数值标签行距系数
      label_offset_fraction = .047, # 柱顶标签相对柱高的向上偏移比例
      top_padding_fraction = .155, # 柱顶标签上方额外留白比例
      y_from_zero = TRUE, # 是否强制纵坐标从零开始
      y_tick_count = 5, # 纵轴自动刻度数量建议值
      y_axis_limits = c(0, NA_real_), # 纵轴显示范围；NA上限表示自动适配数据
      y_axis_breaks = NULL, # 纵轴刻度序列；NULL使用自动刻度
      x_expand = c(.36, .36) # 两根柱左右额外留白宽度
    ),
    panel_C = list( # Panel C：事件捕获柱状图配置
      letter = "C", # 该Panel左上角编号字母
      title = "15-year ASCVD events captured\nin high/very-high risk", # 该Panel标题原文
      y_title = "ASCVD events captured (%)", # 该Panel纵轴标题文本
      x_title = NULL, # 该Panel横轴标题文本
      x_labels = c(`2019` = "2019 SCORE\n(≥5%)", `2025` = "2025 SCORE2\n(≥10%)"), # 横轴2019/2025两组类别标签
      title_size_pt = 16.2, # 该Panel标题字号（pt）
      title_gap_pt = 6, # 该Panel标题到图形主体的间距（pt）
      title_lineheight = 1.03, # 该Panel标题多行文本的行高系数
      axis_title_size_pt = 13.8, # 该Panel坐标轴标题字号（pt）
      axis_tick_size_pt = 13.0, # 该Panel坐标轴刻度字号（pt）
      axis_title_gap_pt = 5, # Panel纵轴标题与刻度文字间距（pt）
      x_tick_gap_pt = 3, # Panel横轴刻度数字与轴线间距（pt）
      panel_margin_pt = c(5, 5, 5, 5), # 该Panel上/右/下/左留白（pt）
      bar_width = .54, # 柱宽（图形坐标单位）
      bar_border_width = .55, # 柱轮廓线宽
      show_ci = TRUE, # 是否绘制置信区间误差线
      errorbar_width = .022, # 误差线端帽宽度
      errorbar_linewidth = .70, # 误差线线宽
      label_mode = "ci", # 柱顶标签模式（ci=百分比与95%CI；n=人数与比例）
      label_digits = 2L, # 柱顶百分比显示小数位数
      label_size_mm = 4.50, # 柱顶数值标签字号（mm）
      label_lineheight = 1.08, # 柱顶双行数值标签行距系数
      label_offset_fraction = .047, # 柱顶标签相对柱高的向上偏移比例
      top_padding_fraction = .155, # 柱顶标签上方额外留白比例
      y_from_zero = TRUE, # 是否强制纵坐标从零开始
      y_tick_count = 5, # 纵轴自动刻度数量建议值
      y_axis_limits = c(0, NA_real_), # 纵轴显示范围；NA上限表示自动适配数据
      y_axis_breaks = NULL, # 纵轴刻度序列；NULL使用自动刻度
      x_expand = c(.36, .36) # 两根柱左右额外留白宽度
    ),
    panel_D = list( # Panel D：干预资格柱状图配置
      letter = "D", # 该Panel左上角编号字母
      title = "Eligibility for LDL-C lowering\nintervention", # 该Panel标题原文
      y_title = "Eligible individuals (%)", # 该Panel纵轴标题文本
      x_title = NULL, # 该Panel横轴标题文本
      x_labels = c(`2019` = "2019 guideline", `2025` = "2025 guideline"), # 横轴2019/2025两组类别标签
      title_size_pt = 16.2, # 该Panel标题字号（pt）
      title_gap_pt = 5, # 该Panel标题到图形主体的间距（pt）
      title_lineheight = 1.03, # 该Panel标题多行文本的行高系数
      axis_title_size_pt = 13.8, # 该Panel坐标轴标题字号（pt）
      axis_tick_size_pt = 13.0, # 该Panel坐标轴刻度字号（pt）
      axis_title_gap_pt = 5, # Panel纵轴标题与刻度文字间距（pt）
      x_tick_gap_pt = 3, # Panel横轴刻度数字与轴线间距（pt）
      panel_margin_pt = c(5, 5, 5, 5), # 该Panel上/右/下/左留白（pt）
      bar_width = .54, # 柱宽（图形坐标单位）
      bar_border_width = .55, # 柱轮廓线宽
      show_ci = FALSE, # 是否绘制置信区间误差线
      errorbar_width = .022, # 误差线端帽宽度
      errorbar_linewidth = .70, # 误差线线宽
      label_mode = "n", # 柱顶标签模式（n=人数/比例，具体显示由show_count控制）
      show_count = FALSE, # FALSE：柱顶只显示百分比；TRUE：显示人数（百分比）
      label_digits = 2L, # 柱顶百分比显示小数位数
      n_digits = 0L, # 柱顶人数显示小数位数
      label_size_mm = 4.25, # 柱顶数值标签字号（mm）
      label_lineheight = 1.08, # 柱顶双行数值标签行距系数
      label_offset_fraction = .040, # 柱顶标签相对柱高的向上偏移比例
      top_padding_fraction = .125, # 柱顶标签上方额外留白比例
      y_from_zero = TRUE, # 是否强制纵坐标从零开始
      y_tick_count = 5, # 纵轴自动刻度数量建议值
      y_axis_limits = c(0, NA_real_), # 纵轴显示范围；NA上限表示自动适配数据
      y_axis_breaks = NULL, # 纵轴刻度序列；NULL使用自动刻度
      x_expand = c(.36, .36) # 两根柱左右额外留白宽度
    )
  ),
  fig2 = list( # Figure 2专用的所有参数
    fig2_export_width = 15.5, # Figure 2导出宽度（英寸）
    fig2_export_height = 10, # Figure 2导出高度（英寸）
    panel_title_size = publication_design$typography$panel_title_pt * publication_design$typography$fig2_export_scale, # Figure 2所有Panel标题字号（pt）
    panel_title_gap_pt = 5, # Figure 2普通Panel标题和图形间距（pt）
    axis_text_size = publication_design$typography$axis_tick_pt * publication_design$typography$fig2_export_scale, # Figure 2各坐标轴刻度字号（pt）
    axis_title_size = publication_design$typography$axis_title_pt * publication_design$typography$fig2_export_scale, # Figure 2各坐标轴标题字号（pt）
    x_tick_gap_pt = c(A = 5, B = 5, C = 5, D = 6), # Figure 2 A/B/C/D各自刻度数字到轴线的间距（pt）
    x_title_gap_pt = c(A = 13, B = 13, C = 13, D = 15), # Figure 2 A/B/C/D各自横轴标题到刻度数字的间距（pt）
    value_label_size = publication_design$typography$fig2_bar_value_mm, # Panel A柱末百分比标签字号（mm）
    count_digits = 1L, # Panel A比例的小数位数
    bar_width = .60, # Panel A治疗转换比例柱宽
    group_axis_limits = c(0, 75), # Panel A分类比例横轴显示范围（%）
    group_axis_breaks = seq(0, 70, 10), # Panel A分类比例横轴刻度（%）
    group_colors = publication_design$colors$transitions, # Panel A六组柱色表
    observed_bar_color = publication_design$colors$observed, # Panel C观察风险柱形填色
    observed_bar_width = .62, # Panel C观察风险柱宽
    observed_bar_border_color = publication_design$colors$bar_edge, # Panel C柱形描边颜色
    bar_border_color = publication_design$colors$bar_edge, # 柱状图轮廓线颜色
    ci_color = publication_design$colors$text, # 误差棒颜色
    ci_linewidth = .60, # 误差棒线宽
    ci_cap_width = .16, # 误差棒端帽宽度
    point_size = 3.0, # 空心圆点大小
    point_stroke = 1.00, # 空心圆轮廓线宽
    point_shape = 21, # 空心圆点形状编号
    connector_color = publication_design$colors$connector, # 配对均值点之间连接线颜色
    connector_width = .58, # 配对均值点之间连接线线宽
    risk_axis_limits = c(0, 14), # 观察风险坐标轴范围（%）
    risk_axis_breaks = seq(0, 14, 2), # 观察风险坐标轴刻度（%）
    reduction_axis_limits = c(0, 2.1), # LDL-C平均降幅坐标范围（mmol/L）
    reduction_axis_breaks = seq(0, 2, .5), # LDL-C平均降幅坐标刻度（mmol/L）
    effect_axis_limits = c(NA_real_, NA_real_), # 2025−2019风险差坐标范围，双NA为自动
    effect_axis_breaks = NULL, # 2025−2019风险差坐标刻度；NULL为自动
    effect_axis_step = 0.50, # 自动风险差横轴刻度间隔（pp）
    effect_axis_tick_digits = 1L, # 风险差横轴刻度显示小数位数
    d_contribution_header = "Contribution (pp; 95% CI)", # Panel D右侧人口加权RD及95%CI列表头
    d_contribution_digits = 3L, # Panel D贡献点估计的小数位数（pp）
    d_contribution_ci_digits = 3L, # Panel D贡献置信区间的小数位数（pp）
    d_contribution_label_layout = "single_line", # Panel D贡献文字排版：单行或双行
    d_contribution_text_size = 4.20, # Panel D贡献点估计及CI文字字号（mm）
    d_contribution_label_lineheight = 1.00, # Panel D贡献双行文字行距系数
    d_contribution_header_size = 4.20, # Panel D贡献列表头字号（mm）
    d_contribution_x_gap_fraction = 0.048, # Panel D森林图横轴至Contribution列的水平间距比例
    d_contribution_right_margin_pt = 50, # Panel D右侧Contribution专用留白（pt）
    d_contribution_header_vjust = -2.00, # Contribution表头的垂直对齐系数（越负越上）
    d_panel_title_gap_pt = 20, # Panel D标题与绘图区之间留白（pt）；避免贡献表头碰撞
    d_contribution_additivity_tolerance = 2e-5, # 贡献之和与Overall RD核对时允许的风险单位误差
    d_rd_color = publication_design$colors$guideline[["2019"]], # Panel D风险差点及区间线的颜色
    dodge_width = .35, # 同一分类内多个系列的水平错位宽度
    row_guide_color = publication_design$colors$panel_guide, # 横向类别辅助线颜色
    row_guide_width = .35, # 横向类别辅助线宽度
    legend_size = publication_design$typography$legend_pt * publication_design$typography$fig2_export_scale, # Figure 2统一图例字号基础值（pt）
    b_legend_position = "right", # Panel B图例位置，right为右侧
    b_legend_ncol = 1L, # Panel B图例列数，1表示上下排列
    b_legend_text_size_pt = 12.4, # Panel B图例标签字号（pt）
    b_legend_key_width_pt = 12, # Panel B图例空心圆图标占宽（pt）
    b_legend_key_height_pt = 19, # Panel B图例空心圆图标占高（pt）
    b_legend_spacing_y_pt = 6, # Panel B两条图例项目间距（pt）
    b_legend_box_spacing_pt = 4, # Panel B绘图区和图例之间距离（pt）
    b_legend_margin_pt = c(0, 0, 0, 2), # Panel B图例边距：上/右/下/左（pt）
    zero_line_color = publication_design$colors$zero_line, # 零效应参考线颜色
    panelC_show_row_guides = FALSE, # 是否显示Panel C的水平参考线
    panel_letters = c(A = "A", B = "B", C = "C", D = "D"), # 各Panel标题的编号文字
    labels = list( # 分类、坐标或图例的英文显示名称集合
      A = "Treatment-transition groups", # Panel A对应的英文标题或轴标签
      B = "Mean minimum required LDL-C reduction", # Panel B对应的英文标题或轴标签
      C = "Observed 15-year ASCVD risk", # Panel C对应的英文标题或轴标签
      D = "2025 − 2019 risk difference" # Panel D对应的英文标题或轴标签
    ),
    x_titles = list( # 各Panel横坐标标题文字集合
      A = "Proportion of individuals (%)", # Panel A对应的英文标题或轴标签
      B = "Mean minimum required LDL-C reduction (mmol/L)", # Panel B对应的英文标题或轴标签
      C = "Observed ASCVD risk (%; 95% CI)", # Panel C对应的英文标题或轴标签
      D = "2025 − 2019 risk difference (pp; 95%CI)" # Panel D对应的英文标题或轴标签
    ),
    guideline_legend_labels = c(`2019` = "2019 guideline", `2025` = "2025 guideline") # 2019和2025图例显示名称
  ),
  fig3 = list( # Figure 3专用配置：新版三面板A/B/C，所有图形参数集中在此
    data_paths = list( # 仅注册Figure 3真正需要的已汇总pooled结果
      hvh = source_paths$fig3_hvh, # Panel A：2019/2025 × S0/SN/ST/SNT最终H/VH及从M升级人数
      counts = source_paths$transition_counts, # Panel B：SN/ST/SNT下六种干预变化组的汇总人数和比例
      effects = source_paths$transition_effects, # Panel C：各干预变化组2025−2019配对RD及95%CI
      overall_effects = source_paths$all_effects, # Panel C：总体2025−2019配对RD及95%CI
      contributions = source_paths$transition_contribution # Panel C：Module 12保存的人口加权贡献及95%CI
    ),
    horizon = report_style$primary_horizon, # 结果随访时点（年，通常15年）
    guidelines = c("2019", "2025"), # Panel A指南先后顺序，不允许自行交换语义
    scenarios = c("SN", "ST", "SNT"), # Panel B/C展示的三种RM策略及顺序
    classification_scenarios = c("SN", "ST", "SNT"), # Panel A仅展示应用RM后的三种策略顺序，不再展示S0
    transition_order = report_style$transition_order, # Panel B六种治疗转换组的纵向顺序
    comparative_order = c( # Panel C仅显示Overall及四种干预变化组；不展示Unchanged intensity
      "overall", # Overall位于最上方
      "newly_intervention", # 新增干预
      "intervention_removed", # 撤销干预
      "treatment_intensified", # 干预增强
      "treatment_deintensified" # 干预减弱
    ),
    transition_labels = report_style$labels$transition, # Panel B/C横向分类对应的英文标签
    overall_label = "Overall", # Panel C总体人群英文标签
    fig3_export_width = 16.0, # Figure 3矢量图导出宽度（英寸）；加宽画幅，为Panel B长分类标签及横向柱形提供空间
    fig3_export_height = 13.0, # Figure 3矢量图导出高度（英寸）；为Panel A轴下双层标签和下排C保留适当高度
    font_family = publication_design$typography$family, # 所有Panel使用的字体族
    rm_colors = publication_design$colors$rm, # S0灰、SN蓝、ST橙、SNT绿；各Panel颜色一致
    panel_prefix_separator = "  ", # Panel字母和标题之间的空格
    layout = list( # 三Panel拼图及全图画幅结构参数
      top_panel_widths = c(0.34, 0.66), # 上排Panel A与B相对宽度：A略收窄、B略放宽；适配A仅三种RM策略
      top_bottom_heights = c(1.00, 1.14), # 上排AB与下排C相对高度；减少下排过高造成的空白
      legend_height_ratio = 0.108, # 图例单行中间行高度略增，留给放大的标签文字
      legend_position = "between_AB_and_C", # 图例固定放在Panel A/B下方、Panel C上方；严格校验对应布局
      top_panel_gap_pt = 24, # Panel A与Panel B之间额外横向间距（pt），由A右侧和B左侧各分配一半，避免文字相碰
      row_panel_gap_pt = 10, # 上排AB与下排C额外纵向间距（pt）；保持标题不遮挡
      figure_margin_pt = c(3, 5, 3, 5) # 整图外边距：上/右/下/左（pt）；减少外围留白
    ),
    panel_A = list( # Panel A：按指南分组的6根并排堆叠柱，浅色虚线段表示RM升级
      letter = "A", # Panel编号
      title = "High/very-high risk classification", # Figure 3A新标题；仅修改显示文本
      x_title = NULL, # 横轴标题：已通过指南分面和策略刻度说明，默认不额外显示
      y_title = "Proportion classified as\nhigh/very-high risk (%)", # 纵轴标题及换行位置
      guideline_labels = c( # Panel A两组指南的分面条显示文字
        `2019` = "2019 guideline", # 2019指南分面标题
        `2025` = "2025 guideline" # 2025指南分面标题
      ),
      strategy_labels = c( # Panel A每组三根柱下方的策略缩写
        SN = "SN", # 临床RM策略
        ST = "ST", # 检测RM策略
        SNT = "SNT" # 联合RM策略
      ),
      title_size_pt = 18.0, # 标题字号（pt）
      title_gap_pt = 7, # 标题下缘至图形主体的距离（pt）
      title_wrap_width = 40L, # Panel A标题按字符宽度换行（仅显示层；原始英文标题不变），适配当前较窄A面板
      title_lineheight = 1.00, # 标题多行的行高倍率
      axis_title_size_pt = 16.4, # Panel A纵轴标题字号（pt）；适度放大
      axis_tick_size_pt = 14.8, # Panel A纵轴刻度数字字号（pt）；适度放大
      x_tick_size_pt = 13.7, # Panel A策略刻度字号（pt）；适度放大
      facet_label_size_pt = 16.4, # 两个指南标签字号（pt）；与Panel A纵轴标题一致，不加粗
      percent_label_size_mm = 4.15, # Panel A柱顶百分比字号（mm）；保证6组数值清晰不相撞
      percent_digits = 1L, # 柱顶最终H/VH比例保留的小数位数
      bar_width = 0.62, # Panel A柱形宽度（离散轴单位）；保持当前柱宽
      x_expand = c(0.65, 0.45), # Panel A每个指南分面的左右横轴扩展量（类别单位）；左侧0.65避免SN柱贴住纵轴
      bar_border_width = 0.38, # 柱形实色部分边框线宽
      bar_border_color = publication_design$colors$bar_edge, # 柱形边框颜色
      upgraded_alpha = 0.45, # RM升级部分的透明度（越小越浅；通过虚线边框辨识）
      upgraded_border_width = 0.65, # RM升级部分虚线边框线宽
      upgraded_border_linetype = "dashed", # RM升级部分边框的线型
      upgraded_border_color = publication_design$colors$bar_edge, # RM升级部分虚线框的颜色
      label_offset_pp = 1.8, # 柱顶百分比文字高于最终柱顶的距离（百分点）
      label_top_padding_pp = 4.0, # 顶部百分比文字与Panel顶部之间额外留白（百分点）
      y_axis_limits = c(0, NA_real_), # 纵轴起止值；NA表示按真实结果及上方标签自动确定
      y_axis_breaks = NULL, # 纵轴固定刻度；NULL表示按y_axis_step自动生成
      y_axis_step = 10, # 纵轴自动刻度间隔（百分点）
      y_axis_tick_digits = 0L, # 纵轴刻度数字小数位数
      show_y_grid = FALSE, # Panel A关闭全部背景网格线，保留坐标轴与刻度
      grid_color = publication_design$colors$panel_guide, # 水平网格线颜色
      grid_width = 0.40, # 水平网格线线宽
      facet_gap_pt = 24, # Panel A两版指南分面之间间距（pt）；适度拉开2019与2025两组柱状图
      strip_gap_pt = 9, # 指南年份分组标签与SN/ST/SNT刻度标签之间的垂直间距（pt）
      x_tick_gap_pt = 4, # 策略刻度文字与横轴轴线之间间距（pt）
      x_title_gap_pt = 7, # 如设置横轴标题，标题与刻度文字的间距（pt）
      y_tick_gap_pt = 4, # 纵轴刻度数字与轴线的间距（pt）
      y_title_gap_pt = 12, # Panel A纵轴标题至纵轴刻度数字的间距（pt）；增大避免过近
      axis_line_width = 0.55, # 坐标轴线宽
      axis_tick_width = 0.40, # 坐标轴刻度线宽
      panel_margin_pt = c(3, 3, 12, 2) # Panel A边距（上/右/下/左，pt）；底部为移至坐标轴下方的年份标签留空间
    ),
    panel_B = list( # Panel B：3种RM策略×6种干预变化组的横向并排条形
      letter = "B", # Panel编号
      title = "Treatment-transition groups", # Panel B标题
      x_title = "Proportion of individuals (%)", # Panel B横轴标题
      y_title = NULL, # Panel B纵轴标题；分类名称已经作为纵轴刻度
      title_size_pt = 18.0, # Panel B标题字号（pt）
      title_gap_pt = 7, # Panel B标题与绘图区顶部之间的间距（pt）
      title_lineheight = 1.00, # 标题行高倍率
      axis_title_size_pt = 16.4, # Panel B横坐标标题字号（pt）；适度放大
      axis_tick_size_pt = 14.8, # Panel B横纵坐标刻度字号（pt）；适度放大
      bar_label_size_mm = 4.50, # Panel B横条末端百分比字号（mm）；仅在show_value_labels=TRUE时使用，当前不显示
      show_value_labels = FALSE, # Panel B是否显示18根横向柱末端百分比；FALSE为完全删除数值文字图层
      percent_digits = 1L, # Panel B数值标签小数位数；仅在show_value_labels=TRUE时生效
      bar_height = 0.15, # 每条横向柱在类别轴上的高度
      bar_border_width = 0.37, # 横向柱形描边线宽
      bar_border_color = publication_design$colors$bar_edge, # 横向柱形边框色
      bar_label_offset_pp = 0.5, # 数值标签与柱形末端的间距（百分点）
      label_right_padding_pp = 5.0, # Panel B显示数值标签时在文字右侧保留的空间（百分点）
      bar_right_padding_pp = 16.0, # Panel B不显示柱末数值时，在最长条形右侧额外留16个百分点；坐标上界按10pp进位
      strategy_row_offsets = c(SN = 0.23, ST = 0, SNT = -0.23), # 同一组内三策略由上到下的偏移量
      y_axis_padding = 0.55, # 分类纵轴上下留白（类别坐标单位）
      x_axis_limits = c(0, NA_real_), # 横轴范围；上限NA为自动
      x_axis_breaks = NULL, # 横轴手动刻度；NULL为自动
      x_axis_step = 10, # 自动生成横轴刻度的间隔（百分点）
      x_axis_digits = 0L, # 横轴刻度小数位数
      show_row_guides = TRUE, # 是否绘制六种干预变化组的水平分隔线
      row_guide_color = publication_design$colors$panel_guide, # 六组参考分隔线颜色
      row_guide_width = 0.34, # 六组参考分隔线线宽
      x_tick_gap_pt = 5, # 横轴刻度数字与轴线之间距离（pt）
      x_title_gap_pt = 7, # 横轴标题与刻度数字之间距离（pt）
      y_tick_gap_pt = 4, # 分类名称与纵轴轴线之间距离（pt）
      axis_line_width = 0.55, # 坐标轴线宽
      axis_tick_width = 0.40, # 坐标轴刻度线宽
      panel_margin_pt = c(3, 28, 6, 4) # Panel B边距（上/右/下/左，pt）；减少多余空白
    ),
    panel_C = list( # Panel C：三策略的2025−2019配对RD森林图与右侧贡献(95%CI)独立数字列
      y_category_labels = c( # Panel C专用精简分类名称：减小纵轴文字占宽，左移森林图绘图区；不更改上游分组
        overall = "Overall", # 总体人群
        newly_intervention = "New intervention", # 新增干预组（简化required）
        intervention_removed = "Intervention removed", # 撤销干预组
        treatment_intensified = "Treatment intensified", # 治疗加强组
        treatment_deintensified = "Treatment de-intensified" # 治疗减弱组（简写，以减少最左侧标签区域）
      ),
      letter = "C", # Panel编号
      title = "Estimated 15-year ASCVD risk difference between guidelines", # Panel C标题
      x_title = "2025 − 2019 risk difference (percentage points)", # 森林图横坐标标题；pp为百分点
      y_title = NULL, # 森林图纵坐标标题；分类名称已经显示
      title_size_pt = 18.0, # Panel C标题字号（pt）
      title_gap_pt = 12, # Panel C标题与数据区域之间距离（pt）
      title_lineheight = 1.00, # 标题行高倍率
      axis_title_size_pt = 16.4, # Panel C横坐标标题字号（pt）；适度放大
      axis_tick_size_pt = 14.8, # Panel C横纵坐标刻度字号（pt）；适度放大
      point_size = 2.85, # RD空心圆点大小
      point_stroke = 0.95, # RD空心圆轮廓线宽
      point_shape = 21L, # RD空心圆点形状
      point_fill_color = publication_design$colors$white, # 空心圆内部填充色
      ci_linewidth = 0.58, # RD置信区间横线线宽
      ci_cap_height = 0.085, # RD置信区间端帽高度
      strategy_row_offsets = c(SN = 0.23, ST = 0, SNT = -0.23), # 同组内三策略RD及文字的上下错位
      y_axis_padding = 0.55, # 分类纵轴上下留白（坐标单位）
      x_axis_limits = c(NA_real_, NA_real_), # RD横轴范围；两个NA表示根据数据自动获取
      x_axis_breaks = NULL, # RD横轴刻度；NULL时采用自动间隔
      x_axis_step = 0.25, # RD横轴自动刻度间隔（百分点）
      x_axis_digits = 2L, # RD横轴数字小数位数
      show_zero_line = TRUE, # 是否绘制2025−2019差值为0的参考线
      zero_line_color = publication_design$colors$zero_line, # 零线颜色
      zero_line_type = "dashed", # 零线线型
      zero_line_width = 0.50, # 零线线宽
      show_row_guides = TRUE, # 是否绘制Panel C五组人群的水平参考线
      row_guide_color = publication_design$colors$panel_guide, # 参考线颜色
      row_guide_width = 0.34, # 参考线线宽
      contribution_header = "Contribution (pp; 95% CI)", # 右侧独立数值列的标题
      contribution_header_size_pt = 16.5, # Contribution列标题字号（pt）；从14增大至16.5
      contribution_text_size_mm = 5.20, # Contribution点估计及95%CI字号（mm）；从4.30放大
      contribution_digits = 3L, # 贡献点估计保留小数位数（pp）
      contribution_ci_digits = 3L, # 贡献区间端点保留小数位数（pp）
      contribution_show_ci = TRUE, # TRUE显示点估计及95%CI；FALSE仅显示点估计
      contribution_label_layout = "single_line", # 数字列显示模式：single_line单行，two_line两行
      contribution_lineheight = 0.95, # 两行Contribution文字的行高倍率
      contribution_plot_widths = c(0.78, 0.22), # Panel C森林图与Contribution数值列宽度比例为72:28，减少数值列右侧空白
      contribution_text_x = 0.02, # 数值与Contribution标题左对齐；统一通过下方左外边距产生留白
      contribution_left_margin_pt = 16, # Contribution列左侧独立留白（pt）；增加与森林图之间的间距
      contribution_title_gap_pt = 8, # Contribution标题与数值区的垂直距离（pt）
      contribution_right_margin_pt = 3, # 数字列最右侧额外留白（pt）
      x_tick_gap_pt = 5, # RD横轴刻度数字与轴线间距（pt）
      x_title_gap_pt = 11, # RD横轴标题与刻度数字间距（pt）
      y_tick_gap_pt = 2, # Panel C类别名称与纵轴之间的间距（pt）；缩小左侧文字占宽
      axis_line_width = 0.55, # 坐标轴线宽
      axis_tick_width = 0.40, # 刻度线线宽
      panel_margin_pt = c(3, 6, 5, 4) # Panel C森林图边距（上/右/下/左，pt）；减少左侧非数据留白
    ),
    legend = list( # Figure 3全宽单行图例：固定展示在Panel A/B下方、Panel C上方
      title = "RM strategy", # 左侧策略颜色分组标题；普通字体，不加粗
      labels = c(SN = "SN", ST = "ST", SNT = "SNT"), # 三种RM策略的实色方块
      component_title = "Classification by RM", # 右侧分类构成分组标题；普通字体，不加粗
      original_label = "Original H/VH", # 实色三色拼合示意：未经RM升级的原始H/VH
      upgraded_labels = c( # 三个虚线方块均保留完整的独立说明，按蓝、橙、绿显示
        SN = "Upgraded from M (SN)", # SN蓝色虚线方块
        ST = "Upgraded from M (ST)", # ST橙色虚线方块
        SNT = "Upgraded from M (SNT)" # SNT绿色虚线方块
      ),
      title_size_mm = 4.60, # 中间图例标题RM strategy/Classification by RM字号（mm），适当增大
      text_size_mm = 4.45, # 单行图例条目字号（mm）；略放大，不改变图例标题字号
      swatch_width = 0.018, # SN/ST/SNT策略实色方块宽度（图例局部归一化坐标）
      swatch_height = 0.46, # 单行图例色块高度（图例局部归一化坐标）
      swatch_border_width = 0.55, # 实色及虚线色块边框线宽
      component_swatch_width = 0.020, # 分类构成示意色块宽度
      component_alpha = 0.45, # RM升级段浅色透明度，与Panel A一致
      component_border_linetype = "dashed", # 升级色块边框线型
      x_title = 0.006, # RM strategy分组标题左端，配合全宽图例整体靠左
      x_swatches = c(SN = 0.105, ST = 0.174, SNT = 0.244), # 三种策略色块横向间距扩大
      component_title_x = 0.300, # Classification by RM标题左端
      original_swatch_x = 0.427, # 三色拼合Original H/VH色块中心
      upgraded_swatch_x = c(SN = 0.533, ST = 0.686, SNT = 0.843), # 三种RM升级色块横向间距扩大，同时防右侧裁切
      marker_text_gap = 0.015, # 色块中心至标签左端的横向位移（归一化坐标）
      y = 0.50, # 所有图例内容固定在同一水平行中心
      legend_margin_pt = c(4, 4, 4, 2) # 图例行外边距：上、右、下、左（pt）
    ),
    contribution_additivity_tolerance = 2e-5 # 贡献值之和与总体RD一致性校验容差（原始风险单位）
  ),
   fig4 = list( # Figure 4三面板：A总体RD、B相对S0变化、C各亚组RD；仅改显示，不改统计量
     fig4_export_width = 13.5, # 矢量图导出宽度（英寸）；防止九个亚组标签拥挤
     fig4_export_height = 8.4, # 两排布局输出高度（英寸）；A/B并排，C独占下排
     font_family = publication_design$typography$family, # 全图统一字体
     layout_top_widths = c(1, 1), # 第一行A与B的相对宽度（默认1:1）
     layout_row_heights = c(1.0, 1.32), # 第一行A/B与第二行C高度比例
     layout_column_gap_pt = 14, # 第一行A/B间水平留白（pt），由左右边距各承担一半
     layout_row_gap_pt = 14, # 第一行与第二行之间的垂直留白（pt）
     panel_heights = c(1.0, 1.32), # 兼容旧字段；实际拼图采用layout_row_heights
     panel_gap_pt = 9, # 兼容旧字段；实际行间距采用layout_row_gap_pt
     panel_margin_pt = c(3, 6, 4, 7), # 每个Panel基本边距：上/右/下/左（pt）
     panel_title_gap_pt = c(A=8, B=8, C=9), # 各Panel标题到数据区距离（pt）
     x_tick_gap_pt = c(A=5, B=5, C=5), # 每个Panel横轴类别到轴线的距离（pt）
     y_tick_gap_pt = c(A=4, B=4, C=4), # 每个Panel纵轴刻度到轴线的距离（pt）
     y_title_gap_pt = c(A=8, B=8, C=8), # 每个Panel纵轴标题与刻度的距离（pt）
     axis_text_size = 14.1, # 横纵轴刻度字号（pt）
     subgroup_text_size = 13.3, # Panel C九个亚组的标签字号（pt）
     axis_title_size = 15.1, # 各Panel纵轴标题字号（pt）
     panel_y_title_size = c(A=15.1, B=12.5, C=15.1), # 单独减小Panel B纵轴标题字号（pt），A/C不变
     panel_title_size = 17.4, # A/B/C面板标题字号（pt）
     point_size = 2.35, # Figure 4 A/B/C空心圆统一缩小；仍先画点再画CI
     point_stroke = 0.70, # Figure 4 A/B/C空心圆轮廓线宽；减细描边
     point_shape = 21L, # 空心圆形状编号
     ci_linewidth = 0.68, # 95%CI线宽
     ci_cap_width = c(A = 0.02, B = 0.02, C = 0.18), # Figure 4 A/B/C各自的CI端帽横向宽度（分别调节；越大端帽越长）
     dodge_width = 0.62, # Panel C四个策略的横向错位宽度
     arr_y_limits = c(0, 3), # 兼容保留：历史ARR纵轴范围；本版Figure 4不使用
     arr_y_breaks = seq(0,3,1), # 兼容保留：历史ARR刻度；本版Figure 4不使用
     rd_y_limits = c(-0.25, 0.42), # A/B/C风险差纵轴默认范围（pp），超出CI时自动扩展
     rd_y_breaks = seq(-0.2, 0.4, 0.2), # A/B/C风险差刻度（pp）
     delta_y_limits = c(-0.38, 0.08), # B相对S0变化的纵轴范围（pp）
     delta_y_breaks = c(-0.3, -0.2, -0.1, 0), # B相对S0变化纵轴刻度
     panel_y_limits = list(
       A = c(-0.08, 0.15),
       B = c(-0.20, 0.01),
       C = c(-0.15, 0.38)
     ), # A/B/C独立纵轴范围（pp）
     panel_y_breaks = list(A=c(-0.05,0,0.05,0.10),B=c(-0.3,-0.2,-0.1,0),C=c(-0.2,-0.1,0,0.1,0.2,0.3)), # A/B/C独立刻度
     zero_line_color = publication_design$colors$zero_line, # 零效应虚线颜色
     zero_line_width = 0.48, # 零效应虚线宽度
     group_divider_color = publication_design$colors$divider, # Panel C亚组大类虚线颜色
     group_divider_width = 0.42, # Panel C亚组大类虚线宽度
     legend_size = 14.2, # 图例标签字号（pt）
     legend_point_size = 2.65, # 底部图例空心圆相应缩小
     legend_key_width_pt = 14, # 图例色点占据的宽度（pt）
     legend_key_height_pt = 15, # 图例色点占据的高度（pt）
     legend_spacing_pt = 10, # 图例条目之间间距（pt）
     legend_margin_pt = c(4,3,2,3), # 底部统一图例的边距（pt）
     legend_ncol = 4L, # 底部四种RM策略单行排列
     favors_label_size = 4.20, # Panel C正负风险差倾向提示文字大小（mm）
     favors_x = 8.35, # Panel C取消Overall后，8个亚组右侧Favors文字横坐标
     favors_top_fraction = 0.05, # Favors 2019文字与上边界距离占纵轴比例
     favors_bottom_fraction = 0.05, # Favors 2025文字与下边界距离占纵轴比例
     panel_letters = c(A="A",B="B",C="C"), # Panel编号文本
     panel_title_wrap_width = c(A=34L, B=100L, C=90L), # 上排标题自动换行宽度（字符），防止并排遮挡
     top_x_labels = c( # 仅A/B横坐标使用的换行短标签；不改变底部完整图例
       S0 = "No RM\n(S0)",
       SN = "Non-laboratory\nRM (SN)",
       ST = "Testing RM\n(ST)",
       SNT = "Combined RM\n(SNT)"
     ),
     top_x_label_size_pt = 11.7, # A/B并排时的横轴标签字号（pt）
     panel_titles = c( # 三个Panel标题（与本次用户截图保持一致）
       A = "Overall 15-year risk difference", # 总体差异
       B = "Change in 2025–2019 RD by RM strategy", # 与S0比较的差异变化
       C = "Subgroup 15-year risk difference" # 九个总体/亚组差异
     ),
     y_titles = c( # A/B/C纵轴标题，单位均为百分点pp
       A = "2025 − 2019 RD (pp)", # Panel A指南间总体RD
       B = "RD with RM − RD without RM (pp)", # Panel B相对于未应用RM的RD变化
       C = "2025 − 2019 RD (pp)" # Panel C各亚组指南间RD
     ),
     legend_labels = c( # 使用新版统一情景名称S0/SN/ST/SNT
       S0 = "No RM (S0)", # 基础SCORE/SCORE2，不应用RM
       SN = "Non-laboratory RM (SN)", # 非实验室风险修饰因素策略
       ST = "Testing RM (ST)", # 检测RM规则
       SNT = "Combined RM (SNT)" # 联合RM规则
     ),
     subgroup_labels = publication_design$labels$short_subgroup, # Panel C显示8个亚组；Overall供Panel A使用
     group_dividers = c(2.5, 4.5, 6.5), # Panel C仅8个亚组：年龄/性别/SES/高血压之间的分隔线
     favors_2019 = "Favors 2019", # 正风险差倾向2019指南
     favors_2025 = "Favors 2025" # 负风险差倾向2025指南
   ),
  supplementary = list( # ==================== 附件图1–3专用参数：不影响正文 Figure 1–4 ====================
    # 一、输出画布与版面比例（PDF/SVG；可独立微调）
    width = 13.0, # 补充图画布宽度（英寸）；适当增加宽度以容纳单行色标标题
    height = 11.0, # 补充图画布高度（英寸）；两面板上下排布，适配放大字体
    panel_height_ratio = c(1, 1), # Panel A与Panel B高度比例
    plot_margin_pt = c(2, 8, 2, 4), # 单个Panel外边距：上、右、下、左（pt）
    title_margin_pt = 10, # Panel标题与表格主体的距离（pt）
    panel_title_size = 19.0, # Panel标题字号（pt）
    base_text_size_pt = 13.5, # 基础字号（pt）
    # 二、各文字元素的字号；geom_text的单位为mm、标题单位为pt
    cell_count_text_size = 4.55, # 单元格第一行：人数/比例或捕获率（mm，加粗）
    cell_ci_text_size = 3.80, # 单元格第二行：风险与95%CI或捕获率CI（mm）
    header_text_size = 4.55, # 顶部年份总标题与行列类别字号（mm）
    label_text_size = 4.55, # 左侧2019表头和风险类别字号（mm）
    region_text_size = 3.80, # 顶部高危分区说明字号（mm）
    footer_text_size = 3.80, # 底部高危分区说明字号（mm）
    legend_text_size = 3.80, # 色标刻度数字字号（mm）
    legend_title_text_size = 3.85, # 色标标题字号（mm），单行显示且独立于刻度
    cell_lineheight = 0.95, # 单元格两行文字各自行距系数
    cell_first_line_offset = 0.15, # 核心估计值相对单元格中心向上偏移（坐标单位）
    cell_second_line_offset = -0.17, # CI相对单元格中心向下偏移（坐标单位）
    # 三、边框、底色、文字对比度（均为附件图自己的设置）
    cell_border_color = "#CFD6DE", # 普通单元格描边颜色
    cell_border_width = 0.32, # 普通单元格线宽（ggplot linewidth单位）
    group_border_color = "#344454", # 高/极高风险四象限分割线颜色
    group_border_width = 1.00, # 高/极高风险分割线线宽
    header_fill = "#F0F3F7", # 顶部分类表头背景
    row_header_fill = "#F3F5F8", # 左侧分类表头背景
    region_fill = "#F6F8FA", # 顶部象限说明带底色
    footer_fill = "#F6F8FA", # 底部象限说明带底色
    empty_fill = "#F5F5F5", # 人数=0或未估计的单元格背景
    text_color = "#1F2933", # 数值与标题的深灰色
    secondary_text_color = "#58636E", # 区域说明与CI的辅助色
    dark_cell_text_color = "#FFFFFF", # 色块较深时的文字颜色
    use_adaptive_text_color = TRUE, # TRUE：按色块亮度自动切换白字/深字
    dark_cell_luminance_threshold = 0.43, # 自动白字的亮度阈值；越大越容易用白字
    # 四、颜色（附件图独立定义，不改变正文Figure 1的颜色）
    risk_limits = c(0, 25), # Panel A所有SN/ST/SNT共用0–25%的色标范围
    risk_breaks = seq(0, 25, 5), # Panel A色标刻度
    risk_colors = c("#FAFCFF", "#DCEAF7", "#95BDE0", "#2166AC"), # Panel A低饱和蓝色色阶
    capture_limits = c(0, 60), # Panel B所有SN/ST/SNT共用0–60%的色标范围
    capture_breaks = seq(0, 60, 10), # Panel B色标刻度
    capture_colors = c("#FEFBFF", "#EADAF5", "#BD91DC", "#773DB1"), # Panel B低饱和紫色色阶
    legend_steps = 120L, # 色标由多少个连续色块构成；越多越平滑
    # 四A、左上角嵌入式横向色标：位于行表头上方、2025列总表头左侧的空白区域
    # Panel A和Panel B分别绘制独立色标；不覆盖表头、不占用矩阵数据单元格。
    legend_position = "upper_left", # 固定为左上角嵌入；颜色含义仍由各Panel的单行标题说明
    legend_x = c(-1.57, 0.27), # 色条横向起止坐标；最大值应小于0.50（矩阵第一列左边缘）
    legend_y = c(6.02, 6.17), # 色条纵向下/上边界；必须位于5.50以上且7.10以下的空白区域
    legend_title_x = -1.57, # 色标单行标题的左对齐坐标
    legend_title_y = 6.58, # 色标标题纵坐标；高于色条且低于Panel标题
    legend_tick_y = 5.77, # 色标刻度文字纵坐标；低于色条，仍位于矩阵上方
    legend_risk_breaks = c(0, 10, 20, 25), # Panel A：蓝色色标刻度，色阶范围仍为0–25%
    legend_capture_breaks = c(0, 20, 40, 60), # Panel B：紫色色标刻度，色阶范围仍为0–60%
    legend_title_lineheight = 1.0, # 保留作为单行标题的行高参数
    # 五、矩阵几何及边距：均由主脚本配置，函数内不再写死
    side_header_x = c(-1.72, -0.69), # 2019总表头矩形的左右边界
    row_label_x = c(-0.68, 0.50), # 五行类别标签区域的左右边界
    row_label_text_x = -0.13, # 行类别文字的横向中心
    side_header_text_x = -1.205, # 2019总表头文字的横向中心
    top_header_y = c(6.57, 6.92), # 顶部2025合并表头纵向边界
    column_header_y = c(5.85, 6.56), # 2025五/四列类别标题纵向边界
    column_label_text_y = 6.17, # 2025类别列名的纵坐标
    region_header_y = c(5.50, 5.84), # 上部四象限说明带纵向边界
    region_header_text_y = 5.67, # 上部四象限说明文字的纵坐标
    footer_y = c(0.11, 0.48), # 下部四象限说明带纵向边界
    footer_text_y = 0.30, # 下部四象限说明文字的纵坐标
    horizontal_split_y = 3.50, # 2019高危/非高危的横向分界
    vertical_split_x = 2.50, # 2025高危/非高危的纵向分界
    matrix_bottom = 0.50, # 风险矩阵最下面一行的下边界
    matrix_top = 5.50, # 风险矩阵最上面一行的上边界
    x_left_padding = 0.07, # 左侧总表头以外的留白
    x_right_padding = 0.20, # 色标右侧留白
    y_limits = c(-0.04, 7.10), # 单Panel纵向绘图区范围
    # 六、是否省略完全无观察对象的2025 VH列
    hide_empty_2025_vh = TRUE, # TRUE：仅当三幅相应数据的VH列总人数=0时隐藏，25格完整CSV仍导出
    # 七、单元格数字格式及特殊状态（绝不裁剪/伪造上游CI）
    cell_count_digits = 1L, # 第一行人数所占总体的百分比保留位数
    cell_value_digits = 1L, # 观察风险与事件捕获比例保留位数
    cell_ci_digits = 1L, # CI端点保留小数位数
    show_less_than_for_tiny_positive = TRUE, # TRUE：附件图中所有非零但按一位小数会显示为0.0%的百分比（人数占比、风险/捕获点估计及CI端点）显示为<0.1%；真实零仍显示0.0%，仅改变文字，不改数值
    warn_unbounded_risk_ci = TRUE, # TRUE：风险95%CI越过0–100%时警告；绝不自动截断数值
    missing_ci_symbol = "—", # CI缺失仅在补充图显示为长破折号；原始CSV状态与数值均不修改
    scenario_labels = c(SN = "SN", ST = "ST", SNT = "SNT"), # 三幅图策略标识
    panel_letters = c(risk = "A", capture = "B"), # 两个Panel的字母
    labels = list( # 所有可编辑的英文标签
      top_header = "2025 risk categories under RM strategies", # 补充图顶部：2025年按RM策略确定的风险类别
      side_header = "2019\nrisk categories\nunder RM\nstrategies", # 补充图左侧：2019年按RM策略确定的风险类别（分行适配窄表头）
      row_categories = c( # 2019行标签
        L = "L (SCORE)", M_not_upgraded = "M ({scenario}−)",
        M_upgraded_H = "H ({scenario}↑)", H = "H (SCORE)", VH = "VH (SCORE)"
      ),
      column_categories = c( # 2025列标签
        L = "L (SCORE2)", M_not_upgraded = "M ({scenario}−)",
        M_upgraded_H = "H ({scenario}↑)", H = "H (SCORE2)", VH = "VH (SCORE2)"
      ),
      top_left = "Neither H/VH", # 2019与2025均非H/VH
      top_right = "Newly H/VH in 2025", # 仅2025为H/VH
      bottom_left = "No longer H/VH in 2025", # 仅2019为H/VH
      bottom_right = "H/VH under both", # 两版指南均为H/VH
      risk_panel = "Risk reclassification and observed %d-year ASCVD risk", # A面板标题
      capture_panel = "Observed %d-year ASCVD event capture", # B面板标题
      risk_legend = "%d-year ASCVD risk (%%)", # A色标标题
      capture_legend = "Event capture (%)", # B色标标题
      empty_risk = "Risk: —" # 无法估计风险时的文字
    )
  )
)
# Figure 3独立来源注册：新版仅读取A/B/C实际需要的5类pooled CSV，旧Panel不再加载。
source_paths$fig3_counts <- figure_layout$fig3$data_paths$counts # Figure 3B：三策略×六类干预转换
source_paths$fig3_effects <- figure_layout$fig3$data_paths$effects # Figure 3C：分组配对风险差
source_paths$fig3_overall_effects <- figure_layout$fig3$data_paths$overall_effects # Figure 3C：总体配对风险差
source_paths$fig3_contributions <- figure_layout$fig3$data_paths$contributions # Figure 3C：人口加权贡献及区间

# =============================================================================
# 1. 只读结果核验、four main and three supplementary figures 绘制及导出
# =============================================================================

# 严格按每图参数显示精度：不能用全局 2 位小数覆盖 Figure 1 示例的 1 位小数。
options(M01_REPORT_OVERWRITE = isTRUE(execution$overwrite)) # 是否允许覆盖本次报告图形文件（绝不影响上游分析）
source(file.path(code_dir, "Function_14_Report_figures_M01.R"), encoding = "UTF-8")
# 14 所需的全部数据校验与辅助函数均在唯一 Function_14 文件中；不加载13文件。
cfg <- rp_load_config(study_config_path, report_style)
style_fp <- rp_style_fingerprint(report_style)
manifest_root <- rp_dir(file.path(report_root, "manifests"))
qc_root <- rp_dir(file.path(report_root, "quality_control"))
figure_root <- rp_dir(file.path(report_root, "figures"))
figure_data_root <- rp_dir(file.path(report_root, "figure_data"))
# 14 是独立报告模块：不读取任何13输出或签名文件，也不依赖13的制表精度。
# 验证当前主脚本参数本身是否符合 07–12 已保存结果的研究约定。
if (isTRUE(execution$strict_internal_style)) {
  rp14_validate_own_style(report_style, figure_layout, execution)
}
# 只读取本次所选择 Figure 必需的已有 pooled CSV，避免单独导出 Figure 2 被其他图缺失文件阻断。
source_paths$fixed_moves <- file.path(analysis_dir,"treatment_transition_matrix_S0_vs_RM","12A_plot_row_transition_summary.csv")
source_paths$fixed_ldl <- file.path(analysis_dir,"treatment_transition_matrix_S0_vs_RM","12A_fixed_LDL_change_MI_mean.csv")
source_paths$fixed_rd <- file.path(analysis_dir,"fixed_S0_treatment_group_RD","12B_plot_fixed_group_RD.csv")
source_paths$all_effects_mi <- file.path(analysis_dir,"strategy_risks_and_effects","per_imputation","all_scenario_risks_and_effects.csv")
active_paths <- rp14_sources_for_figures(source_paths, execution$figures)
check <- rp_preflight(active_paths)
rp_save_csv(check, file.path(manifest_root, "14_input_source_manifest.csv"))
if (any(!check$exists)) stop("14 输入缺失：", paste(check$source[!check$exists], collapse = ", "))
src <- lapply(active_paths, rp_read)
# 所有旧/新情景标签只在 14 的内存中进行一次性规范化，原 pooled CSV 一律只读。
normalized <- rp14_normalize_source_scenarios(src,
  mode = execution$scenario_input_mode,
  reject_mixed = execution$reject_mixed_scenario_generations
)
src <- normalized$data
rp_save_csv(normalized$audit, file.path(qc_root, "14_scenario_code_audit.csv"))
for (k in names(src)) rp_check_ci(src[[k]], k)
if (execution$validate_only) {
  rp_save_csv(
    data.table(check = "source_exists", passed = check$exists, detail = check$source),
    file.path(qc_root, "14_preflight.csv")
  )
  rp_msg("14 validate_only 结束，未输出图像")
} else {
  build <- list(
    `1` = rp14_fig1, `2` = rp14_fig2, `3` = rp14_fig3, `4` = rp14_fig4,
    `5` = function(a, b, c) rp14_supplementary(a, b, c, "SN"),
    `6` = function(a, b, c) rp14_supplementary(a, b, c, "ST"),
    `7` = function(a, b, c) rp14_supplementary(a, b, c, "SNT")
  )
  for (num in execution$figures) {
    id <- as.character(num)
    if (!id %in% names(build)) stop("非法的 Figure 编号：", num)
    rp_msg("14 绘制 ", if (num <= 4) "Main figure " else "Supplementary figure ", if (num <= 4) num else num - 4L)
    obj <- rp_timed(paste("build figure", id), build[[id]](src, report_style, figure_layout))
    prefix <- if (num <= 4) sprintf("Figure_%02d", num) else sprintf("Supplementary_Figure_%02d", num - 4L)
    target_dir <- rp_dir(file.path(figure_root, if (num <= 4) "main_figures" else "supplementary_figures"))
    export_style <- report_style # 导出用样式副本，不影响 13 或其他 Figure
    if (id == "1") {
      export_style$figure$width <- figure_layout$fig1$export_width
      export_style$figure$height <- figure_layout$fig1$export_height
    }
    if (id == "2") {
      export_style$figure$width <- figure_layout$fig2$fig2_export_width
      export_style$figure$height <- figure_layout$fig2$fig2_export_height
    }
    if (id == "3") {
      export_style$figure$width <- figure_layout$fig3$fig3_export_width
      export_style$figure$height <- figure_layout$fig3$fig3_export_height
    }
    if (id == "4") {
      export_style$figure$width <- figure_layout$fig4$fig4_export_width
      export_style$figure$height <- figure_layout$fig4$fig4_export_height
    }
    if (num >= 5) {
      export_style$figure$width <- figure_layout$supplementary$width
      export_style$figure$height <- figure_layout$supplementary$height
    }
    rp_timed(paste("export figure", id), rp14_export(obj$plot, prefix, target_dir, export_style))
    if (execution$export_figure_data) rp14_export_data(obj$plot_data, figure_data_root, prefix)
  }
  saveRDS(
    list(style_fingerprint = style_fp, style = report_style, layout = figure_layout, study_config_path = study_config_path),
    file.path(manifest_root, "14_report_style_signature.rds")
  )
  capture.output(sessionInfo(), file = file.path(manifest_root, "14_sessionInfo.txt"))
  rp_msg("14 图形导出完成：", figure_root)
}


# =============================================================================
# 2. 附件图小样本独立交叉验证（放在绘图和原有QC之后，不影响Figure 1–4）
# 不来源于pooled CI；只取Module07逐插补个体数据与Module10逐插补AJ的Q/U。
# 如果本次仅绘制正文Figure 1–4，默认跳过该专用于补充图的核查。
# =============================================================================
if (isTRUE(supplementary_small_sample_qc$enabled) &&
    any(execution$figures %in% c(5L, 6L, 7L))) {
  rp14_qc_sparse_crossvalidation(
    analysis_dir = analysis_dir,
    study_cfg = cfg,
    opt = supplementary_small_sample_qc,
    scenarios = c('SN','ST','SNT'),
    qc_dir = qc_root
  )
}

