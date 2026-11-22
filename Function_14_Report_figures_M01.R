# =============================================================================
# Function_14_Report_figures_M01.R | Publication layout v3 | 一主脚本一函数文件
# 所有绘图输入仅来自 07-12 已保存的 pooled CSV；末尾独立QC例外：用标准软件核查逐插补AJ，不生成新分析估计。
# 全部研究路径和画图参数均位于 14 主脚本前部。本文件仅放置函数、逐图数据与图层。
# =============================================================================
# 本文件内置 14 专用通用工具函数，因此主脚本只需 source 本文件。
# =============================================================================
# 14 专用内置只读接口（已合并必要辅助函数，无须额外 source）
# 禁止 source 07–12 的主脚本或重拟合模型；独立QC只用标准软件复核AJ，不重做MI/Rubin。
# 所有路径、开关、图形参数都在 14 主脚本前部设置，函数文件不硬编码项目路径。
# =============================================================================
suppressPackageStartupMessages(library(data.table))

rp_msg <- function(...) cat(sprintf('[%s] ',format(Sys.time(),'%H:%M:%S')), ..., '\n', sep='')
rp_dir <- function(p) {dir.create(p,showWarnings=FALSE,recursive=TRUE); invisible(p)}
rp_md5 <- function(p) if(file.exists(p)) unname(as.character(tools::md5sum(p))) else NA_character_
rp_require <- function(p, label='输入') if(!file.exists(p)) stop(label,'不存在：',p,call.=FALSE) else invisible(p)
rp_assert <- function(x, cols, label='数据') {miss<-setdiff(cols,names(x));if(length(miss))stop(label,'缺失字段：',paste(miss,collapse=', '),call.=FALSE); invisible(TRUE)}
rp_guard <- function(path){if(file.exists(path)&&!isTRUE(getOption('M01_REPORT_OVERWRITE',TRUE)))stop('输出文件已存在，overwrite=FALSE：',path);invisible(TRUE)}
rp_save_csv <- function(x, path) {rp_dir(dirname(path));rp_guard(path);data.table::fwrite(x,path,na='',bom=TRUE);invisible(path)}
rp_file <- function(analysis_dir, board, file, layer='pooled') file.path(analysis_dir,board,layer,paste0(file,'.csv'))
rp_read <- function(path, required=NULL, validate_status=FALSE) {
  rp_require(path);d<-data.table::fread(path,encoding='UTF-8')
  if(!is.null(required))rp_assert(d,required,path)
  if(validate_status && 'status'%in%names(d) && any(is.na(d$status)|d$status!='estimable'))
    stop('存在不可估计的 pooled 行：',path,call.=FALSE)
  d
}
rp_optional <- function(path, required=NULL) {
  if(!file.exists(path))return(NULL)
  rp_read(path,required)
}
rp_forbidden_placeholder <- function(x) {
  y<-unlist(lapply(x,as.character),use.names=FALSE)
  if(any(grepl('x[.]xx|n [(]%[)]|FIGURE PLACEHOLDER',y,ignore.case=TRUE)))
    stop('检测到示例占位符，不允许把该文件标记为完成。',call.=FALSE)
  TRUE
}
rp_check_config <- function(cfg, style) {
  rp_assert(as.list(cfg),c('m','times','confidence_level'),'07 study_config')
  if(!all(style$horizons %in% cfg$times))stop('报告时间点不在 07 冻结 cfg$times 中')
  if(!(style$primary_horizon %in% cfg$times))stop('主时间点未由上游分析估计')
  if(!identical(sort(style$scenario_order),sort(c('S0','SN','ST','SNT'))))stop('RM 情景名称未匹配 07')
  invisible(TRUE)
}
rp_style_fingerprint <- function(style) {
  # 将两份主脚本的共同配置做成一个稳定签名，防止图表无意采用不同顺序/精度。
  x<-style[c('horizons','primary_horizon','scenario_order','scope_order','transition_order','digits','labels','colors')]
  p<-tempfile(fileext='.rds');on.exit(unlink(p));saveRDS(x,p,version=2);rp_md5(p)
}
rp_load_config <- function(p,style) {
  rp_require(p,'study_config.rds'); cfg<-readRDS(p);rp_check_config(cfg,style);cfg
}
rp_preflight <- function(paths) {
  data.table(source=names(paths),path=as.character(paths),exists=file.exists(unlist(paths)),
             bytes=as.numeric(file.info(unlist(paths))$size),
             md5=vapply(paths,rp_md5,character(1)))
}
rp_write_manifest <- function(paths,out)rp_save_csv(rp_preflight(paths),out)
rp_scope <- function(d,style,col='scope'){
  rp_assert(d,col);d<-copy(as.data.table(d));d<-d[get(col)%in%style$scope_order];d[,scope_index:=match(get(col),style$scope_order)];setorderv(d,'scope_index');d
}
rp_num <- function(x,digits=2L,group=FALSE) {
  v<-rep('—',length(x));ok<-is.finite(as.numeric(x));val<-as.numeric(x[ok]);
  val[abs(val)<0.5*10^(-digits)]<-0
  v[ok]<-formatC(val,format='f',digits=digits,big.mark=if(group)','else '',decimal.mark='.')
  v
}
rp_n <- function(x)rp_num(x,0L,TRUE)
rp_pct <- function(x,digits=1L,prob=TRUE)rp_num(if(prob)100*x else x,digits)
# 95%CI文本统一使用“下限 to 上限”，数值负号不受影响。
rp_ci <- function(estimate,lcl,ucl,digits=2L,scale=100,separator=" to ") {
  if(length(estimate)!=length(lcl)||length(lcl)!=length(ucl))stop('CI 长度不一致')
  a<-rp_num(estimate*scale,digits);b<-rp_num(lcl*scale,digits);c<-rp_num(ucl*scale,digits)
  ifelse(a=='—'|b=='—'|c=='—','—',paste0(a,' (',b,separator,c,')'))
}
rp_count_pct <- function(n,proportion,digits=1L) {
  a<-rp_n(n);b<-rp_pct(proportion,digits);ifelse(a=='—'|b=='—','—',paste0(a,' (',b,')'))
}
rp_pvalue <- function(x,digits=3L) {
  ifelse(!is.finite(x),'—',ifelse(x<10^(-digits),paste0('<',formatC(10^(-digits),format='f',digits=digits)),rp_num(x,digits)))
}
rp_row <- function(d,filters,required=TRUE,unique=TRUE) {
  z<-as.data.table(d)
  for(nm in names(filters)){if(!nm%in%names(z))stop('来源缺少筛选字段：',nm);v<-filters[[nm]];z<-z[get(nm)%in%v]}
  if(required && nrow(z)==0L)stop('缺少目标结果：',paste(names(filters),filters,sep='=',collapse='; '),call.=FALSE)
  if(unique && nrow(z)>1L)stop('目标结果重复：',paste(names(filters),filters,sep='=',collapse='; '),call.=FALSE)
  z
}
rp_eff <- function(d,scope,scenario,time,estimand,style,scale=100) {
  z<-rp_row(d,list(scope=scope,scenario=scenario,time=time,estimand=estimand))
  rp_ci(z$estimate,z$lcl,z$ucl,if(scale==1)style$digits$rr else if(grepl('risk_natural|risk_2019|risk_2025|observed',estimand))style$digits$risk else style$digits$rd,scale,
    separator=style$ci_interval_separator)
}
rp_group_labels <- function(style)style$labels$scope[style$scope_order]
rp_add_scope <- function(d,style) {
  d<-copy(as.data.table(d));d[,scope_label:=factor(style$labels$scope[scope],levels=style$labels$scope[style$scope_order])];d
}
rp_check_unique <- function(d,keys,label) {
  rp_assert(d,keys,label);if(anyDuplicated(d,by=keys))stop('结果键存在重复：',label,' / ',paste(keys,collapse=', '));TRUE
}
rp_check_ci <- function(d,label) {
  if(is.null(d)||!all(c('estimate','lcl','ucl')%in%names(d)))return(invisible(NULL))
  z<-d[is.finite(estimate)&is.finite(lcl)&is.finite(ucl)]
  if(nrow(z) && any(z$lcl>z$estimate+1e-10 | z$estimate>z$ucl+1e-10))stop('置信区间不包含点估计：',label)
  invisible(TRUE)
}
rp_template <- function(template_dir,number) {
  p<-file.path(template_dir,sprintf('Supplementary_Table_%02d_reference.csv',number))
  rp_require(p);data.table::fread(p,header=FALSE,fill=TRUE,encoding='UTF-8',blank.lines.skip=FALSE)
}
rp_output <- function(tab,number,outdir) {
  rp_save_csv(tab,file.path(outdir,sprintf('Supplementary_Table_%02d.csv',number)))
}
rp_output_main <- function(tab,name,outdir)rp_save_csv(tab,file.path(outdir,paste0(name,'.csv')))
rp_qc_record <- function(table,status,detail) data.table(item=table,status=status,detail=detail)
rp_format_reference <- function(d,style) {
  # 定义性表格内容来自投稿模板，不对文字进行自动医学修订。
  data.table::as.data.table(d)
}
rp_save_docx <- function(named_tables,path,style,heading='Supplementary Tables') {
  if(!requireNamespace('officer',quietly=TRUE)||!requireNamespace('flextable',quietly=TRUE))
    stop('生成 DOCX 需要安装 officer 和 flextable；CSV 已可独立输出。')
  doc<-officer::read_docx();doc<-officer::body_set_default_section(doc,officer::prop_section(page_size=officer::page_size(orient='landscape'),page_margins=officer::page_mar(top=.55,bottom=.55,left=.55,right=.55)));doc<-officer::body_add_par(doc,heading,style='heading 1')
  for(nm in names(named_tables)){
    tab<-as.data.frame(named_tables[[nm]],stringsAsFactors=FALSE)
    doc<-officer::body_add_par(doc,nm,style='heading 2')
    if(!ncol(tab)||!nrow(tab)){doc<-officer::body_add_par(doc,'No estimable rows in saved outputs.');next}
    ft<-flextable::flextable(tab);ft<-flextable::fontsize(ft,size=style$table$font_size,part='all')
    ft<-flextable::font(ft,fontname=style$table$font_family,part='all')
    ft<-flextable::theme_booktabs(ft);ft<-flextable::autofit(ft);ft<-flextable::fit_to_width(ft,max_width=10.5);ft<-flextable::set_table_properties(ft,layout="fixed",opts_word=list(repeat_headers=TRUE,split=FALSE))
    doc<-flextable::body_add_flextable(doc,ft)
  }
  rp_dir(dirname(path));rp_guard(path);print(doc,target=path)
  invisible(path)
}

rp_timed <- function(task,expr){t<-proc.time()[3];status<-'ERROR';rp_msg('START ',task);on.exit({rp_msg(status,' ',task,' | ',round(proc.time()[3]-t,2),' sec')},add=TRUE);z<-force(expr);status<-'OK';z}

if (!requireNamespace('ggplot2', quietly=TRUE)) stop('14 需要 ggplot2')
if (!requireNamespace('patchwork', quietly=TRUE)) stop('14 需要 patchwork')
if (!requireNamespace('scales', quietly=TRUE)) stop('14 需要 scales')


# -----------------------------------------------------------------------------
# 14 独立配置校验、按图选择输入：彻底取消对 13 脚本及签名文件的运行依赖。
# -----------------------------------------------------------------------------
rp14_validate_own_style <- function(s, layout, execution) {
  if (!identical(s$horizons, c(10,15)) || !identical(as.numeric(s$primary_horizon),15))
    stop('14 设置的随访时间与当前已定稿结果不一致；应为 10/15 年，主图 15 年')
  if (!identical(as.character(s$scenario_order),c('S0','SN','ST','SNT')))
    stop('14 的 scenario_order 必须为 S0、SN、ST、SNT')
  if (length(s$scope_order)!=9L || anyDuplicated(s$scope_order))
    stop('14 scope_order 必须包含九个不重复的总体/亚组')
  if (length(s$transition_order)!=6L || anyDuplicated(s$transition_order))
    stop('14 transition_order 必须包含六个不重复转换组')
  if (!all(s$scope_order %in% names(s$labels$scope)) ||
      !all(s$transition_order %in% names(s$labels$transition)))
    stop('14 的英文标签没有覆盖全部展示分组')
  if(!is.character(s$ci_interval_separator) || length(s$ci_interval_separator)!=1L ||
     !identical(s$ci_interval_separator, ' to '))
    stop('95%CI 区间连接词必须为两侧带空格的 to，避免与负号混淆')
  digs <- unlist(s$digits, use.names=TRUE)
  if (!all(is.finite(digs)) || any(digs<0) || any(digs!=floor(digs)))
    stop('14 显示小数位应为非负整数')
  ids <- as.integer(execution$figures)
  if (!length(ids) || anyNA(ids) || any(!ids %in% 1:7) || anyDuplicated(ids))
    stop('14 execution$figures 必须为不重复的 1–7 图号')
  # 出版风格契约：相同视觉语义只允许一个颜色来源；不依赖上游 CSV。
  if(is.null(layout$design$colors) || is.null(layout$design$typography))
    stop('14 缺少 publication_design 统一样式系统')
  dc <- layout$design$colors
  if(!identical(unname(s$colors$guideline_2019),unname(dc$guideline[['2019']])) ||
     !identical(unname(s$colors$guideline_2025),unname(dc$guideline[['2025']])))
    stop('2019/2025 指南颜色不一致；只能在 publication_design$colors$guideline 改色')
  if(!identical(unlist(s$colors$rm[s$scenario_order],use.names=FALSE),
                unname(dc$rm[s$scenario_order])))
    stop('RM 颜色不一致；统一修改 publication_design$colors$rm')
  if(!identical(unname(layout$fig1$guideline_colors),unname(dc$guideline)))
    stop('Figure 1 柱状图的指南色与全局颜色不一致')
  if(!identical(unname(layout$fig1$matrix_colors),unname(dc$risk_gradient)))
    stop('Figure 1蓝色色阶必须与正文统一设计系统一致')
  # 补充图配色已完全独立，不得反向改变正文Figure 1的蓝色色阶。
  for (nm in c('risk_colors','capture_colors')) {
    cols <- layout$supplementary[[nm]]
    if(length(cols)<2L || any(!grepl('^#[A-Fa-f0-9]{6}$', cols)))
      stop('Supplementary ',nm,' 必须由至少2个合法HEX颜色构成')
  }
  if(!identical(layout$fig2$group_colors,dc$transitions))
    stop('Figure 2 六类转换的颜色未与 publication_design 同步')
  cs<-unlist(dc[c('guideline','rm','transitions','risk_gradient','capture_gradient',
                   'observed','text','axis','subtle_text','panel_guide','divider','zero_line',
                   'header_fill','total_fill','empty_fill','border','strong_border','white',
                   'bar_edge','connector')],use.names=FALSE)
  if(any(!grepl('^#[A-Fa-f0-9]{6}$',cs)))stop('出版配色必须为 6 位 HEX，例如 #2166AC')
  if(!all(vapply(c('A','B','C','D'), function(z) {
       o <- layout$fig1[[paste0('panel_',z)]]
       is.list(o) && identical(o$letter,z) && is.character(o$title) &&
       length(o$title)==1L && nzchar(o$title) && is.numeric(o$title_size_pt) &&
       is.finite(o$title_size_pt) && o$title_size_pt>0 &&
       is.numeric(o$title_gap_pt) && is.finite(o$title_gap_pt) && o$title_gap_pt>=0
     }, logical(1L))) ||
      !identical(names(layout$fig2$labels),c('A','B','C','D')) ||
     !identical(names(layout$fig4$panel_titles),c('A','B','C')))
    stop('Figure 1/2/4 的Panel标题命名不完整或顺序不正确')
  if(length(layout$fig1$panel_A$matrix_row_labels)!=5L || length(layout$fig1$panel_A$matrix_column_labels)!=5L)
    stop('Figure 1 风险类别行列标题必须各包含 5 项（包括 Total）')
  if(!all(c('SN','ST','SNT') %in% names(layout$supplementary$scenario_labels)))
    stop('附件图情景标题需要 SN/ST/SNT')
  if(!all(c('L','M_not_upgraded','M_upgraded_H','H','VH') %in%
          names(layout$supplementary$labels$row_categories)) ||
     !all(c('L','M_not_upgraded','M_upgraded_H','H','VH') %in%
          names(layout$supplementary$labels$column_categories)))
    stop('附件图行列类别必须包含完整的 5 个风险类别')
  # 附件图单独的布局/色阶检查：不约束、不改变正文四图。
  supp <- layout$supplementary
  if (!is.logical(supp$hide_empty_2025_vh) || length(supp$hide_empty_2025_vh)!=1L ||
      is.na(supp$hide_empty_2025_vh)) stop('补充图hide_empty_2025_vh必须为TRUE或FALSE')
  for (nm in c('risk_limits','capture_limits')) {
    v <- supp[[nm]]
    if (!is.numeric(v) || length(v)!=2L || any(!is.finite(v)) || v[1]>=v[2])
      stop('附件图色阶范围无效：',nm)
  }
  for (nm in c('cell_count_text_size','cell_ci_text_size','header_text_size',
               'label_text_size','region_text_size','footer_text_size',
               'legend_text_size','group_border_width')) {
    v <- supp[[nm]]
    if (!is.numeric(v) || length(v)!=1L || !is.finite(v) || v<=0)
      stop('附件图参数应为正数：',nm)
  }
  # 嵌入式图例位置及稀疏刻度单独校验，不影响正文图的参数契约。
  for(nm in c('legend_risk_breaks','legend_capture_breaks')) {
    br <- supp[[nm]]
    bounds <- if(nm=='legend_risk_breaks') supp$risk_limits else supp$capture_limits
    if(!is.numeric(br) || length(br)<2L || any(!is.finite(br)) ||
       any(diff(br)<=0) || min(br)<bounds[1] || max(br)>bounds[2])
      stop('补充图左上角色标刻度无效：',nm)
  }
  if (!is.numeric(supp$legend_title_text_size) || length(supp$legend_title_text_size)!=1L ||
      !is.finite(supp$legend_title_text_size) || supp$legend_title_text_size<=0)
    stop('补充图色标标题字号必须是正数')
  if(layout$fig1$panel_A$matrix_fill_max < 100)
    stop('Figure 1 行百分比色阶须覆盖 0–100%，防止 40% 上限导致失真')
  # Figure 1 的跨四行/四列表头几何：对齐数据单元格，不触碰Total。
  a1 <- layout$fig1$panel_A
  for(nm in c('matrix_superheader_x','matrix_superheader_y',
              'matrix_sideheader_x','matrix_sideheader_y')) {
    lim <- a1[[nm]]
    if(!is.numeric(lim) || length(lim)!=2L || any(!is.finite(lim)) || lim[1]>=lim[2])
      stop('Figure 1 合并表头边界无效：',nm)
  }
  tw <- a1$matrix_tile_width; th <- a1$matrix_tile_height
  if(max(abs(a1$matrix_superheader_x-c(2-tw/2,5+tw/2)))>1e-5 ||
     max(abs(a1$matrix_sideheader_y-c(1-th/2,4+th/2)))>1e-5 ||
     abs(a1$matrix_sideheader_x[2]-(1-tw/2))>1e-5 ||
     abs(a1$matrix_superheader_y[1]-(5+th/2))>0.03)
    stop('Figure 1 合并表头没有与4个类别列/行对齐；请核对矩形边界与tile宽高')
  if(a1$matrix_x_limits[1]>a1$matrix_sideheader_x[1] ||
     a1$matrix_y_limits[2]<a1$matrix_superheader_y[2])
    stop('Figure 1 merged header outside matrix display limits')
  # Figure 1 新 D 面板与三列底行排版约定。
  f1 <- layout$fig1
  if(length(f1$layout$bottom_panel_widths)!=3L ||
     any(!is.finite(f1$layout$bottom_panel_widths)) ||
     any(f1$layout$bottom_panel_widths<=0) ||
     !isTRUE(f1$panel_A$show_row_total) || !isTRUE(f1$panel_A$show_col_total))
    stop('Figure 1 应有 B/C/D 三列正宽度，并保留 A 的行列合计')
  for(id in c('B','C','D')){
    o <- f1[[paste0('panel_',id)]]
    need <- c('2019','2025')
    if(!identical(names(o$x_labels),need) ||
       !is.character(o$y_title) || !nzchar(o$y_title) ||
       length(o$panel_margin_pt)!=4L || any(!is.finite(o$panel_margin_pt)) ||
       any(o$panel_margin_pt<0) ||
       !isTRUE(o$y_from_zero) || !o$label_mode%in%c('ci','n'))
      stop('Figure 1 Panel ',id,' 标签、字号或间距配置不完整')
  }
  if(!identical(f1$panel_D$label_mode,'n') || isTRUE(f1$panel_D$show_ci))
    stop('Figure 1 D 只能从 pooled Module 10 读取人数/比例，无 CI 时不得绘制 CI')
  if(!is.logical(f1$panel_D$show_count) ||
     length(f1$panel_D$show_count)!=1L || is.na(f1$panel_D$show_count))
    stop('Figure 1 D show_count 必须为单个 TRUE 或 FALSE')
  if(!is.character(layout$common$panel_prefix_separator) ||
     length(layout$common$panel_prefix_separator)!=1L)
    stop('Panel 字母间距需要单个字符串')
  if(2L %in% ids){
    f2 <- layout$fig2
    if(!is.numeric(f2$d_contribution_additivity_tolerance) ||
       length(f2$d_contribution_additivity_tolerance)!=1L ||
       !is.finite(f2$d_contribution_additivity_tolerance) || f2$d_contribution_additivity_tolerance<=0)
      stop('Figure 2D 加权 RD 一致性容差必须为正数')
    # 右侧贡献及 95%CI 的展示参数必须在主脚本中完整有效。
    if(!is.character(f2$d_contribution_label_layout) ||
       length(f2$d_contribution_label_layout)!=1L ||
       !f2$d_contribution_label_layout %in% c('single_line','two_line'))
      stop('Figure 2D 贡献CI排列只能为 single_line 或 two_line')
    for(nm in c('d_contribution_digits','d_contribution_ci_digits')){
      val<-f2[[nm]]
      if(!is.numeric(val) || length(val)!=1L || !is.finite(val) ||
         val<0 || val!=floor(val))
        stop('Figure 2D ',nm,' 必须为非负整数')
    }
    for(nm in c('d_contribution_text_size','d_contribution_header_size',
                'd_contribution_label_lineheight','d_contribution_right_margin_pt',
                'd_panel_title_gap_pt')){
      val<-f2[[nm]]
      if(!is.numeric(val) || length(val)!=1L || !is.finite(val) || val<=0)
        stop('Figure 2D ',nm,' 必须为正数')
    }
    # Figure 2 横坐标留白与 B 图例位置、字号、行列数校验。
    for(nm in c('x_tick_gap_pt','x_title_gap_pt')){
      val <- f2[[nm]]
      if(!is.numeric(val) || !identical(names(val),c('A','B','C','D')) ||
         any(!is.finite(val)) || any(val<0))
        stop('Figure 2 ',nm,' 应为按 A/B/C/D 命名的非负间距（pt）')
    }
    if(!identical(f2$b_legend_position,'right') ||
       !identical(as.integer(f2$b_legend_ncol),1L))
      stop('Figure 2 B 图例必须位于右侧且按一列纵向排列')
    for(nm in c('b_legend_text_size_pt','b_legend_key_width_pt',
                'b_legend_key_height_pt','b_legend_spacing_y_pt',
                'b_legend_box_spacing_pt')){
      val <- f2[[nm]]
      if(!is.numeric(val) || length(val)!=1L || !is.finite(val) || val<=0)
        stop('Figure 2 ',nm,' 必须为有限正数')
    }
    if(!is.numeric(f2$b_legend_margin_pt) || length(f2$b_legend_margin_pt)!=4L ||
       any(!is.finite(f2$b_legend_margin_pt)) || any(f2$b_legend_margin_pt<0))
      stop('Figure 2 图例边距必须为4个非负数（上/右/下/左，pt）')
    if(!is.numeric(f2$effect_axis_step) || length(f2$effect_axis_step)!=1L ||
       !is.finite(f2$effect_axis_step) || f2$effect_axis_step<=0 ||
       !is.numeric(f2$effect_axis_tick_digits) ||
       length(f2$effect_axis_tick_digits)!=1L ||
       !is.finite(f2$effect_axis_tick_digits) ||
       f2$effect_axis_tick_digits<0 ||
       f2$effect_axis_tick_digits!=floor(f2$effect_axis_tick_digits))
      stop('Figure 2 D 的刻度间距应为正数，显示小数位数应为非负整数')
    if(!is.numeric(f2$d_contribution_x_gap_fraction) ||
       length(f2$d_contribution_x_gap_fraction)!=1L ||
       !is.finite(f2$d_contribution_x_gap_fraction) ||
       f2$d_contribution_x_gap_fraction<0)
      stop('Figure 2D 贡献数值列间距必须为非负数')
  }
  if(3L %in% ids){
    f3 <- layout$fig3
    # Figure 3已由旧版A/B/C/D四面板调整为新版A/B/C三面板。
    # 标题分别放在panel_A、panel_B、panel_C的title字段，不再使用panel_titles。
    # 同时校验Panel编号、标题和顺序，防止旧版配置混入。
    expected_panels <- c('panel_A','panel_B','panel_C')
    if(!all(expected_panels %in% names(f3)) || 'panel_D' %in% names(f3) ||
       !identical(vapply(f3[expected_panels], function(p) p$letter, character(1)),
                  c(panel_A='A',panel_B='B',panel_C='C')) ||
       !all(vapply(f3[expected_panels], function(p) {
         is.character(p$title) && length(p$title)==1L &&
           !is.na(p$title) && nzchar(trimws(p$title))
       }, logical(1))))
      stop('Figure 3必须配置A/B/C三个Panel，且每个Panel有非空标题')
    # 新Figure 3必须是A/B/C三个Panel；分类图包含4个策略，B/C只用3种RM策略。
    if(!identical(f3$guidelines,c('2019','2025')) ||
       !identical(f3$classification_scenarios,c('SN','ST','SNT')) ||
       !identical(f3$scenarios,c('SN','ST','SNT')) ||
       !identical(vapply(f3[c('panel_A','panel_B','panel_C')],`[[`,character(1),'letter'),c(panel_A='A',panel_B='B',panel_C='C')))
      stop('Figure 3新版必须为A/B/C三Panel，Panel A策略为SN/ST/SNT，Panel B/C策略为SN/ST/SNT')
    if(!identical(names(f3$rm_colors),c('S0','SN','ST','SNT')) ||
       !identical(unname(f3$rm_colors),unname(dc$rm)))
      stop('Figure 3四种RM颜色必须与全局S0/SN/ST/SNT颜色一致')
    if(length(f3$transition_order)!=6L || length(f3$comparative_order)!=5L ||
       anyDuplicated(f3$comparative_order) ||
       !identical(f3$comparative_order[1L],'overall') ||
       !setequal(f3$comparative_order,
                 c('overall',setdiff(f3$transition_order,
                                     c('unchanged_untreated','unchanged_intensity')))))
      stop('Figure 3B须保留六种转换组；Panel C须为Overall及排除Unchanged intensity后的四组')
    if(!is.numeric(f3$layout$top_panel_widths) || length(f3$layout$top_panel_widths)!=2L ||
       any(!is.finite(f3$layout$top_panel_widths)) || any(f3$layout$top_panel_widths<=0) ||
       !is.numeric(f3$layout$top_bottom_heights) || length(f3$layout$top_bottom_heights)!=2L ||
       any(!is.finite(f3$layout$top_bottom_heights)) || any(f3$layout$top_bottom_heights<=0) ||
       !is.numeric(f3$layout$legend_height_ratio) || f3$layout$legend_height_ratio<=0 ||
       !is.numeric(f3$fig3_export_width) || f3$fig3_export_width<=0 ||
       !is.numeric(f3$fig3_export_height) || f3$fig3_export_height<=0)
      stop('Figure 3画幅、面板宽高及图例高度参数必须有效')
    for(nm in c('panel_A','panel_B','panel_C')){
      p <- f3[[nm]]
      if(!is.character(p$title) || !nzchar(p$title) ||
         !is.numeric(p$title_size_pt) || p$title_size_pt<=0 ||
         !is.numeric(p$title_gap_pt) || p$title_gap_pt<0 ||
         !is.numeric(p$axis_tick_size_pt) || p$axis_tick_size_pt<=0 ||
         !is.numeric(p$axis_title_size_pt) || p$axis_title_size_pt<=0 ||
         !is.numeric(p$panel_margin_pt) || length(p$panel_margin_pt)!=4L ||
         any(!is.finite(p$panel_margin_pt)) || any(p$panel_margin_pt<0))
        stop('Figure 3 ',nm,' 标题、字号、边距参数无效')
    }
    if(!identical(f3$panel_A$title,'High/very-high risk classification'))
      stop('Figure 3A标题必须为High/very-high risk classification')
    if(!identical(names(f3$panel_A$strategy_labels),f3$classification_scenarios) ||
       !identical(names(f3$panel_A$guideline_labels),f3$guidelines))
      stop('Figure 3A必须为2019/2025各SN/ST/SNT三根柱')
    # Figure 3更新：每个分面的横轴留白与中间单行图例均由主脚本统一配置。
    if(!is.numeric(f3$panel_A$x_expand) || length(f3$panel_A$x_expand)!=2L ||
       any(!is.finite(f3$panel_A$x_expand)) || any(f3$panel_A$x_expand<0) ||
       f3$panel_A$x_expand[1L]<=f3$panel_A$bar_width/2)
      stop('Figure 3A横轴左侧留白必须大于半根柱宽，防止最左侧柱形贴近纵轴')
    lg <- f3$legend
    # 图例必须位于A/B与C之间，且所有标题/色块/标签只用一个Y坐标。
    if(!identical(f3$layout$legend_position,'between_AB_and_C') ||
       !identical(lg$component_title,'Classification by RM') ||
       !identical(names(lg$labels),f3$classification_scenarios) ||
       !identical(names(lg$x_swatches),f3$classification_scenarios) ||
       !identical(names(lg$upgraded_swatch_x),f3$classification_scenarios) ||
       !identical(names(lg$upgraded_labels),f3$classification_scenarios) ||
       !is.numeric(lg$y) || length(lg$y)!=1L || !is.finite(lg$y) ||
       lg$y<=0 || lg$y>=1 ||
       any(!is.finite(c(lg$x_title,lg$component_title_x,lg$original_swatch_x,
         lg$x_swatches,lg$upgraded_swatch_x))) ||
       any(c(lg$x_title,lg$component_title_x,lg$original_swatch_x,
         lg$x_swatches,lg$upgraded_swatch_x)<0) ||
       any(c(lg$x_title,lg$component_title_x,lg$original_swatch_x,
         lg$x_swatches,lg$upgraded_swatch_x)>1) ||
       !all(diff(as.numeric(lg$x_swatches))>0) ||
       !all(diff(as.numeric(lg$upgraded_swatch_x))>0) ||
       any(!nzchar(lg$upgraded_labels)) ||
       !is.numeric(lg$title_size_mm) || lg$title_size_mm<=0 ||
       !is.numeric(lg$text_size_mm) || lg$text_size_mm<=0)
      stop('Figure 3中间单行图例配置无效：需包含SN/ST/SNT和三色RM升级虚线色块')
    # Panel B的数值显示开关及无标签缓冲必须有效，防止回退为旧版拥挤图。
    if(!is.logical(f3$panel_B$show_value_labels) ||
       length(f3$panel_B$show_value_labels)!=1L ||
       is.na(f3$panel_B$show_value_labels) ||
       !is.numeric(f3$panel_B$bar_right_padding_pp) ||
       length(f3$panel_B$bar_right_padding_pp)!=1L ||
       !is.finite(f3$panel_B$bar_right_padding_pp) ||
       f3$panel_B$bar_right_padding_pp<0)
      stop('Figure 3B数值开关或无文字时的右侧留白参数无效')
    if(length(f3$panel_B$strategy_row_offsets)!=3L ||
       !identical(names(f3$panel_B$strategy_row_offsets),f3$scenarios) ||
       length(f3$panel_C$strategy_row_offsets)!=3L ||
       !identical(names(f3$panel_C$strategy_row_offsets),f3$scenarios))
      stop('Figure 3B/C三种RM策略的错位参数必须齐全')
    if(max(abs(f3$panel_B$strategy_row_offsets))+f3$panel_B$bar_height/2>=0.5 ||
       max(abs(f3$panel_C$strategy_row_offsets))+f3$panel_C$ci_cap_height/2>=0.5)
      stop('Figure 3B/C三策略子行之间发生重叠或超出类别边界')
    # Panel C轴标签采用专用英文简写，覆盖显示的五组，不更改底层分组键。
    if(!identical(names(f3$panel_C$y_category_labels),f3$comparative_order) ||
       anyNA(f3$panel_C$y_category_labels) || any(!nzchar(f3$panel_C$y_category_labels)))
      stop('Figure 3C专用纵轴短标签必须完整对应五个显示分组')
    if(!identical(f3$panel_C$contribution_label_layout %in% c('single_line','two_line'),TRUE) ||
       !is.numeric(f3$panel_C$contribution_plot_widths) ||
       length(f3$panel_C$contribution_plot_widths)!=2L ||
       any(!is.finite(f3$panel_C$contribution_plot_widths)) ||
       any(f3$panel_C$contribution_plot_widths<=0) ||
       f3$panel_C$contribution_plot_widths[1L]/sum(f3$panel_C$contribution_plot_widths)<0.5 ||
       f3$panel_C$contribution_plot_widths[1L]/sum(f3$panel_C$contribution_plot_widths)>0.85 ||
       !is.numeric(f3$panel_C$contribution_left_margin_pt) ||
       length(f3$panel_C$contribution_left_margin_pt)!=1L ||
       !is.finite(f3$panel_C$contribution_left_margin_pt) ||
       f3$panel_C$contribution_left_margin_pt<0 ||
       !is.numeric(f3$panel_C$contribution_text_x) ||
       length(f3$panel_C$contribution_text_x)!=1L ||
       !is.finite(f3$panel_C$contribution_text_x) ||
       f3$panel_C$contribution_text_x<0 || f3$panel_C$contribution_text_x>=1)
      stop('Figure 3C森林图区宽度比例须在50%至85%之间，且左侧留白/数值位置须有效')
    if(!is.numeric(f3$contribution_additivity_tolerance) ||
       length(f3$contribution_additivity_tolerance)!=1L ||
       !is.finite(f3$contribution_additivity_tolerance) ||
       f3$contribution_additivity_tolerance<=0)
      stop('Figure 3C贡献汇总校验容差无效')
    for(x in c('panel_A','panel_B','panel_C')){
      a <- f3[[x]]
      lim <- if(x=='panel_A')a$y_axis_limits else a$x_axis_limits
      br <- if(x=='panel_A')a$y_axis_breaks else a$x_axis_breaks
      step <- if(x=='panel_A')a$y_axis_step else a$x_axis_step
      if(!is.numeric(lim) || length(lim)!=2L || any(is.nan(lim)) ||
         !(all(is.na(lim)) && x=='panel_C') && !is.finite(lim[1]) ||
         !is.numeric(step) || length(step)!=1L || !is.finite(step) || step<=0 ||
         (!is.na(lim[2]) && (!is.finite(lim[2]) || lim[2]<=lim[1])) ||
         (!is.null(br) && (!is.numeric(br) || any(!is.finite(br)) ||
            (length(br)>1L && any(diff(br)<=0)))))
        stop('Figure 3 ',x,' 横纵坐标上限或刻度配置无效')
    }
    if(!all(vapply(f3$data_paths,function(x)is.character(x)&&length(x)==1L&&nzchar(x),logical(1))))
      stop('Figure 3有未配置的pooled数据路径')
  }
  if(4L %in% ids){
    # Figure 4双列上排与单列下排的版式参数检查。
    f4 <- layout$fig4
    for(nm in c('layout_top_widths','layout_row_heights')) {
      v <- f4[[nm]]
      if(!is.numeric(v) || length(v)!=2L || any(!is.finite(v)) || any(v<=0))
        stop('Figure 4 ', nm, ' 必须为两个有限正数')
    }
    for(nm in c('layout_column_gap_pt','layout_row_gap_pt','top_x_label_size_pt')) {
      v <- f4[[nm]]
      if(!is.numeric(v) || length(v)!=1L || !is.finite(v) || v<0)
        stop('Figure 4 ', nm, ' 必须为非负数')
    }
    if(!identical(names(f4$top_x_labels),c('S0','SN','ST','SNT')) ||
       any(!nzchar(f4$top_x_labels)) ||
       !identical(names(f4$panel_title_wrap_width),c('A','B','C')) ||
       any(!is.finite(f4$panel_title_wrap_width)) || any(f4$panel_title_wrap_width<10))
      stop('Figure 4 A/B横轴短标签或Panel标题换行参数不完整')
    if(!is.numeric(f4$panel_y_title_size) ||
       !identical(names(f4$panel_y_title_size),c('A','B','C')) ||
       any(!is.finite(f4$panel_y_title_size)) || any(f4$panel_y_title_size<=0))
      stop('Figure 4 A/B/C纵轴标题字号必须为三个正数')
    # Figure 4三Panel CI端帽独立参数校验，避免出现缺失、负值或无效名称。
    caps <- layout$fig4$ci_cap_width
    if(!is.numeric(caps) || !identical(names(caps),c('A','B','C')) ||
       length(caps)!=3L || any(!is.finite(caps)) || any(caps<0))
      stop('Figure 4 ci_cap_width必须是按A/B/C命名的三个非负有限数值')
    if(!identical(unname(layout$fig4$legend_labels[['SN']]),'Non-laboratory RM (SN)'))
      stop('Figure 4 SN展示名称应为Non-laboratory RM (SN)')
    # Figure 4C仅8个亚组，确保分割线与右侧文字位置匹配。
    if(!identical(as.numeric(layout$fig4$group_dividers),c(2.5,4.5,6.5)) ||
       !is.numeric(layout$fig4$favors_x) || length(layout$fig4$favors_x)!=1L ||
       !is.finite(layout$fig4$favors_x) || layout$fig4$favors_x<=8)
      stop('Figure 4C排除Overall后，分隔线及Favors文字坐标需匹配8个亚组')
  }
  sizes <- unlist(layout$design$typography[-1L],use.names=FALSE)
  if(any(!is.finite(sizes)) || any(sizes<=0))
    stop('publication_design 中字体大小和宽图字号倍率必须为正数')
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Scenario-name compatibility; this only renames two IDs in *report memory*.
# It does not change RM membership, outcomes, risk differences, weights, or saved CSV.
# Never combine old and new scenario IDs in one pooled source generation.
# -----------------------------------------------------------------------------
rp14_normalize_source_scenarios <- function(source_data, mode="allow_legacy_rename",
                                           reject_mixed=TRUE) {
  if(!mode %in% c("allow_legacy_rename", "strict_new"))
    stop("14 scenario_input_mode must be allow_legacy_rename or strict_new")
  old <- c("SH","SHT")
  new <- c("SN","SNT")
  canonical <- c("S0","SN","ST","SNT")
  permitted <- c(canonical,old)
  audit <- vector("list",length(source_data)); names(audit)<-names(source_data)
  groups <- character(0)
  for(nm in names(source_data)) {
    d<-source_data[[nm]]
    if(!"scenario"%in%names(d)) {
      audit[[nm]]<-data.table::data.table(source=nm,encoding="no_scenario_column",
         original_scenarios="",normalized_scenarios="",renamed_rows=0L)
      next
    }
    vals<-unique(as.character(d$scenario))
    vals<-vals[!is.na(vals)]
    unknown<-setdiff(vals,permitted)
    if(length(unknown))stop("14 upstream ",nm," has unexpected scenario codes: ",paste(unknown,collapse=", "))
    has_old<-any(vals%in%old)
    has_new<-any(vals%in%new)
    if(has_old && has_new)
      stop("14 upstream ",nm," mixes SH/SHT with SN/SNT; do not merge different run generations")
    typ<-if(has_old)"legacy" else if(has_new)"new" else "S0/ST_only"
    if(typ%in%c("legacy","new") && !nm %in% c("fixed_moves","fixed_ldl","fixed_rd")) groups<-c(groups,typ)
    if(has_old && identical(mode,"strict_new"))
      stop("14 upstream ",nm," uses SH/SHT but report requires SN/SNT; rerun upstream or enable explicit alias mode")
    if(has_old) {
      x<-data.table::copy(data.table::as.data.table(d))
      x[scenario=="SH",scenario:="SN"]
      x[scenario=="SHT",scenario:="SNT"]
      # 12+ pooled report contains one row per outcome/transition key; changing
      # source IDs must not produce any duplicate rows.
      source_data[[nm]]<-x
    }
    audit[[nm]]<-data.table::data.table(source=nm,encoding=typ,
        original_scenarios=paste(sort(vals),collapse="|"),
        normalized_scenarios=paste(sort(unique(as.character(source_data[[nm]]$scenario))),collapse="|"),
        renamed_rows=sum(as.character(d$scenario)%in%old,na.rm=TRUE))
  }
  if(isTRUE(reject_mixed) && length(unique(groups))>1L)
    stop("14 input CSVs contain BOTH old SH/SHT and new SN/SNT generations; ",
         "refresh 07-12 outputs consistently before publication reporting. ",
         "See sources in the input manifest.")
  q<-data.table::rbindlist(audit,use.names=TRUE,fill=TRUE)
  q[,mapping_rule:="SH->SN; SHT->SNT (report memory only)"]
  q[,mixed_generation_check:=if(length(unique(groups))>1L)"MIXED" else "PASS"]
  n<-sum(q$renamed_rows)
  if(n>0L)warning("14: normalized ",n," legacy scenario-code rows in memory: ",
                "SH->SN; SHT->SNT. The upstream CSVs are UNMODIFIED; see 14_scenario_code_audit.csv")
  list(data=source_data,audit=q)
}

rp14_sources_for_figures <- function(paths, figures) {
  required<-list(`1`=c('reclassification','reclassification_direction','hvh_observed','fig1_eligibility'),
    `2`=c('transition_counts','transition_observed','transition_effects','population_mean','all_effects','transition_contribution'),
    `3`=c('fig3_hvh','fig3_counts','fig3_effects','fig3_overall_effects','fig3_contributions','fixed_moves','fixed_ldl','fixed_rd','all_effects'),
    `4`=c('all_effects','all_effects_mi'),
    `5`=c('joint_counts','joint_risks'),`6`=c('joint_counts','joint_risks'),
    `7`=c('joint_counts','joint_risks'))
  if(any(!as.character(figures)%in%names(required)))stop('Figure ids: 1:4 main, 5:7 supplementary')
  keys<-unique(unlist(required[as.character(figures)]));if(any(!keys%in%names(paths)))stop('Missing source registry key')
  paths[keys]
}

# -----------------------------------------------------------------------------
# 模块 00：共享出版样式、字体与数据转换（不做上游统计分析）
# -----------------------------------------------------------------------------
rp14_theme <- function(s, layout, legend='none', margins=NULL) {
  v <- layout$common
  cc <- layout$design$colors
  if (is.null(margins)) margins <- v$figure_margin_pt
  ggplot2::theme_classic(base_size=s$figure$base_size,
                         base_family=s$figure$font_family) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face='bold', size=v$panel_title_size, hjust=0,
                                       margin=ggplot2::margin(b=v$panel_title_gap_pt)),
      plot.subtitle=ggplot2::element_text(size=v$axis_text_size, face='italic',
                                          margin=ggplot2::margin(b=6)),
      axis.title=ggplot2::element_text(size=v$axis_title_size, color=cc$text),
      axis.text=ggplot2::element_text(size=v$axis_text_size, color=cc$text),
      axis.text.x=ggplot2::element_text(angle=0, hjust=.5, vjust=1),
      axis.line=ggplot2::element_line(color=cc$axis, linewidth=.55),
      axis.ticks=ggplot2::element_line(color=cc$axis, linewidth=.4),
      panel.grid=ggplot2::element_blank(),
      legend.position=legend, legend.title=ggplot2::element_blank(),
      legend.text=ggplot2::element_text(size=v$axis_text_size),
      legend.key=ggplot2::element_blank(),
      legend.spacing.x=grid::unit(9,'pt'),
      plot.margin=ggplot2::margin(margins[1], margins[2], margins[3], margins[4]),
      plot.background=ggplot2::element_rect(fill=cc$white, color=NA)
    )
}

rp14_panel <- function(p, letter, title, layout, subtitle=NULL) {
  # 将 Panel 字母和标题合成同一个 title grob，彻底消除 ggplot 独立 tag/title 的左对齐偏差。
  # plot.title.position='plot' 保证标题从该图的左边界开始，而非从坐标轴绘图区开始。
  p + ggplot2::labs(tag=NULL, title=paste0(letter,layout$common$panel_prefix_separator,title), subtitle=subtitle) +
    ggplot2::theme(
      plot.title.position='plot',
      plot.title=ggplot2::element_text(face='bold',size=layout$common$panel_title_size,
        hjust=0,margin=ggplot2::margin(l=0,b=layout$common$panel_title_gap_pt)))
}

rp14_to_pct <- function(d) {
  rp_assert(d, c('estimate','lcl','ucl'))
  z <- data.table::copy(data.table::as.data.table(d))
  z[, `:=`(y=100*estimate, lo=100*lcl, hi=100*ucl)]
  z
}

rp14_finite_ci <- function(d, label, allow_na=FALSE) {
  rp_assert(d, c('estimate','lcl','ucl'), label)
  bad <- !is.finite(d$estimate) | !is.finite(d$lcl) | !is.finite(d$ucl)
  if (any(bad) && !allow_na) stop(label,' 存在无效或缺失 CI，不能静默删行')
  ok <- !bad
  if (any(d$lcl[ok]>d$estimate[ok] | d$estimate[ok]>d$ucl[ok]))
    stop(label,' CI 未包含点估计')
  invisible(TRUE)
}

rp14_scope_factor <- function(x,s) factor(as.character(x), levels=s$scope_order)
rp14_scope_scale <- function(s,layout, expand=ggplot2::expansion(add=c(.55,.55))) {
  labs <- layout$common$subgroup_short_labels
  if (!all(s$scope_order %in% names(labs))) stop('绘图亚组短名称未覆盖 scope_order')
  ggplot2::scale_x_discrete(limits=s$scope_order,
    labels=function(x)unname(labs[as.character(x)]), expand=expand, drop=FALSE)
}

rp14_scope_dividers <- function(layout, enabled=TRUE) {
  if (!isTRUE(enabled)) return(NULL)
  ggplot2::geom_vline(xintercept=layout$common$group_dividers,
    linetype='dashed', color=layout$common$group_divider_color,
    linewidth=layout$common$group_divider_width)
}

rp14_positive_limit <- function(x, pad=.18) {
  v <- as.numeric(x[is.finite(x)])
  if (!length(v)) stop('无法根据置信区间计算纵轴范围')
  c(0, max(0.05,max(v))*(1+pad))
}
rp14_range_limits <- function(lower,upper,pad=.15,always_zero=TRUE,symmetric=FALSE) {
  x<-c(lower,upper);x<-x[is.finite(x)]
  if (!length(x))stop('坐标轴范围不能由空/缺失 CI 推断')
  if (symmetric) {
    limit<-max(0.01,max(abs(x)))*(1+pad)
    return(c(-limit,limit))
  }
  lo<-min(x);hi<-max(x);span<-max(hi-lo,.05)
  if(always_zero){lo<-min(0,lo);hi<-max(0,hi)}
  c(lo-if(lo<0)span*pad else 0,hi+span*pad)
}

rp14_guideline_palette <- function(s) c('2019'=s$colors$guideline_2019,
                                        '2025'=s$colors$guideline_2025)
rp14_rm_palette <- function(s,names_rm) {
  unlist(s$colors$rm[names_rm],use.names=TRUE)
}

rp14_label_compact <- function(x, digits=1) {
  vapply(as.numeric(x), function(a){
    if(!is.finite(a))return('—')
    if(abs(a)>=1000000)return(paste0(rp_num(a/1e6,digits),'m'))
    if(abs(a)>=10000)return(paste0(rp_num(a/1e3,digits),'k'))
    rp_num(a,digits,group=TRUE)
  }, character(1))
}

rp14_est_ci_label <- function(est, lo, hi, digits=2L, suffix='', separator=' to ') {
  paste0(rp_num(est,digits), ' (', rp_num(lo,digits), separator, rp_num(hi,digits), ')', suffix)
}

rp14_forest_text_x <- function(limits, gap_mult=.06, extra_mult=.55) {
  span <- max(diff(limits), 0.05)
  list(x_text = limits[2] + span * gap_mult,
       x_max  = limits[2] + span * (gap_mult + extra_mult))
}

rp14_export <- function(p,filename,dir,s) {
  # 只输出两种真正的矢量格式；PDF和SVG均保留文字/线条的矢量属性。
  # SVG需要svglite包，PDF优先Cairo设备；此处不创建任何位图文件。
  rp_dir(dir)
  allowed <- c('pdf','svg')
  formats <- unique(tolower(trimws(as.character(s$figure$formats))))
  if(!length(formats) || any(!formats %in% allowed))
    stop('14仅支持PDF和SVG矢量输出；请设置figure$formats=c("pdf","svg")')
  for(ext in formats)rp_guard(file.path(dir,paste0(filename,'.',ext)))
  if('pdf'%in%formats){
    device <- if(capabilities('cairo'))grDevices::cairo_pdf else 'pdf'
    ggplot2::ggsave(file.path(dir,paste0(filename,'.pdf')),p,device=device,
      width=s$figure$width,height=s$figure$height,units='in',limitsize=FALSE,
      bg=s$figure$background_color)
  }
  if('svg'%in%formats){
    if(!requireNamespace('svglite',quietly=TRUE))
      stop('SVG矢量图导出需要安装svglite包：install.packages("svglite")')
    ggplot2::ggsave(file.path(dir,paste0(filename,'.svg')),p,
      device=svglite::svglite,width=s$figure$width,height=s$figure$height,
      units='in',limitsize=FALSE,bg=s$figure$background_color)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# FIGURE 1 | A=score-based reclassification; B/C=observed outcomes;
#            D=Module 10 LDL-C intervention eligibility (no model refit).
# All visual parameters live in the MAIN script, figure_layout$fig1.
# -----------------------------------------------------------------------------
rp14_fig1_apply_title <- function(p, layout, panel) {
  p + ggplot2::labs(title=paste0(panel$letter,layout$common$panel_prefix_separator,
                                  panel$title)) +
    ggplot2::theme(
      plot.title.position='plot',
      plot.title=ggplot2::element_text(
        face='bold',size=panel$title_size_pt,lineheight=panel$title_lineheight,
        hjust=0,margin=ggplot2::margin(b=panel$title_gap_pt)))
}

rp14_fig1_matrix <- function(src,s,layout) {
  f1<-layout$fig1
  o<-f1$panel_A
  d<-data.table::copy(data.table::as.data.table(src$reclassification))
  rp_assert(d,c('score_2019_category','score2_2025_category','n','row_proportion'))
  levels<-c('L','M','H','VH')
  grid<-data.table::CJ(score_2019_category=levels,score2_2025_category=levels)
  d<-merge(grid,d,by=c('score_2019_category','score2_2025_category'),all.x=TRUE,sort=FALSE)
  if(nrow(d)!=16L || anyNA(d$n) || any(d$n<0))stop('Figure 1 的 4×4 人数表缺失或包含负数')
  if(data.table::uniqueN(d,by=c('score_2019_category','score2_2025_category'))!=16L)
    stop('Figure 1 reclassification 单元格重复')
  N<-sum(d$n)
  if(!is.finite(N) || N<=0)stop('重分类样本数无效')
  d[,row_total:=sum(n),by=score_2019_category]
  d[,row_pct_check:=ifelse(row_total>0,n/row_total,NA_real_)]
  if(any(abs(d$row_pct_check-d$row_proportion)>.005,na.rm=TRUE))
    stop('重新分类 row_proportion 与人数计算不一致')
  d[,`:=`(x=match(score2_2025_category,levels)+1L,
          y=5L-match(score_2019_category,levels),
          fill_pct=100*row_proportion,
          label=paste0(rp_n(n),'\n(',rp_num(100*row_proportion,o$matrix_proportion_digits),'%)'))]
  d[,text_color:=ifelse(fill_pct>=o$matrix_dark_text_threshold,
    layout$design$colors$white,layout$design$colors$text)]
  # Individual cells use 2019 row denominator; total rows/cols use whole-population N.
  rows<-d[,.(n=sum(n)),by=score_2019_category]
  rows[,`:=`(x=6L,y=5L-match(score_2019_category,levels),
    label=paste0(rp_n(n),'\n(',rp_num(100*n/N,o$matrix_proportion_digits),'%)'))]
  cols<-d[,.(n=sum(n)),by=score2_2025_category]
  cols[,`:=`(x=match(score2_2025_category,levels)+1L,y=0L,
    label=paste0(rp_n(n),'\n(',rp_num(100*n/N,o$matrix_proportion_digits),'%)'))]
  all<-data.table::data.table(x=6,y=0,n=N,
    label=paste0(rp_n(N),'\n(',rp_num(100,o$matrix_proportion_digits),'%)'))
  total_cells<-data.table::rbindlist(list(rows[,.(x,y,label)],cols[,.(x,y,label)],all[,.(x,y,label)]))

  column_heads<-data.table::data.table(x=2:6,y=5,label=unname(o$matrix_column_labels))
  # 旧版左上角仅保留无文字的浅色表头单元格；2019标签移至左侧跨四行合并单元格。
  corner_head<-data.table::data.table(x=1,y=5)
  row_heads<-data.table::data.table(x=1,y=c(4,3,2,1,0),label=unname(o$matrix_row_labels))
  p<-ggplot2::ggplot()+
    # 2025横向合并表头：范围恰好对应四个风险类别列，不覆盖Total。
    ggplot2::annotate('rect',xmin=o$matrix_superheader_x[1],xmax=o$matrix_superheader_x[2],
      ymin=o$matrix_superheader_y[1],ymax=o$matrix_superheader_y[2],
      fill=f1$matrix_header_fill,color=o$matrix_merged_header_border_color,
      linewidth=o$matrix_merged_header_border_width)+
    # 2019纵向合并表头：范围恰好对应四个风险类别行，不覆盖Total。
    ggplot2::annotate('rect',xmin=o$matrix_sideheader_x[1],xmax=o$matrix_sideheader_x[2],
      ymin=o$matrix_sideheader_y[1],ymax=o$matrix_sideheader_y[2],
      fill=f1$matrix_header_fill,color=o$matrix_merged_header_border_color,
      linewidth=o$matrix_merged_header_border_width)+
    ggplot2::geom_tile(data=column_heads,ggplot2::aes(x=x,y=y),
      width=o$matrix_tile_width,height=o$matrix_tile_height,
      fill=f1$matrix_header_fill,color=layout$design$colors$white,
      linewidth=o$matrix_header_border_width)+
    ggplot2::geom_tile(data=row_heads,ggplot2::aes(x=x,y=y),
      width=o$matrix_tile_width,height=o$matrix_tile_height,
      fill=f1$matrix_header_fill,color=layout$design$colors$white,
      linewidth=o$matrix_header_border_width)+
    ggplot2::geom_tile(data=corner_head,ggplot2::aes(x=x,y=y),
      width=o$matrix_tile_width,height=o$matrix_tile_height,
      fill=f1$matrix_header_fill,color=layout$design$colors$white,
      linewidth=o$matrix_header_border_width)+
    ggplot2::geom_tile(data=total_cells,ggplot2::aes(x=x,y=y),
      width=o$matrix_tile_width,height=o$matrix_tile_height,
      fill=f1$matrix_total_fill,color=f1$matrix_total_border,
      linewidth=o$matrix_total_border_width)+
    ggplot2::geom_tile(data=d,ggplot2::aes(x=x,y=y,fill=fill_pct),
      width=o$matrix_tile_width,height=o$matrix_tile_height,
      color=layout$design$colors$white,linewidth=o$matrix_cell_border_width)+
    ggplot2::geom_text(data=d,ggplot2::aes(x=x,y=y,label=label,color=text_color),
      size=o$matrix_text_size_mm,lineheight=o$matrix_cell_lineheight,
      show.legend=FALSE,family=s$figure$font_family)+
    ggplot2::geom_text(data=column_heads,ggplot2::aes(x=x,y=y,label=label),
      size=o$matrix_header_size_mm,lineheight=o$matrix_header_lineheight,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::geom_text(data=row_heads,ggplot2::aes(x=x,y=y,label=label),
      size=o$matrix_header_size_mm,lineheight=o$matrix_header_lineheight,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::annotate('text',
      x=mean(o$matrix_sideheader_x),y=mean(o$matrix_sideheader_y),
      label=o$matrix_side_header_label,size=o$matrix_side_header_size_mm,
      lineheight=o$matrix_side_header_lineheight,angle=o$matrix_sideheader_angle,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::geom_text(data=total_cells,ggplot2::aes(x=x,y=y,label=label),
      size=o$matrix_text_size_mm,lineheight=o$matrix_cell_lineheight,
      family=s$figure$font_family)+
    ggplot2::annotate('text',
      x=mean(o$matrix_superheader_x),y=mean(o$matrix_superheader_y),
      label=o$matrix_column_header,size=o$matrix_top_header_size_mm,
      lineheight=o$matrix_superheader_lineheight,angle=o$matrix_superheader_angle,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::scale_fill_gradientn(colors=f1$matrix_colors,limits=c(0,o$matrix_fill_max),
      oob=scales::squish,guide='none')+
    ggplot2::scale_color_identity()+
    ggplot2::coord_cartesian(xlim=o$matrix_x_limits,ylim=o$matrix_y_limits,
      clip='off',expand=FALSE)+
    ggplot2::theme_void(base_family=s$figure$font_family)+
    ggplot2::theme(legend.position='none',
      plot.margin=do.call(ggplot2::margin,as.list(o$panel_margin_pt)))
  p<-rp14_fig1_apply_title(p,layout,o)
  list(plot=p,body=d,rows=rows,cols=cols,overall=all)
}

rp14_fig1_direction <- function(src,s,layout) {
  f1 <- layout$fig1
  o <- f1$panel_A
  a <- data.table::copy(data.table::as.data.table(src$reclassification_direction))
  rp_assert(a,c('direction','n','proportion'))
  need <- c('lower','unchanged','higher')
  if(anyDuplicated(a$direction) || !setequal(a$direction,need))
    stop('reclassification_direction 必须且只能含 lower/unchanged/higher 三行')
  a <- a[match(need,direction)]
  if(!isTRUE(o$show_direction_counts))
    stop('Figure 1 A 必须显示 lower/unchanged/higher 三项统计')
  if(abs(sum(a$proportion)-1)>.005)stop('重分类方向百分比合计不为 100%')
  a[,`:=`(x=o$summary_col_centers,
    title=unname(o$summary_direction_labels[need]),
    value=paste0(rp_n(n),' (',rp_num(100*proportion,o$summary_proportion_digits),'%)'))]
  max_v <- o$matrix_fill_max
  steps <- as.integer(o$summary_gradient_steps)
  if(!is.finite(steps)||steps<10L)stop('summary_gradient_steps 至少为 10')
  xlo<-o$summary_gradient_x[1]; xhi<-o$summary_gradient_x[2]
  gradient <- data.table::data.table(
    xmin=seq(xlo,xhi,length.out=steps+1L)[-(steps+1L)],
    xmax=seq(xlo,xhi,length.out=steps+1L)[-1L],
    ymin=o$summary_gradient_y[1],ymax=o$summary_gradient_y[2],
    value=seq(0,max_v,length.out=steps))
  ticks <- data.table::data.table(value=seq(0,max_v,length.out=5L))
  ticks[,x:=xlo+(xhi-xlo)*value/max_v]
  p <- ggplot2::ggplot()+
    ggplot2::geom_rect(data=gradient,
      ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=value),
      inherit.aes=FALSE,show.legend=FALSE)+
    ggplot2::scale_fill_gradientn(colors=f1$matrix_colors,
      limits=c(0,max_v),oob=scales::squish,guide='none')+
    ggplot2::annotate('text',x=o$summary_legend_title_x,y=o$summary_legend_title_y,
      label=o$summary_color_title,size=o$summary_text_size_mm,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::geom_text(data=ticks,ggplot2::aes(x=x,label=rp_num(value,0)),
      y=o$summary_legend_y,size=o$summary_legend_tick_size_mm,
      family=s$figure$font_family)+
    ggplot2::geom_text(data=a,ggplot2::aes(x=x,label=title),
      y=o$summary_title_y,size=o$summary_text_size_mm,
      family=s$figure$font_family)+
    ggplot2::geom_text(data=a,ggplot2::aes(x=x,label=value),
      y=o$summary_value_y,size=o$summary_value_size_mm,
      fontface='bold',family=s$figure$font_family)+
    ggplot2::coord_cartesian(xlim=o$summary_x_limits,ylim=o$summary_y_limits,
      expand=FALSE,clip='off')+
    ggplot2::theme_void(base_family=s$figure$font_family)+
    ggplot2::theme(plot.margin=do.call(ggplot2::margin,as.list(o$summary_margin_pt)))
  if(isTRUE(o$summary_divider)) p <- p +
    ggplot2::geom_segment(data=data.table::data.table(x=o$summary_divider_x),
      ggplot2::aes(x=x,xend=x),
      y=o$summary_divider_y[1],yend=o$summary_divider_y[2],
      inherit.aes=FALSE,linewidth=o$summary_divider_width,
      color=f1$summary_divider_color)
  list(plot=p,data=a,gradient=gradient)
}

# Common geometry ONLY: all independent Panel B/C/D labels, sizes and margins are
# taken from figure_layout$fig1$panel_B/$panel_C/$panel_D in the main script.
rp14_fig1_bar <- function(d,s,layout,which) {
  f1 <- layout$fig1
  o <- f1[[paste0('panel_',which)]]
  if(!isTRUE(o$y_from_zero))stop('Figure 1 ',which,' 柱图的 Y 轴必须从 0 开始')
  d <- data.table::copy(data.table::as.data.table(d))
  if(nrow(d)!=2L || anyDuplicated(d$guideline) ||
     !setequal(as.character(d$guideline),c('2019','2025')))
    stop('Figure 1 ',which,' 必须恰好包含 2019/2025 两个估计')
  if(any(!is.finite(d$y) | d$y<0 | d$y>100))
    stop('Figure 1 ',which,' 百分比必须为有限且位于 0–100 的数值')
  if(isTRUE(o$show_ci)){
    if(!all(c('lo','hi')%in%names(d)) ||
       any(!is.finite(d$lo) | !is.finite(d$hi)) ||
       any(d$lo>d$y | d$y>d$hi))
      stop('Figure 1 ',which,' 的 95% CI 不完整')
  } else {
    d[,`:=`(lo=NA_real_,hi=NA_real_)]
  }
  upper<-if(isTRUE(o$show_ci))d$hi else d$y
  base_y<-max(max(upper),1)
  d[,label_y:=upper+base_y*o$label_offset_fraction]
  if(identical(o$label_mode,'ci')){
    d[,bar_label:=paste0(rp_num(y,o$label_digits),'%',
                        '\n(',rp_num(lo,o$label_digits),s$ci_interval_separator,rp_num(hi,o$label_digits),')')]
  } else if(identical(o$label_mode,'n')){
    rp_assert(d,'n','Figure 1 D 人数')
    # Figure 1 D：主脚本show_count控制是否显示人数，比例始终显示。
    # TRUE: 226,855 (86.91%)；FALSE: 86.91%。数据中的n仍保留供QC与CSV导出。
    if(isTRUE(o$show_count)){
      d[,bar_label:=paste0(rp_num(n,o$n_digits,group=TRUE),' (',
                          rp_num(y,o$label_digits),'%)')]
    } else {
      d[,bar_label:=paste0(rp_num(y,o$label_digits),'%')]
    }
  } else stop('Figure 1 ',which,' 柱上标签模式必须为 ci 或 n')
  needed_top<-max(d$label_y)+base_y*o$top_padding_fraction
  lim<-o$y_axis_limits
  if(!is.numeric(lim) || length(lim)!=2L || !is.finite(lim[1]) || lim[1]!=0 ||
     is.nan(lim[2]) || (!is.na(lim[2]) && (!is.finite(lim[2]) || lim[2]<=0)))
    stop('Figure 1 ',which,' y_axis_limits 应为 c(0, NA_real_) 或 c(0, 正数)')
  y_top<-if(is.na(lim[2]))needed_top else lim[2]
  if(!is.finite(y_top) || y_top<needed_top)
    stop('Figure 1 ',which,' Y 轴上界不足以容纳柱顶文字；请增大 y_axis_limits 或使用自动模式')
  y_breaks<-if(is.null(o$y_axis_breaks))pretty(c(0,y_top),n=o$y_tick_count) else o$y_axis_breaks
  if(!is.numeric(y_breaks) || any(!is.finite(y_breaks)) || any(diff(y_breaks)<=0))
    stop('Figure 1 ',which,' Y 轴刻度无效')
  d[,guideline:=factor(as.character(guideline),levels=c('2019','2025'))]
  p <- ggplot2::ggplot(d,ggplot2::aes(x=guideline,y=y,fill=guideline))+
    ggplot2::geom_col(width=o$bar_width,color=f1$bar_border_color,
      linewidth=o$bar_border_width,show.legend=FALSE)
  if(isTRUE(o$show_ci)) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),
      width=o$errorbar_width,linewidth=o$errorbar_linewidth,
      color=f1$bar_border_color)
  }
  p <- p +
    ggplot2::geom_text(ggplot2::aes(y=label_y,label=bar_label),
      size=o$label_size_mm,lineheight=o$label_lineheight,
      color=layout$design$colors$text,vjust=0,
      family=s$figure$font_family)+
    ggplot2::scale_fill_manual(values=f1$guideline_colors,guide='none')+
    ggplot2::scale_x_discrete(labels=unname(o$x_labels[c('2019','2025')]),
      expand=ggplot2::expansion(add=o$x_expand))+
    ggplot2::scale_y_continuous(breaks=y_breaks,
      expand=ggplot2::expansion(mult=c(0,0)))+
    ggplot2::coord_cartesian(ylim=c(0,y_top),clip='off')+
    ggplot2::labs(x=o$x_title,y=o$y_title)+
    rp14_theme(s,layout,'none')+
    ggplot2::theme(
      plot.margin=do.call(ggplot2::margin,as.list(o$panel_margin_pt)),
      axis.title.y=ggplot2::element_text(size=o$axis_title_size_pt,
        margin=ggplot2::margin(r=o$axis_title_gap_pt)),
      axis.title.x=ggplot2::element_text(size=o$axis_title_size_pt,
        margin=ggplot2::margin(t=o$axis_title_gap_pt)),
      axis.text.x=ggplot2::element_text(size=o$axis_tick_size_pt,
        margin=ggplot2::margin(t=o$x_tick_gap_pt)),
      axis.text.y=ggplot2::element_text(size=o$axis_tick_size_pt))
  list(plot=rp14_fig1_apply_title(p,layout,o),data=d)
}

rp14_fig1_hvh <- function(src,s,layout,metric,which) {
  d <- data.table::copy(data.table::as.data.table(src$hvh_observed))
  d <- d[scenario=='S0' & time==s$primary_horizon &
         scope%in%c('2019__overall','2025__overall') & estimand==metric]
  if(nrow(d)!=2L || anyDuplicated(d$scope))
    stop('Figure 1 ',which,' 需要 S0 2019/2025 各一行，且不得重复')
  rp14_finite_ci(d,paste('Figure 1',which))
  d <- rp14_to_pct(d)
  d[,guideline:=ifelse(grepl('^2019',scope),'2019','2025')]
  rp14_fig1_bar(d,s,layout,which)
}

rp14_fig1_eligibility <- function(src,s,layout,cohort_n) {
  # Pooled descriptive output from Module 10, not the per-imputation file.
  # 10 emits MI_mean (count), MI_mean_proportion (count / full cohort N), MI_n.
  # No 95%CI in this source. Do not create a CI using binomial assumptions.
  d <- data.table::copy(data.table::as.data.table(src$fig1_eligibility))
  rp_assert(d,c('scope','guideline','scenario','metric','MI_mean','MI_mean_proportion','MI_n'),
            'Figure 1 D Module 10 intervention requirements')
  d <- d[scope=='overall' & scenario=='S0' & metric=='n_intervention' &
         guideline%in%c('2019','2025')]
  if(nrow(d)!=2L || anyDuplicated(d$guideline) ||
     !setequal(as.character(d$guideline),c('2019','2025')))
    stop('Figure 1 D 必须有 S0 / Overall / n_intervention 的 2019 和 2025 两行')
  if(any(!is.finite(d$MI_mean) | !is.finite(d$MI_mean_proportion) |
         !is.finite(d$MI_n) | d$MI_n<1 | d$MI_mean<0 |
         d$MI_mean_proportion<0 | d$MI_mean_proportion>1))
    stop('Figure 1 D pooled intervention proportion/count invalid')
  if(length(unique(d$MI_n))!=1L) stop('Figure 1 D 2019/2025 使用的 MI 数量不同')
  if(!is.finite(cohort_n) || cohort_n<=0 ||
     any(d$MI_mean>cohort_n) ||
     any(abs(d$MI_mean_proportion-d$MI_mean/cohort_n)>0.005))
    stop('Figure 1 D pooled 干预比例与 Figure 1 A 人群分母不一致')
  d[,`:=`(n=MI_mean,y=100*MI_mean_proportion,
          guideline=as.character(guideline))]
  rp14_fig1_bar(d,s,layout,'D')
}

rp14_fig1 <- function(src,s,layout) {
  a<-rp14_fig1_matrix(src,s,layout)
  dir<-rp14_fig1_direction(src,s,layout)
  b<-rp14_fig1_hvh(src,s,layout,'observed_HVH_risk','B')
  c<-rp14_fig1_hvh(src,s,layout,'observed_HVH_capture','C')
  d<-rp14_fig1_eligibility(src,s,layout,a$overall$n[[1L]])
  # AAA: matrix / BBB: summary strip / CCC: vertical gap / DEF: bottom three charts.
  o <- layout$fig1$layout
  r <- o$summary_height_ratio
  z <- o$matrix_to_bottom_height
  if(length(r)!=2L || any(!is.finite(r)) || any(r<=0) ||
     length(z)!=2L || any(!is.finite(z)) || any(z<=0))
    stop('Figure 1 画幅高度与内部高度比例必须均为两个正数')
  if(!is.finite(o$panel_gap_ratio) || o$panel_gap_ratio<=0 || o$panel_gap_ratio>1)
    stop('Figure 1 panel_gap_ratio 必须为 (0,1]')
  h_matrix<-z[1]*r[1]/sum(r)
  h_summary<-z[1]*r[2]/sum(r)
  h_gap<-z[2]*o$panel_gap_ratio
  p <- patchwork::wrap_plots(
    a$plot,dir$plot,patchwork::plot_spacer(),b$plot,c$plot,d$plot,
    design='AAA\nBBB\nCCC\nDEF',
    widths=o$bottom_panel_widths,
    heights=c(h_matrix,h_summary,h_gap,z[2])
  )
  list(plot=p,plot_data=list(
    reclassification=a$body,
    reclassification_rows=a$rows,
    reclassification_columns=a$cols,
    reclassification_total=a$overall,
    directions=dir$data,summary_gradient=dir$gradient,
    hvh_risk=b$data,hvh_capture=c$data,
    intervention_eligibility=d$data))
}

# -----------------------------------------------------------------------------
# FIGURE 2 | 六种治疗转换（A/B/C/D 均使用完整的六组结果）
# B 使用 mean_D_2019/mean_D_2025：人均所需 LDL-C 降幅；已保存描述统计，无新 CI。
# -----------------------------------------------------------------------------
rp14_transition_order <- function(s) {
  tr<-s$transition_order
  if(length(tr)!=6L||anyDuplicated(tr))stop('转换组参数应恰好含六个互不重复的组')
  tr
}
rp14_transition_yticks <- function(s) {
  t<-rev(rp14_transition_order(s))
  list(breaks=seq_along(t),labels=unname(s$labels$transition[t]))
}
rp14_transition_yaxis <- function(s) {
  tick<-rp14_transition_yticks(s)
  ggplot2::scale_y_continuous(breaks=tick$breaks, labels=tick$labels,
    limits=c(.52,6.48), expand=ggplot2::expansion(mult=0))
}
rp14_transitions_source <- function(src,s) {
  d<-data.table::copy(data.table::as.data.table(src$transition_counts))
  rp_assert(d,c('scope','scenario','transition','n','proportion','T_2019','T_2025','mean_D_2019','mean_D_2025'))
  d<-d[scope=='overall' & scenario=='S0']
  tr<-rp14_transition_order(s)
  if(nrow(d)!=6L || !setequal(d$transition,tr) || anyDuplicated(d$transition))
    stop('Figure 2 的 6 组转换人数/比例不完整或重复')
  d[,`:=`(y_id=match(transition,rev(tr)),percent=100*proportion)]
  d
}


# 旧版 Figure 2 forest/text-column 辅助函数已移除；以以下 v5 示例版为唯一绘图入口。

# -----------------------------------------------------------------------------
# Read-only Module 12 absolute contribution interface for Figure 2D/3D.
# Six exhaustive group population-weighted RDs sum to the Overall paired RD.
# Absolute pp values are signed; unlike the old percentage shares there is no denominator instability.
# -----------------------------------------------------------------------------
rp14_transition_contribution <- function(source, scenario_id, horizon, six_groups,
                                          overall_effect, digits=3L, label='Figure D',
                                          tolerance=2e-5) {
  # Source: treatment_transitions/pooled/risk_difference_contributions.csv.
  # estimate/lcl/ucl are risk-unit pooled values (not percentages).
  # Do NOT use contribution_percentages.csv and do NOT re-estimate weights/RDs/CI.
  x <- data.table::copy(data.table::as.data.table(source))
  rp_assert(x,c('scenario','scope','transition','time','estimand','estimate','lcl','ucl'),
            paste(label,'Module 12 population-weighted RD'))
  overall <- data.table::as.data.table(overall_effect)
  rp_assert(overall,c('estimate','lcl','ucl'),paste(label,'paired Overall RD'))
  if(nrow(overall)!=1L || !all(is.finite(unlist(overall[1L,.(estimate,lcl,ucl)]))))
    stop(label,': expected one finite Overall 2025-2019 paired RD and CI')
  x <- x[scenario==scenario_id & scope=='overall' & time==horizon &
         estimand=='absolute_RD_contribution']
  if(nrow(x)!=length(six_groups) || anyDuplicated(x,by='transition') ||
     !setequal(as.character(x$transition),as.character(six_groups)))
    stop(label,': expected all six mutually exclusive absolute RD contributions for ',scenario_id)
  rp14_finite_ci(x,paste(label,'population-weighted RD'))
  if(any(!is.finite(x$estimate) | !is.finite(x$lcl) | !is.finite(x$ucl)))
    stop(label,': contribution estimate or CI is missing/non-finite')
  if(!is.numeric(tolerance) || length(tolerance)!=1L || !is.finite(tolerance) || tolerance<=0)
    stop(label,': positive contribution additivity tolerance required')
  if(abs(sum(x$estimate)-overall$estimate[[1L]])>tolerance)
    stop(label,': pooled population-weighted RD does not add up to Overall paired RD: ',
         'sum=',signif(sum(x$estimate),7),', overall=',signif(overall$estimate[[1L]],7))
  # The Overall row represents the full-population paired RD, with its own pooled CI.
  total <- data.table::data.table(
    scenario=scenario_id,scope='overall',transition='overall',time=horizon,
    estimand='rd_2025_vs_2019',estimate=overall$estimate[[1L]],
    lcl=overall$lcl[[1L]],ucl=overall$ucl[[1L]])
  z <- data.table::rbindlist(list(x,total),use.names=TRUE,fill=TRUE)
  z[,`:=`(contribution_pp=100*estimate,
          contribution_lcl_pp=100*lcl,
          contribution_ucl_pp=100*ucl,
          contribution_status='estimable')]
  z[,contribution_label:=rp_num(contribution_pp,digits)]
  if(anyNA(z$contribution_label))stop(label,': display conversion failed')
  z
}

# ============================================================================
# Figure 2: exact chart *types* and visual semantics from Figure_02(10).pdf.
# All values read from already pooled CSV; no outcome or MI calculations.
# ============================================================================
rp14_fig2 <- function(src,s,layout) {
  o <- layout$fig2
  tr <- s$transition_order
  cd <- c('overall',setdiff(tr,'unchanged_untreated'))
  cnt <- data.table::copy(src$transition_counts[scope=='overall' & scenario=='S0'])
  rp_check_unique(cnt,'transition','Figure 2 treatment counts')
  if(nrow(cnt)!=6L || !setequal(cnt$transition,tr)) stop('Figure 2: six transition groups missing/duplicated')
  cnt[,transition:=factor(as.character(transition),levels=rev(tr))]
  cnt[,percent:=100*proportion]
  pal <- rp14_guideline_palette(s)
  point_shape <- o$point_shape
  gg_theme <- function(panel_id){
    # 上方所有样式沿用原版；横轴标题与刻度各自有独立可调间距。
    rp14_theme(s,layout,'none') + ggplot2::theme(
      axis.text=ggplot2::element_text(size=o$axis_text_size),
      axis.text.x=ggplot2::element_text(size=o$axis_text_size,
        margin=ggplot2::margin(t=unname(o$x_tick_gap_pt[[panel_id]]))),
      axis.title.x=ggplot2::element_text(size=o$axis_title_size,
        margin=ggplot2::margin(t=unname(o$x_title_gap_pt[[panel_id]]))),
      plot.title=ggplot2::element_text(size=o$panel_title_size,face='bold',lineheight=1.06,margin=ggplot2::margin(b=o$panel_title_gap_pt)),
      panel.grid.major.y=ggplot2::element_blank(),
      panel.grid.major.x=ggplot2::element_blank())
  }
  row_guides <- function(groups){
    ggplot2::geom_hline(yintercept=seq_along(groups),color=o$row_guide_color,
      linewidth=o$row_guide_width)
  }
  f2_axis <- function(vals, limits, breaks, label, step=NULL){
    val <- as.numeric(vals[is.finite(vals)])
    if(!length(val))stop('Figure 2 ',label,': empty axis source')
    if(length(limits)!=2L || any(is.nan(limits)))stop('Figure 2 ',label,': axis limits invalid')
    if(all(is.na(limits))){
      if(length(step)!=1L || !is.finite(step) || step<=0)
        stop('Figure 2 ',label,': both bounds automatic requires positive step')
      lo<-floor(min(c(0,val))/step)*step
      hi<-ceiling(max(c(0,val))/step)*step
      if(hi<=lo)hi<-lo+step
      return(list(limits=c(lo,hi),breaks=seq(lo,hi,by=step)))
    }
    if(any(!is.finite(limits)) || limits[1]>=limits[2])
      stop('Figure 2 ',label,': fixed limits must be finite increasing values')
    lim <- c(min(limits[1],min(val)),max(limits[2],max(val)))
    if(any(lim!=limits)) warning('Figure 2 ',label,': data exceed fixed exemplar bounds; expanded to preserve CI')
    list(limits=lim,breaks=if(any(lim!=limits))pretty(lim,n=5) else breaks)
  }
  labs_tr <- unname(s$labels$transition[rev(tr)])
  labs_cd <- c(overall=unname(s$labels$scope[['overall']]),s$labels$transition)

  # A: six discrete colors and horizontal bars, 1-digit labels.
  pA <- ggplot2::ggplot(cnt,ggplot2::aes(x=percent,y=transition,fill=transition))+
    ggplot2::geom_col(width=o$bar_width,color=o$bar_border_color,linewidth=.32)+
    ggplot2::geom_text(ggplot2::aes(label=paste0(rp_num(percent,o$count_digits),'%')),
      hjust=-.12,size=o$value_label_size,family=s$figure$font_family)+
    ggplot2::scale_fill_manual(values=o$group_colors,guide='none')+
    ggplot2::scale_y_discrete(limits=rev(tr),labels=labs_tr,drop=FALSE)+
    ggplot2::scale_x_continuous(
      limits=c(o$group_axis_limits[1],max(o$group_axis_limits[2],max(cnt$percent)*1.13)),
      breaks=o$group_axis_breaks,expand=ggplot2::expansion(mult=c(0,0)))+
    ggplot2::labs(title=paste0(o$panel_letters[['A']],layout$common$panel_prefix_separator,o$labels$A),x=o$x_titles$A,y=NULL)+gg_theme('A')

  # Panel C: observed risk, horizontal bars with horizontal 95% CI.
  obs <- data.table::copy(src$transition_observed[
    scope=='overall' & scenario=='S0' & time==s$primary_horizon & estimand=='observed_transition_risk'])
  if(nrow(obs)!=6L || anyDuplicated(obs$transition))stop('Figure 2 C: observed risks incomplete')
  rp14_finite_ci(obs,'Figure 2 C observed risks')
  obs[,transition:=factor(as.character(transition),levels=rev(tr))]
  obs[,`:=`(pct=100*estimate,lwr=100*lcl,upr=100*ucl)]
  xb <- f2_axis(c(obs$pct,obs$lwr,obs$upr),o$risk_axis_limits,o$risk_axis_breaks,'C')
  pObserved <- ggplot2::ggplot(obs,ggplot2::aes(x=pct,y=transition))+
    (if(isTRUE(o$panelC_show_row_guides)) row_guides(tr) else NULL)+
    ggplot2::geom_col(fill=o$observed_bar_color,color=o$observed_bar_border_color,linewidth=.25,width=o$observed_bar_width)+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=lwr,xmax=upr),orientation='y',
      width=o$ci_cap_width,linewidth=o$ci_linewidth,color=o$ci_color)+
    ggplot2::scale_y_discrete(limits=rev(tr),labels=labs_tr,drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=xb$breaks,expand=ggplot2::expansion(mult=c(0,0)))+
    ggplot2::coord_cartesian(xlim=xb$limits,clip='off')+
    ggplot2::labs(title=paste0(o$panel_letters[['C']],layout$common$panel_prefix_separator,o$labels$C),x=o$x_titles$C,y=NULL)+gg_theme('C')+
    ggplot2::theme(panel.grid=ggplot2::element_blank())

  # Panel B: required LDL-C reduction, open circles + paired connector.
  pop <- data.table::copy(src$population_mean[
    scope=='overall' & scenario=='S0' & guideline%in%c('2019','2025')])
  dd <- data.table::melt(cnt[,.(transition=as.character(transition),mean_D_2019,mean_D_2025)],
    id.vars='transition',variable.name='guideline',value.name='mean_D')
  dd[,guideline:=sub('mean_D_','',guideline,fixed=TRUE)]
  dd <- data.table::rbindlist(list(dd,pop[,.(transition='overall',guideline=as.character(guideline),mean_D)]))
  dd <- dd[transition%in%cd]
  if(nrow(dd)!=12L || anyDuplicated(dd,by=c('transition','guideline')) || any(!is.finite(dd$mean_D)))
    stop('Figure 2 B: expected 6 groups by 2 guidelines of finite mean_D')
  dd[,transition:=factor(transition,levels=rev(cd))]
  pairs <- data.table::dcast(dd,transition~guideline,value.var='mean_D')
  xc <- f2_axis(dd$mean_D,o$reduction_axis_limits,o$reduction_axis_breaks,'B')
  pIntensity <- ggplot2::ggplot(dd,ggplot2::aes(x=mean_D,y=transition,color=guideline))+
    row_guides(cd)+
    ggplot2::geom_segment(data=pairs,ggplot2::aes(x=`2019`,xend=`2025`,y=transition,yend=transition),
      inherit.aes=FALSE,color=o$connector_color,linewidth=o$connector_width)+
    ggplot2::geom_point(shape=point_shape,fill=layout$design$colors$white,size=o$point_size,stroke=o$point_stroke)+
    ggplot2::scale_color_manual(values=pal,breaks=c('2019','2025'),labels=unname(o$guideline_legend_labels[c('2019','2025')]))+
    ggplot2::scale_y_discrete(limits=rev(cd),labels=unname(labs_cd[rev(cd)]),drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=xc$breaks,expand=ggplot2::expansion(mult=c(.02,.02)))+
    ggplot2::coord_cartesian(xlim=xc$limits,clip='off')+
    ggplot2::labs(title=paste0(o$panel_letters[['B']],layout$common$panel_prefix_separator,o$labels$B),x=o$x_titles$B,y=NULL)+gg_theme('B')

   # Panel D: the pooled paired RD (2025 minus 2019), exactly one row per group.
   # Right column: population-weighted RD (pp), directly from Module 12 pooled results.
   # It is an absolute contribution (P(group)*paired RD), not a contribution ratio.
   eff <- data.table::copy(src$transition_effects[
     scope=='overall' & scenario=='S0' & time==s$primary_horizon &
     estimand=='rd_2025_vs_2019'])
   tot <- data.table::copy(src$all_effects[
     scope=='overall' & scenario=='S0' & time==s$primary_horizon &
     estimand=='rd_2025_vs_2019'])
   if(nrow(tot)!=1L) stop('Figure 2 D: need one S0 Overall paired RD from all_effects')
   tot[,transition:='overall']
   eff <- data.table::rbindlist(list(eff,tot),fill=TRUE)[transition%in%cd]
   if(nrow(eff)!=6L || anyDuplicated(eff,by='transition'))
     stop('Figure 2 D: need five treatment groups plus Overall, each with paired RD')
   rp14_finite_ci(eff,'Figure 2 D paired 2025-2019 RD')
   contribution <- rp14_transition_contribution(src$transition_contribution,'S0',
     s$primary_horizon,tr,tot,o$d_contribution_digits,'Figure 2 D',
     tolerance=o$d_contribution_additivity_tolerance)
   contribution <- contribution[transition%in%cd]
   eff <- merge(eff,contribution[,.(transition,contribution_pp,contribution_lcl_pp,
     contribution_ucl_pp,contribution_label,contribution_status)],
     by='transition',all.x=TRUE,sort=FALSE)
   if(nrow(eff)!=6L || anyNA(eff$contribution_label) ||
      any(!is.finite(eff$contribution_pp) | !is.finite(eff$contribution_lcl_pp) |
          !is.finite(eff$contribution_ucl_pp)))
     stop('Figure 2 D: contribution estimate or 95% CI missing after merging')
   # 仅格式化 Module 12 已有的加权风险差 CI；不会把组内 RD 的 CI 当作贡献 CI。
    # 上下限统一以“ to ”连接，便于与负数符号区分。
   # 保留 contribution_label 及全精度 contribution_*_pp 字段，兼容 Figure 3 和数据导出。
   label_sep <- if(identical(o$d_contribution_label_layout,'two_line')) '\n' else ' '
   eff[,contribution_ci_label:=paste0(
     rp_num(contribution_pp,o$d_contribution_digits),label_sep,'(',
     rp_num(contribution_lcl_pp,o$d_contribution_ci_digits),s$ci_interval_separator,
     rp_num(contribution_ucl_pp,o$d_contribution_ci_digits),')')]
   eff[,`:=`(transition=factor(transition,levels=rev(cd)),
     pct=100*estimate,lwr=100*lcl,upr=100*ucl)]
   xd <- f2_axis(c(eff$lwr,eff$upr,0),o$effect_axis_limits,o$effect_axis_breaks,'D',step=o$effect_axis_step)
   spanD <- diff(xd$limits)
   xCol <- xd$limits[2]+spanD*o$d_contribution_x_gap_fraction
   pD <- ggplot2::ggplot(eff,ggplot2::aes(x=pct,y=transition))+
     row_guides(cd)+ggplot2::geom_vline(xintercept=0,color=o$zero_line_color,linetype='dashed')+
     ggplot2::geom_point(shape=point_shape,fill=layout$design$colors$white,
       color=o$d_rd_color,size=o$point_size,stroke=o$point_stroke)+
     # Draw paired CI last, especially important when its interval is shorter than the marker.
     ggplot2::geom_errorbar(ggplot2::aes(xmin=lwr,xmax=upr),orientation='y',
       width=o$ci_cap_width,linewidth=o$ci_linewidth,color=o$d_rd_color)+
     ggplot2::geom_text(ggplot2::aes(x=xCol,label=contribution_ci_label),
       hjust=0,size=o$d_contribution_text_size,
       lineheight=o$d_contribution_label_lineheight,
       color=layout$design$colors$text,family=s$figure$font_family)+
     ggplot2::annotate('text',x=xCol,y='overall',vjust=o$d_contribution_header_vjust,
       label=o$d_contribution_header,hjust=0,fontface='bold',
       size=o$d_contribution_header_size,family=s$figure$font_family)+
     ggplot2::scale_y_discrete(limits=rev(cd),labels=unname(labs_cd[rev(cd)]),drop=FALSE)+
     ggplot2::scale_x_continuous(breaks=xd$breaks,
       labels=function(x)rp_num(x,o$effect_axis_tick_digits),
       expand=ggplot2::expansion(mult=c(.02,.02)))+
     ggplot2::coord_cartesian(xlim=xd$limits,clip='off')+
     ggplot2::labs(title=paste0(o$panel_letters[['D']],layout$common$panel_prefix_separator,o$labels$D),
       x=o$x_titles$D,y=NULL)+gg_theme('D')+
     ggplot2::theme(
       plot.title=ggplot2::element_text(size=o$panel_title_size,face='bold',lineheight=1.06,
         margin=ggplot2::margin(b=o$d_panel_title_gap_pt)),
       plot.margin=ggplot2::margin(8,o$d_contribution_right_margin_pt,7,10))

   # Figure 2 中只有 Panel B 区分指南年份：两个空心圆图例在 B 最右侧纵向排列。
   # 不使用 patchwork 收集图例，以免被统一放回整张图片的底部。
   pIntensity <- pIntensity +
     ggplot2::guides(color=ggplot2::guide_legend(
       ncol=o$b_legend_ncol,byrow=TRUE,
       override.aes=list(shape=point_shape,fill=layout$design$colors$white,
         size=o$point_size,stroke=o$point_stroke)))+
     ggplot2::theme(
       legend.position=o$b_legend_position,
       legend.direction='vertical',
       legend.justification='center',
       legend.text=ggplot2::element_text(size=o$b_legend_text_size_pt),
       legend.key.width=grid::unit(o$b_legend_key_width_pt,'pt'),
       legend.key.height=grid::unit(o$b_legend_key_height_pt,'pt'),
       legend.spacing.y=grid::unit(o$b_legend_spacing_y_pt,'pt'),
       legend.box.spacing=grid::unit(o$b_legend_box_spacing_pt,'pt'),
       legend.margin=do.call(ggplot2::margin,as.list(o$b_legend_margin_pt)))
   # A 左上、B 右上（自带右侧图例）、C 左下、D 右下；底部不再额外预留图例行。
   p <- patchwork::wrap_plots(pA,pIntensity,pObserved,pD,ncol=2,heights=c(1,1))
  list(plot=p,plot_data=list(counts=cnt,mean_required_reduction=dd,observed=obs,paired_risk_differences=eff))
}

# ============================================================================
# FIGURE 3 | SN/ST/SNT 三面板（A高危分类、B治疗转换、C RD森林图+贡献列）
# Panel C绘图五组但校验原始18行RD和完整6转换组贡献，以免改变统计意义。
# 所有显示参数、原始 CSV 路径、组序和坐标刻度均位于 14 主脚本 fig3。
# ============================================================================
# =============================================================================
# FIGURE 3｜新版正文图：A高危分类（含RM升级段）、B治疗转换、C风险差及贡献（95%CI）。
# 数据全部来自10和12已保存的pooled CSV；不重新进行模型拟合、插补或Bootstrap。
# 所有视觉常数统一由主脚本 figure_layout$fig3$panel_A/B/C、layout、legend控制。
# =============================================================================

rp14_fig3_prepare <- function(src, s, o) {
  k <- as.character(o$scenarios) # Panel B/C有RM的三种策略
  ka <- as.character(o$classification_scenarios) # Panel A仅展示应用RM后的三种策略
  years <- as.character(o$guidelines) # 2019、2025两版指南
  tr <- as.character(o$transition_order) # 六种互斥干预变化组
  cd <- as.character(o$comparative_order) # Panel C显示：Overall及另外四种转换组（共5组）
  full_cd <- c('overall',setdiff(tr,'unchanged_untreated')) # 原有18行RD校验与贡献核对仍使用完整6组
  horizon <- as.numeric(o$horizon) # 分析时间点（年）
  if(!identical(ka,c('SN','ST','SNT')) || !identical(k,c('SN','ST','SNT')) ||
     !identical(years,c('2019','2025'))) stop('Figure 3策略或指南顺序无效')
  if(length(tr)!=6L || anyDuplicated(tr) || length(cd)!=5L ||
     anyDuplicated(cd) || !identical(cd[1L],'overall') ||
     !setequal(cd,c('overall',setdiff(tr,c('unchanged_untreated','unchanged_intensity')))))
    stop('Figure 3B需要六组完整数据；Panel C只显示Overall及排除Unchanged intensity后的4组')

  # 从Module 12六组汇总结果获取整个研究人群分母，三种RM规则下分母应相同。
  cnt <- data.table::copy(data.table::as.data.table(src$fig3_counts))
  rp_assert(cnt,c('scenario','scope','transition','n','proportion'),'Figure 3B transition counts')
  # CSV读取后统一键字段类型，避免data.table::fsetdiff因character/integer不一致而报错。
  cnt[, `:=`(scenario=as.character(scenario), scope=as.character(scope),
             transition=as.character(transition))]
  cnt <- cnt[scope=='overall' & scenario%in%k & transition%in%tr]
  required <- data.table::CJ(scenario=k,transition=tr,sorted=FALSE)
  if(nrow(cnt)!=18L || anyDuplicated(cnt,by=c('scenario','transition')) ||
     nrow(data.table::fsetdiff(required,cnt[,.(scenario,transition)]))>0L)
    stop('Figure 3B必须有SN/ST/SNT×六种转换组的18行pooled汇总')
  if(any(!is.finite(cnt$n)|cnt$n<0|!is.finite(cnt$proportion)|
         cnt$proportion<0|cnt$proportion>1)) stop('Figure 3B人数或比例无效')
  totals <- cnt[,.(N=sum(n),proportion_sum=sum(proportion)),by=scenario]
  if(any(abs(totals$proportion_sum-1)>1e-5))stop('Figure 3B每种RM策略下六组比例合计必须为100%')
  cohort_N <- mean(totals$N)
  if(!is.finite(cohort_N)||cohort_N<=0 || any(abs(totals$N-cohort_N)>max(0.01,1e-5*cohort_N)))
    stop('Figure 3A/B研究总体分母不一致；请核查Module 12的三种RM策略')
  cnt[,`:=`(percent=100*proportion,
            y=as.numeric(setNames(rev(seq_along(tr)),tr)[transition])+
              as.numeric(o$panel_B$strategy_row_offsets[scenario]))]
  if(anyNA(cnt$y))stop('Figure 3B三种RM的行位置无效')

  # A：Module 10的n_hvh是最终H/VH人数，n_upgrade是原M经RM升至H的人数。
  # 注意n_upgrade自身的MI_mean_proportion分母是原始M人数，不能直接当整个人群比例。
  hvh <- data.table::copy(data.table::as.data.table(src$fig3_hvh))
  rp_assert(hvh,c('scope','guideline','scenario','metric','MI_mean','MI_mean_proportion'),
            'Figure 3A Module 10 risk classification')
  # guideline在CSV中可能被fread自动推断为整数2019/2025；集合比较要求两侧类型一致。
  hvh[, `:=`(scope=as.character(scope), guideline=as.character(guideline),
             scenario=as.character(scenario), metric=as.character(metric))]
  hvh <- hvh[scope=='overall' & guideline%in%years & scenario%in%ka &
             metric%in%c('n_hvh','n_upgrade')]
  need_h <- data.table::CJ(guideline=years,scenario=ka,
                          metric=c('n_hvh','n_upgrade'),sorted=FALSE)
  if(nrow(hvh)!=(length(years)*length(ka)*2L) || anyDuplicated(hvh,by=c('guideline','scenario','metric')) ||
     nrow(data.table::fsetdiff(need_h,hvh[,.(guideline,scenario,metric)]))>0L)
    stop('Figure 3A必须有2019/2025×SN/ST/SNT×n_hvh/n_upgrade共12行')
  if(any(!is.finite(hvh$MI_mean)|hvh$MI_mean<0))stop('Figure 3A分类计数无效')
  final_prop <- hvh[metric=='n_hvh',.(guideline,scenario,recorded_fraction=MI_mean_proportion)]
  hvh <- data.table::dcast(hvh,guideline+scenario~metric,value.var='MI_mean')
  hvh <- merge(hvh,final_prop,by=c('guideline','scenario'),sort=FALSE)
  hvh[, `:=`(original_n=n_hvh-n_upgrade,
             original_pct=100*(n_hvh-n_upgrade)/cohort_N,
             upgraded_pct=100*n_upgrade/cohort_N,
             final_pct=100*n_hvh/cohort_N)]
  if(any(hvh$original_n < -1e-5 | hvh$n_upgrade>hvh$n_hvh+1e-5 |
         !is.finite(hvh$recorded_fraction) |
         abs(hvh$recorded_fraction-hvh$n_hvh/cohort_N)>0.005))
    stop('Figure 3A分类数、总体分母及pooled比例不一致')
  # 三种RM策略下，原始H/VH人数应保持一致；最终H/VH应等于原始H/VH加升级人数。
  if(any(hvh[,max(abs(original_n-original_n[1L])),by=guideline]$V1>0.05))
    stop('Figure 3A：三种RM策略的原始H/VH人数不一致；需排查上游分类')
  hvh[,`:=`(guideline_label=factor(unname(o$panel_A$guideline_labels[guideline]),
                                   levels=unname(o$panel_A$guideline_labels[years])),
             scenario=factor(scenario,levels=ka),
             upgraded_center=original_pct+upgraded_pct/2)]
  data.table::setorder(hvh,guideline_label,scenario)

  # C：不同RM规则下的指南间配对RD和贡献均直接读取Module 12的pooled估计及95%CI。
  eff <- data.table::copy(data.table::as.data.table(src$fig3_effects))
  tot <- data.table::copy(data.table::as.data.table(src$fig3_overall_effects))
  rp_assert(eff,c('scenario','scope','transition','time','estimand','estimate','lcl','ucl'),
            'Figure 3C treatment effects')
  rp_assert(tot,c('scenario','scope','time','estimand','estimate','lcl','ucl'),
            'Figure 3C overall effects')
  # 统一Panel C的分组键类型，避免不同模块CSV的自动类型推断影响合并与校验。
  eff[, `:=`(scenario=as.character(scenario),scope=as.character(scope),
             transition=as.character(transition),estimand=as.character(estimand))]
  tot[, `:=`(scenario=as.character(scenario),scope=as.character(scope),
             estimand=as.character(estimand))]
  # 先保持与旧版相同的5种非Unchanged untreated转换组（3策略×5组=15行），
  # 加入Overall后完整校验原18行，再仅删除绘图中的Unchanged intensity一组。
  eff <- eff[scope=='overall' & scenario%in%k & transition%in%setdiff(full_cd,'overall') &
               time==horizon & estimand=='rd_2025_vs_2019']
  tot <- tot[scope=='overall' & scenario%in%k & time==horizon &
               estimand=='rd_2025_vs_2019']
  if(nrow(tot)!=3L || anyDuplicated(tot,by='scenario'))
    stop('Figure 3C必须包含三种RM规则的总体配对RD')
  tot[,transition:='overall']
  eff <- data.table::rbindlist(list(eff,tot),use.names=TRUE,fill=TRUE)
  required_C <- data.table::CJ(scenario=k,transition=full_cd,sorted=FALSE)
  if(nrow(eff)!=18L || anyDuplicated(eff,by=c('scenario','transition')) ||
     nrow(data.table::fsetdiff(required_C,eff[,.(scenario,transition)]))>0L)
    stop('Figure 3C绘图前必须完整校验三策略×六组（含Overall）18行RD')
  group_n <- cnt[,.(scenario,transition,group_n=n)]
  group_n <- data.table::rbindlist(list(group_n,
    totals[,.(scenario,transition='overall',group_n=N)]),use.names=TRUE)
  eff <- merge(eff,group_n,by=c('scenario','transition'),all.x=TRUE,sort=FALSE)
  eff[,estimable:=is.finite(estimate)&is.finite(lcl)&is.finite(ucl)]
  if(any(is.na(eff$group_n) | (eff$group_n>0 & !eff$estimable)))
    stop('Figure 3C非空干预变化组的RD或CI不得缺失')
  rp14_finite_ci(eff[estimable==TRUE],'Figure 3C paired RD')
  contributions <- data.table::rbindlist(lapply(k,function(sc){
    rp14_transition_contribution(src$fig3_contributions,sc,horizon,tr,tot[scenario==sc],
      digits=o$panel_C$contribution_digits,label=paste0('Figure 3C ',sc),
      tolerance=o$contribution_additivity_tolerance)
  }),use.names=TRUE,fill=TRUE)
  contributions <- contributions[transition%in%full_cd] # 完整六组先合并，避免静默丢掉原有一组
  eff <- merge(eff,contributions[,.(scenario,transition,contribution_pp,
    contribution_lcl_pp,contribution_ucl_pp,contribution_status)],
    by=c('scenario','transition'),all.x=TRUE,sort=FALSE)
  if(nrow(eff)!=18L || any(!is.finite(eff$contribution_pp) |
      !is.finite(eff$contribution_lcl_pp) | !is.finite(eff$contribution_ucl_pp)))
    stop('Figure 3C完整18行贡献数值或95%CI缺失')
  # 只筛选画图用的五组×三策略=15行；Overall估计及六组贡献可加性校验均已完成。
  eff <- eff[transition%in%cd]
  if(nrow(eff)!=15L || anyDuplicated(eff,by=c('scenario','transition')) ||
     nrow(data.table::fsetdiff(data.table::CJ(scenario=k,transition=cd,sorted=FALSE),
                               eff[,.(scenario,transition)]))>0L)
    stop('Figure 3C显示结果应恰好包含五组×三种RM策略共15行')
  copt <- o$panel_C
  if(isTRUE(copt$contribution_show_ci)){
    line_sep <- if(identical(copt$contribution_label_layout,'two_line'))'\n' else ' '
    eff[,contribution_print:=paste0(rp_num(contribution_pp,copt$contribution_digits),
      line_sep,'(',rp_num(contribution_lcl_pp,copt$contribution_ci_digits),
      s$ci_interval_separator,rp_num(contribution_ucl_pp,copt$contribution_ci_digits),')')]
  } else eff[,contribution_print:=rp_num(contribution_pp,copt$contribution_digits)]
  eff[,`:=`(pct=100*estimate,lo=100*lcl,hi=100*ucl,
            y=as.numeric(setNames(rev(seq_along(cd)),cd)[transition])+
              as.numeric(copt$strategy_row_offsets[scenario]))]
  if(anyNA(eff$y))stop('Figure 3C三种RM策略的行偏移位置无效')
  list(hvh=hvh,counts=cnt,risk_differences=eff,cohort_n=cohort_N)
}

rp14_fig3_legacy <- function(src,s,layout, classification_only=FALSE){
  o <- layout$fig3
  a <- o$panel_A; b <- o$panel_B; copt <- o$panel_C
  d <- rp14_fig3_prepare(src,s,o)
  palette <- o$rm_colors

  # 根据真实数据和预设间隔构造坐标轴，不截断点估计、置信区间或柱顶标签。
  make_axis <- function(values,lim,breaks,step,tag,include_zero=TRUE){
    z <- as.numeric(values[is.finite(values)])
    if(!length(z))stop('Figure 3',tag,'没有可绘制值')
    if(length(lim)!=2L || any(is.nan(lim)) || !is.finite(step) || step<=0)
      stop('Figure 3',tag,'坐标轴配置无效')
    if(all(is.na(lim))){
      lo <- floor(min(c(z,if(include_zero)0))/step)*step
      hi <- ceiling(max(c(z,if(include_zero)0))/step)*step
    } else if(is.finite(lim[1L]) && is.na(lim[2L])){
      lo <- lim[1L]
      if(min(z)<lo-1e-9)stop('Figure 3',tag,'数据低于固定坐标轴下限')
      hi <- ceiling(max(z)/step)*step
    } else if(all(is.finite(lim)) && lim[2L]>lim[1L]){
      lo <- lim[1L];hi <- lim[2L]
      if(min(z)<lo-1e-9 || max(z)>hi+1e-9)
        stop('Figure 3',tag,'数据超出固定坐标轴；请调整主脚本参数')
    } else stop('Figure 3',tag,'坐标轴上下限无效')
    if(hi<=lo)hi<-lo+step
    zbreaks <- if(is.null(breaks))seq(lo,hi,by=step) else breaks
    list(limits=c(lo,hi),breaks=zbreaks)
  }
  panel_title <- function(p){
    # Panel A较窄时自动将标题显示为两行；不改变主脚本冻结的英文标题原文。
    label <- p$title
    if(identical(p$letter,'A') && !is.null(p$title_wrap_width)) {
      wrapped <- strwrap(label, width=as.integer(p$title_wrap_width))
      label <- paste(wrapped,collapse='\n')
    }
    paste0(p$letter,o$panel_prefix_separator,label)
  }
  # 每个Panel的标题、坐标字号、刻度间距与绘图边距彼此独立。
  theme_one <- function(p,kind){
    m <- p$panel_margin_pt
    if(kind=='A')m[2] <- m[2]+o$layout$top_panel_gap_pt/2
    if(kind=='B')m[4] <- m[4]+o$layout$top_panel_gap_pt/2
    if(kind%in%c('A','B'))m[3] <- m[3]+o$layout$row_panel_gap_pt/2
    if(kind=='C')m[1] <- m[1]+o$layout$row_panel_gap_pt/2
    ggplot2::theme_classic(base_size=s$figure$base_size,base_family=o$font_family)+
      ggplot2::theme(
        plot.title.position='plot',
        plot.title=ggplot2::element_text(face='bold',size=p$title_size_pt,
          lineheight=p$title_lineheight,hjust=0,margin=ggplot2::margin(b=p$title_gap_pt)),
        axis.title.x=ggplot2::element_text(size=p$axis_title_size_pt,
          margin=ggplot2::margin(t=p$x_title_gap_pt)),
        axis.title.y=ggplot2::element_text(size=p$axis_title_size_pt,
          margin=ggplot2::margin(r=if(kind=='A')p$y_title_gap_pt else 0)),
        axis.text.x=ggplot2::element_text(size=if(kind=='A')p$x_tick_size_pt else p$axis_tick_size_pt,
          margin=ggplot2::margin(t=p$x_tick_gap_pt)),
        axis.text.y=ggplot2::element_text(size=p$axis_tick_size_pt,
          margin=ggplot2::margin(r=p$y_tick_gap_pt)),
        axis.line=ggplot2::element_line(color=layout$design$colors$axis,
          linewidth=p$axis_line_width),
        axis.ticks=ggplot2::element_line(linewidth=p$axis_tick_width),
        legend.position='none',
        plot.margin=do.call(ggplot2::margin,as.list(m)),
        plot.background=ggplot2::element_rect(fill='white',color=NA))
  }
  # 将类别从上往下按主脚本顺序绘制；三种RM策略在同一行按各自的偏移排列。
  y_axis <- function(groups,p,labels_override=NULL){
    # 默认沿用Panel B治疗转换的完整名称；C允许传入主脚本中更紧凑的专用短标签。
    nm <- if(is.null(labels_override)) {
      ifelse(groups=='overall',o$overall_label,unname(o$transition_labels[groups]))
    } else {
      unname(labels_override[groups])
    }
    if(anyNA(nm))stop('Figure 3纵轴类别标签缺失')
    ggplot2::scale_y_continuous(
      breaks=rev(seq_along(groups)),labels=nm,
      limits=c(1-p$y_axis_padding,length(groups)+p$y_axis_padding),
      expand=ggplot2::expansion(mult=0))
  }

  # A：两版指南各三根独立柱；柱底原始H/VH，浅色虚线段为M经RM升级。
  aa <- make_axis(c(d$hvh$final_pct+a$label_offset_pp+a$label_top_padding_pp,0),
    a$y_axis_limits,a$y_axis_breaks,a$y_axis_step,'A')
  pA <- ggplot2::ggplot(d$hvh,ggplot2::aes(x=scenario,fill=scenario))+
    (if(isTRUE(a$show_y_grid))ggplot2::geom_hline(yintercept=aa$breaks,
       color=a$grid_color,linewidth=a$grid_width)else NULL)+
    ggplot2::geom_col(ggplot2::aes(y=original_pct),width=a$bar_width,
      color=a$bar_border_color,linewidth=a$bar_border_width,show.legend=FALSE)+
    ggplot2::geom_tile(data=d$hvh[upgraded_pct>0],
      ggplot2::aes(y=upgraded_center,height=upgraded_pct),width=a$bar_width,
      color=a$upgraded_border_color,linetype=a$upgraded_border_linetype,
      linewidth=a$upgraded_border_width,alpha=a$upgraded_alpha,show.legend=FALSE)+
    ggplot2::geom_text(ggplot2::aes(y=final_pct+a$label_offset_pp,
      label=paste0(rp_num(final_pct,a$percent_digits),'%')),
      size=a$percent_label_size_mm,vjust=0,family=o$font_family,
      color=layout$design$colors$text,show.legend=FALSE)+
    ggplot2::scale_fill_manual(values=palette,guide='none',drop=FALSE)+
    ggplot2::scale_x_discrete(limits=o$classification_scenarios,
      labels=unname(a$strategy_labels[o$classification_scenarios]),drop=FALSE,
      expand=ggplot2::expansion(add=a$x_expand))+
    ggplot2::scale_y_continuous(breaks=aa$breaks,
      labels=function(x)rp_num(x,a$y_axis_tick_digits),
      expand=ggplot2::expansion(mult=0))+
    ggplot2::facet_wrap(~guideline_label,nrow=1,strip.position='bottom')+
    ggplot2::coord_cartesian(ylim=aa$limits,clip='off')+
    ggplot2::labs(title=panel_title(a),x=a$x_title,y=a$y_title)+theme_one(a,'A')+
    ggplot2::theme(
      panel.spacing.x=grid::unit(a$facet_gap_pt,'pt'),
      # strip.placement='outside'：指南标签排在SN/ST/SNT刻度文字下方。
      # 若不设置，ggplot2会把年份条放在坐标轴与策略标签之间，造成标签上下颠倒。
      strip.placement='outside',
      strip.background=ggplot2::element_blank(),
      strip.text.x.bottom=ggplot2::element_text(face='plain',size=a$facet_label_size_pt,
        margin=ggplot2::margin(t=a$strip_gap_pt,b=0)))

  if (classification_only) return(list(plot=pA, data=d$hvh))

  # B：保留六组×三RM共18根横向柱；默认彻底删除柱末百分比数值，减轻拥挤。
  # 关闭数值后不再使用旧的5pp文字预留空间，而是仅留很小的bar_right_padding_pp。
  b_label_space <- if(isTRUE(b$show_value_labels))
    b$bar_label_offset_pp+b$label_right_padding_pp else b$bar_right_padding_pp
  ab <- make_axis(c(d$counts$percent+b_label_space,0),
    b$x_axis_limits,b$x_axis_breaks,b$x_axis_step,'B')
  pB <- ggplot2::ggplot(d$counts)+
    (if(isTRUE(b$show_row_guides))ggplot2::geom_hline(yintercept=seq_along(o$transition_order),
      color=b$row_guide_color,linewidth=b$row_guide_width)else NULL)+
    ggplot2::geom_rect(ggplot2::aes(xmin=0,xmax=percent,
      ymin=y-b$bar_height/2,ymax=y+b$bar_height/2,fill=scenario),
      color=b$bar_border_color,linewidth=b$bar_border_width)+
    # TRUE才添加数值文字图层；FALSE时不会绘制任何18组柱末百分比。
    (if(isTRUE(b$show_value_labels))ggplot2::geom_text(
       ggplot2::aes(x=percent+b$bar_label_offset_pp,y=y,
         label=paste0(rp_num(percent,b$percent_digits),'%')),
       hjust=0,size=b$bar_label_size_mm,family=o$font_family)else NULL)+
    ggplot2::scale_fill_manual(values=palette,guide='none')+
    y_axis(o$transition_order,b)+
    ggplot2::scale_x_continuous(breaks=ab$breaks,
      labels=function(x)rp_num(x,b$x_axis_digits),expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(xlim=ab$limits,clip='off')+
    ggplot2::labs(title=panel_title(b),x=b$x_title,y=b$y_title)+theme_one(b,'B')

  # C：独立森林图区与文字区，避免通过过大margin强行容纳贡献数值导致右侧留白。
  ac <- make_axis(c(d$risk_differences$lo[d$risk_differences$estimable],
                    d$risk_differences$hi[d$risk_differences$estimable],0),
    copt$x_axis_limits,copt$x_axis_breaks,copt$x_axis_step,'C')
  pCforest <- ggplot2::ggplot(d$risk_differences)+
    (if(isTRUE(copt$show_row_guides))ggplot2::geom_hline(yintercept=seq_along(o$comparative_order),
       color=copt$row_guide_color,linewidth=copt$row_guide_width)else NULL)+
    (if(isTRUE(copt$show_zero_line))ggplot2::geom_vline(xintercept=0,
       color=copt$zero_line_color,linetype=copt$zero_line_type,
       linewidth=copt$zero_line_width)else NULL)+
    ggplot2::geom_point(data=d$risk_differences[estimable==TRUE],
      ggplot2::aes(x=pct,y=y,color=scenario),shape=copt$point_shape,
      fill=copt$point_fill_color,size=copt$point_size,stroke=copt$point_stroke)+
    ggplot2::geom_errorbar(data=d$risk_differences[estimable==TRUE],
      ggplot2::aes(xmin=lo,xmax=hi,y=y,color=scenario),orientation='y',
      width=copt$ci_cap_height,linewidth=copt$ci_linewidth)+
    ggplot2::geom_text(data=d$risk_differences[estimable==FALSE],
      ggplot2::aes(x=0,y=y),label='—',inherit.aes=FALSE,
      size=copt$point_size,color=layout$design$colors$text)+
    ggplot2::scale_color_manual(values=palette,guide='none')+
    y_axis(o$comparative_order,copt,labels_override=copt$y_category_labels)+
    ggplot2::scale_x_continuous(breaks=ac$breaks,
      labels=function(x)rp_num(x,copt$x_axis_digits),
      expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(xlim=ac$limits,clip='off')+
    ggplot2::labs(title=panel_title(copt),x=copt$x_title,y=copt$y_title)+theme_one(copt,'C')

  # C的贡献列每行保持与森林图同一个y坐标；95%CI端点均直接使用已保存值。
  # 标题与森林图独立占位，避免贡献标题压在Panel C总标题上。
  right_margin <- copt$panel_margin_pt
  right_margin[1] <- right_margin[1]+o$layout$row_panel_gap_pt/2
  right_margin[2] <- copt$contribution_right_margin_pt
  right_margin[4] <- copt$contribution_left_margin_pt # 空白留在数值与标题之前，二者左对齐
  pCcontrib <- ggplot2::ggplot(d$risk_differences)+
    ggplot2::geom_text(ggplot2::aes(x=copt$contribution_text_x,y=y,
      label=contribution_print),hjust=0,
      size=copt$contribution_text_size_mm,
      lineheight=copt$contribution_lineheight,family=o$font_family,
      color=layout$design$colors$text)+
    ggplot2::scale_x_continuous(limits=c(0,1),expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_y_continuous(limits=c(1-copt$y_axis_padding,
      length(o$comparative_order)+copt$y_axis_padding),
      expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(clip='off')+
    ggplot2::labs(title=copt$contribution_header)+
    ggplot2::theme_void(base_family=o$font_family)+
    ggplot2::theme(
      plot.title.position='plot',
      plot.title=ggplot2::element_text(face='bold',size=copt$contribution_header_size_pt,
        hjust=0,margin=ggplot2::margin(b=copt$contribution_title_gap_pt)),
      plot.margin=do.call(ggplot2::margin,as.list(right_margin)))
  lower <- patchwork::wrap_plots(pCforest,pCcontrib,ncol=2,
    widths=copt$contribution_plot_widths)

  # 图例单行：按“策略色块3个 + 原始H/VH三色拼块 + RM升级虚线块3个”展示。
  # 放置在A/B和C之间；全部文字与色块共享同一个Y坐标，避免换行。
  # 颜色与上游结果完全不变；所有位置和字号均由主脚本统一控制。
  leg <- o$legend
  footer <- ggplot2::ggplot() +
    ggplot2::annotate('text',x=leg$x_title,y=leg$y,label=leg$title,
      hjust=0,fontface='plain',size=leg$title_size_mm,family=o$font_family,
      color=layout$design$colors$text)
  for(sc in o$classification_scenarios){
    xpos <- unname(leg$x_swatches[[sc]])
    footer <- footer + ggplot2::annotate('rect',
      xmin=xpos-leg$swatch_width/2,xmax=xpos+leg$swatch_width/2,
      ymin=leg$y-leg$swatch_height/2,
      ymax=leg$y+leg$swatch_height/2,
      fill=palette[[sc]],color=layout$design$colors$bar_edge,
      linewidth=leg$swatch_border_width)+
      ggplot2::annotate('text',x=xpos+leg$marker_text_gap,y=leg$y,
        label=leg$labels[[sc]],hjust=0,size=leg$text_size_mm,
        family=o$font_family,color=layout$design$colors$text)
  }
  footer <- footer + ggplot2::annotate('text',
      x=leg$component_title_x,y=leg$y,
      label=leg$component_title,hjust=0,fontface='plain',
      size=leg$title_size_mm,family=o$font_family,
      color=layout$design$colors$text)
  # 原始H/VH可以属于任何RM策略，使用三色实色拼接示意，避免只显示蓝色而误导。
  for(i in seq_along(o$classification_scenarios)){
    sc <- o$classification_scenarios[[i]]
    n_sc <- length(o$classification_scenarios)
    xmin <- leg$original_swatch_x-leg$component_swatch_width/2+
      (i-1)*leg$component_swatch_width/n_sc
    footer <- footer + ggplot2::annotate('rect',
      xmin=xmin,xmax=xmin+leg$component_swatch_width/n_sc,
      ymin=leg$y-leg$swatch_height/2,
      ymax=leg$y+leg$swatch_height/2,
      fill=palette[[sc]],color=NA)
  }
  footer <- footer +
    ggplot2::annotate('rect',
      xmin=leg$original_swatch_x-leg$component_swatch_width/2,
      xmax=leg$original_swatch_x+leg$component_swatch_width/2,
      ymin=leg$y-leg$swatch_height/2,
      ymax=leg$y+leg$swatch_height/2,
      fill=NA,color=layout$design$colors$bar_edge,
      linewidth=leg$swatch_border_width)+
    ggplot2::annotate('text',
      x=leg$original_swatch_x+leg$marker_text_gap,y=leg$y,
      label=leg$original_label,hjust=0,size=leg$text_size_mm,
      family=o$font_family,color=layout$design$colors$text)
  # RM升级段：必须用蓝/橙/绿三种浅色虚线色块，且分别有文字标签。
  for(sc in o$classification_scenarios){
    xpos <- unname(leg$upgraded_swatch_x[[sc]])
    footer <- footer + ggplot2::annotate('rect',
      xmin=xpos-leg$component_swatch_width/2,
      xmax=xpos+leg$component_swatch_width/2,
      ymin=leg$y-leg$swatch_height/2,
      ymax=leg$y+leg$swatch_height/2,
      fill=scales::alpha(palette[[sc]],leg$component_alpha),
      color=layout$design$colors$bar_edge,
      linetype=leg$component_border_linetype,
      linewidth=leg$swatch_border_width)+
      ggplot2::annotate('text',
        x=xpos+leg$marker_text_gap,y=leg$y,
        label=leg$upgraded_labels[[sc]],hjust=0,size=leg$text_size_mm,
        family=o$font_family,color=layout$design$colors$text)
  }
  footer <- footer +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE,clip='off')+
    ggplot2::theme_void(base_family=o$font_family)+
    ggplot2::theme(plot.margin=do.call(ggplot2::margin,
      as.list(leg$legend_margin_pt)))

  # 最终图按上至下顺序拼接：Panel A+B -> 单行图例 -> Panel C。
  # 图例使用独立的完整画幅宽度行，而非附着于Panel A、B或C某一面板。
  # 使用wrap_elements隔离A/B轴布局，避免B长分组标签压缩A绘图宽度。
  upper <- patchwork::wrap_plots(patchwork::wrap_elements(full=pA),
    patchwork::wrap_elements(full=pB),ncol=2,widths=o$layout$top_panel_widths)
  height_top <- unname(o$layout$top_bottom_heights[[1L]])
  height_bottom <- unname(o$layout$top_bottom_heights[[2L]])
  height_legend <- height_top*o$layout$legend_height_ratio
  # 独立封装图例行：使图例以整张画布左边界为基准，不受Panel C轴标签列对齐影响。
  footer_full <- patchwork::wrap_elements(full=footer)
  full <- patchwork::wrap_plots(upper,footer_full,lower,ncol=1,
    heights=c(height_top,height_legend,height_bottom)) +
    patchwork::plot_annotation(theme=ggplot2::theme(
      plot.margin=do.call(ggplot2::margin,as.list(o$layout$figure_margin_pt))))
  list(plot=full,plot_data=list(
    high_very_high_classification=d$hvh,
    transition_counts=d$counts,
    estimated_risk_difference=d$risk_differences))
}

# ============================================================================
# Figure 4: formerly Figure 3; four outline colors and IDENTICAL hollow circles.
# ============================================================================
rp14_fig4_legacy <- function(src,s,layout) {
  # Figure 4：A=四策略总体指南间RD；B=三种RM相对S0的配对RD变化；C=8个亚组四策略RD（不显示Overall）。
  # 所有CI均取自Module 11已合并输出，B仅做严格符号转换，不重新估计不确定性。
  o <- layout$fig4
  pal <- unlist(s$colors$rm[s$scenario_order],use.names=TRUE)
  all <- data.table::copy(data.table::as.data.table(src$all_effects))
  rp_assert(all,c('scenario','scope','time','estimand','estimate','lcl','ucl'),'Figure 4 paired RDs')
  all[,`:=`(scenario=as.character(scenario),scope=as.character(scope),
             estimand=as.character(estimand))]
  rd <- all[time==s$primary_horizon & scope%in%s$scope_order &
            scenario%in%s$scenario_order & estimand=='rd_2025_vs_2019']
  rp_check_unique(rd,c('scope','scenario','estimand'),'Figure 4 RD')
  if(nrow(rd)!=length(s$scope_order)*length(s$scenario_order))
    stop('Figure 4 A/C: requires 9 subgroups x 4 strategies with paired 2025-2019 RD')
  rp14_finite_ci(rd,'Figure 4 paired RD')
  rd[,`:=`(pct=100*estimate,lo=100*lcl,hi=100*ucl,
            scenario=factor(scenario,levels=s$scenario_order),
            scope=factor(scope,levels=s$scope_order))]
  overall <- rd[scope=='overall']
  if(nrow(overall)!=4L)stop('Figure 4A: four overall strategy RDs required')

  # Module 11的additional_reduction_difference_2025_vs_2019为 RD_S0 - RD_strategy。
  # Panel B要求Δ=RD_strategy-RD_S0，因此对已经pool的点估计及95%CI同步取负并交换上下界。
  inc <- data.table::copy(data.table::as.data.table(src$rm_incremental))
  rp_assert(inc,c('scenario','scope','time','estimand','estimate','lcl','ucl'),
            'Figure 4B paired incremental RD')
  inc[,`:=`(scenario=as.character(scenario),scope=as.character(scope),
             estimand=as.character(estimand))]
  inc <- inc[time==s$primary_horizon & scope=='overall' &
    scenario%in%c('SN','ST','SNT') &
    estimand=='additional_reduction_difference_2025_vs_2019']
  if('rm_comparison'%in%names(inc)){
    # 兼容历史SH_vs_S0/SHT_vs_S0编码，但不修改pooled CSV。
    inc[,rm_comparison:=as.character(rm_comparison)]
    inc[rm_comparison=='SH_vs_S0',rm_comparison:='SN_vs_S0']
    inc[rm_comparison=='SHT_vs_S0',rm_comparison:='SNT_vs_S0']
    inc <- inc[rm_comparison%in%c('SN_vs_S0','ST_vs_S0','SNT_vs_S0')]
  }
  if(nrow(inc)!=3L || anyDuplicated(inc,by='scenario'))
    stop('Figure 4B: needs SN/ST/SNT versus S0, three independently pooled paired comparisons')
  rp14_finite_ci(inc,'Figure 4B incremental CI')
  inc[,`:=`(pct=-100*estimate,lo=-100*ucl,hi=-100*lcl,
             scenario=factor(scenario,levels=s$scenario_order))]
  if(any(inc$lo>inc$pct | inc$pct>inc$hi))stop('Figure 4B: reversed paired CI invalid')

  # 固定轴范围不足时自动扩大，只改变显示区间，不裁切真实CI。
  axis_one <- function(q,id) {
    wanted <- o$panel_y_limits[[id]]
    lim <- c(min(wanted[1],q$lo,0),max(wanted[2],q$hi,0))
    if(any(lim!=wanted)) warning('Figure 4 ',id,'：置信区间超过预设轴范围，已扩大显示范围')
    list(limits=lim,breaks=o$panel_y_breaks[[id]])
  }
  # 独立panel版式：标题/坐标/间距均来自14主脚本fig4，减少重复硬编码。
  panel_theme <- function(id){
    mm <- o$panel_margin_pt
    # 上排A/B各承担一半列间隔；上排共同承担行间距。
    if(id=='A') mm[2] <- mm[2]+o$layout_column_gap_pt/2
    if(id=='B') mm[4] <- mm[4]+o$layout_column_gap_pt/2
    if(id %in% c('A','B')) mm[3] <- mm[3]+o$layout_row_gap_pt/2
    if(id=='C') mm[1] <- mm[1]+o$layout_row_gap_pt/2
    rp14_theme(s,layout,'none',margins=mm) + ggplot2::theme(
      plot.title.position='plot',
      plot.title=ggplot2::element_text(size=o$panel_title_size,face='bold',lineheight=1.06,
        margin=ggplot2::margin(b=unname(o$panel_title_gap_pt[[id]]))),
      axis.title.y=ggplot2::element_text(size=unname(o$panel_y_title_size[[id]]),
        margin=ggplot2::margin(r=unname(o$y_title_gap_pt[[id]]))),
      axis.text.y=ggplot2::element_text(size=o$axis_text_size,
        margin=ggplot2::margin(r=unname(o$y_tick_gap_pt[[id]]))),
      axis.text.x=ggplot2::element_text(size=if(id=='C')o$subgroup_text_size else o$top_x_label_size_pt,
        margin=ggplot2::margin(t=unname(o$x_tick_gap_pt[[id]]))),
      axis.line=ggplot2::element_line(color=layout$design$colors$axis),
      legend.position='none')
  }
  legend_scale <- function(){
    ggplot2::scale_color_manual(values=pal,breaks=s$scenario_order,
      labels=unname(o$legend_labels[s$scenario_order]),drop=FALSE)
  }
  plot_small <- function(q,id,keys){
    ax <- axis_one(q,id)
    q[,scenario:=factor(as.character(scenario),levels=keys)]
    ggplot2::ggplot(q,ggplot2::aes(x=scenario,y=pct,color=scenario))+
      ggplot2::geom_hline(yintercept=0,color=o$zero_line_color,
        linetype='dashed',linewidth=o$zero_line_width)+
      # Figure 4 A/B：先画空心圆，再画95%CI；CI线条位于圆圈上层。
      ggplot2::geom_point(shape=o$point_shape,fill=layout$design$colors$white,
        size=o$point_size,stroke=o$point_stroke)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),
        linewidth=o$ci_linewidth,width=unname(o$ci_cap_width[[id]]),show.legend=FALSE)+
      legend_scale()+
      ggplot2::scale_x_discrete(limits=keys,
        labels=unname(o$top_x_labels[keys]),drop=FALSE,
        expand=ggplot2::expansion(add=c(.30,.30)))+
      ggplot2::scale_y_continuous(breaks=ax$breaks,
        expand=ggplot2::expansion(mult=0))+
      ggplot2::coord_cartesian(ylim=ax$limits,clip='off')+
      ggplot2::labs(title=paste0(o$panel_letters[[id]],'  ',paste(strwrap(o$panel_titles[[id]],width=unname(o$panel_title_wrap_width[[id]])),collapse='\n')),
        x=NULL,y=o$y_titles[[id]])+panel_theme(id)
  }
  pA <- plot_small(data.table::copy(overall),'A',s$scenario_order)
  pB <- plot_small(data.table::copy(inc),'B',c('SN','ST','SNT'))

  # Panel C单独排除Overall；Panel A仍使用原始总体RD，统计结果不变。
  subgroup_scopes <- setdiff(s$scope_order, 'overall')
  rd_c <- data.table::copy(rd[as.character(scope) %in% subgroup_scopes])
  if(nrow(rd_c)!=length(subgroup_scopes)*length(s$scenario_order) ||
     anyDuplicated(rd_c,by=c('scope','scenario')))
    stop('Figure 4C: requires 8 subgroups x 4 strategies (Overall excluded)')
  rd_c[,scope:=factor(as.character(scope),levels=subgroup_scopes)]
  axc <- axis_one(rd_c,'C')
  dodge <- ggplot2::position_dodge(width=o$dodge_width)
  pC <- ggplot2::ggplot(rd_c,ggplot2::aes(x=scope,y=pct,color=scenario,group=scenario))+
    ggplot2::geom_vline(xintercept=o$group_dividers,linetype='dashed',
      color=o$group_divider_color,linewidth=o$group_divider_width)+
    ggplot2::geom_hline(yintercept=0,linetype='dashed',
      color=o$zero_line_color,linewidth=o$zero_line_width)+
    # Figure 4 C：先画空心圆，再画95%CI；两层共用同一个dodge位置。
    ggplot2::geom_point(position=dodge,shape=o$point_shape,
      fill=layout$design$colors$white,size=o$point_size,stroke=o$point_stroke)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lo,ymax=hi),position=dodge,
      linewidth=o$ci_linewidth,width=unname(o$ci_cap_width[['C']]),show.legend=FALSE)+
    legend_scale()+
    ggplot2::scale_x_discrete(limits=subgroup_scopes,
      labels=unname(o$subgroup_labels[subgroup_scopes]),drop=FALSE)+
    ggplot2::scale_y_continuous(breaks=axc$breaks,
      expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(ylim=axc$limits,clip='off')+
    ggplot2::labs(title=paste0('C  ',o$panel_titles[['C']]),x=NULL,y=o$y_titles[['C']])+
    panel_theme('C')+
    ggplot2::annotate('text',x=o$favors_x,
      y=axc$limits[2]-o$favors_top_fraction*diff(axc$limits),
      label=o$favors_2019,hjust=1,vjust=1,fontface='italic',
      size=o$favors_label_size,family=o$font_family)+
    ggplot2::annotate('text',x=o$favors_x,
      y=axc$limits[1]+o$favors_bottom_fraction*diff(axc$limits),
      label=o$favors_2025,hjust=1,vjust=0,fontface='italic',
      size=o$favors_label_size,family=o$font_family)
  # 只从Panel C提取一个图例，保证A/B没有重复图例且最终仍单行。
  pC <- pC + ggplot2::guides(color=ggplot2::guide_legend(
    nrow=1,ncol=o$legend_ncol,byrow=TRUE,
    override.aes=list(shape=o$point_shape,fill=layout$design$colors$white,
      size=o$legend_point_size,stroke=o$point_stroke)))+
    ggplot2::theme(legend.position='bottom',
      legend.text=ggplot2::element_text(size=o$legend_size),
      legend.title=ggplot2::element_blank(),
      legend.key.width=grid::unit(o$legend_key_width_pt,'pt'),
      legend.key.height=grid::unit(o$legend_key_height_pt,'pt'),
      legend.spacing.x=grid::unit(o$legend_spacing_pt,'pt'),
      legend.margin=do.call(ggplot2::margin,as.list(o$legend_margin_pt)))
  # 第一行A/B并列，第二行C跨整行；A/B维持各自独立的纵轴范围。
  top <- patchwork::wrap_plots(pA,pB,ncol=2,widths=o$layout_top_widths)
  whole <- patchwork::wrap_plots(top,pC,ncol=1,heights=o$layout_row_heights)
  list(plot=whole,plot_data=list(overall_risk_differences=overall,
    changes_vs_no_RM=inc,subgroup_risk_differences=rd_c))
}

# ============================================================================
# Supplementary Figure 1–3：SN / ST / SNT 风险重分类及观察事件捕获
# 仅从 Module 10/11 的 pooled CSV 读取；25格全精度数据完整保留。
# 所有画幅、颜色、字号、边框、定位与文字参数均在主脚本 supplementary 中。
# 仅改变附件三图；正文 Figure 1–4 函数不修改。
# ============================================================================
rp14_supplementary <- function(src, s, layout, scenario) {
  o <- layout$supplementary
  k <- scenario
  lev <- c('L', 'M_not_upgraded', 'M_upgraded_H', 'H', 'VH')
  counts <- data.table::copy(src$joint_counts[
    guideline == 'joint' & scope == 'overall' & scenario == k])
  d <- data.table::copy(src$joint_risks[
    scope == 'overall' & scenario == k & time == s$primary_horizon])
  rp_assert(counts, c('cell_id','origin_2019','origin_2025','n','proportion'), 'joint counts')
  rp_assert(d, c('cell_id','estimand','estimate','lcl','ucl'), 'joint risks')
  if(nrow(counts)!=25L || anyDuplicated(counts$cell_id) ||
     !setequal(counts$origin_2019, lev) || !setequal(counts$origin_2025, lev))
    stop('Supplementary: 25 independent cells required')
  if(abs(sum(counts$proportion)-1)>1e-6)
    stop('Supplementary: joint proportions do not sum to 1')
  risk <- d[estimand=='joint_category_risk']
  cap <- d[estimand=='joint_category_capture']
  if(nrow(risk)!=25L || nrow(cap)!=25L ||
     anyDuplicated(risk$cell_id) || anyDuplicated(cap$cell_id))
    stop('Supplementary: incomplete joint risk/capture pooled rows')
  a <- merge(counts, risk[,.(cell_id,estimate,lcl,ucl,status)], by='cell_id')
  b <- merge(counts, cap[,.(cell_id,estimate,lcl,ucl,status)], by='cell_id')
  if(nrow(a)!=25L || nrow(b)!=25L)stop('Supplementary: count/risk merge incomplete')
  if(any(is.finite(b$estimate)) && abs(sum(b$estimate,na.rm=TRUE)-1)>1e-5)
    stop('Supplementary: normalized event capture must total 1')
  if(any(a$n!=b$n))stop('Supplementary: count mismatch between panels')

  # 显示列数自动随是否省略全空VH列切换；不删除原始25个单元格及导出CSV。
  vh_empty <- all(counts[origin_2025=='VH', n] == 0)
  visible_levels <- if(isTRUE(o$hide_empty_2025_vh) && vh_empty) lev[1:4] else lev
  if(isTRUE(o$hide_empty_2025_vh) && !vh_empty)
    warning('Supplementary ',k,': 2025 VH列存在参与者，自动保留，防止隐藏真实数据')
  ncols <- length(visible_levels)
  right <- ncols + 0.5
  v_split <- o$vertical_split_x
  # 色标改为每个Panel左上角嵌入，不再为右侧图例保留额外横向空间。
  if(!identical(o$legend_position, 'upper_left'))
    stop('补充图legend_position必须为upper_left')
  legend_x <- o$legend_x
  if(length(legend_x)!=2L || any(!is.finite(legend_x)) ||
     legend_x[1]>=legend_x[2] || legend_x[1]<o$side_header_x[1] || legend_x[2]>=0.5)
    stop('补充图左上角色标横向坐标必须在侧表头左缘与数据矩阵左缘之间')
  if(length(o$legend_y)!=2L || any(!is.finite(o$legend_y)) ||
     o$legend_y[1]>=o$legend_y[2] || o$legend_y[1]<=o$matrix_top ||
     o$legend_y[2]>=o$y_limits[2] ||
     !(o$legend_tick_y>o$matrix_top && o$legend_tick_y<o$legend_y[1]) ||
     !(o$legend_title_y>o$legend_y[2] && o$legend_title_y<o$y_limits[2]))
    stop('补充图左上角色标Y坐标与矩阵或标题区域冲突')
  xlim <- c(o$side_header_x[1] - o$x_left_padding,
            right + o$x_right_padding)
  # 所有标签与分区均基于显示列数动态定位，避免隐藏VH列后错误留白。
  right_center <- (v_split + right)/2

  for(nm in c('a','b')) {
    z <- get(nm)
    z[,`:=`(x=match(origin_2025,lev), y=6-match(origin_2019,lev),
             pct=100*estimate, lo=100*lcl, hi=100*ucl)]
    if(anyNA(z$x) || anyNA(z$y))stop('Supplementary: invalid category labels')
    assign(nm,z)
  }
  if(isTRUE(o$warn_unbounded_risk_ci)) {
    bad <- a[n>0 & (is.finite(lo) & lo<0 | is.finite(hi) & hi>100)]
    if(nrow(bad)) warning('Supplementary ',k,': ',nrow(bad),
      '个单元格风险95%CI越过[0,100]%，仅提示核查上游区间算法，图中原样保留')
  }

  # 附件图专用百分比格式：只改变标签，完整精度结果/CSV均不变。
  # 规则：真实零=0.0；0<百分比<0.05时显示<0.1（按位数自动计算）；
  # 负值（包括负CI下限）照常显示，不进行裁剪；缺失CI仍用“—”。
  fmt_pct_number <- function(v, digits) {
    ans <- rp_num(v,digits)
    if(isTRUE(o$show_less_than_for_tiny_positive)) {
      tiny <- is.finite(v) & v>0 & v < 0.5*10^(-digits)
      ans[tiny] <- paste0('<',formatC(10^(-digits),format='f',digits=digits))
    }
    ans
  }
  fmt_pct <- function(v, digits) paste0(fmt_pct_number(v,digits),'%')
  ci_lab <- function(x, low, high) {
    ans <- rep(paste0('(',o$missing_ci_symbol,')'),length(x))
    ok <- is.finite(x) & is.finite(low) & is.finite(high)
    ans[ok] <- paste0('(',fmt_pct_number(low[ok],o$cell_ci_digits),
                      s$ci_interval_separator,
                      fmt_pct_number(high[ok],o$cell_ci_digits),')')
    ans
  }
  a[,`:=`(cell_line1=ifelse(n>0,
      paste0(rp_n(n),' (',fmt_pct(100*proportion,o$cell_count_digits),')'),'—'),
    cell_line2=ifelse(n==0,'',ifelse(is.finite(pct),
      paste0(fmt_pct(pct,o$cell_value_digits),' ',ci_lab(pct,lo,hi)),o$labels$empty_risk)))]
  b[,`:=`(cell_line1=ifelse(n>0 & is.finite(pct),fmt_pct(pct,o$cell_value_digits),'—'),
    cell_line2=ifelse(n==0,'',ifelse(is.finite(pct),ci_lab(pct,lo,hi),o$missing_ci_symbol)))]

  # 深色热图自动使用白字：通过真实渐变色的RGB亮度确定，而非只用数值大小。
  cell_colors <- function(z, colors, limits) {
    out <- rep(o$text_color,nrow(z))
    if(!isTRUE(o$use_adaptive_text_color))return(out)
    mapper <- scales::col_numeric(palette=colors,domain=limits)
    ok <- is.finite(z$pct) & z$n>0
    if(any(ok)) {
      hex <- mapper(pmin(limits[2],pmax(limits[1],z$pct[ok])))
      rgb <- grDevices::col2rgb(hex)/255
      lum <- as.numeric(c(0.2126,0.7152,0.0722) %*% rgb)
      out[ok] <- ifelse(lum < o$dark_cell_luminance_threshold,
                         o$dark_cell_text_color,o$text_color)
    }
    out
  }

  display_scenario <- unname(o$scenario_labels[[k]])
  lab19 <- unname(gsub('{scenario}',display_scenario,o$labels$row_categories,fixed=TRUE))
  lab25 <- unname(gsub('{scenario}',display_scenario,o$labels$column_categories,fixed=TRUE))
  header_color <- o$header_fill
  border <- o$cell_border_color
  draw_panel <- function(z, type) {
    is_risk <- identical(type,'risk')
    # 左上角色标仅保留四个可辨识刻度，色阶的实际上下限仍保持不变。
    ticks <- if(is_risk)o$legend_risk_breaks else o$legend_capture_breaks
    lim <- if(is_risk)o$risk_limits else o$capture_limits
    cols <- if(is_risk)o$risk_colors else o$capture_colors
    if(any(z$pct<lim[1] | z$pct>lim[2],na.rm=TRUE))
      warning('Supplementary ',k,' ',type,': 色阶外的真实数值仅压缩颜色，不改变估计及CI')
    z <- data.table::copy(z[origin_2025 %in% visible_levels])
    z[,`:=`(text_color=cell_colors(z,cols,lim), fill_value=ifelse(n>0,pct,NA_real_))]
    step <- as.integer(o$legend_steps)
    gg <- data.table::data.table(
      xmin=seq(legend_x[1],legend_x[2],length.out=step+1L)[- (step+1L)],
      xmax=seq(legend_x[1],legend_x[2],length.out=step+1L)[-1L],
      ymin=o$legend_y[1],ymax=o$legend_y[2],
      value=seq(lim[1],lim[2],length.out=step))
    tick_data <- data.table::data.table(
      x=legend_x[1]+(ticks-lim[1])/diff(lim)*diff(legend_x),
      y=o$legend_tick_y,label=as.character(ticks))
    hdr <- data.table::data.table(x=seq_len(ncols),label=lab25[seq_len(ncols)])
    rw <- data.table::data.table(y=5:1,label=lab19)
    p <- ggplot2::ggplot() +
      # 热图单元格；先画背景色，再分别绘制第一行、第二行数字。
      ggplot2::geom_tile(data=z,ggplot2::aes(x=x,y=y,fill=fill_value),
        width=1,height=1,color=border,linewidth=o$cell_border_width)+
      ggplot2::geom_text(data=z,ggplot2::aes(x=x,y=y+o$cell_first_line_offset,
        label=cell_line1,color=text_color),size=o$cell_count_text_size,
        fontface='bold',lineheight=o$cell_lineheight,family=s$figure$font_family)+
      ggplot2::geom_text(data=z,ggplot2::aes(x=x,y=y+o$cell_second_line_offset,
        label=cell_line2,color=text_color),size=o$cell_ci_text_size,
        lineheight=o$cell_lineheight,family=s$figure$font_family)+
      # 顶部2025年合并表头、类别列名及两个高危分区说明条。
      ggplot2::annotate('rect',xmin=.5,xmax=right,
        ymin=o$top_header_y[1],ymax=o$top_header_y[2],
        fill=header_color,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('text',x=(.5+right)/2,y=mean(o$top_header_y),
        label=o$labels$top_header,fontface='bold',size=o$header_text_size,
        family=s$figure$font_family)+
      ggplot2::geom_rect(data=hdr,
        ggplot2::aes(xmin=x-.5,xmax=x+.5,
          ymin=o$column_header_y[1],ymax=o$column_header_y[2]),
        fill=header_color,color=border,linewidth=o$cell_border_width)+
      ggplot2::geom_text(data=hdr,
        ggplot2::aes(x=x,y=o$column_label_text_y,label=label),
        size=o$header_text_size,fontface='bold',family=s$figure$font_family)+
      ggplot2::annotate('rect',xmin=.5,xmax=v_split,
        ymin=o$region_header_y[1],ymax=o$region_header_y[2],
        fill=o$region_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('rect',xmin=v_split,xmax=right,
        ymin=o$region_header_y[1],ymax=o$region_header_y[2],
        fill=o$region_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('text',x=(.5+v_split)/2,y=o$region_header_text_y,
        label=o$labels$top_left,size=o$region_text_size,fontface='italic',
        color=o$secondary_text_color,family=s$figure$font_family)+
      ggplot2::annotate('text',x=right_center,y=o$region_header_text_y,
        label=o$labels$top_right,size=o$region_text_size,fontface='italic',
        color=o$secondary_text_color,family=s$figure$font_family)+
      # 左边2019合并表头、5个风险类别标签。
      ggplot2::geom_rect(data=rw,ggplot2::aes(
        xmin=o$row_label_x[1],xmax=o$row_label_x[2],ymin=y-.5,ymax=y+.5),
        fill=o$row_header_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::geom_text(data=rw,ggplot2::aes(x=o$row_label_text_x,y=y,label=label),
        fontface='bold',size=o$label_text_size,family=s$figure$font_family)+
      ggplot2::annotate('rect',xmin=o$side_header_x[1],xmax=o$side_header_x[2],
        ymin=o$matrix_bottom,ymax=o$matrix_top,
        fill=o$row_header_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('text',x=o$side_header_text_x,y=3,
        label=o$labels$side_header,fontface='bold',size=o$label_text_size,
        lineheight=1.0,family=s$figure$font_family)+
      # 底部分区标题随显示列宽自动居中。
      ggplot2::annotate('rect',xmin=.5,xmax=v_split,
        ymin=o$footer_y[1],ymax=o$footer_y[2],
        fill=o$footer_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('rect',xmin=v_split,xmax=right,
        ymin=o$footer_y[1],ymax=o$footer_y[2],
        fill=o$footer_fill,color=border,linewidth=o$cell_border_width)+
      ggplot2::annotate('text',x=(.5+v_split)/2,y=o$footer_text_y,
        label=o$labels$bottom_left,size=o$footer_text_size,fontface='italic',
        color=o$secondary_text_color,family=s$figure$font_family)+
      ggplot2::annotate('text',x=right_center,y=o$footer_text_y,
        label=o$labels$bottom_right,size=o$footer_text_size,fontface='italic',
        color=o$secondary_text_color,family=s$figure$font_family)+
      # 两条粗线明确四象限；数据区不额外上色以免干扰风险热图。
      ggplot2::annotate('segment',x=v_split,xend=v_split,
        y=o$footer_y[1],yend=o$column_header_y[2],
        linewidth=o$group_border_width,color=o$group_border_color)+
      ggplot2::annotate('segment',x=.5,xend=right,
        y=o$horizontal_split_y,yend=o$horizontal_split_y,
        linewidth=o$group_border_width,color=o$group_border_color)+
      # 左上角嵌入式色条/刻度；Panel A蓝色、Panel B紫色，三种RM使用相同范围。
      ggplot2::geom_rect(data=gg,ggplot2::aes(
        xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=value),
        inherit.aes=FALSE,color=NA)+
      ggplot2::geom_text(data=tick_data,ggplot2::aes(x=x,y=y,label=label),
        size=o$legend_text_size,color=o$secondary_text_color,
        family=s$figure$font_family)+
      ggplot2::annotate('text',x=o$legend_title_x,
        y=o$legend_title_y,
        label=if(is_risk)sprintf(o$labels$risk_legend,s$primary_horizon)
              else o$labels$capture_legend,
        hjust=0,fontface='bold',size=o$legend_title_text_size,
        lineheight=o$legend_title_lineheight,family=s$figure$font_family)+
      ggplot2::scale_fill_gradientn(colors=cols,limits=lim,
        oob=scales::squish,na.value=o$empty_fill,guide='none')+
      ggplot2::scale_color_identity(guide='none')+
      ggplot2::coord_cartesian(xlim=xlim,ylim=o$y_limits,clip='off',expand=FALSE)+
      ggplot2::theme_void(base_family=s$figure$font_family,
                          base_size=o$base_text_size_pt)+
      ggplot2::theme(plot.title=ggplot2::element_text(
        size=o$panel_title_size,face='bold',hjust=0,
        margin=ggplot2::margin(b=o$title_margin_pt)),
        plot.title.position='plot',
        plot.margin=do.call(ggplot2::margin,as.list(o$plot_margin_pt)))+
      ggplot2::labs(title=paste0(
        o$panel_letters[[type]],layout$common$panel_prefix_separator,
        if(is_risk)sprintf(o$labels$risk_panel,s$primary_horizon)
        else sprintf(o$labels$capture_panel,s$primary_horizon)))
    p
  }
  p <- patchwork::wrap_plots(draw_panel(a,'risk'),draw_panel(b,'capture'),
                             ncol=1,heights=o$panel_height_ratio)
  # 完整25格CSV，无论是否隐藏VH列、是否以破折号显示缺失CI，均保留真实上游值。
  list(plot=p,plot_data=list(participants_and_risk=a,
                              normalized_event_capture=b))
}

# 保留供追踪的数据 CSV，绝不把含舍入标注的版本作为下游分析来源。
rp14_export_data <- function(data,dir,name) {
  rp_dir(dir)
  for(k in names(data))rp_save_csv(data[[k]],file.path(dir,paste0(name,'_',k,'.csv')))
  invisible(TRUE)
}
# =============================================================================
# 补充图独立统计核查（仅 QC；不参与任何绘图或论文结果计算）
# 目的：对联合分类中小样本单元格，将上游逐插补 AJ 估计 Q/U 与独立软件
#       cmprsk::cuminc()（CIF及方差）和 survival::survfit()（多状态AJ CIF）比较。
# 读取：Module 07 完成插补后的分组数据 + Module 10 逐插补估计CSV。
# 约束：不调用上游 m09_aj()；不重做Rubin合并；不改上游任何文件。
# =============================================================================

# 标准实现1：cmprsk 的竞争风险累积发生函数及其方差（原始概率尺度）。
rp14_qc_cmprsk <- function(followup, status, horizon) {
  # group固定为“All”，避免小组内部再引入第二套分层定义。
  obj <- cmprsk::cuminc(ftime=followup, fstatus=status,
                        group=rep('All',length(status)), cencode=0)
  keys <- grep(' 1$', names(obj), value=TRUE) # 事件编码1=ASCVD
  if(length(keys)==0L) {
    if(!any(status==1L & followup<=horizon))
      return(list(risk=0, se=0, note='no ASCVD events by horizon'))
    stop('cmprsk未返回ASCVD原因1的CIF，检查结局编码')
  }
  if(length(keys)!=1L)stop('cmprsk返回多个ASCVD曲线，不能自动选择')
  z<-obj[[keys]]
  ix<-which(z$time<=horizon+1e-10)
  if(!length(ix))return(list(risk=0,se=0,note='no earlier CIF time'))
  j<-tail(ix,1L)
  v<-z$var[j]
  if(!is.finite(z$est[j]) || !is.finite(v) || v<0)
    stop('cmprsk风险或方差非有限/负值')
  list(risk=as.numeric(z$est[j]),se=sqrt(as.numeric(v)),note='ok')
}

# 标准实现2：survival 多状态Aalen–Johansen，读取ASCVD吸收状态的占据概率。
rp14_qc_survival <- function(followup,status,horizon) {
  dat<-data.frame(time=as.numeric(followup),
    outcome=factor(as.integer(status),levels=c(0L,1L,2L),
                   labels=c('censor','ASCVD','competing')))
  fit<-survival::survfit(survival::Surv(time,outcome)~1,data=dat)
  if(is.null(fit$time) || !length(fit$time))return(0)
  ps<-fit$pstate
  if(is.null(dim(ps)))ps<-matrix(ps,nrow=length(fit$time))
  # 不同survival版本的状态名称可能存在fit$states而非pstate列名中。
  if(is.null(colnames(ps)) && !is.null(fit$states) && length(fit$states)==ncol(ps))
    colnames(ps)<-as.character(fit$states)
  cols<-which(colnames(ps)=='ASCVD')
  if(length(cols)!=1L)stop('survival多状态结果无法确定ASCVD列')
  ix<-which(fit$time<=horizon+1e-10)
  if(!length(ix))return(0)
  as.numeric(ps[tail(ix,1L),cols])
}

# 重建Module 10的分组逻辑（仅用于交叉核查，不修改个体数据）。
rp14_qc_origin <- function(x,year,scenario) {
  base_col<-if(year=='2019')'score_2019_category' else 'score2_2025_category'
  # 兼容正式分析的旧代码SH/SHT及后续显示层SN/SNT命名。
  aliases<-switch(scenario,SN=c('SN','SH'),ST='ST',SNT=c('SNT','SHT'))
  upgrades<-paste0('upgrade_',year,'_',aliases)
  found<-upgrades[upgrades%in%names(x)]
  if(length(found)!=1L)stop('缺少或重复的RM升级列：',year,'/',scenario,
                            '；候选=',paste(upgrades,collapse=','))
  if(!base_col%in%names(x))stop('缺少原始风险分类：',base_col)
  base<-as.character(x[[base_col]])
  up<-as.integer(x[[found]])
  if(anyNA(base)||anyNA(up)||any(!base%in%c('L','M','H','VH'))||any(!up%in%0:1))
    stop('原始风险类别或RM升级标记缺失/不合法')
  if(any(up==1L & base!='M'))stop('发现非M人群被升级，不符合Module 10定义')
  ifelse(base=='M',ifelse(up==1L,'M_upgraded_H','M_not_upgraded'),base)
}

rp14_qc_sparse_crossvalidation <- function(analysis_dir, study_cfg, opt, scenarios, qc_dir) {
  if(!isTRUE(opt$enabled))return(invisible(NULL))
  if(!requireNamespace('cmprsk',quietly=TRUE))
    stop('小样本交叉核查需要标准R包cmprsk，请先install.packages("cmprsk")')
  if(!requireNamespace('survival',quietly=TRUE))
    stop('小样本交叉核查需要标准R包survival，请先install.packages("survival")')
  if(length(opt$max_cell_n)!=1L || !is.finite(opt$max_cell_n) || opt$max_cell_n<2)
    stop('max_cell_n必须为不小于2的整数')
  if(length(opt$point_tolerance)!=1L || !is.finite(opt$point_tolerance) || opt$point_tolerance<=0)
    stop('point_tolerance必须是正数（概率尺度）')
  if(length(opt$se_relative_warning)!=1L || !is.finite(opt$se_relative_warning) || opt$se_relative_warning<=0)
    stop('se_relative_warning必须是正数')
  horizon<-as.numeric(opt$horizon)
  if(length(horizon)!=1L || !horizon%in%study_cfg$times)
    stop('小样本QC时间点不属于Module07已计算的研究时间点')
  mcount<-as.integer(study_cfg$m)
  if(!length(mcount)||!is.finite(mcount)||mcount<2L)stop('无效多重插补次数')

  risk_file<-file.path(analysis_dir,'observed_risk_and_event_capture',
    'per_imputation','joint_origin_category_risks_and_event_capture.csv')
  if(!file.exists(risk_file))stop('缺少逐插补原始观察风险：',risk_file,
    '；请先完成Module10统一导出，不可用pooled CI反推个体方差')
  mi_risk<-data.table::fread(risk_file)
  rp_assert(mi_risk,c('imputation_id','scenario','scope','guideline','cell_id',
    'origin_2019','origin_2025','time','estimand','Q','U','cell_n','cell_status'),
    'Module10逐插补AJ结果')
  mi_risk<-mi_risk[scope=='overall' & guideline=='joint' &
    estimand=='joint_category_risk' & time==horizon]
  if(!nrow(mi_risk))stop('逐插补风险CSV缺少总体15年联合分类结果')
  rows<-list();ii<-0L
  lev<-c('L','M_not_upgraded','M_upgraded_H','H','VH')
  for(m in seq_len(mcount)) {
    path<-file.path(analysis_dir,'strategy_assignments',sprintf('m%02d.rds',m))
    if(!file.exists(path))stop('缺少Module07逐插补完整数据：',path)
    x<-data.table::as.data.table(readRDS(path))
    rp_assert(x,c('participant_id','followup_years','event_status'),paste('Module07 MI',m))
    if(anyDuplicated(x$participant_id)||anyNA(x$participant_id))
      stop('MI ',m,' participant_id重复或缺失')
    tt<-as.numeric(x$followup_years);st<-as.integer(x$event_status)
    if(any(!is.finite(tt)|tt<=0) || anyNA(st) || any(!st%in%0:2))
      stop('MI ',m,' 生存时间或竞争结局编码无效')
    for(k in scenarios){
      candidates<-switch(k,SN=c('SN','SH'),ST='ST',SNT=c('SNT','SHT'))
      d<-mi_risk[imputation_id==m & scenario%in%candidates]
      if(nrow(d)!=25L || anyDuplicated(d$cell_id) ||
         !setequal(as.character(d$cell_id),as.vector(outer(lev,lev,paste,sep='__'))))
        stop('MI ',m,'/',k,'：逐插补联合分类必须恰好包含25个唯一单元格')
      a<-rp14_qc_origin(x,'2019',k)
      b<-rp14_qc_origin(x,'2025',k)
      keys<-paste(a,b,sep='__')
      groups<-split(seq_along(keys),keys,drop=TRUE)
      # 逐格先核对人数；选取小样本n<=阈值（包含n=0/1用于审计）。
      for(j in seq_len(nrow(d))){
        rec<-d[j]
        idx<-groups[[as.character(rec$cell_id)]]
        if(is.null(idx))idx<-integer()
        n<-length(idx)
        if(n != as.integer(rec$cell_n))stop('MI ',m,'/',k,'/',rec$cell_id,
          '：Module07重建分类人数不等于Module10原始单元格人数')
        if(n>opt$max_cell_n)next
        ii<-ii+1L
        n1<-sum(st[idx]==1L & tt[idx]<=horizon)
        n2<-sum(st[idx]==2L & tt[idx]<=horizon)
        cmrisk<-cmse<-svrisk<-NA_real_
        msg<-''
        if(n<2L) {
          msg<-if(n==0L)'empty_cell' else 'n_below_2'
        }else{
          ref1<-tryCatch(rp14_qc_cmprsk(tt[idx],st[idx],horizon),
            error=function(e)e)
          ref2<-tryCatch(rp14_qc_survival(tt[idx],st[idx],horizon),
            error=function(e)e)
          if(inherits(ref1,'error')||inherits(ref2,'error')){
            msg<-paste(c(if(inherits(ref1,'error'))paste0('cmprsk: ',conditionMessage(ref1)),
                         if(inherits(ref2,'error'))paste0('survival: ',conditionMessage(ref2))),
                       collapse='; ')
          }else{
            cmrisk<-ref1$risk;cmse<-ref1$se;svrisk<-ref2
            msg<-if(n1==0L)'no_ASCVD_events_by_horizon' else 'compared'
          }
        }
        q<-as.numeric(rec$Q);u<-as.numeric(rec$U)
        upstream_se<-if(is.finite(u)&&u>=0)sqrt(u) else NA_real_
        dcm<-if(is.finite(cmrisk)&&is.finite(q))abs(cmrisk-q) else NA_real_
        dsv<-if(is.finite(svrisk)&&is.finite(q))abs(svrisk-q) else NA_real_
        # 软件之间的方差方法可能并非有限样本严格等价：SE差异只作预警。
        se_rel<-if(is.finite(cmse)&&is.finite(upstream_se))
          abs(cmse-upstream_se)/max(cmse,upstream_se,1e-10) else NA_real_
        rows[[ii]]<-data.table::data.table(
          scenario=k,imputation_id=m,time=horizon,
          origin_2019=as.character(rec$origin_2019),
          origin_2025=as.character(rec$origin_2025),cell_id=as.character(rec$cell_id),
          n=n,n_ascvd_by_horizon=n1,n_competing_death_by_horizon=n2,
          upstream_cell_status=as.character(rec$cell_status),
          upstream_risk=q,upstream_variance=u,upstream_se=upstream_se,
          cmprsk_risk=cmrisk,cmprsk_se=cmse,survival_risk=svrisk,
          abs_diff_cmprsk=dcm,abs_diff_survival=dsv,
          point_match=if(is.finite(dcm)&&is.finite(dsv))
            dcm<=opt$point_tolerance && dsv<=opt$point_tolerance else NA,
          se_relative_gap=se_rel,
          se_review=if(is.finite(se_rel))se_rel>opt$se_relative_warning else NA,
          qc_state=msg)
      }
    }
    rm(x);gc(FALSE)
    rp_msg('补充图稀疏单元格交叉验证：MI ',m,'/',mcount,' 已处理')
  }
  if(!length(rows))stop('没有选中小样本单元格；请调整 max_cell_n')
  all<-data.table::rbindlist(rows,use.names=TRUE,fill=TRUE)
  data.table::setorder(all,scenario,imputation_id,n,cell_id)
  summary<-all[,.(
    checked_cells=.N,
    software_comparisons=sum(!is.na(point_match)),
    point_matched=sum(point_match==TRUE,na.rm=TRUE),
    point_mismatched=sum(point_match==FALSE,na.rm=TRUE),
    reference_failed=sum(!is.na(qc_state)&grepl('cmprsk:|survival:',qc_state)),
    zero_ascvd_cells=sum(n>=2 & n_ascvd_by_horizon==0),
    insufficient_n_cells=sum(n<2),
    se_review_flags=sum(se_review==TRUE,na.rm=TRUE),
    max_abs_risk_gap_cmprsk=if(all(is.na(abs_diff_cmprsk)))NA_real_ else max(abs_diff_cmprsk,na.rm=TRUE),
    max_abs_risk_gap_survival=if(all(is.na(abs_diff_survival)))NA_real_ else max(abs_diff_survival,na.rm=TRUE)
  ),by=scenario]
  rp_save_csv(all,file.path(qc_dir,'14_supplementary_sparse_cell_crossvalidation.csv'))
  rp_save_csv(summary,file.path(qc_dir,'14_supplementary_sparse_cell_summary.csv'))
  info<-c(
    '目的：检验Module10逐插补Aalen-Johansen观察风险的点估计及标准误，不重做MI/Rubin合并。',
    sprintf('样本量标准：每个MI、每个SN/ST/SNT联合分类单元格 n <= %d；包括无事件/极小n作为审计行。',as.integer(opt$max_cell_n)),
    sprintf('时间点：%g年；风险单位：0–1（不是%%）。',horizon),
    sprintf('点估计绝对差阈值：%.8g；标准误相对差异预警阈值：%.3g。',opt$point_tolerance,opt$se_relative_warning),
    paste0('cmprsk版本：',as.character(utils::packageVersion('cmprsk'))),
    paste0('survival版本：',as.character(utils::packageVersion('survival'))),
    '核查软件：cmprsk::cuminc 用于竞争风险CIF/方差；survival::survfit 的多状态AJ用于CIF点估计。',
    'upstream_risk 与 upstream_variance 直接来自Module10 per_imputation原始Q和U。',
    'point_match要求同时与两种软件的CIF点估计一致。',
    'se_relative_gap仅用于标记差异；不同的有限样本方差实现不能仅据差异大小判断错误。',
    '风险为0、无事件、人数不足2或软件返回错误时保留审计记录，不伪造CI。',
    '此核查不改变主图、补充图、上游pooled或per_imputation文件。'
  )
  path<-file.path(qc_dir,'14_supplementary_sparse_cell_methodology.txt')
  rp_guard(path);writeLines(info,path,useBytes=TRUE)
  fail<-sum(all$point_match==FALSE,na.rm=TRUE)
  failed_ref<-sum(grepl('cmprsk:|survival:',all$qc_state))
  if(fail || failed_ref)
    warning('补充图小样本交叉验证发现 ',fail,' 个风险点估计不匹配、',failed_ref,
            ' 个标准软件无法比较的单元格；请检查QC CSV。')
  nse<-sum(all$se_review==TRUE,na.rm=TRUE)
  if(nse)warning('补充图小样本QC：',nse,' 个单元格的标准误相对差异超过阈值；请检查影响函数实现与随访/删失模式（不代表自动判定错误）')
  rp_msg('补充图小样本交叉验证完成：',nrow(all),' 条逐插补单元格；结果保存到 ',qc_dir)
  invisible(list(detail=all,summary=summary))
}

# Figure 4 contribution point estimates: form the weighted contrast within each MI.
rp14_subgroup_contributions <- function(per_mi, horizon=15L) {
  rd <- data.table::as.data.table(per_mi)[time==horizon & estimand=="rd_2025_vs_2019"]
  keys <- c("imputation_id","scenario")
  overall <- rd[scope=="overall", .(imputation_id,scenario,N=n_standardized,Q_overall=Q)]
  groups <- c(age_lt_60="Age",age_ge_60="Age",female="Sex",male="Sex",
              high_ses="SES",low_ses="SES",no_hypertension="Hypertension",
              hypertension="Hypertension")
  sub <- merge(rd[scope %in% names(groups)],overall,by=keys)
  if (nrow(sub)!=8L*nrow(overall) || !nrow(sub) ||
      anyDuplicated(sub,by=c(keys,"scope")) ||
      any(!is.finite(sub$n_standardized) | !is.finite(sub$N) | sub$N<=0))
    stop("Figure 4 incomplete subgroup contribution inputs")
  sub[, `:=`(dimension=unname(groups[scope]), contribution_pp=100*n_standardized/N*Q)]
  check <- sub[, .(sum_n=sum(n_standardized),N=N[1],residual_pp=sum(contribution_pp)-100*Q_overall[1]),
               by=.(imputation_id,scenario,dimension)]
  if (any(check$sum_n!=check$N | abs(check$residual_pp)>1e-6))
    stop("Figure 4 subgroup contributions do not reproduce overall RD")
  list(summary=sub[, .(n=mean(n_standardized),contribution_pp=mean(contribution_pp)),
                   by=.(scope,scenario,dimension)],qa=check)
}



rp14_fig3 <- function(src,s,layout) {
  a <- rp14_fig3_legacy(src,s,layout,classification_only=TRUE)
  groups <- s$transition_order
  labels <- setNames(c("Unchanged untreated","New intervention","Intervention removed",
                      "Treatment intensified","Unchanged intensity","Treatment de-intensified"),groups)
  palette <- c(S0="#888888",SN="#2166AC",ST="#E66101",SNT="#1B7837")
  dodge <- ggplot2::position_dodge(width=.65)
  group_axis <- function() ggplot2::scale_y_discrete(labels=labels)
  theme <- ggplot2::theme_bw(base_size=9)
  b <- data.table::copy(src$fixed_moves)
  b[, group:=factor(row_group,rev(groups))]
  c <- data.table::copy(src$fixed_ldl[fixed_group!=groups[1]])
  c[, group:=factor(fixed_group,rev(groups[-1]))]
  d <- data.table::copy(src$all_effects[scope=="overall" & time==15 & estimand=="rd_2025_vs_2019"])
  e <- data.table::copy(src$fixed_rd[time==15 & fixed_group %in% c("G2","G3","G4","G6")])
  e[, group:=factor(groups[fixed_order],rev(groups[c(2,3,4,6)]))]
  stopifnot(nrow(b)==18L,nrow(c)==30L,nrow(d)==4L,nrow(e)==16L)
  bplot <- ggplot2::ggplot(b,ggplot2::aes(pct_moved,group,fill=scenario))+
    ggplot2::geom_col(position=dodge)+group_axis()+ggplot2::scale_fill_manual(values=palette)+
    ggplot2::labs(title="B. Treatment-decision reclassification",x="Reclassified within S0 group (%)",y=NULL,fill=NULL)+theme
  cplot <- ggplot2::ggplot(c,ggplot2::aes(mean_delta_mmol_l,group,fill=scenario))+
    ggplot2::geom_col(position=dodge)+ggplot2::facet_wrap(~guideline)+group_axis()+
    ggplot2::scale_fill_manual(values=palette)+ggplot2::labs(title="C. Change in required LDL-C reduction",
      x="RM minus S0 (mmol/L)",y=NULL,fill=NULL)+theme
  forest <- function(data, y, title) {
    data[, `:=`(effect_pp=100*estimate,low_pp=100*lcl,high_pp=100*ucl)]
    ggplot2::ggplot(data,ggplot2::aes(effect_pp,.data[[y]],colour=scenario))+
      ggplot2::geom_vline(xintercept=0,linetype=2)+
      ggplot2::geom_errorbar(ggplot2::aes(xmin=low_pp,xmax=high_pp),orientation="y",width=.14,position=dodge)+
      ggplot2::geom_point(position=dodge)+ggplot2::scale_colour_manual(values=palette)+
      ggplot2::labs(title=title,x="2025 minus 2019 (pp)",y=NULL,colour=NULL)+theme
  }
  d[, strategy:=factor(scenario,rev(names(palette)))]
  bottom <- patchwork::wrap_plots(cplot,forest(d,"strategy","D. Overall 15-year RD"),
    forest(e,"group","E. Fixed S0 group 15-year RD")+group_axis(),nrow=1,widths=c(1.35,.8,1.2))
  top <- patchwork::wrap_plots(a$plot,bplot,nrow=1)
  list(plot=patchwork::wrap_plots(top,bottom,ncol=1,guides="collect") &
         ggplot2::theme(legend.position="bottom"),
       plot_data=list(classification=a$data,reclassification=b,ldl_change=c,overall_rd=d,fixed_group_rd=e))
}

rp14_fig4 <- function(src,s,layout) {
  scopes <- c("age_lt_60","age_ge_60","female","male","high_ses","low_ses","no_hypertension","hypertension")
  labels <- setNames(c("<60 years","≥60 years","Female","Male","High SES","Low SES","No HTN","HTN"),scopes)
  palette <- c(S0="#888888",SN="#2166AC",ST="#E66101",SNT="#1B7837")
  contribution <- rp14_subgroup_contributions(src$all_effects_mi)
  effects <- src$all_effects[time==15 & scope %in% scopes]
  get_effect <- function(metric,prefix) {
    x <- effects[estimand==metric, .(scope,scenario,estimate=100*estimate,lcl=100*lcl,ucl=100*ucl)]
    if(nrow(x)!=32L || anyDuplicated(x,by=c("scope","scenario"))) stop("Incomplete Figure 4 ",metric)
    data.table::setnames(x,c("estimate","lcl","ucl"),paste0(prefix,c("_est","_lo","_hi")))
    x
  }
  data <- Reduce(function(a,b) merge(a,b,by=c("scope","scenario")),
    list(contribution$summary,get_effect("rd_2025_vs_natural","nc25"),
         get_effect("rd_2019_vs_natural","nc19"),get_effect("rd_2025_vs_2019","rd")))
  stopifnot(nrow(data)==32L)
  data[, y:=37-((match(scope,scopes)-1)*4+match(scenario,names(palette)))-
         (match(scope,scopes)-1)%/%2]
  data[, label:=ifelse(scenario=="S0",paste(labels[scope],scenario,sep=" · "),scenario)]
  data[, n_label:=ifelse(scenario=="S0",formatC(as.integer(round(n)),format="d",big.mark=","),"")]
  for(prefix in c("nc25","nc19","rd"))
    data[, (prefix):=sprintf("%.2f (%.2f, %.2f)",get(paste0(prefix,"_est")),
                             get(paste0(prefix,"_lo")),get(paste0(prefix,"_hi")))]
  data[, contribution:=sprintf("%.2f",contribution_pp)]
  y_limits <- c(min(data$y)-.6,max(data$y)+.6)
  column <- function(field,title) {
    ggplot2::ggplot(data,ggplot2::aes(x=0,y=y,label=.data[[field]]))+
      ggplot2::geom_text(size=2.6,hjust=0)+ggplot2::scale_y_continuous(limits=y_limits,expand=c(0,0))+
      ggplot2::coord_cartesian(xlim=c(0,1),clip="off")+ggplot2::labs(title=title)+
      ggplot2::theme_void()+ggplot2::theme(plot.title=ggplot2::element_text(size=9,face="bold"),
        axis.text.x=ggplot2::element_text(colour="transparent",size=8),axis.title.x=ggplot2::element_blank())
  }
  forest <- ggplot2::ggplot(data,ggplot2::aes(rd_est,y,colour=scenario))+
    ggplot2::geom_vline(xintercept=0,linetype=2)+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=rd_lo,xmax=rd_hi),width=.2,orientation="y")+
    ggplot2::geom_point(size=1.5)+ggplot2::scale_colour_manual(values=palette,guide="none")+
    ggplot2::scale_y_continuous(limits=y_limits,expand=c(0,0),breaks=NULL)+
    ggplot2::labs(title="Forest plot",x="Favors 2025   |   Favors 2019",y=NULL)+ggplot2::theme_bw(base_size=8)
  plot <- patchwork::wrap_plots(column("label","Subgroup"),column("n_label","N"),
    column("nc25","2025 − NC (pp)"),column("nc19","2019 − NC (pp)"),
    column("rd","2025 − 2019 (pp)"),forest,column("contribution","Contribution (pp)"),
    nrow=1,widths=c(1.35,.65,1.35,1.35,1.45,2.15,1.15))
  list(plot=plot,plot_data=list(subgroup_risk_differences=data,contribution_qa=contribution$qa))
}


