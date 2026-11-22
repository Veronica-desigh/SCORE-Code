# =============================================================================
# M01 | RM 与 S0 置信区间审计 v2.2（修复Q/U映射并增加合并前自检；独立诊断，基于 2026-11-10 提供的 07–12 主/函数脚本）
# =============================================================================
# 目的：解释 2025–2019 风险差在应用 RM 后部分 95% CI 极窄的原因。
#       核查：配对协方差、5 类解析 IF 来源、MI 内/间方差、治疗路径趋同、
#             09/11 已导出结果、RM 与 S0 增量差值的完整四状态配对方差。
# 只读输入：07 strategy_assignments/mXX.rds，08 risk_prediction/bundle_mXX.rds，
#           study_config.rds；09/11 导出 CSV 如存在，则执行一致性验证。
# 默认只读缓存。可选独立验证在内存重拟合，不修改07–12结果。
# M01_AUDIT_GRADIENT_DIRECTIONS=2：中心差分核对完整IF（两个步长）。
# M01_AUDIT_BOOT_REPS=20：配对重抽样试跑；正式比较建议增加重复次数。
# M01_AUDIT_REFIT_MI=1；M01_AUDIT_SPLINE=1 同批比较两原因别模型的LDL样条。
# 可选验证为选定MI、最长时点、S0/SNT的插补内审计，不估计最终MI区间覆盖率。
# 所有诊断仅写入：analysis_results/rm_ci_variance_audit_v2_2/。
#
# 使用方法（优先）：source(file.choose(), encoding="UTF-8")，选择本文件。
# 如果代码和结果仍在原 07–12 主脚本所示 D 盘默认目录，无需设置环境变量。
# 如项目搬迁，可选手动设置：
# Sys.setenv(M01_PROJECT_ROOT="D:/.../04_Code")
# Sys.setenv(M01_CODE_DIR="D:/.../LJW_code/Main_code/M01_ASCVD")
# Sys.setenv(M01_ANALYSIS_OUTPUT="D:/.../LJW_outcome/data_intervention/M01_ASCVD/analysis_results")
# Sys.setenv(M01_AUDIT_SCOPES="overall,age_ge_60,low_ses")
#
# 情景：内部 S0 / SH / ST / SHT；论文对应 S0 / SN / ST / SNT。
# 单位：计算时为概率（0–1），写作字段 *_pp 时才换算为百分点。
# 重要限制：这属于当前固定 RM 决策规则/结局模型下的抽样及 MI 推断诊断，
#           不包含模型错设、RM 定义结构性不确定性、随机干预的 Monte Carlo 误差。
# =============================================================================

options(stringsAsFactors = FALSE)
audit_time_start <- Sys.time()
.audit_norm <- function(p) {
  if (!length(p) || is.na(p[1L]) || !nzchar(p[1L])) return("")
  normalizePath(p[1L], winslash="/", mustWork=FALSE)
}
.audit_nonempty <- function(x) unique(x[!is.na(x) & nzchar(trimws(x))])
.audit_env <- function(name) Sys.getenv(name, unset="")
.audit_exist <- function(d,f) file.exists(file.path(d,f))
.audit_label <- c(S0="S0",SH="SN",ST="ST",SHT="SNT")
.audit_required_fns <- paste0("Function_", c("07_Guideline_strategies", "08_Counterfactual_risk",
                                          "09_Main_analysis"), "_M01.R")
.audit_root_default <- "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code"
.audit_this_file <- function() {
  f <- unlist(lapply(sys.frames(),function(x) {
    v <- x$ofile
    if (!is.null(v) && length(v)==1L && file.exists(v)) v else character()
  }),use.names=FALSE)
  if(length(f)) tail(f,1L) else ""
}
.audit_script_path <- .audit_this_file()
.audit_script_dir <- if(nzchar(.audit_script_path)) dirname(.audit_norm(.audit_script_path)) else getwd()
.audit_work_dir <- .audit_norm(getwd())
.audit_root_from_code <- function(d){
  if(!nzchar(d)) return(character())
  x<-gsub("\\\\","/",sub("/+$","",d))
  suffix <- "/LJW_code/Main_code/M01_ASCVD"
  if(endsWith(tolower(x),tolower(suffix))) substr(x,1L,nchar(x)-nchar(suffix)) else character()
}
.audit_parse_root <- function(file) {
  if (!file.exists(file)) return(character())
  x <- tryCatch(readLines(file,warn=FALSE,encoding="UTF-8"),error=function(e) character())
  # 匹配 07–12 主脚本实际的 Sys.getenv('M01_PROJECT_ROOT','D:/...')
  rex <- "Sys[.]getenv[[:space:]]*\\([[:space:]]*['\"]M01_PROJECT_ROOT['\"][[:space:]]*,[[:space:]]*['\"]([^'\"]+)['\"]"
  a <- regmatches(x,regexec(rex,x,perl=TRUE))
  .audit_nonempty(vapply(a,function(z)if(length(z)>1L)z[2L]else "",character(1)))
}
.audit_near_dirs <- .audit_nonempty(c(.audit_script_dir,.audit_work_dir,
  dirname(.audit_script_dir),dirname(.audit_work_dir),
  file.path(.audit_script_dir,"LJW_code/Main_code/M01_ASCVD"),
  file.path(.audit_work_dir,"LJW_code/Main_code/M01_ASCVD")))
.audit_main_names <- paste0(c("07_Guideline_strategies", "08_Counterfactual_risk",
                              "09_Main_analysis", "10_RM_classification",
                              "11_RM_effects", "12_Treatment_transitions"),"_M01.R")
.audit_main_paths <- unique(unlist(lapply(.audit_near_dirs,function(d)
  file.path(d,.audit_main_names)),use.names=FALSE))
.audit_main_paths <- .audit_main_paths[file.exists(.audit_main_paths)]
.audit_main_roots <- .audit_nonempty(unlist(lapply(.audit_main_paths,.audit_parse_root),use.names=FALSE))
.audit_root_candidates <- .audit_nonempty(c(
  .audit_env("M01_PROJECT_ROOT"),.audit_root_from_code(.audit_env("M01_CODE_DIR")),
  .audit_root_from_code(.audit_script_dir),.audit_root_from_code(.audit_work_dir),
  .audit_main_roots,.audit_root_default))
.audit_code_dirs <- .audit_nonempty(c(.audit_env("M01_CODE_DIR"),
  .audit_script_dir,.audit_work_dir,.audit_near_dirs,
  unlist(lapply(.audit_root_candidates,function(r)file.path(r,"LJW_code/Main_code/M01_ASCVD")),use.names=FALSE)))
.audit_code_table <- data.frame(
  path=.audit_code_dirs,
  folder_exists=vapply(.audit_code_dirs,dir.exists,logical(1)),
  n_functions=vapply(.audit_code_dirs,function(x)sum(file.exists(file.path(x,.audit_required_fns))),integer(1)),
  stringsAsFactors=FALSE)
.audit_good_code <- .audit_code_table$path[.audit_code_table$n_functions==length(.audit_required_fns)]

cat("\n================ M01 RM CI 审计 v2.2 | 路径检查 ================\n")
cat("当前工作目录：", .audit_work_dir,"\n审计脚本目录：",.audit_script_dir,"\n",sep="")
cat("M01_PROJECT_ROOT：", if(nzchar(.audit_env("M01_PROJECT_ROOT"))) .audit_env("M01_PROJECT_ROOT") else "<未设置>","\n",sep="")
cat("M01_CODE_DIR：",if(nzchar(.audit_env("M01_CODE_DIR"))) .audit_env("M01_CODE_DIR") else "<未设置>","\n",sep="")
cat("M01_ANALYSIS_OUTPUT：",if(nzchar(.audit_env("M01_ANALYSIS_OUTPUT"))) .audit_env("M01_ANALYSIS_OUTPUT") else "<未设置>","\n",sep="")

.audit_output_dirs <- .audit_nonempty(c(.audit_env("M01_ANALYSIS_OUTPUT"),
  unlist(lapply(.audit_root_candidates,function(r)
    file.path(r,"LJW_outcome/data_intervention/M01_ASCVD/analysis_results")),use.names=FALSE),
  file.path(.audit_script_dir,"analysis_results"),
  file.path(.audit_work_dir,"analysis_results"),
  file.path(dirname(.audit_script_dir),"analysis_results")))
.audit_out_status <- list(); .audit_good_out <- ""; .audit_config <- NULL
for (d in .audit_output_dirs) {
  cp<-file.path(d,"study_config.rds")
  r<-list(path=d,folder=dir.exists(d),config=file.exists(cp),m=NA_integer_,
          missing07="",missing08="",reason="")
  if(r$config){
    cfgx<-tryCatch(readRDS(cp),error=function(e)e)
    if (inherits(cfgx, "error")) {
      r$reason <- paste0("无法读取配置：", conditionMessage(cfgx))
    } else if (is.null(cfgx$m) || length(cfgx$m)!=1L ||
               !is.finite(suppressWarnings(as.numeric(cfgx$m))) ||
               as.integer(cfgx$m)<2L) {
      r$reason <- "study_config.rds 中 m 无效"
    } else {
      r$m<-as.integer(cfgx$m)
      ix<-seq_len(r$m)
      f07<-file.path(d,"strategy_assignments",sprintf("m%02d.rds",ix))
      f08<-file.path(d,"risk_prediction",sprintf("bundle_m%02d.rds",ix))
      r$missing07<-paste(basename(f07[!file.exists(f07)]),collapse=", ")
      r$missing08<-paste(basename(f08[!file.exists(f08)]),collapse=", ")
      if(!nzchar(r$missing07)&&!nzchar(r$missing08)) {
        if(!nzchar(.audit_good_out)){
          .audit_good_out<-d; .audit_config<-cfgx
        }
      }else r$reason<-"至少一套 MI 缺少07/08产物"
    }
  } else r$reason<-"study_config.rds 不存在"
  .audit_out_status[[length(.audit_out_status)+1L]] <- r
}
.audit_out_table<-do.call(rbind,lapply(.audit_out_status,function(r)
  data.frame(path=r$path,folder_exists=r$folder,has_config=r$config,
             mi_m=r$m,missing_07=r$missing07,missing_08=r$missing08,
             problem=r$reason,stringsAsFactors=FALSE)))
.audit_diagnostic <- function() {
  cat("\n[详细诊断：代码目录候选]\n")
  for(i in seq_len(nrow(.audit_code_table))) {
    d<-.audit_code_table$path[i]
    cat(sprintf("  %02d | %s | 目录=%s | 07/08/09函数=%d/3\n",i,d,
      if(dir.exists(d))"存在"else"不存在",.audit_code_table$n_functions[i]))
    if(dir.exists(d)) for(f in .audit_required_fns)
      cat(if(file.exists(file.path(d,f)))"     [有] "else"     [缺] ",f,"\n",sep="")
  }
  cat("\n[详细诊断：结果目录候选]\n")
  for(i in seq_len(nrow(.audit_out_table))) {
    r<-.audit_out_table[i,]
    cat(sprintf("  %02d | %s | 目录=%s | config=%s | MI=%s\n",i,r$path,
      r$folder_exists,r$has_config,if(is.na(r$mi_m))"?"else as.character(r$mi_m)))
    if(nzchar(r$missing_07))cat("       缺失07：",r$missing_07,"\n",sep="")
    if(nzchar(r$missing_08))cat("       缺失08：",r$missing_08,"\n",sep="")
    if(nzchar(r$problem))cat("       原因：",r$problem,"\n",sep="")
  }
  cat("\n[处理方法]\n")
  cat("1. 如果文件在原 D 盘目录，请确认原来的 07 和 08 已完整运行。\n")
  cat('2. 项目迁移后执行：Sys.setenv(M01_PROJECT_ROOT="D:/你的项目/04_Code")\n')
  cat('3. 或分别执行：Sys.setenv(M01_CODE_DIR="含 Function_07~09 的目录")\n')
  cat('               Sys.setenv(M01_ANALYSIS_OUTPUT="含 study_config.rds 的 analysis_results 目录")\n')
  cat('4. 可用 file.choose() 选中 Function_09：\n')
  cat('   Sys.setenv(M01_CODE_DIR=dirname(normalizePath(file.choose(),winslash="/")))\n')
  cat("5. 如果只上传了 R 脚本、没有07/08的真实RDS数据，无法执行统计审计。\n")
}
if(!length(.audit_good_code) || !nzchar(.audit_good_out)) {
  .audit_diagnostic()
  stop("未定位到完整的原07–09函数与07/08 MI结果；详细原因见上方。",call.=FALSE)
}
.audit_code_dir <- .audit_norm(.audit_good_code[1L]); .audit_out_dir <- .audit_norm(.audit_good_out)
.audit_root_auto <- .audit_root_from_code(.audit_code_dir)
.audit_root <- if(length(.audit_root_auto)) .audit_root_auto[1L] else
  if(nzchar(.audit_env("M01_PROJECT_ROOT"))) .audit_env("M01_PROJECT_ROOT") else .audit_root_default
cat("[定位成功] 原函数目录：",.audit_code_dir,"\n[定位成功] 分析结果目录：",.audit_out_dir,
    "\n[定位成功] MI套数：",.audit_config$m,"\n",sep="")
if(nzchar(.audit_env("M01_CODE_DIR")) && .audit_norm(.audit_env("M01_CODE_DIR"))!=.audit_code_dir)
  cat("[提示] 先前M01_CODE_DIR不可用，已自动改用有效路径。\n")
if(nzchar(.audit_env("M01_ANALYSIS_OUTPUT")) && .audit_norm(.audit_env("M01_ANALYSIS_OUTPUT"))!=.audit_out_dir)
  cat("[提示] 先前M01_ANALYSIS_OUTPUT不可用，已自动改用有效路径。\n")

# 当前R会话的环境变量调整用于加载原函数，原脚本和已有结果不会被写入。
Sys.setenv(M01_CODE_DIR=.audit_code_dir,M01_ANALYSIS_OUTPUT=.audit_out_dir,
           M01_PROJECT_ROOT=.audit_root)
if(!requireNamespace("data.table",quietly=TRUE))stop("请先安装 data.table：install.packages('data.table')")
if(!requireNamespace("Matrix",quietly=TRUE))stop("请先安装 Matrix：install.packages('Matrix')")
if(!requireNamespace("survival",quietly=TRUE))stop("请先安装 survival：install.packages('survival')")
suppressPackageStartupMessages(source(file.path(.audit_code_dir,"Function_09_Main_analysis_M01.R"),
                               local=FALSE,encoding="UTF-8"))
.audit_cfg <- m07_config()
.audit_cfg$out_dir <- .audit_out_dir
.audit_cfg$code_dir <- .audit_code_dir
if (!identical(as.integer(.audit_cfg$m),as.integer(.audit_config$m)))
  stop("主函数 m07_config 与自动检查到的 MI 套数不一致")
.audit_scope_text <- .audit_env("M01_AUDIT_SCOPES")
.audit_scopes <- if (nzchar(.audit_scope_text)) {
  .audit_nonempty(trimws(strsplit(.audit_scope_text, ",", fixed=TRUE)[[1L]]))
} else {
  c("overall", "age_ge_60")
}
.audit_scenarios <- m07_rms()
if(!identical(.audit_scenarios,c("S0","SH","ST","SHT")))
  stop("当前07内部情景已改变，不应继续使用原审计定义：",paste(.audit_scenarios,collapse=","))
.audit_components<-c("beta_ASCVD","beta_competing_death", "baseline_ASCVD",
                     "baseline_competing_death","empirical_standardization")
.audit_tolerance <- 1e-8

# 注意：审计不在源文件中调用 run、source(主脚本)、m08_predict_mi() 或重建缓存。
.audit_comparison <- function(a,b,relative_tolerance=1e-6,absolute_tolerance=1e-13){
  if(length(a)!=1L||length(b)!=1L||!is.finite(a)||!is.finite(b))return(FALSE)
  abs(a-b)<=absolute_tolerance+relative_tolerance*max(abs(a),abs(b))
}
.audit_component_IF <- function(b,state,key){
  eng<-b$engine
  s<-b$states[[state]]
  kk<-match(key,s$keys)
  if(is.na(kk))stop("08缓存缺少 target：",state,"/",key)
  ix<-b$targets$index[[key]]
  if(!length(ix))stop("scope为空：",key)
  n<-length(ix);N<-eng$n;H<-length(eng$times)
  ans<-setNames(lapply(seq_along(.audit_components),function(k)matrix(0,N,H)),.audit_components)
  for(h in seq_len(H)) {
    ans[[1]][,h]<-as.numeric(eng$source1$iid_beta%*%s$B1[kk,,h])
    ans[[2]][,h]<-as.numeric(eng$source2$iid_beta%*%s$B2[kk,,h])
    ans[[3]][,h]<-m01_weighted_baseline_iid(eng$source1,s$G1[kk,,h])
    ans[[4]][,h]<-m01_weighted_baseline_iid(eng$source2,s$G2[kk,,h])
    ans[[5]][ix,h]<-(s$risk[ix,h]-s$estimate[kk,h])/n
  }
  pred<-m08_prediction(b,state,key)
  recon<-Reduce(`+`,ans)
  err<-max(abs(recon-pred$iid))
  if(!is.finite(err)||err>.audit_tolerance)
    stop("原08的 IF 无法由五项解析分量重构：",state," / ",key,"；max=",signif(err,5))
  list(parts=ans,pred=pred,recon_error=err)
}
.audit_crossmatrix <- function(p25,p19,h){
  mat<-do.call(cbind,lapply(.audit_components,function(v)
    p25$parts[[v]][,h]-p19$parts[[v]][,h]))
  colnames(mat)<-.audit_components
  crossprod(mat)
}
.audit_record_incremental <- function(ref,cur,idx,scenario,scope,time,dfcom){
  IF<-cur$iid-ref$iid
  Q<-cur$Q-ref$Q
  data.table::data.table(imputation_id=idx,scope=scope,scenario_internal=scenario,
    scenario_paper=.audit_label[[scenario]],time=time,comparison="RM_minus_S0_RD",
    Q=Q,U=sum(IF^2),dfcom=dfcom,
    covariance_RD_RM_S0=sum(cur$iid*ref$iid),
    corr_RD_RM_S0=if(sum(cur$iid^2)*sum(ref$iid^2)>0)
      sum(cur$iid*ref$iid)/sqrt(sum(cur$iid^2)*sum(ref$iid^2))else NA_real_)
}
# 明确列名映射：perMI 保存 Q_RD/U_RD；增量表 incMI 保存 Q/U。
# 不再在 data.table 的表达式里隐式引用不存在的 Q、U。
.audit_pooled <- function(tab,scope_cols,q_col="Q",u_col="U") {
  d <- data.table::as.data.table(tab)
  if (!nrow(d)) return(data.table::data.table())
  needed <- unique(c(scope_cols,"imputation_id",q_col,u_col,"dfcom"))
  missing <- setdiff(needed,names(d))
  if (length(missing)) {
    stop("Rubin合并缺少字段：",paste(missing,collapse=", "),
         "；现有列：",paste(names(d),collapse=", "),call.=FALSE)
  }
  if (anyDuplicated(d,by=c(scope_cols,"imputation_id"))) {
    stop("Rubin合并遇到重复的scope/scenario/time/MI键，不能合并。",call.=FALSE)
  }
  d[, {
    qq <- as.numeric(.SD[[q_col]])
    uu <- as.numeric(.SD[[u_col]])
    dd <- as.numeric(.SD[["dfcom"]])
    mm <- as.integer(.SD[["imputation_id"]])
    if (length(mm)!=.audit_cfg$m || anyNA(mm) ||
        !setequal(mm,seq_len(.audit_cfg$m))) {
      stop("Rubin合并遇到不完整/不匹配的MI编号：",
           paste(mm,collapse=","),"；预期1–",.audit_cfg$m,call.=FALSE)
    }
    if (any(!is.finite(qq)) || any(!is.finite(uu)) || any(uu<0) ||
        any(!is.finite(dd)) || any(dd<=0)) {
      stop("Rubin合并的Q、U或dfcom含缺失/非法值。",call.=FALSE)
    }
    z <- m09_rubin(qq,uu,min(dd),.audit_cfg$confidence_level)
    list(Q_pooled=z$estimate,estimate_pp=100*z$estimate,
      se_pp=100*z$se,CI_lower_pp=100*z$lcl,CI_upper_pp=100*z$ucl,
      within_variance=z$within_variance,between_variance=z$between_variance,
      total_variance=z$total_variance,
      fraction_between=if(z$total_variance>0)
        (1+1/z$M)*z$between_variance/z$total_variance else NA_real_,
      t_df=z$df,M=z$M,
      CI_rounds_equal_at_2dp=(round(100*z$lcl,2)==round(100*z$ucl,2)))
  }, by=scope_cols, .SDcols=c("imputation_id",q_col,u_col,"dfcom")]
}

# 在加载大型07/08 RDS前，用极小模拟数据独立检验两种列名布局。
.audit_pooled_selftest <- function() {
  mm <- seq_len(.audit_cfg$m)
  q <- 0.001 + mm*1e-6
  u <- rep(1e-7,length(mm))
  d <- data.table::data.table(imputation_id=mm,scope="selftest",time=.audit_cfg$times[1L],
                              Q_RD=q,U_RD=u,dfcom=rep(1000,length(mm)))
  expect <- m09_rubin(q,u,1000,.audit_cfg$confidence_level)
  got_rd <- .audit_pooled(d,c("scope","time"),q_col="Q_RD",u_col="U_RD")
  d[,`:=`(Q=Q_RD,U=U_RD)]
  got_delta <- .audit_pooled(d,c("scope","time"),q_col="Q",u_col="U")
  for (ans in list(got_rd,got_delta)) {
    if(nrow(ans)!=1L ||
       !isTRUE(all.equal(ans$Q_pooled,expect$estimate,tolerance=1e-12)) ||
       !isTRUE(all.equal(ans$total_variance,expect$total_variance,tolerance=1e-12))) {
      stop("Rubin合并列映射自检失败：不应继续读取分析RDS。",call.=FALSE)
    }
  }
  cat("[自检通过] Rubin合并已验证 Q_RD/U_RD 与 Q/U 两种字段结构。\n")
  invisible(TRUE)
}
.audit_pooled_selftest()

.audit_integer <- function(name, default, minimum=0L) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset=as.character(default))))
  if (!is.finite(value) || value < minimum || value > .Machine$integer.max || value != trunc(value))
    stop('Invalid integer setting: ', name)
  as.integer(value)
}

# Independent point predictions. Use survival's weighted Cox fits, basehaz and
# predict methods, then integrate on the same grid as 08 without its IF code.
.audit_refit_predict <- function(x, w, cuts, times, scopes, scenarios, chunk, spline=FALSE) {
  keep <- which(w > 0)
  fit_data <- x[keep]
  fit_data[, audit_weight := w[keep]]
  fit <- function(event, ldl) {
    terms <- sub('^ldl_ascvd$', ldl, m08_outcome_terms())
    if (spline) {
      basis <- splines::ns(x[[ldl]],df=3)
      ldl_knots <- attr(basis,'knots')
      ldl_boundary <- attr(basis,'Boundary.knots')
      terms[terms==ldl] <- paste0('splines::ns(',ldl,', knots=ldl_knots, Boundary.knots=ldl_boundary)')
    }
    formula <- reformulate(terms, paste0('survival::Surv(followup_years, ',event,')'))
    withCallingHandlers(
      survival::coxph(formula, data=fit_data, weights=audit_weight, ties='breslow',
                      robust=FALSE, singular.ok=FALSE, x=TRUE, y=TRUE),
      warning=function(c) stop(conditionMessage(c),call.=FALSE))
  }
  asc <- fit('ascvd_event','ldl_ascvd')
  death <- fit('competing_death_event','ldl_death')
  cumulative <- function(model) {
    h <- survival::basehaz(model, centered=TRUE)
    # Match 08’s tolerance for event times just above an interior grid cut.
    # At the last cut, 08 excludes events strictly after the analysis horizon.
    c(0,h$hazard)[1L+findInterval(pmin(cuts+1e-12,max(cuts)),h$time)]
  }
  A <- cumulative(asc); D <- cumulative(death)
  dA <- diff(c(0,A)); prevA <- head(c(0,A),-1L); prevD <- head(c(0,D),-1L)
  horizon <- match(times,cuts)
  if(anyNA(horizon)) stop('Audit horizons missing from prediction grid')
  scope_matrix <- vapply(scopes,function(i) w[keep]*as.numeric(i[keep]),numeric(length(keep)))
  denom <- colSums(scope_matrix)
  if(any(denom<=0)) stop('An audited subgroup is empty in this refit')
  death_hr <- exp(predict(death,newdata=fit_data,type='lp',reference='sample'))
  result <- list()
  for (k in scenarios) {
    sums <- array(0,c(length(scopes),length(times),2L))
    for (start in seq.int(1L,length(keep),by=chunk)) {
      ix <- seq.int(start,min(start+chunk-1L,length(keep)))
      nd <- data.table::copy(fit_data[ix])
      lp <- lapply(c('2019','2025'),function(year) {
        nd[, ldl_ascvd := x[[paste0('Lcf_',year,'_',k)]][keep[ix]]]
        exp(predict(asc,newdata=nd,type='lp',reference='sample'))
      })
      risk <- matrix(0,length(ix),2L)
      for (j in seq_len(max(horizon))) {
        for (y in 1:2) risk[,y] <- risk[,y]+dA[j]*lp[[y]]*exp(-prevA[j]*lp[[y]]-prevD[j]*death_hr[ix])
        h <- match(j,horizon)
        if(!is.na(h)) sums[,h,] <- sums[,h,]+crossprod(scope_matrix[ix,,drop=FALSE],risk)
      }
    }
    for (h in seq_along(times)) result[[length(result)+1L]] <- data.table::data.table(
      scope=names(scopes),time=times[h],scenario_internal=k,
      risk19=sums[,h,1L]/denom,risk25=sums[,h,2L]/denom)
  }
  out <- data.table::rbindlist(result)
  if(any(!is.finite(out$risk19) | !is.finite(out$risk25) |
         out$risk19<0 | out$risk19>1 | out$risk25<0 | out$risk25>1))
    stop('Independent grid CIF outside [0,1]')
  out[, RD := risk25-risk19]
  out
}

.audit_independent <- function(x, bundle, m, cfg, scopes, scenarios, reps, directions, spline) {
  x <- m08_prepare_model(x)
  n <- nrow(x)
  keys <- c('scope','time','scenario_internal')
  times <- max(cfg$times)
  h <- match(times,bundle$engine$times)
  groups <- m07_scope_spec(x)[scopes]
  evaluate <- function(w, alternative=FALSE) .audit_refit_predict(
    x,w,bundle$engine$cuts,times,groups,scenarios,cfg$prediction_chunk,alternative)
  cat(sprintf("[独立审计] MI %d：重新拟合原模型并核对风险预测\n",m))
  baseline <- evaluate(rep(1,n))
  gradients <- lapply(seq_len(nrow(baseline)),function(i) {
    row <- baseline[i]
    key <- paste0('scope__',row$scope)
    a <- m08_prediction(bundle,paste0('2019_',row$scenario_internal),key)
    b <- m08_prediction(bundle,paste0('2025_',row$scenario_internal),key)
    baseline[i,`:=`(cached_2019=a$estimate[h],cached_2025=b$estimate[h],
                    analytic_SE=sqrt(sum((b$iid[,h]-a$iid[,h])^2)))]
    cbind(risk19=a$iid[,h],risk25=b$iid[,h],RD=b$iid[,h]-a$iid[,h])
  })
  ref <- match(paste(baseline$scope,baseline$time,'S0'),
               paste(baseline$scope,baseline$time,baseline$scenario_internal))
  baseline[, analytic_SE_ratio_vs_S0 := ifelse(analytic_SE[ref]>0,analytic_SE/analytic_SE[ref],NA_real_)]
  baseline[, analytic_SE_increment := vapply(seq_len(.N),function(i)
    sqrt(sum((gradients[[i]][,'RD']-gradients[[ref[i]]][,'RD'])^2)),numeric(1L))]
  baseline[, `:=`(imputation_id=m, max_abs_risk_error=pmax(abs(risk19-cached_2019),abs(risk25-cached_2025)))]
  baseline[, status := ifelse(max_abs_risk_error<=1e-8,'agree','mismatch')]
  if(any(baseline$status!='agree')) warning('Independent point predictions differ from cached risks. Inspect the baseline comparison before interpreting gradients.')
  derivative <- list()
  for (direction in seq_len(directions)) {
    cat(sprintf("[独立审计] MI %d 权重扰动 %d/%d\n",m,direction,directions))
    z <- sample(c(-1,1), n, replace=TRUE)
    analytic <- t(vapply(gradients,function(g)as.numeric(crossprod(z,g)),numeric(3L)))
    for (step in c(0.01,0.001)) {
      plus <- evaluate(1+step*z); minus <- evaluate(1-step*z)
      stopifnot(identical(plus[,..keys],baseline[,..keys]),identical(minus[,..keys],baseline[,..keys]))
      for (j in seq_along(c('risk19','risk25','RD'))) {
        metric <- c('risk19','risk25','RD')[j]
        numeric <- (plus[[metric]]-minus[[metric]])/(2*step)
        error <- abs(numeric-analytic[,j])
        derivative[[length(derivative)+1L]] <- cbind(baseline[,..keys], data.table::data.table(
          imputation_id=m,direction=direction,step=step,metric=metric,analytic=analytic[,j],
          finite_difference=numeric,absolute_error=error,
          status=ifelse(baseline$status=='agree' & error<=1e-9+0.01*pmax(abs(numeric),abs(analytic[,j])),'agree','review')))
      }
    }
  }
  draws <- failures <- list()
  for (i in seq_len(reps)) {
    if(i==1L || i%%10L==0L) cat(sprintf('[独立审计] MI %d bootstrap %d/%d\n',m,i,reps))
    weights <- tabulate(sample.int(n,n,replace=TRUE),nbins=n)
    for (model in if(spline)c('original','LDL_spline_ASCVD_and_death')else 'original') {
      fit <- tryCatch(evaluate(weights,model!='original'),error=identity)
      if(inherits(fit,'error')) {
        failures[[length(failures)+1L]] <- data.table::data.table(imputation_id=m,replicate=i,model=model,message=conditionMessage(fit))
      } else {
        fit[, `:=`(imputation_id=m,replicate=i,model=model)]
        draws[[length(draws)+1L]] <- fit
      }
    }
  }
  draws <- data.table::rbindlist(draws)
  summary <- data.table::data.table()
  if(nrow(draws)) {
    # Fit success is replicate-wide. Use the same successful replicates across
    # model specifications as well, so model differences are paired.
    matched <- draws[,.(models=data.table::uniqueN(model)),by=replicate][models==(1L+as.integer(spline)),replicate]
    paired <- draws[replicate%in%matched]
    if(nrow(paired)) {
      ref <- paired[scenario_internal=='S0',.(replicate,model,scope,time,S0_RD=RD)]
      paired <- merge(paired,ref,by=c('replicate','model','scope','time'))
      paired[, RM_minus_S0_RD := RD-S0_RD]
      summary <- paired[, .(successful=.N, requested=reps,
        bootstrap_U19=var(risk19),bootstrap_U25=var(risk25),bootstrap_cov=cov(risk19,risk25),
        bootstrap_SE_RD=sd(RD),bootstrap_SE_S0=sd(S0_RD),bootstrap_SE_increment=sd(RM_minus_S0_RD),
        bootstrap_mean_RD=mean(RD),
        bootstrap_percentile_lcl=unname(quantile(RD,.025)),
        bootstrap_percentile_ucl=unname(quantile(RD,.975))),by=.(imputation_id,model,scope,time,scenario_internal)]
      summary <- merge(summary,baseline[,c(keys,'analytic_SE','analytic_SE_ratio_vs_S0','analytic_SE_increment'),with=FALSE],by=keys)
      summary[, `:=`(bootstrap_RM_over_S0=ifelse(bootstrap_SE_S0>0,bootstrap_SE_RD/bootstrap_SE_S0,NA_real_),
        bootstrap_over_analytic_SE=ifelse(model=='original' & analytic_SE>0,bootstrap_SE_RD/analytic_SE,NA_real_),
        bootstrap_increment_over_analytic_SE=ifelse(model=='original' & analytic_SE_increment>0,
          bootstrap_SE_increment/analytic_SE_increment,NA_real_),
        approximate_SE_relative_MC_error=ifelse(successful>1L,1/sqrt(2*(successful-1L)),NA_real_),
        bootstrap_status=ifelse(successful<100L | successful/requested<0.9,'pilot_or_insufficient','review_estimates'))]
    }
  }
  if(reps>0L && !nrow(summary)) warning('No matched bootstrap estimates. Inspect failures.csv.')
  list(baseline=baseline,gradient=data.table::rbindlist(derivative),draws=draws,
       failures=data.table::rbindlist(failures),summary=summary)
}

.audit_boot_reps <- .audit_integer('M01_AUDIT_BOOT_REPS',0L)
.audit_gradient_directions <- .audit_integer('M01_AUDIT_GRADIENT_DIRECTIONS',0L)
.audit_refit_mi <- .audit_integer('M01_AUDIT_REFIT_MI',1L,1L)
.audit_refit_enabled <- .audit_boot_reps>0L || .audit_gradient_directions>0L
.audit_spline <- identical(.audit_env('M01_AUDIT_SPLINE'),'1')
if(.audit_refit_enabled && .audit_refit_mi>.audit_cfg$m) stop('Selected refit MI is not available')
if(.audit_boot_reps==1L) stop('Set bootstrap repeats to 0 or at least 2')
if(.audit_spline && .audit_boot_reps<2L) stop('Spline sensitivity needs at least two bootstrap repeats')
.audit_refit_results <- NULL

audit_per_mi<-list();audit_cross<-list();audit_transition<-list();audit_incremental<-list()
ni<-nc<-nt<-nz<-0L
cat("\n================ 开始逐套 MI 的只读审计 ================\n")
for (m in seq_len(.audit_cfg$m)) {
  cat(sprintf("\n[审计] MI %02d/%02d：读取08预测bundle和07个体策略\n",m,.audit_cfg$m))
  # 原 m08_bundle() 会重新计算输入/代码签名，并拒绝过期缓存；不会重拟合。
  b<-tryCatch(m08_bundle(.audit_cfg,m),error=function(e)e)
  if(inherits(b,"error"))stop("MI ",m," 的08缓存验证失败：",conditionMessage(b),
    "\n需要先检查代码文件是否与08缓存对应；不要绕过签名复用旧缓存。",call.=FALSE)
  x<-m07_read(.audit_cfg,m)
  if(b$engine$n!=nrow(x))stop("MI",m,"：07/08参与者总数不一致")
  if (!identical(m07_hash_object(as.character(x$participant_id)),b$engine$participant_id_md5))
    stop("MI",m,"：07与08中参与者顺序/ID不一致，不能计算配对IF")
  if (!identical(as.numeric(b$engine$times),as.numeric(.audit_cfg$times)))
    stop("MI",m,"：08缓存分析时点与07配置不一致")
  available_scopes<-names(m07_scope_spec(x))
  unknown<-setdiff(.audit_scopes,available_scopes)
  if(length(unknown))stop("不存在的亚组：",paste(unknown,collapse=", "),
    "；可选：",paste(available_scopes,collapse=", "))
  for(scope in .audit_scopes) {
    key<-paste0("scope__",scope); ix<-b$targets$index[[key]]
    if(!length(ix))stop("08缓存缺少亚组：",scope)
    if(!identical(sort(as.integer(ix)),sort(which(m07_scope_spec(x)[[scope]]))))
      stop("MI",m,"：08缓存的亚组成员与07不一致：",scope)
    cat("[审计] MI",m," / ",scope," / n=",length(ix),"\n",sep="")
    ref19<-.audit_component_IF(b,"2019_S0",key)
    ref25<-.audit_component_IF(b,"2025_S0",key)
    refIF<-ref25$pred$iid-ref19$pred$iid
    refQ<-ref25$pred$estimate-ref19$pred$estimate
    for(k in .audit_scenarios){
      if (k=="S0") {
        p19 <- ref19
        p25 <- ref25
      } else {
        p19 <- .audit_component_IF(b,paste0("2019_",k),key)
        p25 <- .audit_component_IF(b,paste0("2025_",k),key)
      }
      for(h in seq_along(.audit_cfg$times)){
        t<-.audit_cfg$times[h]
        v19<-p19$pred$iid[,h];v25<-p25$pred$iid[,h];vD<-v25-v19
        U19<-sum(v19^2);U25<-sum(v25^2);cov19_25<-sum(v19*v25);Urd<-sum(vD^2)
        Uindependent<-U19+U25
        if(!.audit_comparison(Urd,Uindependent-2*cov19_25))
          stop("协方差恒等式失败：MI",m,"/",scope,"/",k,"/",t)
        C<-.audit_crossmatrix(p25,p19,h)
        if(!.audit_comparison(sum(C),Urd))
          stop("五项IF分量方差重构失败：MI",m,"/",scope,"/",k,"/",t)
        Qrd<-p25$pred$estimate[h]-p19$pred$estimate[h]
        # 与原09 m09_one_effect直接复核（使用实际Q/U，不是仅仅比较自己写的函数）。
        pp<-setNames(list(p19$pred,p25$pred),c(paste0("2019_",k),paste0("2025_",k)))
        weights<-setNames(c(1,-1),c(paste0("2025_",k),paste0("2019_",k)))
        original<-m09_one_effect(pp,weights,"rd_2025_vs_2019",scope,k,m,
                                 b$engine$complete_df)
        if(!.audit_comparison(original$Q[h],Qrd)||!.audit_comparison(original$U[h],Urd))
          stop("独立计算与原09的Q/U不一致：MI",m,"/",scope,"/",k,"/",t)
        pair_risk <- b$states[[paste0("2025_",k)]]$risk[ix,h]-
                     b$states[[paste0("2019_",k)]]$risk[ix,h]
        same_ldl<-abs(x[[paste0("Lcf_2025_",k)]][ix]-
                      x[[paste0("Lcf_2019_",k)]][ix])<1e-12
        if(!.audit_comparison(mean(pair_risk),Qrd))
          stop("08每人风险差平均值无法重构总体RD：MI",m,"/",scope,"/",k,"/",t)
        # 相同反事实 LDL 对应相同模型风险（并不意味着其总 IF 分量逐人均为零）。
        max_same_diff<-if(any(same_ldl))max(abs(pair_risk[same_ldl]))else 0
        if(max_same_diff>1e-10)
          stop("相同反事实LDL却预测出不同风险：MI",m,"/",scope,"/",k,"/",t)
        ni<-ni+1L
        audit_per_mi[[ni]]<-data.table::data.table(
          imputation_id=m,scope=scope,scenario_internal=k,scenario_paper=.audit_label[[k]],
          time=t,n_population=length(ix),Q_2019=p19$pred$estimate[h],Q_2025=p25$pred$estimate[h],
          Q_RD=Qrd,RD_pp=100*Qrd,
          U_2019=U19,U_2025=U25,Cov_2019_2025=cov19_25,U_RD=Urd,
          U_independent_assumption=Uindependent,
          covariance_cancellation_fraction=if(Uindependent>0)2*cov19_25/Uindependent else NA_real_,
          variance_ratio_paired_over_independent=if(Uindependent>0)Urd/Uindependent else NA_real_,
          corr_IF_2019_2025=if(U19>0 && U25>0)cov19_25/sqrt(U19*U25) else NA_real_,
          SE_RD_pp_within_MI=100*sqrt(Urd),
          fraction_same_counterfactual_LDL=mean(same_ldl),
          fraction_same_individual_risk=mean(abs(pair_risk)<1e-12),
          max_risk_difference_when_same_LDL=max_same_diff,
          U_main_effect_components=sum(diag(C)),
          U_cross_component_terms=Urd-sum(diag(C)),
          IF_reconstruction_error=max(p19$recon_error,p25$recon_error),
          original_m09_Q_error=Qrd-original$Q[h],
          original_m09_U_error=Urd-original$U[h],
          dfcom=b$engine$complete_df)
        # 每个分量的对角方差 + 10对带符号的交叉项；可验证总体方差来源。
        for(j in seq_along(.audit_components)){
          nc<-nc+1L
          audit_cross[[nc]]<-data.table::data.table(
            imputation_id=m,scope=scope,scenario_internal=k,time=t,
            kind="variance_component",term_1=.audit_components[j],term_2=.audit_components[j],
            variance_contribution=C[j,j])
          if(j<length(.audit_components))for(l in seq.int(j+1L,length(.audit_components))){
            nc<-nc+1L
            audit_cross[[nc]]<-data.table::data.table(
              imputation_id=m,scope=scope,scenario_internal=k,time=t,
              kind="twice_covariance",term_1=.audit_components[j],term_2=.audit_components[l],
              variance_contribution=2*C[j,l])
          }
        }
        # 依真实07 transition_ 情景分组：可量化“unchanged intensity”对精度的解释力度。
        transition<-as.character(x[[paste0("transition_",k)]][ix])
        if(anyNA(transition)||!setequal(unique(transition),unique(transition[!is.na(transition)])))
          stop("07治疗转换存在NA：MI",m,"/",scope,"/",k)
        for(tr in m07_transitions()){
          is_tr<-transition==tr;n_tr<-sum(is_tr)
          nt<-nt+1L
          audit_transition[[nt]]<-data.table::data.table(
            imputation_id=m,scope=scope,scenario_internal=k,scenario_paper=.audit_label[[k]],
            time=t,transition=tr,n=n_tr,proportion=n_tr/length(ix),
            share_same_counterfactual_LDL=if(n_tr)mean(same_ldl[is_tr]) else NA_real_,
            RD_within_transition_pp=if(n_tr)100*mean(pair_risk[is_tr]) else NA_real_,
            contribution_to_scope_RD_pp=if(n_tr)100*sum(pair_risk[is_tr])/length(ix) else 0)
        }
        if(k!="S0"){
          cur<-list(Q=Qrd,iid=vD); ref<-list(Q=refQ[h],iid=refIF[,h])
          nz<-nz+1L
          audit_incremental[[nz]]<-.audit_record_incremental(ref,cur,m,k,scope,t,b$engine$complete_df)
        }
      }
      if(k!="S0"){rm(p19,p25);gc(FALSE)}
    }
    rm(ref19,ref25,refIF,refQ);gc(FALSE)
  }
  if(.audit_refit_enabled && m==.audit_refit_mi) {
    set.seed(.audit_integer('M01_AUDIT_REFIT_SEED',20261110L,1L))
    .audit_refit_results <- tryCatch(
      .audit_independent(x,b,m,.audit_cfg,.audit_scopes,c('S0','SHT'),
        .audit_boot_reps,.audit_gradient_directions,.audit_spline),
      error=function(e) {
        warning('Independent verification failed: ',conditionMessage(e),call.=FALSE)
        list(failures=data.table::data.table(imputation_id=m,stage='independent_verification',
                                            message=conditionMessage(e)))
      })
  }
  rm(x,b);gc(FALSE)
}

# ---------------------------- 合并与诊断输出 -------------------------------
perMI<-data.table::rbindlist(audit_per_mi,fill=TRUE)
parts<-data.table::rbindlist(audit_cross,fill=TRUE)
tr<-data.table::rbindlist(audit_transition,fill=TRUE)
incMI<-data.table::rbindlist(audit_incremental,fill=TRUE)
if(!nrow(perMI))stop("没有生成任何审计记录")
pooled<-.audit_pooled(perMI,c("scope","scenario_internal","scenario_paper","time"),
                      q_col="Q_RD",u_col="U_RD")
extra<-perMI[,.(mean_IF_corr=mean(corr_IF_2019_2025),
  mean_cancellation_fraction=mean(covariance_cancellation_fraction),
  mean_variance_ratio=mean(variance_ratio_paired_over_independent),
  mean_identical_LDL_fraction=mean(fraction_same_counterfactual_LDL),
  mean_identical_risk_fraction=mean(fraction_same_individual_risk),
  mean_U19=mean(U_2019),mean_U25=mean(U_2025),
  mean_Cov19_25=mean(Cov_2019_2025)),
  by=.(scope,scenario_internal,scenario_paper,time)]
pooled<-merge(pooled,extra,by=c("scope","scenario_internal","scenario_paper","time"))
refSE<-pooled[scenario_internal=="S0",.(scope,time,S0_SE_pp=se_pp,S0_CI_width_pp=CI_upper_pp-CI_lower_pp)]
pooled<-merge(pooled,refSE,by=c("scope","time"),all.x=TRUE)
pooled[,`:=`(SE_ratio_RM_over_S0=ifelse(S0_SE_pp>0,se_pp/S0_SE_pp,NA_real_),
  CI_width_ratio_RM_over_S0=ifelse(S0_CI_width_pp>0,
    (CI_upper_pp-CI_lower_pp)/S0_CI_width_pp,NA_real_))]
# Positive contributions explain how much of the S0 minus RM variance is lost.
variance_ref <- pooled[scenario_internal=='S0',.(scope,time,
  marginal_S0=mean_U19+mean_U25,covariance_S0=2*mean_Cov19_25,
  within_S0=within_variance,between_S0=between_variance,total_S0=total_variance)]
variance_change <- merge(pooled[scenario_internal!='S0'],variance_ref,by=c('scope','time'))
variance_change[, `:=`(
  marginal_variance_reduction_pp2=1e4*(marginal_S0-mean_U19-mean_U25),
  covariance_cancellation_gain_pp2=1e4*(2*mean_Cov19_25-covariance_S0),
  within_variance_reduction_pp2=1e4*(within_S0-within_variance),
  MI_between_contribution_pp2=1e4*(1+1/.audit_cfg$m)*(between_S0-between_variance),
  total_variance_reduction_pp2=1e4*(total_S0-total_variance))]
if(any(abs(variance_change[,marginal_variance_reduction_pp2+covariance_cancellation_gain_pp2-
                            within_variance_reduction_pp2])>1e-8) ||
   any(abs(variance_change[,within_variance_reduction_pp2+MI_between_contribution_pp2-
                            total_variance_reduction_pp2])>1e-8)) stop('Variance change decomposition failed')
incPooled<-.audit_pooled(incMI,c("scope","scenario_internal","scenario_paper","time","comparison"),
                         q_col="Q",u_col="U")

# 检查组别加权RD能重建总体RD（与原12的质量控制保持一致）。
trQA<-tr[,.(sum_group_RD_pp=sum(contribution_to_scope_RD_pp),
  summed_n=sum(n),sum_group_proportion=sum(proportion)),
  by=.(imputation_id,scope,scenario_internal,time)]
trQA<-merge(trQA,perMI[,.(imputation_id,scope,scenario_internal,time,RD_pp,n_population)],
            by=c("imputation_id","scope","scenario_internal","time"))
trQA[,point_error_pp:=sum_group_RD_pp-RD_pp]
trQA[,passed:=abs(point_error_pp)<1e-7 &
  summed_n==n_population & abs(sum_group_proportion-1)<1e-10]
if(any(!trQA$passed))stop("六组治疗转换分解无法重构总体风险差")

# 源 09/11 已导出CSV的独立交叉核对：缺少文件明确报告 'missing'，不能假装通过。
.read_saved <- function(path){
  if(!file.exists(path))return(NULL)
  tryCatch(data.table::fread(path,encoding="UTF-8"),error=function(e){
    warning("无法读取已有CSV：",path," / ",conditionMessage(e));NULL})
}
.compare_rows <- function(source,path,a,b,join,expected_fields,actual_fields) {
  if(is.null(b))return(data.table::data.table(source=source,path=path,status="file_missing",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details="文件缺失，未执行对照"))
  miss<-setdiff(c(join,actual_fields),names(b))
  if(length(miss))return(data.table::data.table(source=source,path=path,status="columns_missing",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details=paste(miss,collapse=",")))
  merged<-merge(a,b,by=join,all.x=TRUE,suffixes=c("_audit","_saved"))
  if(!nrow(merged))return(data.table::data.table(source=source,path=path,status="no_matched_rows",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details="没有匹配到要核对的亚组/情景/时间点"))
  v1<-expected_fields;v2<-actual_fields
  # merge() 会给双方同名字段自动加 _audit / _saved 后缀。
  eQ<-abs(merged[[paste0(v1[1],"_audit")]]-merged[[paste0(v2[1],"_saved")]])
  eU<-abs(merged[[paste0(v1[2],"_audit")]]-merged[[paste0(v2[2],"_saved")]])
  qscale<-pmax(abs(merged[[paste0(v1[1],"_audit")]]),
               abs(merged[[paste0(v2[1],"_saved")]]))
  uscale<-pmax(abs(merged[[paste0(v1[2],"_audit")]]),
               abs(merged[[paste0(v2[2],"_saved")]]))
  bad<-(!is.finite(eQ)|!is.finite(eU)|
          eQ>pmax(1e-11,1e-7*qscale)|eU>pmax(1e-15,1e-5*uscale))
  data.table::data.table(source=source,path=path,
    status=if(any(bad))"mismatch"else"matched",
    checked_rows=nrow(merged),failed_rows=sum(bad),
    max_abs_Q_error=if(all(is.na(eQ)))NA_real_ else max(eQ,na.rm=TRUE),
    max_abs_U_error=if(all(is.na(eU)))NA_real_ else max(eU,na.rm=TRUE),
    details=if(any(bad))"Q/U与原09/11导出值不一致"else"逐MI Q/U一致")
}
# 实际 Function_09 / Function_11 用 m09_emit() 写入以下路径。
a09<-file.path(.audit_out_dir,"strategy_risks_and_effects/per_imputation/strategy_risks_and_risk_differences.csv")
a11<-file.path(.audit_out_dir,"strategy_risks_and_effects/per_imputation/all_scenario_risks_and_effects.csv")
aInc<-file.path(.audit_out_dir,"risk_modifier_incremental_effects/per_imputation/incremental_benefits_and_guideline_difference.csv")
s09<-.read_saved(a09);s11<-.read_saved(a11);sInc<-.read_saved(aInc)
expected<-perMI[,.(imputation_id,scope,scenario=scenario_internal,time,estimand="rd_2025_vs_2019",Q=Q_RD,U=U_RD)]
filtered09<-if(is.null(s09))NULL else s09[scenario=="S0" & estimand=="rd_2025_vs_2019"]
filtered11<-if(is.null(s11))NULL else s11[estimand=="rd_2025_vs_2019"]
comp09<-.compare_rows("Function_09",a09,expected[scenario=="S0"],filtered09,
  c("imputation_id","scope","scenario","time","estimand"),c("Q","U"),c("Q","U"))
comp11<-.compare_rows("Function_11",a11,expected,filtered11,
  c("imputation_id","scope","scenario","time","estimand"),c("Q","U"),c("Q","U"))
# Function_11 中 additional_reduction_difference_2025_vs_2019
# 与本审计定义的 [RD(RM)-RD(S0)] 正好互为相反数；其U应完全相等。
expectedInc<-incMI[,.(imputation_id,scope,scenario=scenario_internal,time,
  estimand="additional_reduction_difference_2025_vs_2019",rm_comparison=paste0(scenario_internal,"_vs_S0"),
  Q=-Q,U=U)]
filteredInc<-if(is.null(sInc))NULL else sInc[estimand=="additional_reduction_difference_2025_vs_2019" &
    rm_comparison %in% paste0(c("SH","ST","SHT"),"_vs_S0")]
compInc<-.compare_rows("Function_11_incremental",aInc,expectedInc,filteredInc,
  c("imputation_id","scope","scenario","time","estimand","rm_comparison"),
  c("Q","U"),c("Q","U"))
savedQA<-data.table::rbindlist(list(comp09,comp11,compInc),fill=TRUE)

# 第二层对照：检查原09/11已经合并的点估计、95%CI及Rubin U/B/T。
# 与逐MI检查不同，这一层能识别"单套MI正确、合并结果或旧CSV不一致"。
.compare_pooled_rows <- function(source,path,a,b,join){
  if(is.null(b))return(data.table::data.table(source=source,path=path,status="file_missing",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details="合并CSV不存在，未对照"))
  cols<-c("estimate","lcl","ucl","within_variance","between_variance","total_variance")
  missing<-setdiff(c(join,cols),names(b))
  if(length(missing))return(data.table::data.table(source=source,path=path,status="columns_missing",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details=paste(missing,collapse=",")))
  x<-merge(a,b,by=join,all.x=TRUE,suffixes=c("_audit","_saved"))
  if(!nrow(x))return(data.table::data.table(source=source,path=path,status="no_matched_rows",
    checked_rows=0L,failed_rows=NA_integer_,max_abs_Q_error=NA_real_,max_abs_U_error=NA_real_,
    details="没有与审计条件对应的合并行"))
  diffs<-lapply(cols,function(nm){
    va<-x[[paste0(nm,"_audit")]]; vb<-x[[paste0(nm,"_saved")]]
    d<-abs(va-vb)
    # 风险本身容许数值约1e-11；方差在很小尺度下仍要求相对一致。
    tol <- if (nm %in% c("estimate","lcl","ucl")) {
      pmax(1e-11,1e-7*pmax(abs(va),abs(vb)))
    } else {
      pmax(1e-15,1e-5*pmax(abs(va),abs(vb)))
    }
    list(error=d,bad=(!is.finite(d)|!is.finite(tol)|d>tol))
  })
  bad<-Reduce(`|`,lapply(diffs,`[[`,"bad"))
  maxdiff<-function(v){v<-v[is.finite(v)];if(length(v))max(v)else NA_real_}
  data.table::data.table(source=source,path=path,status=if(any(bad))"mismatch"else"matched",
    checked_rows=nrow(x),failed_rows=sum(bad),
    max_abs_Q_error=maxdiff(diffs[[1]]$error),max_abs_U_error=maxdiff(diffs[[4]]$error),
    details=if(any(bad))"已合并点估计/CI或Rubin方差不一致"else"已合并估计和Rubin方差均一致")
}
pExpected<-pooled[,.(scope,scenario=scenario_internal,time,estimand="rd_2025_vs_2019",
  estimate=Q_pooled,lcl=CI_lower_pp/100,ucl=CI_upper_pp/100,
  within_variance,between_variance,total_variance)]
pExpectedInc<-incPooled[,.(scope,scenario=scenario_internal,time,
  estimand="additional_reduction_difference_2025_vs_2019",
  rm_comparison=paste0(scenario_internal,"_vs_S0"),
  estimate=-Q_pooled,lcl=-CI_upper_pp/100,ucl=-CI_lower_pp/100,
  within_variance,between_variance,total_variance)]
p09Path<-file.path(.audit_out_dir,"strategy_risks_and_effects/pooled/strategy_risks_and_risk_differences.csv")
p11Path<-file.path(.audit_out_dir,"strategy_risks_and_effects/pooled/all_scenario_risks_and_effects.csv")
pIncPath<-file.path(.audit_out_dir,"risk_modifier_incremental_effects/pooled/incremental_benefits_and_guideline_difference.csv")
p09<-.read_saved(p09Path);p11<-.read_saved(p11Path);pInc<-.read_saved(pIncPath)
if(!is.null(p09))p09<-p09[scenario=="S0" & estimand=="rd_2025_vs_2019"]
if(!is.null(p11))p11<-p11[estimand=="rd_2025_vs_2019"]
if(!is.null(pInc))pInc<-pInc[estimand=="additional_reduction_difference_2025_vs_2019" &
  rm_comparison %in% paste0(c("SH","ST","SHT"),"_vs_S0")]
savedQA<-data.table::rbindlist(list(savedQA,
  .compare_pooled_rows("Function_09_pooled",p09Path,pExpected[scenario=="S0"],p09,
                       c("scope","scenario","time","estimand")),
  .compare_pooled_rows("Function_11_pooled",p11Path,pExpected,p11,
                       c("scope","scenario","time","estimand")),
  .compare_pooled_rows("Function_11_incremental_pooled",pIncPath,pExpectedInc,pInc,
                       c("scope","scenario","time","estimand","rm_comparison"))),fill=TRUE)
# 路径信息+输入签名概述，以便区分数据问题和程序错误。
meta<-data.table::data.table(item=c("analysis_dir","code_dir","config","M","times","scopes",
  "scenarios","function_07_MD5","function_08_MD5","function_09_MD5","audit_date"),
  value=c(.audit_out_dir,.audit_code_dir,file.path(.audit_out_dir,"study_config.rds"),
    as.character(.audit_cfg$m),paste(.audit_cfg$times,collapse=","),paste(.audit_scopes,collapse=","),
    paste(.audit_scenarios,collapse=","),
    unname(tools::md5sum(file.path(.audit_code_dir,.audit_required_fns))),
    format(Sys.time(),"%Y-%m-%d %H:%M:%S")))

# 全部写入一个全新目录；旧 07–12 产物保持不动。
.audit_dest <- file.path(.audit_out_dir,"rm_ci_variance_audit_v2_2")
dir.create(.audit_dest,recursive=TRUE,showWarnings=FALSE)
.write <- function(x,f) data.table::fwrite(x,file.path(.audit_dest,f),bom=TRUE,na="NA")
.write(perMI,"01_per_MI_paired_variance.csv")
.write(pooled,"02_pooled_paired_variance.csv")
.write(parts,"03_per_MI_IF_components_covariance.csv")
.write(tr,"04_treatment_transition_mechanism.csv")
.write(trQA,"05_transition_reconstruction_QA.csv")
.write(incMI,"06_per_MI_RM_minus_S0_RD.csv")
.write(incPooled,"07_pooled_RM_minus_S0_RD.csv")
.write(savedQA,"08_saved_09_11_Q_U_checks.csv")
.write(meta,"09_input_paths_signatures.csv")
.write(variance_change,"11_variance_reduction_vs_S0.csv")
if(!is.null(.audit_refit_results)) {
  refit_dir <- tempfile(pattern=sprintf("independent_MI%02d_B%d_D%d_spline%d_seed%d_",
    .audit_refit_mi,.audit_boot_reps,.audit_gradient_directions,as.integer(.audit_spline),
    .audit_integer("M01_AUDIT_REFIT_SEED",20261110L,1L)),tmpdir=.audit_dest)
  dir.create(refit_dir,recursive=TRUE,showWarnings=FALSE)
  for (name in names(.audit_refit_results)) {
    result <- .audit_refit_results[[name]]
    if(ncol(result)) data.table::fwrite(result,file.path(refit_dir,paste0(name,".csv")))
  }
}

report<-c(
  "M01 RM风险差 CI 审计 v2.2（只读输出）",
  paste0("生成日期：",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("代码目录：",.audit_code_dir),paste0("分析目录：",.audit_out_dir),
  paste0("插补套数：",.audit_cfg$m,"；分析时点：",paste(.audit_cfg$times,collapse=","),"年"),
  "",
  "1. 此结果仅属于固定Cox模型、固定RM升级规则和固定干预LDL随机分位数下的统计推断。",
  "2. 请查看02_pooled_paired_variance.csv的SE_ratio_RM_over_S0（小于1意味着RM较S0更精确），",
  "   mean_IF_corr、mean_cancellation_fraction、mean_identical_LDL_fraction和fraction_between。",
  "3. 查看03_per_MI_IF_components_covariance.csv：五项分量的独立方差与两两协方差必须合计为U_RD。",
  "4. 查看04_treatment_transition_mechanism.csv：Unchanged intensity及其他分组的比例和RD贡献。",
  "5. 查看06/07：RM_minus_S0_RD = RD(RM)-RD(S0)，与原11中additional_reduction_difference的符号相反。",
  "6. 08_saved_09_11_Q_U_checks.csv如显示mismatch，不能使用审计结果解释原论文的CI。",
  "7. 如02的CI_rounds_equal_at_2dp=TRUE，论文显示(0.00–0.00)不代表标准误为零。",
  "8. 本审计的IF分量重构验证了代数一致性，不等于独立证明Cox/Breslow解析梯度正确。",
  "9. 11_variance_reduction_vs_S0.csv将总方差减少拆为边际方差、协方差抵消和MI间方差。",
  "   可选独立验证见independent_MI目录，比较完整IF数值导数与配对重拟合SE。",
  "10. RM测量误差、RM决策规则不确定性、模型错设和Monte Carlo误差不在该常规CI之内。",
  "11. v2.2修复了perMI(Q_RD/U_RD)与incMI(Q/U)字段映射，并在读取大RDS前进行合并自检。",
  "",
  "[已导出CSV对照检查]"
)
for(i in seq_len(nrow(savedQA)))report<-c(report,paste0("  ",savedQA$source[i],": ",savedQA$status[i],
 "；matched_rows=",savedQA$checked_rows[i],"；failed=",savedQA$failed_rows[i]))
for(scope_name in .audit_scopes)for(time_point in .audit_cfg$times){
  z<-pooled[pooled$scope==scope_name & pooled$time==time_point]
  if(nrow(z)==0L) next
  report<-c(report,"",paste0("[",scope_name," / ",time_point,"年]"))
  for(k in .audit_scenarios){
    a<-z[scenario_internal==k]
    if(!nrow(a))next
    report<-c(report,sprintf("  %-3s RD=%+.6f pp；95%% CI=(%+.6f,%+.6f) pp；SE=%.7f pp；IF corr=%.7f；同LDL比例=%.2f%%；SE相对S0=%.4f",
       .audit_label[[k]],a$estimate_pp,a$CI_lower_pp,a$CI_upper_pp,a$se_pp,
       a$mean_IF_corr,100*a$mean_identical_LDL_fraction,a$SE_ratio_RM_over_S0))
  }
}
report<-c(report,"",paste0("完成耗时（秒）：",round(as.numeric(difftime(Sys.time(),audit_time_start,units="secs")),1)))
writeLines(report,file.path(.audit_dest,"10_audit_readme_and_summary.txt"),useBytes=TRUE)

cat("\n================ 审计完成 ================\n")
print(pooled[,.(scope,scenario_paper,time,estimate_pp,CI_lower_pp,CI_upper_pp,
                 se_pp,SE_ratio_RM_over_S0,mean_IF_corr,mean_identical_LDL_fraction,
                 fraction_between)])
cat("\n09/11结果一致性：\n")
print(savedQA[,.(source,status,checked_rows,failed_rows,max_abs_Q_error,max_abs_U_error)])
if(any(savedQA$status=="mismatch"))
 warning("发现已导出09/11 Q或U与审计不一致；请检查对应版本及缓存签名。",call.=FALSE)
cat("\n诊断文件目录：",.audit_dest,"\n",sep="")
if(.audit_refit_enabled) cat("独立验证目录：",refit_dir,"\n",sep="") else cat("本次仅完成缓存审计。独立重拟合与Bootstrap可通过脚本顶部的环境变量启用。\n")
