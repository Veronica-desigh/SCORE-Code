# 线程与进度工具（所有 S01 Function 文件自包含；主脚本可分别设置 workers=3L）。
s01_check_workers <- function(workers, tasks) {
  val <- suppressWarnings(as.numeric(workers))
  n <- suppressWarnings(as.integer(val))
  if(length(val)!=1L || is.na(n) || !is.finite(val) || n<1L || val!=n)
    stop("workers 必须为正整数，例如 workers <- 3L")
  if(length(tasks)!=1L || !is.finite(tasks) || tasks<1L)
    stop("并行任务数必须为正整数")
  min(n,as.integer(tasks))
}
s01_print_progress <- function(stage, done, total, started, detail="") {
  elapsed <- proc.time()[["elapsed"]]-started
  remaining <- if(done>0L && done<total) elapsed*(total-done)/done else if(done>=total) 0 else NA_real_
  message(sprintf("[%s] S01 %-10s | %d/%d (%5.1f%%) | 耗时 %.1fs (%.2fmin) | 预计剩余 %s%s",
    format(Sys.time(),"%H:%M:%S"), stage,done,total,100*done/total,
    elapsed,elapsed/60,if(is.finite(remaining))sprintf("%.1fs",remaining) else "--",
    if(nzchar(detail)) paste0(" | ",detail) else ""))
}

# S01 共用计算。每份插补保留原队列、评分、结局与反事实 LDL-C 随机分位数。
# 策略只改变原始中危者的 RM 触发；最低降幅、竞争风险和配对 IF 调用 M01。
s01_years <- c("2019", "2025")
s01_alias <- c(S0="S0", SN="SH", ST="ST", SNT="SHT")
s01_scenario <- function(id, parent, rule="original", arg="", reference=parent) {
  list(id=id, parent=parent, rule=rule, arg=arg, reference=reference)
}
s01_originals <- function() lapply(names(s01_alias),function(id) s01_scenario(id,id))

# 固定16项RM ID及Word LOO/OAT的编号；禁止根据CSV当前排列自动重新编号。
s01_rm_order <- function() c("deprivation","obesity","inactive","stress","family",
  "inflammatory","psychiatric","hiv","af","sleep_apnoea","fatty_liver",
  "south_asian","premature_menopause","hdp","lpa","crp")
s01_proxy_exclusions <- function(rule,year) {
  if(!rule %in% c("remove2","strict")) return(character())
  basic <- c("family","stress")
  if(rule=="remove2") return(basic)
  c(basic,"sleep_apnoea",if(year=="2019") "fatty_liver" else "crp")
}
s01_check_dictionary <- function(dict) {
  stopifnot(all(c("domain","type","include_2019","include_2025") %in% names(dict)))
  if(!setequal(dict$domain,s01_rm_order()) || anyDuplicated(dict$domain))
    stop("RM_domains_v1.8必须含完整16项固定ID；不允许静默新增或漏项")
  if(!setequal(dict$domain[dict$type=="T"],c("lpa","crp")))
    stop("ST实验室RM必须为lpa/crp")
  for(year in s01_years) {
    actual <- dict$domain[dict[[paste0("include_",year)]]==1L]
    required <- if(year=="2019") c("family","stress","sleep_apnoea","fatty_liver")
      else c("family","stress","sleep_apnoea","crp")
    if(!all(required %in% actual)) stop("严格删除所需RM未列入",year,"年定义: ",paste(setdiff(required,actual),collapse=","))
  }
  invisible(TRUE)
}
s01_items <- function(dict,year,parent) {
  types <- switch(parent,S0=character(),SN="H",ST="T",SNT=c("H","T"),stop("未知原始策略: ",parent))
  dict$domain[dict[[paste0("include_",year)]]==1L & dict$type %in% types]
}

s01_flags <- function(x, items) {
  m <- as.matrix(as.data.frame(x)[,paste0("rm_",items),drop=FALSE])
  storage.mode(m) <- "integer"
  if (anyNA(m) || any(!m %in% 0:1)) stop("完成插补 RM 必须为 0/1")
  colnames(m) <- items
  m
}
s01_rule <- function(x,dict,sc,year,definitions=list()) {
  items <- s01_items(dict,year,sc$parent)
  base <- as.character(x[[if(year=="2019") "score_2019_category" else "score2_2025_category"]])
  score <- x[[if(year=="2019") "score_2019_percent" else "score2_2025_percent"]]
  moderate <- base=="M"
  window <- moderate
  if(sc$rule=="common") items <- intersect(items,s01_items(dict,setdiff(s01_years,year),sc$parent))
  if(sc$rule %in% c("remove2","strict"))
    items <- setdiff(items,s01_proxy_exclusions(sc$rule,year))
  if(sc$rule=="drop") items <- setdiff(items,sc$arg)
  # OAT中当年不适用的RM返回NULL，绝不可用S0或SNT替代。
  if(sc$rule=="one" && !sc$arg %in% items) return(NULL)
  if(sc$rule=="one") items <- sc$arg
  flags <- if(length(items)) s01_flags(x,items) else matrix(FALSE,nrow(x),0)
  if(sc$rule=="hmr") {
    bounds <- if(year=="2019") c(3,4,5) else c(6,8,10)
    if(!sc$arg %in% c("A","B")) stop("HMR类型必须为A或B")
    window <- moderate & score>=bounds[match(sc$arg,c("A","B"))] & score<bounds[3]
  }
  if(sc$rule=="p90" && "deprivation" %in% items) {
    cut <- definitions$p90
    if(length(cut)!=1L || !is.finite(cut)) stop("缺少从固定基线队列核实的Townsend P90")
    flags[,"deprivation"] <- as.integer(x$townsend_deprivation_index>=cut)
  }
  if(sc$rule=="crp" && "crp" %in% items) {
    value <- as.numeric(as.character(x$rm_source_crp))
    if(any(flags[,"crp"]==1L & !is.finite(value)))
      stop("CRP-RES: hs-CRP已判阳性但缺少可判定的数值，不能把未知当阴性")
    flags[,"crp"] <- as.integer(flags[,"crp"]==1L & is.finite(value) & value<=10)
  }
  # 新STRICT只改变RM集合；既不调用旧的definitions$strict，也不收紧保留RM定义。
  if(sc$rule=="observed") {
    f <- definitions$observed
    if(!is.function(f)) stop("缺少OBS观察确认函数")
    updated <- f(x,year,flags)
    if(!identical(dim(updated),dim(flags)) || !identical(colnames(updated),colnames(flags)) ||
       anyNA(updated) || any(!updated %in% 0:1) || any(updated>flags)) stop("OBS确认矩阵不合法")
    flags <- updated
  }
  domain_flags <- NULL
  hit <- rowSums(flags)>0
  if(sc$rule=="domain") {
    map <- definitions$domains
    if(is.null(map) || !all(dict$domain %in% names(map)) ||
       any(!map[dict$domain] %in% paste0("D",1:6))) stop("需要完整D1-D6领域映射")
    domain_flags <- vapply(paste0("D",1:6),function(d)
      rowSums(flags[,map[items]==d,drop=FALSE])>0,logical(nrow(x)))
    hit <- rowSums(domain_flags)>=as.integer(sc$arg)
  }
  if(anyNA(hit) || anyNA(window)) stop("RM触发或中危窗口存在缺失")
  list(moderate=moderate,window=window,upgrade=window&hit,items=items,
       score_level=base,level=ifelse(window&hit,"H",base),domain_flags=domain_flags)
}

# 始终使用M01已经冻结的干预目标、随机分位数和LDL-C处理函数。
s01_classify <- function(x,dict,sc,year,cfg,definitions) {
  z <- s01_rule(x,dict,sc,year,definitions)
  if(is.null(z)) return(NULL)
  z$ldl <- m07_one(x$ldl_cholesterol,z$level,x$z_common,cfg)
  z
}

# 每个分类人数使用明确分母。相对父情景的损失以父情景升级人数为分母。
s01_count_metric_names <- function() c(
  "score_based_moderate","score_based_high_very_high","rm_eligible_window","rm_upgraded",
  "upgrades_lost_vs_parent","upgrades_lost_vs_snt","upgrades_lost_vs_snt_initial_moderate",
  "final_high_very_high","intervention_eligible","new_intervention_vs_s0",
  "ldlc_requirement_changed_vs_parent","mean_min_ldlc_reduction",
  "mean_min_ldlc_reduction_change_vs_parent")
s01_counts <- function(z,parent,s0,snt) {
  n <- length(z$level); moderate <- sum(z$moderate)
  if(length(snt$upgrade)!=n || length(parent$upgrade)!=n) stop("比较策略人数不一致")
  lost_snt <- sum(snt$upgrade & !z$upgrade)
  numerators <- c(moderate,sum(z$score_level%in%c("H","VH")),sum(z$window),
    sum(z$upgrade),sum(parent$upgrade&!z$upgrade),lost_snt,lost_snt,
    sum(z$level%in%c("H","VH")),sum(z$ldl$A==1),
    sum(z$ldl$A==1&s0$ldl$A==0),sum(abs(z$ldl$D-parent$ldl$D)>1e-10),
    NA_real_,NA_real_)
  denoms <- c(n,n,moderate,moderate,sum(parent$upgrade),sum(snt$upgrade),moderate,
    n,n,n,n,n,n)
  values <- ifelse(denoms>0,numerators/denoms,NA_real_)
  values[12:13] <- c(mean(z$ldl$D),mean(z$ldl$D-parent$ldl$D))
  data.frame(metric=s01_count_metric_names(),n=numerators,denominator=denoms,value=values)
}

s01_transition <- function(a,b) {
  six <- c("unchanged_untreated","newly_intervention","intervention_removed",
           "treatment_intensified","unchanged_intensity","treatment_deintensified")
  tr <- m07_derive_transition(a$ldl$A,b$ldl$A,a$ldl$D,b$ldl$D)
  n6 <- tabulate(match(tr,six),length(six))
  n3 <- c(sum(n6[c(2,4)]),sum(n6[c(1,5)]),sum(n6[c(3,6)]))
  stopifnot(sum(n6)==length(tr),sum(n3)==length(tr))
  data.frame(group=c(six,"intensified","unchanged","deintensified"),
    grouping=rep(c(6L,3L),c(6L,3L)),n=c(n6,n3),denominator=length(tr),
    proportion=c(n6,n3)/length(tr))
}
s01_paths <- function(root) {
  list(code=file.path(root,"LJW_code","Main_code","M01_ASCVD"),
       clean=file.path(root,"LJW_outcome","data_clean","M01_ASCVD"),
       main=file.path(root,"LJW_outcome","data_intervention","M01_ASCVD","analysis_results"),
       out=Sys.getenv("S01_OUTPUT",file.path(root,"LJW_outcome","data_intervention",
          "S01_different_RM_strategies","analysis_results")))
}
s01_setup <- function(root) {
  p <- s01_paths(root)
  out <- tolower(normalizePath(p$out,winslash="/",mustWork=FALSE))
  main <- tolower(normalizePath(dirname(p$main),winslash="/",mustWork=TRUE))
  clean <- tolower(normalizePath(p$clean,winslash="/",mustWork=TRUE))
  if(any(out %in% c(main,clean)) || any(startsWith(out,paste0(c(main,clean),"/")))) {
    stop("S01_OUTPUT 必须与 M01 结果目录隔离")
  }
  Sys.setenv(M01_CODE_DIR=p$code,M01_CLEAN_DIR=p$clean,M01_ANALYSIS_OUTPUT=p$main)
  source(file.path(p$code,"Function_09_Main_analysis_M01.R"),encoding="UTF-8")
  p
}
s01_data <- function(cfg,p,i,ids=NULL) {
  f <- file.path(cfg$clean_dir,"MI_M01","MI_completed",sprintf("MI_%02d_analysis_M01.rds",i))
  a <- readRDS(f)
  master <- readRDS(file.path(p$main,"strategy_assignments","observed_score_master.rds"))
  if(is.list(master) && !is.data.frame(master)) stop("请检查 M01 observed_score_master 的版本")
  if(!is.null(ids)) a <- a[match(ids,as.character(a$participant_id)),]
  m07_strategy_one(a,master,cfg,i)
}
s01_predict <- function(engine,x,z,idx,cfg) {
  x[["ldl_2019_S01"]] <- z$ldl$Lcf
  targets <- list(index=list(overall=idx),meta=data.table::data.table(
    target_key="overall",target_type="scope",scenario="S01"))
  state <- m01_strategy_state(engine,x,targets,"2019_S01",cfg)
  m01_target_prediction(engine,state,targets,"overall")
}
s01_effects <- function(pred,sc,i,dfcom,scope="overall",weight=1) {
  key <- function(id) paste(id,c("2025","2019"),sep=":")
  d <- key(sc$id)
  terms <- list(risk_2019=setNames(1,d[2]),risk_2025=setNames(1,d[1]),
    rd_2025_vs_2019=setNames(c(1,-1),d),
    delta_rd_vs_s0=setNames(c(1,-1,-1,1),c(d,key("S0"))),
    delta_rd_vs_snt=setNames(c(1,-1,-1,1),c(d,key("SNT"))),
    gamma_rd_vs_parent=setNames(c(1,-1,-1,1),c(d,key(sc$parent))),
    gamma_rd_vs_reference=setNames(c(1,-1,-1,1),c(d,key(sc$reference))))
  rows <- lapply(names(terms),function(nm) {
    ww <- terms[[nm]]
    # 年份特有OAT仅输出适用年份的风险；跨年份RD和两种delta均不存在。
    if(!all(names(ww)%in%names(pred))) return(NULL)
    ref <- switch(nm,delta_rd_vs_s0="S0",delta_rd_vs_snt="SNT",
      gamma_rd_vs_parent=sc$parent,gamma_rd_vs_reference=sc$reference,"")
    ans <- m09_one_effect(pred,ww,nm,scope,sc$id,i,dfcom,reference=ref)
    ans$Q <- ans$Q*weight
    ans$U <- ans$U*weight^2
    ans
  })
  data.table::rbindlist(rows,fill=TRUE)
}

# 正式批次按插补和方向保存，整批完成才发布文件。单侧 ONE 无配对 RD。
s01_run <- function(root,family,scenarios,definitions=list(),mechanisms=FALSE,workers=3L) {
  p <- s01_setup(root)
  cfg <- m07_config()
  cfg$code_dir <- p$code
  cfg$clean_dir <- p$clean
  cfg$out_dir <- p$main
  cfg$times <- 15
  dict <- m07_rm_definitions()
  s01_check_dictionary(dict)
  sample_n <- as.integer(Sys.getenv("S01_TEST_N","0"))
  ids <- NULL
  if(sample_n>0L) {
    if(!nzchar(Sys.getenv("S01_OUTPUT"))) stop("小样本测试必须指定独立 S01_OUTPUT")
    master <- readRDS(file.path(p$main,"strategy_assignments","observed_score_master.rds"))
    set.seed(20261107)
    ids <- sort(sample(as.character(master$participant_id),sample_n))
    cfg$m <- min(cfg$m,as.integer(Sys.getenv("S01_TEST_M","2")))
  }
  if(is.function(scenarios)) scenarios <- scenarios(dict)
  if(sample_n==0L) {
    audit_path <- file.path(p$out,"quality_control","S01_M01_input_manifest.csv")
    if(!file.exists(audit_path)) stop("请先运行 00 校验 M01 父情景")
    input_manifest <- m07_manifest(cfg)
    audited <- data.table::fread(audit_path)
    if(!identical(input_manifest$md5,audited$md5)) stop("M01 插补文件变更，请重跑 00")
  }
  all_sc <- c(s01_originals(),scenarios)
  all_id <- vapply(all_sc,`[[`,character(1),"id")
  if(anyDuplicated(all_id) || any(!vapply(all_sc,`[[`,character(1),"parent") %in% names(s01_alias)))
    stop("情景ID重复或父策略不属于S0/SN/ST/SNT")
  dir <- file.path(p$out,"per_imputation",family)
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  workers <- s01_check_workers(workers, cfg$m)
  started <- proc.time()[["elapsed"]]
  message(sprintf("[%s] S01 %s 开始 | %d 份插补 | %d 个PSOCK工作进程 | PID %d",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), family, cfg$m, workers, Sys.getpid()))
  on.exit(message(sprintf("[%s] S01 %s 结束/中止 | 总耗时 %.1f 秒 (%.2f 分钟)",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), family,
    proc.time()[["elapsed"]]-started, (proc.time()[["elapsed"]]-started)/60)),add=TRUE)
  # Windows 使用 PSOCK。每个 worker 独立读写一份插补，避免复制整张队列给其他进程。
  if(workers>1L) {
    code <- Sys.getenv("S01_CODE_DIR",file.path(root,"LJW_code","Sensitivity_code","S01_different_RM_strategies"))
    # 每个工作进程只加载与当前主脚本一一对应的 Function_XX 文件。
    helper <- getOption("S01_CURRENT_FUNCTION", NULL)
    if(!is.character(helper) || length(helper)!=1L || !grepl("^Function_[0-9]{2}_.+_S01\\.R$",helper))
      stop("请在当前主脚本 source 前设置 options(S01_CURRENT_FUNCTION=basename(function_file))")
    cl <- parallel::makePSOCKcluster(workers,outfile="",setup_timeout=30)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(root,dir_code,helper) {
      Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
      assign("code",dir_code,envir=.GlobalEnv)
      source(file.path(dir_code,helper),encoding="UTF-8")
      s01_setup(root)
      data.table::setDTthreads(1L)
      NULL
    },root,normalizePath(code,winslash="/",mustWork=TRUE),helper)
  }
  # 每批最多 workers 份；工作进程实时打印情景进度，主进程打印完成率/耗时/预计剩余时间。
  timing <- list()
  batches <- split(seq_len(cfg$m),ceiling(seq_len(cfg$m)/workers))
  for(batch in batches) {
    args <- list(p=p,cfg=cfg,ids=ids,dict=dict,all_sc=all_sc,definitions=definitions,
      mechanisms=mechanisms,dir=dir,family=family)
    batch_start <- proc.time()[["elapsed"]]
    rows <- if(workers==1L) do.call(lapply,c(list(batch,s01_run_imputation),args)) else
      do.call(parallel::parLapplyLB,c(list(cl,batch,s01_run_imputation),args))
    timing <- c(timing,rows)
    m07_csv(data.table::rbindlist(timing),file.path(p$out,"quality_control",paste0("S01_",family,"_timing.csv")))
    s01_print_progress(family,length(timing),cfg$m,started,
      detail=sprintf("本批 %.1f 秒; worker=%d",proc.time()[["elapsed"]]-batch_start,workers))
  }
  invisible(dir)
}

s01_progress <- function(family,i,label) {
  message(sprintf("[%s] S01 %s | 插补 %02d | %s",format(Sys.time(),"%H:%M:%S"),family,i,label))
}

# 单份插补在一个进程内完成拟合、策略计算和结果写入。耗时包含读取和保存。
s01_run_imputation <- function(i,p,cfg,ids,dict,all_sc,definitions,mechanisms,dir,family) {
  started <- proc.time()[["elapsed"]]
  s01_progress(family,i,"开始")
  x <- s01_data(cfg,p,i,ids)
  model <- m08_prepare_model(x)
  fit <- m08_fit(model)
  engine <- m01_make_grid_engine(list(csc_fit=fit,complete_data_df=fit$df),model,cfg)
  idx <- seq_len(nrow(x))
  counts <- transitions <- effects <- domains <- assignments <- pred <- assigned <- list()
  # 第一步：同一插补内生成全部情景。随后才计算固定S0/SNT参照，确保没有前向引用缺失。
  for(j in seq_along(all_sc)) {
    sc <- all_sc[[j]]
    start_sc <- proc.time()[["elapsed"]]
    z <- setNames(lapply(s01_years,function(y) s01_classify(x,dict,sc,y,cfg,definitions)),s01_years)
    assigned[[sc$id]] <- z
    for(year in s01_years) {
      a <- z[[year]]
      if(is.null(a)) next
      k <- paste(sc$id,year,sep=":")
      if(sc$rule=="original") {
        key <- paste(year,s01_alias[[sc$id]],sep="_")
        stopifnot(identical(as.integer(a$upgrade),x[[paste0("upgrade_",key)]]),
          isTRUE(all.equal(a$ldl$Lcf,x[[paste0("Lcf_",key)]])))
      }
      if(!is.null(a$domain_flags)) {
        flags <- a$domain_flags[a$moderate,,drop=FALSE]
        nn <- c(colSums(flags),colSums(!flags),sum(rowSums(flags)==1),sum(rowSums(flags)>=2))
        domains[[k]] <- data.frame(scenario_id=sc$id,guideline=year,
          metric=c(paste0("D",1:6,"_positive"),paste0("D",1:6,"_negative"),"exactly_1","at_least_2"),
          n=nn,denominator=nrow(flags),proportion=if(nrow(flags))nn/nrow(flags) else NA_real_)
      }
      assignments[[k]] <- data.frame(participant_id=x$participant_id,
        score_based_moderate=a$moderate,rm_eligible_window=a$window,rm_upgrade=a$upgrade,
        final_risk=a$level,required_ldlc_reduction=a$ldl$D,intervention=a$ldl$A)
      attr(assignments[[k]],"rm_items") <- a$items
      pred[[k]] <- s01_predict(engine,x,a,idx,cfg)
    }
    s01_progress(family,i,sprintf("情景 %d/%d %s 完成 | %.1f 秒",j,length(all_sc),sc$id,
      proc.time()[["elapsed"]]-start_sc))
  }
  # 纯删除情景的个体级嵌套检验；2019实验室strict、2019 CRP-RES必须等于原父策略。
  if(all(c("SNT_REMOVE2","SNT_STRICT","SN_REMOVE2","SN_STRICT","ST_STRICT")%in%names(assigned))) {
    for(year in s01_years) {
      for(prefix in c("SN","SNT")) {
        original <- assigned[[prefix]][[year]]$upgrade
        remove2 <- assigned[[paste0(prefix,"_REMOVE2")]][[year]]$upgrade
        strict <- assigned[[paste0(prefix,"_STRICT")]][[year]]$upgrade
        if(any(remove2 & !original) || any(strict & !remove2))
          stop("RM集合嵌套失败: ",year,"/",prefix," 原始⊇REMOVE2⊇STRICT")
      }
      if(any(assigned$ST_STRICT[[year]]$upgrade & !assigned$ST[[year]]$upgrade))
        stop("严格删除不能增加ST升级人数")
    }
    for(id in c("ST_STRICT","ST_HSCRP_GT2_LE10","SNT_HSCRP_GT2_LE10")) {
      parent <- if(grepl("^ST_",id))"ST" else "SNT"
      a <- assigned[[id]][["2019"]]
      b <- assigned[[parent]][["2019"]]
      if(is.null(a) || any(a$upgrade!=b$upgrade) || any(abs(a$ldl$D-b$ldl$D)>1e-10))
        stop("2019应当不受ST-STRICT或CRP-RES改变: ",id)
    }
  }
  # 第二步：计算共同原始SNT参照与跨年配对差（含与S0、SNT的两类delta）。
  for(sc in all_sc) {
    z <- assigned[[sc$id]]
    for(year in s01_years) {
      a <- z[[year]]; k <- paste(sc$id,year,sep=":")
      if(is.null(a)) {
        counts[[k]] <- data.frame(scenario_id=sc$id,guideline=year,metric="not_applicable",
          n=NA_real_,denominator=NA_real_,value=NA_real_)
        next
      }
      parent <- assigned[[sc$parent]][[year]]
      counts[[k]] <- cbind(scenario_id=sc$id,guideline=year,
        s01_counts(a,parent,assigned$S0[[year]],assigned$SNT[[year]]))
      if(length(a$upgrade) != nrow(x)) stop("升级人数向量长度错误")
    }
    if(all(vapply(z,Negate(is.null),logical(1))))
      transitions[[sc$id]] <- cbind(scenario_id=sc$id,s01_transition(z[["2019"]],z[["2025"]]))
    effects[[sc$id]] <- s01_effects(pred,sc,i,engine$complete_df)
  }
  mechanism <- if(mechanisms) s01_mechanisms(x,engine,assigned,all_sc,cfg,i) else NULL
  m07_atomic(list(imputation_id=i,scenarios=all_sc,counts=counts,transitions=transitions,
    mechanism=mechanism,effects=effects,domains=domains,assignments=assignments,p90=definitions$p90),
    file.path(dir,sprintf("S01_%s_imp%02d.rds",family,i)))
  seconds <- proc.time()[["elapsed"]]-started
  s01_progress(family,i,sprintf("完成 | %.1f 秒",seconds))
  data.frame(direction=family,imputation_id=i,pid=Sys.getpid(),seconds=seconds)
}

# 正式规则由研究团队冻结后提供。函数接收原版触发矩阵，返回同维度收紧矩阵。
s01_definitions <- function(code) {
  p <- Sys.getenv("S01_DEFINITIONS",file.path(code,"config","definitions_approved_S01.R"))
  if(!file.exists(p)) return(list())
  e <- new.env(parent=globalenv())
  sys.source(p,envir=e)
  if(!exists("definitions",e,inherits=FALSE) || !is.list(e$definitions)) stop("定义脚本必须提供 definitions 列表")
  e$definitions
}

s01_observed_base <- function(root) {
  readRDS(file.path(s01_paths(root)$clean,"data_analysis_base_incomplete_M01.rds"))
}

# 05 在插补前按 M01 操作定义生成逐人的 *_observed。NA 是原始资料未能确认。
s01_observed_rule <- function(root) {
  base <- s01_observed_base(root)
  id <- as.character(base$participant_id)
  if(anyNA(id) || anyDuplicated(id)) stop("M01 原始证据 ID 缺失或重复")
  columns <- grep("^rm_.*_observed$",names(base),value=TRUE)
  observed <- as.data.frame(base)[,columns,drop=FALSE]
  rm(base) # PSOCK 只传输原始 RM 证据与 ID。
  function(x,year,flags) {
    # 05 的家族史与应激保留代理字段名，07 以 family、stress 登记。
    items <- colnames(flags)
    items[items=="family"] <- "family_proxy"
    items[items=="stress"] <- "stressor_proxy"
    cols <- paste0("rm_",items,"_observed")
    ix <- match(as.character(x$participant_id),id)
    if(anyNA(ix) || !all(cols %in% names(observed))) stop("M01 原始确认字段或 ID 不完整")
    evidence <- as.matrix(observed[ix,cols,drop=FALSE])
    if(any(!is.na(evidence) & !evidence %in% 0:1)) stop("M01 原始确认字段必须为 0/1/NA")
    flags[] <- as.integer(flags==1 & !is.na(evidence) & evidence==1)
    flags
  }
}

# 从主分析固定队列的原始观测计算一次 P90。先复算 P80 验证人口与分位数口径。
s01_townsend_p90 <- function(root) {
  base <- s01_observed_base(root)
  q <- stats::quantile(base$townsend_deprivation_index,c(.8,.9),na.rm=TRUE,type=7)
  saved <- unique(base$rm_townsend_q80)
  if(length(saved)!=1L || !isTRUE(all.equal(unname(q[1]),as.numeric(saved),tolerance=1e-10))) {
    stop("原始 Townsend P80 无法复现 M01 冻结值，请核对人口分布")
  }
  unname(q[2])
}

# 核验已保存 M01 预测与正式表的对应关系，同时登记插补文件指纹。
# 每份插补独立重放 M01 原始S0/SN/ST/SNT风险，主进程按插补编号确定性合并。
s01_validate_one <- function(i,cfg) {
  start <- proc.time()[["elapsed"]]
  message(sprintf("[%s] S01 VALIDATE | 插补 %02d/%02d 开始",format(Sys.time(),"%H:%M:%S"),i,cfg$m))
  b <- m08_bundle(cfg,i)
  rows <- lapply(names(s01_alias),function(id) {
    k <- s01_alias[[id]]
    pd <- setNames(lapply(s01_years,function(y)
      m08_prediction(b,paste(y,k,sep="_"),"scope__overall")),
      paste(id,s01_years,sep=":"))
    s01_effects(pd,s01_scenario(id,id),i,b$engine$complete_df)
  })
  seconds <- proc.time()[["elapsed"]]-start
  message(sprintf("[%s] S01 VALIDATE | 插补 %02d 完成 | %.1f 秒",
    format(Sys.time(),"%H:%M:%S"),i,seconds))
  list(effects=data.table::rbindlist(rows,fill=TRUE),
       timing=data.frame(direction="VALIDATE",imputation_id=i,pid=Sys.getpid(),seconds=seconds))
}

s01_validate <- function(root,workers=3L) {
  started <- proc.time()[["elapsed"]]
  on.exit(message(sprintf("[%s] S01 VALIDATE 结束/中止 | 总耗时 %.1f 秒",
    format(Sys.time(),"%H:%M:%S"),proc.time()[["elapsed"]]-started)),add=TRUE)
  p <- s01_setup(root)
  cfg <- m07_config()
  cfg$code_dir <- p$code
  cfg$clean_dir <- p$clean
  cfg$out_dir <- p$main
  manifest <- m07_manifest(cfg)
  cached <- readRDS(file.path(p$main,"strategy_assignments","input_manifest.rds"))
  stopifnot(identical(as.character(manifest$md5),as.character(cached$md5)))
  a <- m07_read(cfg,1L)
  n <- nrow(a)
  master <- readRDS(file.path(p$main,"strategy_assignments","observed_score_master.rds"))
  if(n!=nrow(master) || !setequal(as.character(a$participant_id),as.character(master$participant_id)))
    stop("M01已保存评分主表与07分析队列人数/ID不一致：",n)
  workers <- s01_check_workers(workers,cfg$m)
  message(sprintf("[%s] S01 VALIDATE 开始 | %d份插补 | %d个PSOCK工作进程",
    format(Sys.time(),"%H:%M:%S"),cfg$m,workers))
  if(workers>1L){
    code <- Sys.getenv("S01_CODE_DIR",file.path(root,"LJW_code","Sensitivity_code","S01_different_RM_strategies"))
    helper <- getOption("S01_CURRENT_FUNCTION","Function_00_Validate_inputs_S01.R")
    cl <- parallel::makePSOCKcluster(workers,outfile="",setup_timeout=30)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(root,code,helper){
      Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
      assign("code",code,envir=.GlobalEnv)
      source(file.path(code,helper),encoding="UTF-8")
      s01_setup(root)
      data.table::setDTthreads(1L)
      NULL
    },root,normalizePath(code,winslash="/",mustWork=TRUE),helper)
  }
  results <- vector("list",cfg$m)
  batches <- split(seq_len(cfg$m),ceiling(seq_len(cfg$m)/workers))
  for(batch in batches){
    rows <- if(workers==1L) lapply(batch,s01_validate_one,cfg=cfg) else
      parallel::parLapplyLB(cl,batch,s01_validate_one,cfg=cfg)
    for(k in seq_along(batch)) results[[batch[[k]]]] <- rows[[k]]
    completed <- sum(!vapply(results,is.null,logical(1)))
    s01_print_progress("VALIDATE",completed,cfg$m,started)
    times <- data.table::rbindlist(lapply(results,function(z) if(is.null(z))NULL else z$timing),fill=TRUE)
    m07_csv(times,file.path(p$out,"quality_control","S01_VALIDATE_timing.csv"))
  }
  z <- m09_pool(data.table::rbindlist(lapply(results,`[[`,"effects"),fill=TRUE),cfg)
  z <- z[time==15 & estimand %in% c("risk_2019","risk_2025","rd_2025_vs_2019")]
  ref <- data.table::fread(file.path(p$main,"strategy_risks_and_effects","pooled",
    "all_scenario_risks_and_effects.csv"))
  ref <- ref[scope=="overall" & time==15 & estimand %in% c("risk_2019","risk_2025","rd_2025_vs_2019")]
  ref$scenario <- names(s01_alias)[match(ref$scenario,s01_alias)]
  audit <- merge(z,ref[,.(scenario,estimand,m01_estimate=estimate,m01_lcl=lcl,m01_ucl=ucl)],
    by=c("scenario","estimand"))
  if(nrow(audit)!=12L || any(abs(audit$estimate-audit$m01_estimate)>1e-8) ||
     any(abs(audit$lcl-audit$m01_lcl)>1e-8) || any(abs(audit$ucl-audit$m01_ucl)>1e-8)) {
    stop("M01 父情景与主分析正式汇总表不一致")
  }
  m07_csv(manifest,file.path(p$out,"quality_control","S01_M01_input_manifest.csv"))
  m07_csv(audit,file.path(p$out,"quality_control","S01_M01_reference_replay.csv"))
  m07_csv(data.frame(n=n,moderate_2019=sum(a$score_2019_category=="M"),
    moderate_2025=sum(a$score2_2025_category=="M")),file.path(p$out,"quality_control","S01_M01_cohort_check.csv"))
  invisible(audit)
}

# 仅原始资料确认阳性且符合原版定义的 RM 可触发。与相应原版比较升级损失、新增干预、降幅及配对风险差。
s01_obs_scenarios <- function(dict) {
  lapply(c("SN","ST","SNT"),function(p) s01_scenario(paste0(p,"_OBSERVED_CONFIRMED"),p,"observed"))
}
