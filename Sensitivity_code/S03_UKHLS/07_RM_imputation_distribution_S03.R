# S03 07：仅对六项风险修饰因素（RM）缺失进行填补和分布统计
# 注：04纳排、05评分及06干预变化均已固定；07不改变其人数、ID、评分和结果。
# ========================= 前端可调参数 =========================
THREADS <- 5L                  # 默认并行线程（并行独立插补链）
M <- 5L                       # 多重插补套数
MAXIT <- 5L                   # 每套插补迭代次数
SEED <- 20261112L              # 固定随机种子
TEST_N <- 0L                   # 必须与05/06测试参数一致；0=正式全量
MIN_DONORS <- 10L              # 拟合样本至少具有的非缺失观测数
WAVE3_ACTIVITY <- "not_identifiable"  # Wave3无同定义成人体力活动问项，不外推
# ==================================================================
M <- as.integer(Sys.getenv("S03_M", as.character(M)))
MAXIT <- as.integer(Sys.getenv("S03_MAXIT", as.character(MAXIT)))
threads <- as.integer(Sys.getenv("S03_THREADS", as.character(THREADS)))
test_n <- as.integer(Sys.getenv("S03_ANALYSIS_TEST_N", as.character(TEST_N)))
stopifnot(threads >= 1L, test_n >= 0L, M >= 1L, MAXIT >= 1L, MIN_DONORS >= 2L)
if (!identical(WAVE3_ACTIVITY, "not_identifiable")) stop("Wave 3体力活动代理不可无依据跨波次填补")
arg <- grep("^--file=", commandArgs(), value=TRUE)
fr <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
start <- if (length(fr)) dirname(tail(fr,1)[[1]]) else if(length(arg)) dirname(sub("^--file=","",arg[1])) else getwd()
root <- normalizePath(start, winslash="/",mustWork=TRUE)
root_override <- Sys.getenv("S03_PROJECT_ROOT", "")
if (nzchar(root_override)) {
  root <- normalizePath(root_override, winslash = "/", mustWork = TRUE)
}

while(!dir.exists(file.path(root,"LJW_code","Main_code","M01_ASCVD"))) {
  up <- dirname(root);if(identical(up,root))stop("找不到项目根目录");root <- up
}
clean <- file.path(root,"LJW_outcome","data_clean","S03_UKHLS")
out <- file.path(root,"LJW_outcome","data_intervention","S03_UKHLS")
if(test_n>0L)out<-file.path(out,paste0("test_",test_n))
if (Sys.getenv("S03_SMOKING_MODE", "same_wave") == "wave2_prior") out <- file.path(out,"wave2_smoking_sensitivity")
for(p in c("rm_imputation","rm_imputation/per_imputation","rm_imputation/pooled","quality_control","logs"))
  dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
t0<-proc.time()[["elapsed"]]
logfile<-file.path(out,"logs","07_RM_imputation_distribution_S03_runtime.log")
log<-function(...) {txt<-sprintf("[%s +%.1fs] %s",format(Sys.time(),"%H:%M:%S"),proc.time()[["elapsed"]]-t0,paste0(...,collapse=""));message(txt);cat(txt,"\n",file=logfile,append=TRUE)}
write_out <- function(d, path) s03_write_csv(d, path)

# Excel-compatible UTF-8 CSV with BOM. Never rewrite source data in GBK.
s03_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  bytes <- readBin(tmp, what = "raw", n = file.info(tmp)$size)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(bytes, con)
  invisible(path)
}

RM <- c("rm_south_asian","rm_obesity_central","rm_hscrp_gt2",
        "rm_social_deprivation","rm_physical_inactivity_sport_proxy",
        "rm_chronic_inflammatory_arthritis_proxy")
# 连续源变量优先插补；其他RM直接插补二元状态。观测RM始终保持不变。
TARGET <- c("bmi","waist_cm","hscrp", "rm_south_asian",
            "rm_social_deprivation","rm_physical_inactivity_sport_proxy",
            "rm_chronic_inflammatory_arthritis_proxy")
BIN <- c("rm_south_asian","rm_social_deprivation","rm_physical_inactivity_sport_proxy",
         "rm_chronic_inflammatory_arthritis_proxy")
# Fixed observed clinical predictors: never imputed, never changed.
PREDICT <- c("age","sex","wave","current_smoking","tc","hdl","ldl","sbp","dbp","egfr")
# Conditional MICE: other RM source variables also predict missing RM information.
# Avoid deterministic descendants (binary obesity/hs-CRP) when their source variables are used.
USE_OTHER_RM_PREDICTORS <- TRUE
MIN_BINARY_EVENTS <- 5L   # each binary outcome class must have >= this many observations

rm_distribution <- function(d, id, stage) {
  do.call(rbind,lapply(RM,function(v) {
    x<-d[[v]];n<-nrow(d);known<-sum(!is.na(x));pos<-sum(x %in% 1);neg<-sum(x %in% 0)
    data.frame(imputation_id=id,stage=stage,variable=v,n_total=n,
               n_positive=pos,n_negative=neg,n_observed=known,n_missing=n-known,
               positive_pct_observed=if(known)100*pos/known else NA_real_,
               positive_pct_total=if(n)100*pos/n else NA_real_,
               missing_pct=if(n)100*(n-known)/n else NA_real_)
  }))
}
# 利用已知逻辑确定性判定肥胖：BMI>=30 或已知腰围达到性别阈值 => 阳性；
# 两项均可排除后才判阴性。保持04已观察的判定值。
resolve_obesity <- function(d, original = NULL) {
  bw<-ifelse(d$sex==1,102,88)
  yes <- (d$bmi >=30) | (d$waist_cm >= bw)
  no <- !is.na(d$bmi) & !is.na(d$waist_cm) & d$bmi<30 & d$waist_cm<bw
  derived <- ifelse(yes %in% TRUE,1,ifelse(no %in% TRUE,0,NA_real_))
  if (is.null(original)) return(derived)
  ifelse(!is.na(original),original,derived)
}

# 每个链独立运行，用于PSOCK并行，自动保留拟合诊断和事件记录。
run_chain <- function(task, data, method, predictor, where, maxit) {
  fit <- mice::mice(data,m=1L,maxit=maxit,method=method,predictorMatrix=predictor,
                    where=where,seed=task$seed,printFlag=FALSE,
                    remove.constant=TRUE,remove.collinear=TRUE)
  list(id=task$id,values=mice::complete(fit,1),trace_mean=fit$chainMean,
       trace_var=fit$chainVar,events=fit$loggedEvents)
}

if(!requireNamespace("mice",quietly=TRUE)) stop("07需要安装 mice 包")
sample_n <- as.integer(Sys.getenv("S03_SAMPLE_N", "0"))
if (is.na(sample_n) || sample_n < 0L) stop("S03_SAMPLE_N invalid")
if (sample_n > 0L) clean <- file.path(clean, paste0("smoke_test_", sample_n))
smoking_mode <- Sys.getenv("S03_SMOKING_MODE", "same_wave")
if (!smoking_mode %in% c("same_wave","wave2_prior")) stop("Unsupported smoking mode")
if (smoking_mode == "wave2_prior") clean <- file.path(clean,"wave2_smoking_sensitivity")
cohort_path<-file.path(clean,"04_eligible_cohort_S03.csv")
score_path<-file.path(out,"complete_case","05_score_strategy_participants_S03.csv")
if(!file.exists(cohort_path)||!file.exists(score_path))stop("需要先运行04和05")
full<-read.csv(cohort_path,na.strings="",check.names=FALSE)
score<-read.csv(score_path,na.strings="",check.names=FALSE)
if(!all(c("pidp","wave",RM,TARGET,PREDICT)%in%names(full)))stop("04队列缺少RM/源变量/预测变量")
if(anyDuplicated(full$pidp)||anyDuplicated(score$pidp))stop("参与者ID重复")
# 以05评分结果为唯一分析名单；避免07测试与05测试分母错位。
index<-match(score$pidp,full$pidp)
if(anyNA(index)||length(index)!=nrow(score))stop("05 ID无法一一映射04")
d<-full[index,,drop=FALSE]
stopifnot(identical(as.numeric(d$pidp),as.numeric(score$pidp)))
log("开始07；固定队列n=",nrow(d),"；线程=",threads,"；M=",M,"；每套迭代=",MAXIT)
pre<-rm_distribution(d,0L,"observed")
write_out(pre,file.path(out,"quality_control","07_RM_missingness_before_S03.csv"))

# 将wave3体力活动系统性缺失显式标记；这些值不能从Wave2问项无条件外推。
structural<-d$wave==3 & is.na(d$rm_physical_inactivity_sport_proxy)
# 根据已观察信息更新肥胖确定性可判定项（不覆盖04已知状态）
old_rm<-d[RM]
d$rm_obesity_central<-resolve_obesity(d,old_rm$rm_obesity_central)
post_deterministic<-rm_distribution(d,0L,"after_deterministic")
write_out(post_deterministic,file.path(out,"quality_control","07_RM_missingness_after_deterministic_S03.csv"))

# MICE建模数据：所有评分变量仅作预测，不对其插补。
vars<-unique(c(PREDICT,TARGET))
model<-d[vars]
for(v in BIN)model[[v]]<-factor(model[[v]],levels=0:1)
for(v in c("sex","wave","current_smoking"))model[[v]]<-factor(model[[v]])
methods<-setNames(rep("",ncol(model)),names(model))
predmat<-matrix(0L,ncol(model),ncol(model),dimnames=list(names(model),names(model)))
where<-is.na(model)
where[,PREDICT]<-FALSE
if(anyNA(model[PREDICT])) {
  # 只允许完整的固定预测字段进入模型，避免不插补的预测变量造成拟合失败。
  dropped<-PREDICT[vapply(model[PREDICT],anyNA,logical(1))]
  if(length(dropped))log("固定预测变量存在缺失，排除预测矩阵：",paste(dropped,collapse=","))
} else dropped<-character()
valid_preds<-setdiff(PREDICT,dropped)
# 只在同一波次确有阳性/阴性观测支持时，才插补该波次的RM。
status<-do.call(rbind,lapply(TARGET,function(v){
  tab<-do.call(rbind,lapply(sort(unique(d$wave)),function(w){
    ix<-d$wave==w
    obs<-model[[v]][ix];obs<-obs[!is.na(obs)]
    k<-length(obs);distinct<-length(unique(obs))
    binary_target <- v %in% BIN
    class_ok <- if(binary_target) {
      counts <- table(factor(as.character(obs), levels=c("0","1")))
      all(counts >= MIN_BINARY_EVENTS)
    } else TRUE
    ok <- k>=MIN_DONORS && distinct>=2L && class_ok
    if(v=="rm_physical_inactivity_sport_proxy"&&w==3)ok<-FALSE
    data.frame(variable=v,wave=w,n_observed=k,n_missing=sum(is.na(model[[v]][ix])),
               distinct_observed=distinct,imputation_supported=ok,
               reason=if(ok)"modelled" else if(v=="rm_physical_inactivity_sport_proxy"&&w==3)"wave3_no_comparable_adult_item" else "insufficient_same_wave_observed")
  }))
  tab
}))
write_out(status,file.path(out,"quality_control","07_RM_imputation_support_by_wave_S03.csv"))
active<-character()
for(v in TARGET){
  if(!anyNA(model[[v]]))next
  if(!any(status$variable==v & status$imputation_supported))next
  active<-c(active,v)
  methods[v]<-if(v %in% BIN)"logreg" else "pmm"
  for(w in sort(unique(d$wave))){
    if(!status$imputation_supported[status$variable==v & status$wave==w])where[d$wave==w,v]<-FALSE
  }
  predmat[v,valid_preds] <- 1L
}
# Cross-RM predictors are accepted only if all their missing cells will be imputed
# (or the predictor is fully observed). Otherwise a structural NA could enter the
# conditional design matrix; those columns are excluded rather than faked as 0.
rm_predictors <- character()
if (USE_OTHER_RM_PREDICTORS && length(active)) {
  eligible_predictor <- function(v) {
    if (!v %in% names(model)) return(FALSE)
    if (any(is.na(model[[v]]) & !where[,v])) return(FALSE)
    if (length(unique(stats::na.omit(model[[v]]))) < 2L) return(FALSE)
    TRUE
  }
  rm_predictors <- TARGET[vapply(TARGET, eligible_predictor, logical(1))]
  for (v in active) predmat[v, setdiff(rm_predictors, v)] <- 1L
}
# Enforce no unavailable (not modelled) predictors. Fixed baseline predictors remain fixed.
for(v in active) {
  for (other in TARGET) {
    if (other != v && (other %in% setdiff(TARGET,active) ||
        any(is.na(model[[other]]) & !where[,other]))) predmat[v,other] <- 0L
  }
}
log("RM cross-predictors: ", if(length(rm_predictors)) paste(rm_predictors,collapse=", ") else "none eligible")
# Protect observed BMI, waist and hs-CRP values. Each variable's observed cells
# are not imputed. Their binary RM descendants are derived only after each MI chain.
# 预测器如为单一取值，也排除，避免共线性和常量问题。
usable<-vapply(model,function(x)length(unique(na.omit(x)))>1L,logical(1))
predmat[,!usable]<-0L
diag(predmat)<-0L
write_out(data.frame(variable=names(methods),method=unname(methods),predmat,check.names=FALSE),
          file.path(out,"quality_control","07_RM_predictor_matrix_S03.csv"))
write_out(data.frame(target=rep(rownames(predmat),each=ncol(predmat)),
                     predictor=rep(colnames(predmat),times=nrow(predmat)),
                     used=as.integer(as.vector(t(predmat))),
                     stringsAsFactors=FALSE),
          file.path(out,"quality_control","07_RM_predictor_edges_S03.csv"))
log("RM需插补字段：",if(length(active))paste(active,collapse=", ") else "无可支持的目标")
log("Wave3体力活动结构性缺失：",sum(structural),"人；不会将其按阴性填入")

seeds<-SEED+seq_len(M)
tasks<-lapply(seq_len(M),function(i)list(id=i,seed=seeds[i]))
if(length(active)){
  nworkers<-min(threads,M)
  if(nworkers>1L){
    cl<-parallel::makePSOCKcluster(nworkers,setup_timeout=30)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterEvalQ(cl,requireNamespace("mice",quietly=TRUE))
  }
  chains<-vector("list",M)
  for(ids in split(seq_len(M),ceiling(seq_len(M)/max(1L,nworkers)))){
    if(nworkers>1L)chains[ids]<-parallel::parLapplyLB(cl,tasks[ids],run_chain,
      data=model,method=methods,predictor=predmat,where=where,maxit=MAXIT)
    else chains[ids]<-lapply(tasks[ids],run_chain,data=model,method=methods,predictor=predmat,where=where,maxit=MAXIT)
    log("插补链完成：",max(ids),"/",M,"；已用时",round(proc.time()[["elapsed"]]-t0,1),"秒")
  }
} else chains<-lapply(tasks,function(task)list(id=task$id,values=model,trace_mean=NULL,trace_var=NULL,events=NULL))

completed<-vector("list",M)
checks<-vector("list",M)
distributions<-vector("list",M)
trace_tables<-list(); event_tables<-list()
for(i in seq_len(M)){
  x<-d
  vals<-chains[[i]]$values
  for(v in active){
    original<-x[[v]]
    filled<-vals[[v]]
    if(is.factor(filled))filled<-as.numeric(as.character(filled))
    replace<-is.na(original) & where[,v]
    x[[v]][replace]<-filled[replace]
  }
  # 被插补的BMI/腰围、hsCRP必须先恢复RM，保留原有确定的RM状态。
  x$rm_obesity_central<-resolve_obesity(x,old_rm$rm_obesity_central)
  predicted_hscrp<-ifelse(is.na(x$hscrp),NA_real_,as.numeric(x$hscrp>2))
  x$rm_hscrp_gt2<-ifelse(!is.na(old_rm$rm_hscrp_gt2),old_rm$rm_hscrp_gt2,predicted_hscrp)
  # Continuous RM source measurements are also immutable where observed.
  source_mismatch <- vapply(c("bmi","waist_cm","hscrp"),function(v){
    observed <- !is.na(d[[v]])
    sum(is.na(x[[v]][observed]) | x[[v]][observed] != d[[v]][observed])
  },integer(1))
  if(any(source_mismatch>0L))stop("Observed RM source measurements changed: ",
    paste(names(source_mismatch)[source_mismatch>0L],collapse=","))
  # 不覆盖所有原本已观察到的六项RM值。
  mismatch<-vapply(RM,function(v){ix<-!is.na(old_rm[[v]]);sum(is.na(x[[v]][ix])|x[[v]][ix]!=old_rm[[v]][ix])},integer(1))
  if(any(mismatch))stop("RM原观测值被改变：",paste(names(mismatch)[mismatch>0],collapse=","))
  if(anyNA(x$pidp)||!identical(as.numeric(x$pidp),as.numeric(d$pidp)))stop("插补后ID变化")
  if(any(!is.na(x$rm_physical_inactivity_sport_proxy[structural])))stop("Wave3结构性缺失被填补")
  diff_fixed<-vapply(PREDICT,function(v) !identical(x[[v]],d[[v]]),logical(1))
  if(any(diff_fixed))stop("纳排或SCORE预测信息被修改")
  completed[[i]]<-x
  distributions[[i]]<-rm_distribution(x,i,"post_imputation")
  checks[[i]]<-do.call(rbind,lapply(RM,function(v){
    orig<-old_rm[[v]];current<-x[[v]]
    data.frame(imputation_id=i,variable=v,n_observed_changed=0L,
               n_newly_filled=sum(is.na(orig)&!is.na(current)),
               n_remaining_missing=sum(is.na(current)),n_structural_unfilled=if(v=="rm_physical_inactivity_sport_proxy")sum(is.na(current)&structural) else 0L)
  }))
  if(!is.null(chains[[i]]$trace_mean)){
    tmp<-as.data.frame.table(chains[[i]]$trace_mean,responseName="chain_mean")
    names(tmp)[1:3]<-c("variable","iteration","chain")
    tmp$chain_variance<-as.vector(chains[[i]]$trace_var);tmp$imputation_id<-i
    trace_tables[[length(trace_tables)+1L]]<-tmp
  }
  if(!is.null(chains[[i]]$events) && nrow(chains[[i]]$events))
    event_tables[[length(event_tables)+1L]]<-data.frame(imputation_id=i,chains[[i]]$events)
  log("输出检查",i,"/",M)
}
names(completed)<-sprintf("imp%02d",seq_len(M))
saveRDS(completed,file.path(out,"rm_imputation","07_RM_imputed_cohorts_S03.rds"))
write_out(do.call(rbind,distributions),file.path(out,"rm_imputation","per_imputation","07_RM_distribution_per_imputation_S03.csv"))
write_out(do.call(rbind,checks),file.path(out,"quality_control","07_RM_imputation_checks_S03.csv"))
if(length(trace_tables))write_out(do.call(rbind,trace_tables),file.path(out,"quality_control","07_RM_chain_trace_S03.csv"))
if(length(event_tables))write_out(do.call(rbind,event_tables),file.path(out,"quality_control","07_RM_logged_events_S03.csv"))
pooled<-aggregate(do.call(rbind,distributions)[c("n_positive","n_negative","n_observed","n_missing","positive_pct_observed","positive_pct_total","missing_pct")],
  by=do.call(rbind,distributions)["variable"],FUN=function(z)if(all(is.na(z)))NA_real_ else mean(z,na.rm=TRUE))
pooled$n_total<-nrow(d);pooled$m<-M
# Summarize between-imputation spread rather than hiding uncertainty behind MI means.
all_dist <- do.call(rbind,distributions)
spread <- do.call(rbind,lapply(split(all_dist,all_dist$variable),function(z)
  data.frame(variable=z$variable[1L],m=nrow(z),
    positive_n_sd=stats::sd(z$n_positive),
    positive_n_min=min(z$n_positive),positive_n_max=max(z$n_positive),
    positive_pct_sd=stats::sd(z$positive_pct_total),
    remaining_missing_max=max(z$n_missing))))
write_out(spread,file.path(out,"quality_control","07_RM_between_imputation_spread_S03.csv"))
write_out(pooled,file.path(out,"rm_imputation","pooled","07_RM_distribution_mean_S03.csv"))
write_out(data.frame(imputation_id=seq_len(M),seed=seeds,maxit=MAXIT,threads=threads),
          file.path(out,"quality_control","07_RM_run_settings_S03.csv"))
log("07完成；固定分母=",nrow(d),"；总耗时=",round(proc.time()[["elapsed"]]-t0,1),"秒")
