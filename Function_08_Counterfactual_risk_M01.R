# MODULE 08 | Cox + grid CIF + analytical influence functions
source(file.path(Sys.getenv("M01_CODE_DIR", "."),"Function_07_Guideline_strategies_M01.R"),local=FALSE)
if(!requireNamespace("survival",quietly=TRUE))stop("Install survival")

m08_outcome_terms <- function() c("age_centered","sex","ethnicity","region","ses_class","rm_deprivation","current_smoker","rm_inactive","alcohol_score","diet_score","rm_obesity","systolic_blood_pressure","diastolic_blood_pressure","antihypertensive_medication_baseline","ldl_ascvd","hdl_cholesterol","triglycerides","egfr_2021","family_history_cvd","chronic_liver_disease","chronic_respiratory_disease","rm_inflammatory","rm_psychiatric","rm_af","rm_lpa","rm_crp")
m08_prepare_model <- function(x) {
  x<-data.table::copy(data.table::as.data.table(x))
  if(!"ethnicity"%in%names(x) && "rm_source_south_asian"%in%names(x))x[,ethnicity:=rm_source_south_asian]
  if(!"rm_inactive"%in%names(x) && "rm_source_inactive"%in%names(x))x[,rm_inactive:=rm_source_inactive]
  if(!"ascvd_event"%in%names(x))x[,ascvd_event:=as.integer(event_status==1L)]
  if(!"competing_death_event"%in%names(x))x[,competing_death_event:=as.integer(event_status==2L)]
  if(!"current_smoker"%in%names(x))x[,current_smoker:=as.integer(as.character(smoking_status)=="current")]
  x[,age_centered:=as.numeric(age)-mean(age)]
  x[,ldl_ascvd:=as.numeric(ldl_cholesterol)]
  x[,ldl_death:=as.numeric(ldl_cholesterol)]
  m07_assert(x,c(m08_outcome_terms(),"ldl_death","followup_years","event_status","participant_id"),"08 unified outcome model")
  for(v in c("sex","region","ses_class","ethnicity")){
    levels<-switch(v,sex=c("female","male"),region=c("England","Wales","Scotland"),ses_class=c("low","high"),ethnicity=c("0","1"))
    x[[v]]<-factor(as.character(x[[v]]),levels=levels)
    if(anyNA(x[[v]]))stop("08 invalid/missing model factor ",v)
  }
  continuous<-c("age_centered","systolic_blood_pressure","diastolic_blood_pressure","ldl_ascvd","ldl_death","hdl_cholesterol","triglycerides","egfr_2021","followup_years")
  binary<-c("rm_deprivation","current_smoker","rm_inactive","alcohol_score","diet_score","rm_obesity","antihypertensive_medication_baseline","family_history_cvd","chronic_liver_disease","chronic_respiratory_disease","rm_inflammatory","rm_psychiatric","rm_af","rm_lpa","rm_crp")
  for(v in continuous) {x[[v]]<-as.numeric(as.character(x[[v]]));if(any(!is.finite(x[[v]])))stop("08 missing/nonfinite: ",v)}
  for(v in binary) {x[[v]]<-as.integer(as.character(x[[v]]));if(anyNA(x[[v]])||any(!x[[v]]%in%0:1))stop("08 invalid 0/1 covariate: ",v)}
  x[,event_status:=as.integer(event_status)]
  if(anyNA(x$event_status)||any(!x$event_status %in% 0:2)||any(x$followup_years<=0))stop("08 event coding invalid")
  x[,ascvd_event:=as.integer(event_status==1L)];x[,competing_death_event:=as.integer(event_status==2L)]
  x
}
m08_fit <- function(x) {
  if(!requireNamespace("survival",quietly=TRUE))stop("Install survival")
  vars<-m08_outcome_terms()
  f1<-stats::as.formula(paste0("survival::Surv(followup_years,ascvd_event) ~ ",paste(vars,collapse="+")))
  f2<-stats::as.formula(paste0("survival::Surv(followup_years,competing_death_event) ~ ",paste(replace(vars,which(vars=="ldl_ascvd"),"ldl_death"),collapse="+")))
  a<-m07_timed("prediction","ASCVD Cox fit",survival::coxph(f1,data=x,x=TRUE,y=TRUE,model=FALSE,ties="breslow",singular.ok=FALSE))
  b<-m07_timed("prediction","competing death Cox fit",survival::coxph(f2,data=x,x=TRUE,y=TRUE,model=FALSE,ties="breslow",singular.ok=FALSE))
  if(anyNA(stats::coef(a))||anyNA(stats::coef(b)))stop("Cox model rank deficient")
  list(models=list(a,b),n=nrow(x),df=max(1,nrow(x)-max(length(coef(a)),length(coef(b)))-1),formulas=c(deparse(f1),deparse(f2)))
}
m08_targets <- function(x) {
  scopes<-m07_scope_spec(x)
  scenarios<-m07_rms();idx<-list();rows<-list()
  add<-function(key,ii,typ,scope,scenario="shared",transition_id=NA_character_) {
    if(length(ii)==0L)return(invisible(NULL))
    idx[[key]]<<-as.integer(ii)
    rows[[length(rows)+1L]]<<-data.table(target_key=key,target_type=typ,scope=scope,scenario=scenario,transition_id=transition_id,
                                        n_standardized=length(ii))
  }
  for(scope in names(scopes))add(paste0("scope__",scope),which(scopes[[scope]]),"scope",scope)
  for(k in scenarios) for(scope in names(scopes)) for(tr in m07_transitions()){
    ii<-which(scopes[[scope]] & x[[paste0("transition_",k)]]==tr)
    add(paste0("transition__",k,"__",scope,"__",tr),ii,"transition",scope,k,tr)
  }
  list(index=idx,meta=data.table::rbindlist(rows))
}
m01_bin_sum <- function(x, bin, nbin) {
  if (is.null(dim(x))) x <- matrix(x,ncol=1L)
  out <- matrix(0,nbin,ncol(x))
  if (nrow(x)) {
    z <- rowsum(x,bin,reorder=FALSE)
    out[as.integer(rownames(z)),] <- z
  }
  out
}

m01_event_bins <- function(time,cuts) {
  j <- findInterval(time,cuts)
  on <- j>0L & abs(time-cuts[pmax(j,1L)])<1e-12
  j[!on] <- j[!on]+1L
  if (any(j<1L | j>length(cuts))) stop("事件超出积分网格")
  j
}

m01_cox_source <- function(fit,cuts,expected_time,expected_status) {
  # 输入一个非分层、无权重、右删失Breslow Cox；输出流式IF所需充分量。
  # 风险集从完整随访计算；只截断积分时点，不重新删失/重拟合Cox。
  if (!identical(fit$method,"breslow") || !is.null(fit$strata) ||
      ncol(fit$y)!=2L || anyNA(coef(fit)) || any(!is.finite(coef(fit))))
    stop("网格IF要求非分层、满秩、右删失Breslow Cox")
  if (!is.null(fit$weights) && any(fit$weights!=1)) stop("当前引擎不接受加权Cox")
  time <- as.numeric(fit$y[,1L]); status <- as.integer(fit$y[,2L])
  if (!identical(time,as.numeric(expected_time)) || !identical(status,as.integer(expected_status)))
    stop("模型训练顺序/结局与MI数据不一致")
  x <- as.matrix(fit$x); beta <- coef(fit)
  if (ncol(x)!=length(beta) || !identical(colnames(x),names(beta))) stop("Cox矩阵列不一致")
  risk <- exp(fit$linear.predictors)
  if (any(!is.finite(risk)) || any(risk<=0)) stop("Cox相对危险无效")
  et <- sort(unique(time[status==1L & time<=max(cuts)]))
  if (!length(et)) stop("报告期内该原因没有事件，不能推断")
  ord <- order(time); first <- match(et,time[ord])
  # 同时点的事件和删失均处于该时点风险集中；累计从该时点第一行开始。
  s0 <- rev(cumsum(rev(risk[ord])))[first]
  d <- tabulate(match(time[status==1L & time<=max(cuts)],et),length(et))
  hazard <- d/s0
  means_z <- matrix(0,length(et),ncol(x))
  for (p in seq_len(ncol(x)))
    means_z[,p] <- rev(cumsum(rev(risk[ord]*x[ord,p])))[first]/s0-fit$means[p]
  bin <- m01_event_bins(et,cuts)
  iid_beta <- as.matrix(stats::residuals(fit,type="dfbeta"))
  if (!identical(dim(iid_beta),dim(x)) || any(!is.finite(iid_beta))) stop("无效beta IF")
  event_position <- integer(length(time))
  own <- which(status==1L & time<=max(cuts))
  event_position[own] <- match(time[own],et)
  list(x=x,means=fit$means,beta=beta,risk=risk,iid_beta=iid_beta,
    event_times=et,s0=s0,hazard=hazard,bin=bin,event_position=event_position,
    exit_position=findInterval(time,et),
    dA=as.numeric(m01_bin_sum(hazard,bin,length(cuts))),
    C=m01_bin_sum(means_z*hazard,bin,length(cuts)))
}

m01_make_grid_engine <- function(context,x,cfg) {
  step <- cfg$grid_days/365.25
  regular <- if(step<=max(cfg$times))seq(step,max(cfg$times),by=step)else numeric()
  cuts <- sort(unique(c(regular,cfg$times)))
  f <- context$csc_fit$models
  a <- m01_cox_source(f[[1L]],cuts,x$followup_years,x$ascvd_event)
  b <- m01_cox_source(f[[2L]],cuts,x$followup_years,x$competing_death_event)
  if (!"ldl_ascvd" %in% names(a$beta)) stop("缺少可干预的LDL主效应")
  if (any(grepl(":",names(a$beta),fixed=TRUE))) stop("含交互项时需扩展LDL设计矩阵导数")
  list(source1=a,source2=b,cuts=cuts,times=cfg$times,n=nrow(x),
    report_index=match(cfg$times,cuts),participant_id_md5=m07_hash_object(as.character(x$participant_id)),
    complete_df=context$complete_data_df,grid_days=cfg$grid_days)
}

m01_strategy_state <- function(engine,x,targets,strategy,cfg) {
  # 输入模型充分量及一个策略；分块累加梯度，不创建n×target×事件时间IF张量。
  a <- engine$source1; b <- engine$source2
  scenario <- sub("^[0-9]+_","",strategy)
  meta <- targets$meta
  # 使用显式局部副本，避免data.table列名与函数变量同名时选中所有情景。
  if (strategy!="natural") {
    selected_scenario <- scenario
    meta <- targets$meta[target_type=="scope" | get("scenario")==selected_scenario]
  }
  keys <- meta$target_key; index <- targets$index[keys]
  sizes <- lengths(index)
  W <- Matrix::sparseMatrix(i=unlist(index,use.names=FALSE),
    j=rep(seq_along(keys),sizes),x=rep(1/sizes,sizes),dims=c(engine$n,length(keys)))
  J <- length(engine$cuts); H <- length(engine$times); K <- length(keys)
  G1 <- G2 <- array(0,c(K,J,H))
  B1 <- array(0,c(K,length(a$beta),H)); B2 <- array(0,c(K,length(b$beta),H))
  risk <- matrix(0,engine$n,H)
  h1 <- c(0,head(cumsum(a$dA),-1L)); h2 <- c(0,head(cumsum(b$dA),-1L))
  delta <- if(strategy=="natural")rep(0,engine$n) else x[[paste0("ldl_",strategy)]]-x$ldl_cholesterol
  ldl_column <- match("ldl_ascvd",names(a$beta))
  chunk_starts<-seq.int(1L,engine$n,by=cfg$prediction_chunk)
  for (chunk_i in seq_along(chunk_starts)) {
    start<-chunk_starts[chunk_i]
    if(chunk_i==1L || chunk_i%%5L==0L || chunk_i==length(chunk_starts))
      m07_log("08",strategy," | subjects ",min(engine$n,start+cfg$prediction_chunk-1L),"/",engine$n,
        " (chunk ",chunk_i,"/",length(chunk_starts),")")
    ix <- seq.int(start,min(engine$n,start+cfg$prediction_chunk-1L))
    r1 <- a$risk[ix]*exp(a$beta[ldl_column]*delta[ix]); r2 <- b$risk[ix]
    z1 <- sweep(a$x[ix,,drop=FALSE],2L,a$means,"-")
    z1[,ldl_column] <- z1[,ldl_column]+delta[ix]
    z2 <- sweep(b$x[ix,,drop=FALSE],2L,b$means,"-")
    wc <- W[ix,,drop=FALSE]
    S <- exp(-tcrossprod(r1,h1)-tcrossprod(r2,h2))
    F <- matrix(0,length(ix),J)
    now <- u1 <- u2 <- numeric(length(ix))
    for (j in seq_len(J)) {
      increment <- S[,j]*r1*a$dA[j]
      now <- now+increment; F[,j] <- now
      u1 <- u1+increment*(1-r1*h1[j])
      u2 <- u2-increment*r2*h2[j]
      h <- match(j,engine$report_index)
      if (!is.na(h)) {
        risk[ix,h] <- now
        B1[,,h] <- B1[,,h]+as.matrix(Matrix::crossprod(wc,z1*u1))
        B2[,,h] <- B2[,,h]+as.matrix(Matrix::crossprod(wc,z2*u2))
      }
    }
    for (h in seq_len(H)) {
      jj <- seq_len(engine$report_index[h])
      tail <- risk[ix,h]-F[,jj,drop=FALSE]
      G1[,jj,h] <- G1[,jj,h]+as.matrix(Matrix::crossprod(wc,(S[,jj,drop=FALSE]-tail)*r1))
      G2[,jj,h] <- G2[,jj,h]-as.matrix(Matrix::crossprod(wc,tail*r2))
    }
  }
  if (any(!is.finite(risk)) || any(risk<0 | risk>1))
    stop("网格CIF超出[0,1]，请缩小grid_days；禁止裁切掩盖积分问题")
  # Breslow累计基线的beta依赖：C[j]=sum_bin mean(Z at risk)*dA。
  for (h in seq_len(H)) {
    B1[,,h] <- B1[,,h]-G1[,,h]%*%a$C
    B2[,,h] <- B2[,,h]-G2[,,h]%*%b$C
  }
  estimate <- as.matrix(Matrix::crossprod(W,risk))
  list(strategy=strategy,keys=keys,estimate=estimate,risk=risk,
       G1=G1,G2=G2,B1=B1,B2=B2,times=engine$times)
}

m01_weighted_baseline_iid <- function(source,g) {
  # 精确处理每个源事件的风险集，而不是把参与者退出时间舍入到网格边界。
  w <- g[source$bin]
  cumulative <- c(0,cumsum(w*source$hazard/source$s0))
  out <- -source$risk*cumulative[source$exit_position+1L]
  own <- which(source$event_position>0L)
  out[own] <- out[own]+w[source$event_position[own]]/source$s0[source$event_position[own]]
  out
}

m01_target_prediction <- function(engine,state,targets,key) {
  k <- match(key,state$keys)
  if (is.na(k)) stop("策略缓存没有该目标：",state$strategy,"/",key)
  ix <- targets$index[[key]]; n <- length(ix)
  if (!n) stop("空目标不能生成IF")
  H <- length(engine$times); iid <- matrix(0,engine$n,H)
  for (h in seq_len(H)) {
    iid[,h] <- as.numeric(engine$source1$iid_beta%*%state$B1[k,,h])+
      as.numeric(engine$source2$iid_beta%*%state$B2[k,,h])+
      m01_weighted_baseline_iid(engine$source1,state$G1[k,,h])+
      m01_weighted_baseline_iid(engine$source2,state$G2[k,,h])
    iid[ix,h] <- iid[ix,h]+(state$risk[ix,h]-state$estimate[k,h])/n
  }
  if (any(!is.finite(iid))) stop("非有限网格IF")
  if (max(abs(colSums(iid)))>1e-6)
    stop("IF未中心化；检查Cox收敛、风险集或ID顺序，禁止继续合并")
  list(estimate=as.numeric(state$estimate[k,]),iid=iid,times=engine$times,n_standardized=n)
}


m08_predict_mi <- function(m,cfg) {
  out<-m07_dir(cfg,"08");targetpath<-file.path(out,sprintf("bundle_m%02d.rds",m))
  p07<-file.path(m07_dir(cfg,"07"),sprintf("m%02d.rds",m))
  if(!file.exists(p07))stop("07 result missing: MI",m)
  signature<-m08_signature(cfg,m)
  cached<-file.path(out,sprintf("signature_m%02d.rds",m))
  if(file.exists(targetpath)&&file.exists(cached)&&identical(readRDS(cached),signature)) {
    m07_log("08", "MI",m," cache validated; skipping")
    return(data.table(imputation_id=m,cache_reused=TRUE))
  }
  x<-readRDS(p07)
  model_data<-m08_prepare_model(x)
  model_path<-file.path(out,sprintf("model_m%02d.rds",m))
  model_sig<-paste(m07_hash_object(model_data[,c("participant_id","event_status","followup_years",m08_outcome_terms(),"ldl_death"),with=FALSE]),paste(m08_outcome_terms(),collapse=","),
    m07_md5(file.path(Sys.getenv("M01_CODE_DIR"),"Function_08_Counterfactual_risk_M01.R")),as.character(utils::packageVersion("survival")),sep="|")
  sigmodel<-paste0(model_path,".signature.rds")
  if(file.exists(model_path)&&file.exists(sigmodel)&&identical(readRDS(sigmodel),model_sig)) {
    fit<-readRDS(model_path);m07_log("08","MI",m," model cache hit")
  } else {
    fit<-m07_timed("08",paste0("MI",m," Cox fits"),m08_fit(model_data))
    m07_atomic(fit,model_path);m07_atomic(model_sig,sigmodel)
  }
  ctx<-list(csc_fit=fit,complete_data_df=fit$df)
  eng<-m07_timed("08",paste0("MI",m," grid IF preparations"),m01_make_grid_engine(ctx,model_data,cfg))
  targets<-m08_targets(x)
  rm(model_data);gc(FALSE)
  states<-list()
  # Scenario columns in 07 are Lcf_year_scenario; old engine expects ldl_year_scenario
  for(state in m07_states()){
    m07_log("08","MI",m," counterfactual ",state," (",length(states)+1L,"/",length(m07_states()),")")
    col<-if(state!="natural")paste0("ldl_",state) else "natural"
    if(state!="natural")x[[col]]<-x[[paste0("Lcf_",state)]]
        state_path<-file.path(out,'state_checkpoints',sprintf('mi%02d_%s.rds',m,state))
    state_sig<-paste(signature,state,sep='|');sig_path<-paste0(state_path,'.sig')
    if(file.exists(state_path)&&file.exists(sig_path)&&identical(readRDS(sig_path),state_sig)){
      states[[state]]<-readRDS(state_path);m07_log('08','MI',m,' ',state,' checkpoint reused')
    }else{
      states[[state]]<-m07_timed('08',paste0('MI',m,' ',state,' CIF and IF gradients'),m01_strategy_state(eng,x,targets,state,cfg))
      m07_atomic(states[[state]],state_path);m07_atomic(state_sig,sig_path)
    }
    if(state!="natural")x[[col]]<-NULL
  }
  eng$source1$x<-NULL;eng$source2$x<-NULL
  m07_atomic(list(engine=eng,targets=targets,states=states,mi=m,signature=signature),targetpath)
  m07_atomic(signature,cached)
  data.table(imputation_id=m,cache_reused=FALSE,n=nrow(x),model_p=length(coef(fit$models[[1L]])),n_targets=nrow(targets$meta))
}
m08_prediction <- function(bundle,state,key){
  if(!state%in%names(bundle$states))stop("Unknown state: ",state)
  m01_target_prediction(bundle$engine,bundle$states[[state]],bundle$targets,key)
}

# Verify the full prediction dependency on every downstream read.
m08_signature <- function(cfg,m){
 p07<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
 paste(m07_md5(p07),m07_md5(file.path(cfg$code_dir,'Function_08_Counterfactual_risk_M01.R')),m07_md5(file.path(cfg$code_dir,'Function_07_Guideline_strategies_M01.R')),cfg$grid_days,paste(cfg$times,collapse=','),cfg$prediction_chunk,sep='|')
}
m08_bundle <- function(cfg,m){
 p<-file.path(m07_dir(cfg,'08'),sprintf('bundle_m%02d.rds',m));if(!file.exists(p))stop('Missing prediction bundle; run prediction stage')
 b<-readRDS(p);if(!identical(b$signature,m08_signature(cfg,m)))stop('Stale prediction bundle: rerun prediction stage for MI ',m);b
}



