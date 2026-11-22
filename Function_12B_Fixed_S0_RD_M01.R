# MODULE 12B | Fixed S0 six-group 2025-2019 ASCVD RD and RM-minus-S0 changes.
# Reuses frozen stage 07 allocations, stage 08 Cox fits and Cox/IF engine.
# NO Cox refit, NO new MI, NO replacement of stage 08 bundle/checkpoints.
# Independently rebuilds fixed-S0 prediction targets to preserve paired full-n IF.
source(file.path(Sys.getenv('M01_CODE_DIR','.'),
                 'Function_09_Main_analysis_M01.R'),local=FALSE)
source(file.path(Sys.getenv('M01_CODE_DIR','.'),
                 'Function_12A_Transition_matrix_M01.R'),local=FALSE)
if(!requireNamespace('Matrix',quietly=TRUE))stop('Install Matrix')
suppressPackageStartupMessages(library(data.table))

m12b_version <- function() '12B_fixed_S0_paired_RD_v1.1_training_data_QA'
m12b_board <- function(cfg){
  p<-file.path(cfg$out_dir,'fixed_S0_treatment_group_RD')
  dir.create(p,recursive=TRUE,showWarnings=FALSE);p
}
m12b_signature <- function(cfg,m){
  p07<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
  p08<-file.path(m07_dir(cfg,'08'),sprintf('bundle_m%02d.rds',m))
  model<-file.path(m07_dir(cfg,'08'),sprintf('model_m%02d.rds',m))
  code<-file.path(cfg$code_dir,c('Function_07_Guideline_strategies_M01.R',
    'Function_08_Counterfactual_risk_M01.R','Function_09_Main_analysis_M01.R',
    'Function_12A_Transition_matrix_M01.R','Function_12B_Fixed_S0_RD_M01.R'))
  dep<-c(p07,p08,model,file.path(cfg$out_dir,'study_config.rds'),code)
  if(any(!file.exists(dep)))stop('12B missing dependencies: ',paste(dep[!file.exists(dep)],collapse=', '))
  list(version=m12b_version(),
       dependency_md5=setNames(vapply(dep,m07_md5,character(1L)),basename(dep)),
       reported_times=cfg$times,grid_days=cfg$grid_days,
       chunk_size=cfg$prediction_chunk)
}
m12b_validate_input <- function(x,m){
  x<-as.data.table(x);m07_unique(x)
  m07_assert(x,c('participant_id','transition_S0','ldl_cholesterol',
    as.vector(outer(c('2019','2025'),m07_rms(),function(y,k)paste0('Lcf_',y,'_',k)))),
    paste0('12B MI ',m))
  s<-as.character(x$transition_S0)
  if(anyNA(s)||any(!s%in%m12a_levels()))stop('12B incomplete S0 fixed group')
  if('imputation_id'%in%names(x) && any(x$imputation_id!=m))
    stop('12B MI indexing mismatch')
  invisible(TRUE)
}
m12b_fixed_targets <- function(x){
  fixed<-as.character(x$transition_S0)
  keys<-c('scope__overall',paste0('scope__fixed_S0__',m12a_levels()))
  ix<-c(list(seq_len(nrow(x))),lapply(m12a_levels(),function(g)which(fixed==g)))
  names(ix)<-keys
  sizes<-lengths(ix)
  # m01_strategy_state requires all populated targets; keep empty groups
  # in the metadata for output but exclude from numerical prediction.
  available<-which(sizes>0L)
  idx<-lapply(ix[available],as.integer)
  meta<-data.table(target_key=names(idx),target_type='scope',
                   scope=names(idx),scenario='shared',
                   transition_id=NA_character_,n_standardized=lengths(idx))
  list(targets=list(index=idx,meta=meta),all_keys=keys,sizes=sizes)
}
m12b_load_engine <- function(x, bundle, cfg, m) {
  # Stage 08 first prepares *every* model variable, including the two
  # cause-specific outcomes. Compare with precisely that training input.
  p <- file.path(m07_dir(cfg, '08'), sprintf('model_m%02d.rds', m))
  fit <- readRDS(p)
  if (!is.list(fit$models) || length(fit$models) != 2L)
    stop('12B invalid original stage-08 Cox fit, MI ', m)

  model_data <- m08_prepare_model(x)
  if (nrow(model_data) != nrow(x) ||
      !identical(as.character(model_data$participant_id),
                 as.character(x$participant_id)))
    stop('12B stage-08 preprocessing unexpectedly changed participant order, MI ', m)

  engine <- bundle$engine
  if (!identical(as.numeric(engine$times), as.numeric(cfg$times)) ||
      !identical(as.integer(engine$n), as.integer(nrow(model_data))) ||
      !identical(as.integer(engine$grid_days), as.integer(cfg$grid_days)))
    stop('12B original stage-08 prediction engine configuration mismatch, MI ', m)
  if (!identical(engine$participant_id_md5,
                 m07_hash_object(as.character(model_data$participant_id))))
    stop('12B stage-08 prediction participant ID/order mismatch, MI ', m)

  for (j in seq_len(2L)) {
    mod <- fit$models[[j]]
    original <- engine[[paste0('source', j)]]
    design <- as.matrix(mod$x)
    beta <- stats::coef(mod)
    expected_time <- as.numeric(model_data$followup_years)
    expected_status <- if (j == 1L) {
      as.integer(model_data$ascvd_event)
    } else {
      as.integer(model_data$competing_death_event)
    }
    if (is.null(mod$y) || nrow(mod$y) != nrow(model_data) ||
        nrow(design) != nrow(model_data))
      stop('12B saved Cox fit has unexpected training size: MI ', m,
           ', cause ', j)
    observed_time <- as.numeric(mod$y[, 1L])
    observed_status <- as.integer(mod$y[, 2L])

    time_mismatch <- which(is.na(observed_time) != is.na(expected_time) |
                           (!is.na(observed_time) & !is.na(expected_time) &
                            observed_time != expected_time))
    status_mismatch <- which(is.na(observed_status) != is.na(expected_status) |
                             (!is.na(observed_status) & !is.na(expected_status) &
                              observed_status != expected_status))
    if (length(time_mismatch) || length(status_mismatch)) {
      first_time <- if (length(time_mismatch)) time_mismatch[1L] else NA_integer_
      first_status <- if (length(status_mismatch)) status_mismatch[1L] else NA_integer_
      stop('12B stage-08 Cox training data mismatch after m08_prepare_model: MI ',
           m, ', cause ', j,
           '; time mismatches=', length(time_mismatch),
           ' (first row=', first_time, ')',
           '; event mismatches=', length(status_mismatch),
           ' (first row=', first_status, '). ',
           'Check that 07 MI, 08 model and 08 bundle belong to the same run.')
    }
    if (is.null(design) || ncol(design) != length(beta) ||
        !identical(colnames(design), names(beta)) ||
        !identical(names(beta), names(original$beta)))
      stop('12B saved Cox design columns differ from stage-08 engine: MI ', m,
           ', cause ', j)
    if (anyNA(design) || any(!is.finite(design)) ||
        anyNA(beta) || any(!is.finite(beta)) ||
        length(mod$means) != length(original$means) ||
        max(abs(as.numeric(beta) - as.numeric(original$beta))) > 1e-10 ||
        max(abs(as.numeric(mod$means) - as.numeric(original$means))) > 1e-10 ||
        max(abs(exp(mod$linear.predictors) - original$risk)) > 1e-8)
      stop('12B saved Cox model and frozen prediction engine disagree: MI ', m,
           ', cause ', j)
    engine[[paste0('source', j)]]$x <- design
  }
  rm(model_data)
  m07_log('12B', 'MI ', m,
          ' frozen 08 Cox training time/status/design validated (no model refit)')
  engine
}

m12b_state_fixed <- function(state,engine,x,bundle,fx,cfg,m) {
  col<-paste0('ldl_',state)
  lcf<-paste0('Lcf_',state)
  if(!lcf%in%names(x))stop('12B missing counterfactual LDL: ',lcf)
  if(col%in%names(x))stop('12B counterfactual temporary column already exists: ',col)
  data.table::set(x,j=col,value=as.numeric(x[[lcf]]))
  on.exit(data.table::set(x,j=col,value=NULL),add=TRUE)
  m07_log('12B','MI ',m,' START prediction ',state,
          ' | fixed targets=',nrow(fx$targets$meta))
  st<-m07_timed('12B',paste0('MI ',m,' fixed target CIF+IF ',state),
                m01_strategy_state(engine,x,fx$targets,state,cfg))
  keys<-fx$targets$meta$target_key
  ans<-setNames(lapply(keys,function(key){
    m01_target_prediction(engine,st,fx$targets,key)
  }),keys)
  # A crucial audit: same overall target must reproduce the saved stage-08
  # point estimate and COMPLETE influence function, not just an SE.
  ref<-m08_prediction(bundle,state,'scope__overall')
  observed<-ans[['scope__overall']]
  diff_q<-max(abs(observed$estimate-ref$estimate))
  diff_if<-max(abs(observed$iid-ref$iid))
  if(!is.finite(diff_q)||!is.finite(diff_if)||
     diff_q>cfg$qa_tolerance||diff_if>cfg$qa_tolerance)
    stop('12B overall prediction/IF failed stage 08 comparison: ',state,
         ' | point diff=',diff_q,' IF diff=',diff_if)
  m07_log('12B','MI ',m,' ',state,' validated against stage 08; max IF error=',
          format(diff_if,digits=3))
  list(predictions=ans,
       qa=data.table(imputation_id=m,state=state,
                     max_abs_Q_error=diff_q,max_abs_IF_error=diff_if,
                     passed=TRUE))
}
m12b_scenario_checkpoint <- function(m,public,x,bundle,engine,fx,cfg,signature){
  internal<-if(public=='S0') 'S0' else unname(m12a_scenarios()[[public]])
  cp<-file.path(m12b_board(cfg),'checkpoints',
                sprintf('mi%02d_%s_RD.rds',m,public))
  sp<-paste0(cp,'.signature.rds')
  sign<-list(source=signature,scenario=public,internal=internal,
             groups=m12a_levels())
  if(file.exists(cp)&&file.exists(sp)&&identical(readRDS(sp),sign)){
    z<-readRDS(cp)
    if(!isTRUE(all(z$qa$passed)))stop('12B corrupted checkpoint QA ',public)
    m07_log('12B','MI ',m,' ',public,' verified scenario checkpoint reused')
    return(z)
  }
  years<-c('2019','2025'); result<-vector('list',2L)
  for(j in seq_along(years)){
    state<-paste0(years[j],'_',internal)
    result[[j]]<-m12b_state_fixed(state,engine,x,bundle,fx,cfg,m)
  }
  keys<-fx$all_keys
  rds<-vector('list',length(keys));names(rds)<-keys
  for(key in keys){
    a<-result[[1L]]$predictions[[key]]
    b<-result[[2L]]$predictions[[key]]
    if(is.null(a) && is.null(b)){
      rds[[key]]<-list(Q=rep(NA_real_,length(cfg$times)),
                       iid=NULL,n=0L)
    } else if(is.null(a)||is.null(b)){
      stop('12B one-sided fixed prediction ',public,'/',key)
    } else {
      if(!identical(dim(a$iid),dim(b$iid)) ||
         length(a$estimate)!=length(cfg$times))stop('12B unpaired prediction')
      rds[[key]]<-list(Q=b$estimate-a$estimate,
                       iid=b$iid-a$iid,n=a$n_standardized)
    }
  }
  z<-list(imputation_id=m,scenario=public,internal=internal,
          rds=rds,qa=rbindlist(lapply(result,`[[`,'qa')))
  if(!all(z$qa$passed))stop('12B scenario model QA failed')
  m07_atomic(z,cp);m07_atomic(sign,sp)
  m07_log('12B','MI ',m,' ',public,' checkpoint committed')
  z
}
m12b_metric <- function(m,public,internal,key,metric,Q,iid,cfg,n_group,n_total,dfcom){
  if(is.null(iid)){
    u<-rep(NA_real_,length(cfg$times))
  }else{
    if(!identical(dim(iid),c(as.integer(n_total),as.integer(length(cfg$times)))) ||
       any(!is.finite(iid)))stop('12B invalid influence-function shape')
    if(max(abs(colSums(iid)))>cfg$qa_tolerance)
      stop('12B IF not centred: ',metric,' ',public,' ',key)
    u<-colSums(iid^2)
  }
  code<-if(key=='scope__overall')'Overall'else{
    paste0('G',match(sub('^scope__fixed_S0__','',key),m12a_levels()))
  }
  label<-if(code=='Overall')'Overall' else m12a_labels()$label[match(code,m12a_labels()$code)]
  data.table(imputation_id=m,scenario=public,scenario_internal=internal,
    scope='overall',fixed_group=code,fixed_order=if(code=='Overall')7L else as.integer(sub('G','',code)),
    fixed_label=label,estimand=metric,pool_scale='natural',
    time=cfg$times,Q=as.numeric(Q),U=as.numeric(u),dfcom=dfcom,
    n_standardized=n_group,n_total=n_total,
    group_weight=n_group/n_total,
    cell_status=if(n_group==0L)'empty_fixed_group'else 'estimable')
}
m12b_one <- function(m,cfg,signature){
  x<-as.data.table(m07_read(cfg,m))
  m12b_validate_input(x,m)
  bundle<-m08_bundle(cfg,m)
  fx<-m12b_fixed_targets(x)
  engine<-NULL
  identity<-m07_hash_object(as.character(x$participant_id))
  fixed_hash<-m07_hash_object(as.character(x$transition_S0))
  count<-fx$sizes
  result<-list(); qa<-list(); model_audit<-list()
  ref<-NULL
  groups<-fx$all_keys
  n<-nrow(x)
  for(public in c('S0',names(m12a_scenarios()))){
    scp<-file.path(m12b_board(cfg),'checkpoints',
                   sprintf('mi%02d_%s_RD.rds',m,public))
    sigp<-paste0(scp,'.signature.rds')
    wants<-list(source=signature,scenario=public,
                internal=if(public=='S0')'S0'else unname(m12a_scenarios()[[public]]),
                groups=m12a_levels())
    if(file.exists(scp)&&file.exists(sigp)&&identical(readRDS(sigp),wants)){
      sc<-readRDS(scp)
      if(!all(sc$qa$passed))stop('12B stale invalid model QA')
      m07_log('12B','MI ',m,' ',public,' scenario checkpoint reused')
    }else{
      if(is.null(engine)) engine<-m12b_load_engine(x,bundle,cfg,m)
      # Use the same workhorse for S0 and all RM strategies.
      sc<-m12b_scenario_checkpoint(m,public,x,bundle,engine,fx,cfg,signature)
    }
    model_audit[[length(model_audit)+1L]]<-sc$qa
    if(public=='S0')ref<-sc
    for(key in groups){
      now<-sc$rds[[key]];base<-ref$rds[[key]]
      ng<-count[[key]];p<-ng/n
      if(ng==0L){
        # Empty fixed group: conditional risk is not estimable, but
        # its unconditional whole-cohort contribution is identically zero.
        rd<-delta<-rep(NA_real_,length(cfg$times))
        if_rd<-if_delta<-NULL
        weighted<-weighted_delta<-rep(0,length(cfg$times))
        if_weighted<-if_weighted_delta<-matrix(0,n,length(cfg$times))
      }else{
        rd<-now$Q;delta<-now$Q-base$Q
        if_rd<-now$iid; if_delta<-now$iid-base$iid
        ip<-numeric(n)
        if(key!='scope__overall'){
          member<-which(as.character(x$transition_S0)==
                        sub('^scope__fixed_S0__','',key))
          ip<-(as.numeric(seq_len(n)%in%member)-p)/n
        }
        weighted<-p*rd;weighted_delta<-p*delta
        if_weighted<-p*if_rd+tcrossprod(ip,rd)
        if_weighted_delta<-p*if_delta+tcrossprod(ip,delta)
      }
      rows<-list(
        m12b_metric(m,public,sc$internal,key,'fixed_rd_2025_vs_2019',rd,if_rd,cfg,ng,n,bundle$engine$complete_df),
        m12b_metric(m,public,sc$internal,key,'fixed_delta_rd_vs_S0',delta,if_delta,cfg,ng,n,bundle$engine$complete_df),
        m12b_metric(m,public,sc$internal,key,'weighted_rd_contribution',weighted,if_weighted,cfg,ng,n,bundle$engine$complete_df),
        m12b_metric(m,public,sc$internal,key,'weighted_delta_rd_contribution',weighted_delta,if_weighted_delta,cfg,ng,n,bundle$engine$complete_df))
      result[[length(result)+1L]]<-rbindlist(rows)
    }
    # Exact pooled-risk decomposition is checked within each MI on the full
    # n-vector IF, separately for each time and RM scenario.
    all_rd<-sc$rds[['scope__overall']]
    overall0<-ref$rds[['scope__overall']]
    whole_delta<-all_rd$Q-overall0$Q
    whole_if_delta<-all_rd$iid-overall0$iid
    weighted_sum<-rep(0,length(cfg$times))
    delta_sum<-rep(0,length(cfg$times))
    if_sum<-matrix(0,n,length(cfg$times))
    if_delta_sum<-matrix(0,n,length(cfg$times))
    for(key in groups[-1L]){
      ng<-count[[key]]
      if(ng==0L)next
      p<-ng/n
      g<-sub('^scope__fixed_S0__','',key)
      membership<-as.numeric(as.character(x$transition_S0)==g)
      ip<-(membership-p)/n
      current<-sc$rds[[key]];baseline<-ref$rds[[key]]
      weighted_sum<-weighted_sum+p*current$Q
      delta_sum<-delta_sum+p*(current$Q-baseline$Q)
      if_sum<-if_sum+p*current$iid+tcrossprod(ip,current$Q)
      if_delta_sum<-if_delta_sum+p*(current$iid-baseline$iid)+
        tcrossprod(ip,current$Q-baseline$Q)
    }
    q_err<-max(abs(weighted_sum-all_rd$Q))
    delta_err<-max(abs(delta_sum-whole_delta))
    if_err<-max(abs(if_sum-all_rd$iid))
    delta_if_err<-max(abs(if_delta_sum-whole_if_delta))
    pass<-all(is.finite(c(q_err,delta_err,if_err,delta_if_err))) &&
      all(c(q_err,delta_err,if_err,delta_if_err)<cfg$qa_tolerance)
    qa[[length(qa)+1L]]<-data.table(imputation_id=m,scenario=public,
      max_abs_RD_additivity_error=q_err,
      max_abs_delta_RD_additivity_error=delta_err,
      max_abs_IF_additivity_error=if_err,
      max_abs_delta_IF_additivity_error=delta_if_err,
      n=n,nonempty_groups=sum(count[-1L]>0),
      passed=pass)
    if(!pass)stop('12B fixed-group contribution reconstruction failed MI ',m,
                  ' scenario ',public,' Q=',q_err,' IF=',if_err)
    m07_log('12B','MI ',m,' ',public,' RD/IF reconstruction PASSED | ',
            'max IF err=',format(if_err,digits=3))
    if(public!='S0')rm(sc)
    gc(FALSE)
  }
  list(imputation_id=m,rows=rbindlist(result),qa=rbindlist(qa),
       model_audit=rbindlist(model_audit),
       id_signature=identity,fixed_signature=fixed_hash,
       sizes=data.table(imputation_id=m,fixed_group=c(paste0('G',1:6),'Overall'),
        fixed_label=c(m12a_labels()$label,'Overall'),
        n=c(unname(count[-1L]),n),n_total=n))
}
m12b_task <- function(m,cfg){
  root<-m12b_board(cfg)
  cp<-file.path(root,'checkpoints',sprintf('mi%02d_final.rds',m))
  sigp<-paste0(cp,'.signature.rds')
  sig<-m12b_signature(cfg,m)
  if(file.exists(cp)&&file.exists(sigp)&&identical(readRDS(sigp),sig)){
    z<-readRDS(cp)
    if(!all(z$qa$passed)||!all(z$model_audit$passed))stop('12B final cache QA failed')
    m07_log('12B','MI ',m,' verified final cache reused')
    return(z)
  }
  z<-m07_timed('12B',sprintf('MI %02d fixed S0 groups all RM strategies',m),
               m12b_one(m,cfg,sig))
  m07_atomic(z,cp);m07_atomic(sig,sigp)
  z
}
m12b_plot <- function(pooled,metric,cfg){
  z<-copy(pooled[estimand==metric])
  z[,`:=`(scenario_order=match(scenario,c('S0','SN','ST','SNT')),
           rd_pp=100*estimate,lcl_pp=100*lcl,ucl_pp=100*ucl,
           se_pp=100*se,zero_line=0,
           scenario_label=ifelse(scenario=='S0','No RM (S0)',
             ifelse(scenario=='SN','Non-laboratory RM (SN)',
             ifelse(scenario=='ST','Laboratory RM (ST)','Combined RM (SNT)'))),
           x_axis_label='Risk difference (percentage points)',
           risk_difference_definition=if(grepl('delta',metric,fixed=TRUE))
             'RD with RM minus RD without RM' else '2025 minus 2019',
           followup_years=time)]
  setorderv(z,c('fixed_order','scenario_order','time'))
  z
}
m12b_export <- function(all,cfg){
  if(length(all)!=cfg$m)stop('12B MI count mismatch')
  if(length(unique(vapply(all,`[[`,character(1),'id_signature')))!=1L ||
     length(unique(vapply(all,`[[`,character(1),'fixed_signature')))!=1L)
    stop('12B participant IDs or fixed S0 group membership differ across imputations')
  if(!all(vapply(all,function(z)all(z$qa$passed)&&all(z$model_audit$passed),logical(1))))
    stop('12B model / decomposition QA failed')
  root<-m12b_board(cfg)
  d<-rbindlist(lapply(all,`[[`,'rows'))
  if(anyDuplicated(d,by=c('imputation_id','scenario','fixed_group','time','estimand')))
    stop('12B duplicated MI/target/estimand')
  expected<-cfg$m*4L*7L*length(cfg$times)*4L
  if(nrow(d)!=expected)stop('12B incomplete outputs: expected ',expected,' rows, got ',nrow(d))
  pool<-m09_pool(d,cfg,keys=c('scenario','scenario_internal','scope',
     'fixed_group','fixed_order','fixed_label','estimand','pool_scale','time'))
  sizes<-rbindlist(lapply(all,`[[`,'sizes'))
  size_pool<-sizes[,.(mean_n=mean(n),mean_n_total=mean(n_total),m=.N),
                   by=.(fixed_group,fixed_label)]
  pool<-merge(pool,size_pool,by=c('fixed_group','fixed_label'),all.x=TRUE)
  pool[mean_n==0,`:=`(status='empty_fixed_group',
    estimate=NA_real_,se=NA_real_,lcl=NA_real_,ucl=NA_real_)]
  # By definition weighted empty-group contributions are zero and remain estimable.
  pool[mean_n==0 & estimand%in%c('weighted_rd_contribution',
             'weighted_delta_rd_contribution'),`:=`(
             status='structural_zero',estimate=0,se=0,lcl=0,ucl=0)]
  setorderv(pool,c('estimand','fixed_order','scenario','time'))
  m07_csv(d,file.path(root,'12B_estimands_per_imputation.csv'))
  m07_csv(pool,file.path(root,'12B_pooled_RD_effects.csv'))
  m07_csv(sizes,file.path(root,'12B_fixed_S0_group_sizes_per_imputation.csv'))
  m07_csv(size_pool,file.path(root,'12B_fixed_S0_group_sizes_MI_mean.csv'))
  metrics<-c(fixed_rd_2025_vs_2019='12B_plot_fixed_group_RD.csv',
             fixed_delta_rd_vs_S0='12B_plot_fixed_group_RD_change.csv',
             weighted_rd_contribution='12B_plot_weighted_RD_contribution.csv',
             weighted_delta_rd_contribution='12B_plot_weighted_RD_contribution_change.csv')
  for(metric in names(metrics)){
    dplot<-m12b_plot(pool,metric,cfg)
    if(metric%in%c('fixed_delta_rd_vs_S0','weighted_delta_rd_contribution'))
      dplot<-dplot[scenario!='S0']
    m07_csv(dplot,file.path(root,metrics[[metric]]))
  }
  m07_csv(rbindlist(lapply(all,`[[`,'qa')),
          file.path(root,'12B_RD_and_IF_additivity_QA.csv'))
  m07_csv(rbindlist(lapply(all,`[[`,'model_audit')),
          file.path(root,'12B_frozen_08_prediction_QA.csv'))
  deps<-rbindlist(lapply(all,function(z){
    s<-m12b_signature(cfg,z$imputation_id)
    data.table(imputation_id=z$imputation_id,version=s$version,
               dependency=names(s$dependency_md5),md5=unname(s$dependency_md5))
  }))
  m07_csv(deps,file.path(root,'12B_input_manifest.csv'))
  # A second lightweight check AFTER pooling on estimates and their weighted
  # contributions. Inference is already formed from paired MI-level IFs.
  verify<-pool[estimand%in%c('fixed_rd_2025_vs_2019',
                              'fixed_delta_rd_vs_S0',
                              'weighted_rd_contribution',
                              'weighted_delta_rd_contribution')]
  core<-verify[fixed_group=='Overall',.(scenario,time,estimand,overall=estimate)]
  weighted<-verify[fixed_group!='Overall' &
    estimand%in%c('weighted_rd_contribution','weighted_delta_rd_contribution'),
    .(weighted_sum=sum(estimate)),by=.(scenario,time,estimand)]
  weighted[,estimand:=ifelse(estimand=='weighted_rd_contribution',
                            'fixed_rd_2025_vs_2019','fixed_delta_rd_vs_S0')]
  check<-merge(weighted,core,by=c('scenario','time','estimand'))
  check[,`:=`(abs_error=abs(weighted_sum-overall),
              passed=abs(weighted_sum-overall)<cfg$qa_tolerance)]
  if(nrow(check)!=4L*length(cfg$times)*2L || any(!check$passed))
    stop('12B pooled RD contributions not additive')
  m07_csv(check,file.path(root,'12B_pooled_RD_additivity_QA.csv'))
  m07_log('12B','EXPORT DONE | 4 scenarios, fixed G1-G6 + Overall, '
          ,paste(cfg$times,collapse='/'),' years | ',cfg$m,' MI | output=',root)
  invisible(list(pooled=pool,qa=check))
}
