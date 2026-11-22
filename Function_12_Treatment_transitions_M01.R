# MODULE 12 | six-category transitions, group risks, contribution decomposition
source(file.path(Sys.getenv("M01_CODE_DIR","."),"Function_11_RM_effects_M01.R"),local=FALSE)

m12_one <- function(m,cfg) {
  x<-m07_read(cfg,m);bundle<-m08_bundle(cfg,m);nall<-nrow(x)
  scopes<-m07_scope_spec(x);descs<-list();effects<-list();contrib<-list();obs<-list();qa<-list()
  i<-j<-kout<-kobs<-kqa<-1L
  for(k in m07_rms()) for(scope in names(scopes)){
    m07_log("12","MI",m," | scenario=",k," scope=",scope," (6 groups)")
    scope_ix<-which(scopes[[scope]]);n<-length(scope_ix)
    names(scope_ix)<-NULL
    if(n==0L)stop("Empty scope")
    if_sums<-setNames(lapply(cfg$times,function(t)numeric(nall)),as.character(cfg$times))
    for(group in m07_transitions()){
      ii<-scope_ix[x[[paste0("transition_",k)]][scope_ix]==group];ng<-length(ii)
      row<-data.table(imputation_id=m,scenario=k,scope=scope,transition=group,n=ng,n_scope=n,proportion=ng/n,
        mean_baseline_ldl=if(ng)mean(x$ldl_cholesterol[ii]) else NA_real_)
      for(year in c("2019","2025")){
        D<-x[[paste0("D_",year,"_",k)]][ii]
        row[[paste0("T_",year)]]<-sum(D)
        row[[paste0("mean_D_",year)]]<-if(ng)mean(D)else NA_real_
      }
      descs[[i]]<-row;i<-i+1L
            if(!ng){
        for(t in cfg$times){
          r<-m09_aj_row(NA_real_,rep(NA_real_,nall),m,t,k,scope,"observed_transition_risk");r[,transition:=group];obs[[kobs]]<-r;kobs<-kobs+1L
          for(metric in c("risk_natural","risk_2019","risk_2025","rd_2019_vs_natural","rd_2025_vs_natural","rd_2025_vs_2019")){
            r<-m09_aj_row(NA_real_,rep(NA_real_,nall),m,t,k,scope,metric);r[,transition:=group];effects[[j]]<-r;j<-j+1L
          }
          r<-m09_aj_row(0,numeric(nall),m,t,k,scope,"absolute_RD_contribution");r[,transition:=group];contrib[[kout]]<-r;kout<-kout+1L
        };next
      }
      key<-paste0("transition__",k,"__",scope,"__",group)
      pd<-m09_pred(bundle,scope,k,transitions=group)
      if(is.null(pd))stop("08 transition target missing ",key)
      eff<-m09_contrasts(pd,key,k,m,bundle$engine$complete_df,extra=FALSE)
      scope_value<-scope;eff[,`:=`(scope=scope_value,transition=group)]
      effects[[j]]<-eff;j<-j+1L
      for(t in cfg$times){
        o<-m09_aj_embed(x,ii,t)
        obs[[kobs]]<-m09_aj_row(o$Q,o$iid,m,t,k,scope,"observed_transition_risk")
        obs[[kobs]][,transition:=group];kobs<-kobs+1L
      }
      for(h in seq_along(cfg$times)){
        # Contribution of group: P(group|scope)*RD(group), including uncertainty of group size.
        g25<-pd[[paste0("2025_",k)]];g19<-pd[[paste0("2019_",k)]]
        rd<-g25$estimate[h]-g19$estimate[h]
        ifrd<-g25$iid[,h]-g19$iid[,h]
        p<-ng/n; ip<-numeric(nall);ip[scope_ix]<-(as.integer(scope_ix%in%ii)-p)/n
        ic<-p*ifrd+rd*ip
        tm<-as.character(cfg$times[h]);if_sums[[tm]]<-if(is.null(if_sums[[tm]]))ic else if_sums[[tm]]+ic
        contrib[[kout]]<-data.table(imputation_id=m,scenario=k,scope=scope,transition=group,time=cfg$times[h],
          subgroup_var="",reference="",comparison="",estimand="absolute_RD_contribution",pool_scale="natural",
          Q=p*rd,U=sum(ic^2),dfcom=bundle$engine$complete_df,n_standardized=n)
        kout<-kout+1L
      }
    }
    # All six categories reconstruct overall/subgroup scope RD at both Q and full-n IF level.
    for(h in seq_along(cfg$times)){
      subrows<-contrib[vapply(contrib,function(z)!is.null(z)&&z$scenario[1]==k&&z$scope[1]==scope&&z$time[1]==cfg$times[h],logical(1))]
      total<-sum(vapply(subrows,function(z)z$Q,numeric(1)))
      whole<-m09_pred(bundle,scope,k)
      rdwhole<-whole[[paste0("2025_",k)]]$estimate[h]-whole[[paste0("2019_",k)]]$estimate[h]
      whole_if<-whole[[paste0("2025_",k)]]$iid[,h]-whole[[paste0("2019_",k)]]$iid[,h]
      iferror<-max(abs(if_sums[[as.character(cfg$times[h])]]-whole_if))
      qa[[kqa]]<-data.table(imputation_id=m,scenario=k,scope=scope,time=cfg$times[h],
        sum_contributions=total,overall_RD=rdwhole,difference=total-rdwhole,
        max_abs_IF_reconstruction_error=iferror,
        passed=abs(total-rdwhole)<cfg$qa_tolerance && iferror<cfg$qa_tolerance)
      if(!qa[[kqa]]$passed)stop("Contribution sum != overall RD: ",k," ",scope)
      kqa<-kqa+1L
    }
  }
  list(description=rbindlist(descs,fill=TRUE),effects=rbindlist(effects,fill=TRUE),
       contributions=rbindlist(contrib,fill=TRUE),observed=rbindlist(obs,fill=TRUE),qa=rbindlist(qa,fill=TRUE))
}


m12_task <- function(m,cfg)m07_cached(cfg,'12',m,c(file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m)),file.path(m07_dir(cfg,'08'),sprintf('bundle_m%02d.rds',m))),function(m)m12_one(m,cfg))
m12_export <- function(z,cfg){
 m09_emit(z,cfg,'treatment_transitions',c(effects='transition_risks_and_effects',contributions='risk_difference_contributions',observed='observed_transition_risks'))
 d<-rbindlist(lapply(z,`[[`,'description'));m12_validate_transition_counts(d)
  keys<-c('scenario','scope','transition');vals<-setdiff(names(d),c(keys,'imputation_id'))
 m07_export(cfg,'treatment_transitions','transition_counts_and_required_reductions',d,m10_descriptive_pool(d,keys,vals))
 ref<-d[scenario=='S0',.(imputation_id,scope,transition,reference_proportion=proportion)]
 changes<-merge(d,ref,by=c('imputation_id','scope','transition'));changes[,change_proportion_vs_S0:=proportion-reference_proportion]
 m07_export(cfg,'treatment_transitions','transition_proportion_changes',changes,m10_descriptive_pool(changes,keys,'change_proportion_vs_S0'))
 c<-m09_pool(rbindlist(lapply(z,`[[`,'contributions')),cfg)
 # Use contribution totals, not the unweighted sum of group-specific RDs.
 overall<-c[,.(total_RD=sum(estimate)),by=.(scenario,scope,time)]
 c<-merge(c,overall,by=c('scenario','scope','time'));c[,contribution_percent:=ifelse(abs(total_RD)>1e-8,100*estimate/total_RD,NA_real_)];c[,percentage_status:=ifelse(abs(total_RD)>1e-8,'descriptive','unstable_near_zero_denominator')]
 m07_result(cfg,'treatment_transitions','contribution_percentages',c)
 qa<-rbindlist(lapply(z,`[[`,'qa'));if(!all(qa$passed))stop('Contribution QA failed');m07_result(cfg,'analysis_quality_control','contribution_additivity',qa)
 files<-list.files(cfg$out_dir,pattern='[.]csv$',recursive=TRUE,full.names=TRUE)
 m07_result(cfg,'analysis_quality_control','output_manifest',data.table(result_path=files,bytes=file.info(files)$size,md5=vapply(files,m07_md5,character(1))))
 m12_mean_reduction_exports(cfg)
}

# Descriptive postprocessing only; no MI, model refit or CI arithmetic in reports.
m12_mean_reduction_exports <- function(cfg){
 p<-file.path(cfg$out_dir,'treatment_transitions/per_imputation/transition_counts_and_required_reductions.csv')
 d<-fread(p);a<-melt(d,id.vars=c('imputation_id','scope','scenario','transition'),measure.vars=c('mean_D_2019','mean_D_2025'),variable.name='guideline',value.name='mean_D')
 a[,guideline:=sub('mean_D_','',guideline,fixed=TRUE)]
 dif<-merge(a[guideline=='2019'],a[guideline=='2025'],by=c('imputation_id','scope','scenario','transition'),suffixes=c('_19','_25'))
 dif<-dif[,.(imputation_id,scope,scenario,transition,guideline='2025-2019',mean_D=mean_D_25-mean_D_19)]
 a<-rbindlist(list(a,dif),use.names=TRUE)
 ref<-a[scenario=='S0',.(imputation_id,scope,transition,guideline,reference_mean_D=mean_D)]
 a<-merge(a,ref,by=c('imputation_id','scope','transition','guideline'));a[,change_vs_S0:=mean_D-reference_mean_D]
 keys<-c('scope','scenario','transition','guideline')
 # Do not silently discard an MI with an empty transition category.
 pool<-function(q,keys)q[,.(mean_D=mean(mean_D),change_vs_S0=mean(change_vs_S0),m=.N,m_valid=sum(is.finite(mean_D))),by=keys]
 m07_export(cfg,'treatment_transitions','transition_mean_required_reduction_contrasts',a,pool(a,keys))
 q<-fread(file.path(cfg$out_dir,'intervention_requirements/per_imputation/classification_and_requirements.csv'))
 q<-q[metric=='mean_D',.(imputation_id,scope,scenario,guideline,mean_D=value)]
 dd<-merge(q[guideline=='2019'],q[guideline=='2025'],by=c('imputation_id','scope','scenario'),suffixes=c('_19','_25'))
 q<-rbindlist(list(q,dd[,.(imputation_id,scope,scenario,guideline='2025-2019',mean_D=mean_D_25-mean_D_19)]))
 ref<-q[scenario=='S0',.(imputation_id,scope,guideline,reference_mean_D=mean_D)]
 q<-merge(q,ref,by=c('imputation_id','scope','guideline'));q[,change_vs_S0:=mean_D-reference_mean_D]
 m07_export(cfg,'intervention_requirements','population_mean_required_reduction_contrasts',q,pool(q,c('scope','scenario','guideline')))
}



# Counts-only integrity QA: preserves treatment categories and all calculations.
m12_validate_transition_counts <- function(d){
  keys<-c('imputation_id','scenario','scope','transition')
  if(anyDuplicated(d,by=keys))stop('12: duplicate transition per MI/scope/scenario')
  q<-d[,.(cell_count=.N,denoms=uniqueN(n_scope),total_n=sum(n),
    group_n=n_scope[1L],total_proportion=sum(proportion),
    min_n=min(n)),by=.(imputation_id,scenario,scope)]
  if(any(q$cell_count!=length(m07_transitions())|q$denoms!=1L|
    abs(q$total_n-q$group_n)>1e-6|abs(q$total_proportion-1)>1e-9|
    q$min_n<0))stop('12: six treatment categories do not exhaust each scope')
  invisible(TRUE)
}
