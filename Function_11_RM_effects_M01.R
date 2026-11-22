# MODULE 11 | RM strategy effects and heterogeneity tests
source(file.path(Sys.getenv("M01_CODE_DIR","."),"Function_10_RM_classification_M01.R"),local=FALSE)


m11_task <- function(m,cfg){
 m07_cached(cfg,'11',m,c(file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m)),file.path(m07_dir(cfg,'08'),sprintf('bundle_m%02d.rds',m))),function(m)m11_one(m,cfg))
}
m11_one <- function(m,cfg){
 b<-m08_bundle(cfg,m);df<-b$engine$complete_df;scopes<-names(m07_scope_spec(m07_read(cfg,m)))
 effects<-hets<-inc<-ih<-list();spec<-m07_subgroups()
 comparisons<-list(SH_vs_S0=c('S0','SH'),ST_vs_S0=c('S0','ST'),SHT_vs_S0=c('S0','SHT'),SHT_vs_SH=c('SH','SHT'))
 # Process one comparison at a time, keeping only the needed subgroup IF contrasts.
 for(scope in scopes){
   m07_log('11','MI',m,' scope ',scope)
   pd<-setNames(lapply(m07_states(),function(s)m08_prediction(b,s,paste0('scope__',scope))),m07_states())
   for(k in m07_rms())effects[[length(effects)+1L]]<-m09_contrasts(pd,scope,k,m,df)
   for(cmp in names(comparisons)){
     a<-comparisons[[cmp]][1];z<-comparisons[[cmp]][2]
     weights<-list(additional_reduction_2019=setNames(c(1,-1),paste0('2019_',c(a,z))),additional_reduction_2025=setNames(c(1,-1),paste0('2025_',c(a,z))),additional_reduction_difference_2025_vs_2019=setNames(c(1,-1,-1,1),c(paste0('2025_',c(a,z)),paste0('2019_',c(a,z)))))
     for(nm in names(weights)){r<-m09_one_effect(pd,weights[[nm]],nm,scope,z,m,df);r[,rm_comparison:=cmp];inc[[length(inc)+1L]]<-r}
   }
 }
 for(k in m07_rms())hets[[k]]<-m09_het(b,k,m,df)
 for(cmp in names(comparisons))for(j in seq_len(nrow(spec))){
   a<-comparisons[[cmp]][1];z<-comparisons[[cmp]][2];v<-spec[j]
   calc<-function(scope){pd<-setNames(lapply(m07_states(),function(s)m08_prediction(b,s,paste0('scope__',scope))),m07_states());w<-setNames(c(1,-1,-1,1),c(paste0('2025_',c(a,z)),paste0('2019_',c(a,z))));list(Q=Reduce(`+`,Map(function(s,w)pd[[s]]$estimate*w,names(w),w)),iid=Reduce(`+`,Map(function(s,w)pd[[s]]$iid*w,names(w),w)))}
   c1<-calc(v$comparison);r1<-calc(v$reference)
   ih[[length(ih)+1L]]<-data.table(imputation_id=m,time=cfg$times,scenario=z,rm_comparison=cmp,scope=paste0('heterogeneity_',v$subgroup_var),subgroup_var=v$subgroup_var,reference=v$reference,comparison=v$comparison,estimand='incremental_benefit_heterogeneity',pool_scale='natural',Q=c1$Q-r1$Q,U=colSums((c1$iid-r1$iid)^2),dfcom=df)
 }
 list(effects=rbindlist(effects,fill=TRUE),heterogeneity=rbindlist(hets,fill=TRUE),incremental=rbindlist(inc,fill=TRUE),incremental_heterogeneity=rbindlist(ih,fill=TRUE))
}

# Export gate only: does not recalculate effect estimates or inferential quantities.
m11_validate_outputs <- function(z,cfg){
  if(length(z)!=cfg$m)stop('11: completed MI objects differ from cfg$m')
  for(j in seq_along(z)){
    for(k in c('effects','heterogeneity','incremental','incremental_heterogeneity')){
      d<-z[[j]][[k]]
      if(is.null(d)||!nrow(d))stop('11: missing ',k,' for MI ',j)
      if('imputation_id'%in%names(d)&&any(d$imputation_id!=j))
        stop('11: MI indexing mismatch: ',k,' / ',j)
    }
    if(!setequal(unique(z[[j]]$effects$scenario),m07_rms()))
      stop('11: incomplete RM scenarios for MI ',j)
  }
  invisible(TRUE)
}
