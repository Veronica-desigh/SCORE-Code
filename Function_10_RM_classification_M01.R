# MODULE 10 | SINGLE classification / observed partition computation (2026-11-04)
# One 07 RDS read, one H/T group index, one per-MI cache. The formerly separate
# m07_cached("10") and observed_partition_mXX caches are NOT read.
# All saved filenames stay compatible with 11-14; no models are refitted here.
source(file.path(Sys.getenv('M01_CODE_DIR','.'),'Function_09_Main_analysis_M01.R'),local=FALSE)

m10_mean <- function(tab,keys){
  if(!nrow(tab))return(tab)
  if('value'%in%names(tab))tab[,.(MI_mean=mean(value),
    MI_mean_proportion=if(all(is.na(proportion)))NA_real_ else mean(proportion,na.rm=TRUE),
    MI_n=uniqueN(imputation_id)),by=keys]
  else if('n'%in%names(tab))tab[,.(MI_mean_n=mean(n),MI_mean_proportion=mean(percent),
    MI_n=uniqueN(imputation_id)),by=keys]
  else stop('Cannot summarize')
}
m10_descriptive_pool <- function(d,keys,values){
  d[,lapply(.SD,function(v)if(all(is.na(v)))NA_real_ else mean(v)),by=keys,.SDcols=values]
}

# The only definition of the four original-M H/T membership sets.
m10_ht_groups <- function(x,scopes){
  result<-list()
  for(sc in names(scopes))for(year in c('2019','2025')){
    ix<-which(scopes[[sc]])
    base<-as.character(x[[if(year=='2019')'score_2019_category' else 'score2_2025_category']])
    pop<-ix[base[ix]=='M']
    H<-as.integer(as.character(m07_rm_binary(x,year,'H')))
    T<-as.integer(as.character(m07_rm_binary(x,year,'T')))
    if(anyNA(H[pop])||anyNA(T[pop])||any(!H[pop]%in%0:1)||any(!T[pop]%in%0:1))
      stop('10: invalid H/T binary status for ',year,' ',sc)
    meta<-CJ(h=0:1,t=0:1)
    meta[,`:=`(scope=sc,guideline=year,scenario='S0',cell_id=paste0('H',h,'T',t))]
    members<-lapply(seq_len(nrow(meta)),function(j)pop[H[pop]==meta$h[j]&T[pop]==meta$t[j]])
    if(anyDuplicated(unlist(members))||!setequal(unlist(members),pop))
      stop('10: H/T is not an exhaustive, disjoint original-M partition: ',year,' ',sc)
    result[[paste(sc,year,sep='__')]]<-list(population=pop,members=members,meta=meta)
  }
  result
}

# The older overlap.csv is still exported, but derived from the SAME partition
# membership and not computed independently. Maintains all downstream columns.
m10_overlap_from_counts <- function(ct){
  z<-copy(ct[guideline%in%c('2019','2025') & scenario=='S0' & !is.na(h) & !is.na(t),
    .(imputation_id,scope,guideline,h,t,n,denominator,percent=proportion)])
  if(!nrow(z)||anyDuplicated(z,by=c('imputation_id','scope','guideline','h','t')))
    stop('10: invalid H/T count keys')
  z[,`:=`(positive_n=sum(n[h==1L|t==1L]),cells=.N),by=.(imputation_id,scope,guideline)]
  if(any(z$cells!=4L))stop('10: H/T must have exactly four cells')
  z[,share_SHT_upgraded:=ifelse(h==0L&t==0L|positive_n<=0,NA_real_,n/positive_n)]
  z[,c('positive_n','cells'):=NULL]
  setcolorder(z,c('imputation_id','scope','guideline','h','t','n','denominator','percent','share_SHT_upgraded'))
  z
}

# Classification/requirements and upgraded risk only. H/T counts/capture and
# high-risk capture are computed centrally, not recomputed inside this function.
m10_core <- function(m,cfg,x,scopes){
  N<-nrow(x);desc<-risk<-list()
  for(sc in names(scopes)){
    ix<-which(scopes[[sc]]);n<-length(ix)
    for(year in c('2019','2025')){
      base<-x[[if(year=='2019')'score_2019_category' else 'score2_2025_category']]
      base_n<-sum(base[ix]=='M')
      s0_n<-sum(x[[paste0('A_',year,'_S0')]][ix])
      for(k in m07_rms()){
        nm<-paste(year,k,sep='_')
        cat<-x[[paste0('risk_',nm)]]
        up<-which(scopes[[sc]] & x[[paste0('upgrade_',nm)]]==1L)
        hv<-which(scopes[[sc]] & cat%in%c('H','VH'))
        a<-x[[paste0('A_',nm)]]
        added<-sum(a[ix]==1 & x[[paste0('A_',year,'_S0')]][ix]==0)
        count<-c(n_upgrade=length(up),n_hvh=length(hv),n_intervention=sum(a[ix]),
          n_additional=added,relative_additional_vs_S0=if(s0_n>0)added/s0_n else NA_real_,
          sum_D=sum(x[[paste0('D_',nm)]][ix]),mean_D=mean(x[[paste0('D_',nm)]][ix]))
        for(metric in names(count)){
          denom<-if(metric=='n_upgrade')base_n else if(metric%in%c('n_hvh','n_intervention'))n else NA_real_
          desc[[length(desc)+1L]]<-data.table(imputation_id=m,scope=sc,guideline=year,scenario=k,
            metric=metric,value=as.numeric(count[[metric]]),denominator=denom,
            proportion=if(is.finite(denom)&&denom>0)as.numeric(count[[metric]])/denom else NA_real_)
        }
        for(tm in cfg$times){
          if(length(up)){
            z<-m09_aj_embed(x,up,tm)
            r<-m09_aj_row(z$Q,z$iid,m,tm,k,paste(year,sc,sep='__'),'observed_upgraded_risk')
          }else r<-m09_aj_row(NA_real_,rep(NA_real_,N),m,tm,k,
            paste(year,sc,sep='__'),'observed_upgraded_risk')
          risk[[length(risk)+1L]]<-r
        }
      }
    }
  }
  cross<-x[,.(n=.N),by=.(score_2019_category,score2_2025_category)]
  cross[,`:=`(imputation_id=m,denominator=N,percent=n/N)]
  list(description=rbindlist(desc,fill=TRUE),observed=rbindlist(risk,fill=TRUE),
       reclassification=cross)
}

# Both H/VH and its complement are included in the denominator:
# sum_j P(j|scope)*AJ_j(t). This is the SAME event-burden normalization as H/T.
# Paired S0 differences use the paired influence contributions from the same run.
m10_hvh_outputs <- function(m,cfg,x,scopes){
  if(!identical(as.character(m07_rms()),c('S0','SH','ST','SHT')))
    stop('10: H/VH paired comparisons require scenario ordering S0, SH, ST, SHT')
  observed<-paired<-list()
  for(sc in names(scopes))for(year in c('2019','2025')){
    ix<-which(scopes[[sc]])
    for(tm in cfg$times){
      base<-NULL
      for(k in m07_rms()){
        high<-ix[as.character(x[[paste0('risk_',year,'_',k)]][ix])%in%c('H','VH')]
        rest<-ix[!ix%in%high]
        z<-m10_partition(x,ix,list(high,rest),tm)
        values<-list(risk=z$risks[[1]]$Q,capture=z$capture[1])
        influences<-list(risk=z$risks[[1]]$iid,capture=z$capture_if[,1L])
        for(metric in names(values)){
          r<-m09_aj_row(values[[metric]],influences[[metric]],m,tm,k,
            paste(year,sc,sep='__'),paste0('observed_HVH_',metric))
          observed[[length(observed)+1L]]<-r
          if(k!='S0'){
            q<-m09_aj_row(values[[metric]]-base$values[[metric]],
              influences[[metric]]-base$influences[[metric]],m,tm,k,sc,
              paste0('observed_HVH_',metric,'_change_vs_S0'))
            q[,guideline:=year]
            paired[[length(paired)+1L]]<-q
          }
        }
        if(k=='S0')base<-list(values=values,influences=influences)
      }
    }
  }
  list(observed=rbindlist(observed,fill=TRUE),paired=rbindlist(paired,fill=TRUE))
}

# Raw RM domains are retained; any H/T flags are counted FROM the authoritative
# four-cell partition rather than recomputed from separate H/T indexing.
m10_domains <- function(m,x,scopes,ht){
  domains<-list();defs<-m07_rm_definitions()
  for(sc in names(scopes))for(year in c('2019','2025')){
    group<-ht[[paste(sc,year,sep='__')]]
    mi<-group$population
    H<-vapply(group$members,length,integer(1L))
    pairs<-group$meta
    total<-length(mi)
    anycounts<-c(any_H=sum(H[pairs$h==1L]),any_T=sum(H[pairs$t==1L]),
      any_HT=sum(H[pairs$h==1L|pairs$t==1L]))
    cols<-c(paste0('rm_',defs[get(paste0('include_',year))==1,domain]),names(anycounts))
    for(v in cols){
      if(v%in%names(anycounts))one<-as.numeric(anycounts[[v]])
      else{
        val<-as.integer(as.character(x[[v]][mi]))
        if(anyNA(val)||any(!val%in%0:1))stop('10: invalid RM domain ',v,': ',year,' ',sc)
        one<-sum(val==1L)
      }
      for(level in 0:1){
        ct<-if(level==1L)one else total-one
        domains[[length(domains)+1L]]<-data.table(imputation_id=m,scope=sc,
          guideline=year,domain=v,level=level,n=ct,denominator=total,
          proportion=if(total)ct/total else NA_real_)
      }
    }
  }
  rbindlist(domains,fill=TRUE)
}

m10_reclassification_grid <- function(cross,m,N){
  grid<-CJ(score_2019_category=c('L','M','H','VH'),
           score2_2025_category=c('L','M','H','VH'))
  out<-merge(grid,cross,by=c('score_2019_category','score2_2025_category'),all.x=TRUE)
  out[is.na(n),`:=`(n=0L,imputation_id=m,denominator=N,percent=0)]
  out[,row_denominator:=sum(n),by=score_2019_category]
  out[,column_denominator:=sum(n),by=score2_2025_category]
  out[,row_proportion:=ifelse(row_denominator>0,n/row_denominator,NA_real_)]
  out
}

m10_reduction_changes <- function(d){
 q<-d[metric=='sum_D',.(imputation_id,scope,guideline,scenario,total_required_reduction=value)]
 ref<-q[scenario=='S0',.(imputation_id,scope,guideline,reference_total=total_required_reduction)]
 q<-merge(q,ref,by=c('imputation_id','scope','guideline'));q[,change_vs_S0:=total_required_reduction-reference_total]
 a<-q[guideline=='2019',.(imputation_id,scope,scenario,total_2019=total_required_reduction,change_2019=change_vs_S0)]
 b<-q[guideline=='2025',.(imputation_id,scope,scenario,total_2025=total_required_reduction,change_2025=change_vs_S0)]
 v<-merge(a,b,by=c('imputation_id','scope','scenario'));v[,`:=`(guideline_difference=total_2025-total_2019,change_guideline_difference_vs_S0=change_2025-change_2019)]
 list(totals=q,differences=v)
}


m10_partition <- function(x, population, members, time) {
  N<-nrow(x); n<-length(population)
  if(anyDuplicated(unlist(members)) || !setequal(unlist(members),population))stop('Partition must be exhaustive and disjoint')
  K<-length(members); a<-numeric(K); ia<-matrix(0,N,K); risks<-vector('list',K)
  for(j in seq_len(K)){
    ix<-members[[j]]; nj<-length(ix); p<-if(n)nj/n else NA_real_
    z<-if(nj)m09_aj_embed(x,ix,time)else list(Q=NA_real_,iid=rep(NA_real_,N))
    risks[[j]]<-z
    if(!nj)next
    ip<-numeric(N);ip[population]<-(as.integer(population%in%ix)-p)/n
    a[j]<-p*z$Q;ia[,j]<-p*z$iid+z$Q*ip
  }
  A<-sum(a); iA<-rowSums(ia)
  cap<-if(is.finite(A)&&A>0)a/A else rep(NA_real_,K)
  ic<-if(is.finite(A)&&A>0)ia/A-outer(iA,a/A^2)else matrix(NA_real_,N,K)
  valid<-is.finite(A)&&A>0
  if(valid && (abs(sum(cap)-1)>1e-10 || max(abs(rowSums(ic)))>1e-9))stop('Capture point/IF additivity failed')
  list(risks=risks,capture=cap,capture_if=ic,burden=a,denominator=A,
       passed=valid,if_sum_error=if(valid)max(abs(rowSums(ic)))else NA_real_)
}
m10_origin_category <- function(x,year,scenario){
  b<-as.character(x[[if(year=='2019')'score_2019_category'else 'score2_2025_category']])
  up<-x[[paste0('upgrade_',year,'_',scenario)]]
  if(anyNA(b)||any(!b%in%c('L','M','H','VH'))||anyNA(up))stop('Invalid risk category/upgrade')
  if(any(up==1L & b!='M'))stop('Only original M can upgrade')
  expected<-ifelse(up==1L,'H',b)
  if(!identical(expected,as.character(x[[paste0('risk_',year,'_',scenario)]])))stop('Origin category disagrees with final risk category')
  ifelse(b=='M',ifelse(up==1L,'M_upgraded_H','M_not_upgraded'),b)
}

m10_v5_one <- function(m,cfg,x,ht){
  N<-nrow(x);risk<-counts<-qa<-labels<-list()
  add_partition<-function(pop,members,meta,type){
    stopifnot(nrow(meta)==length(members))
    c0<-copy(meta);c0[,`:=`(imputation_id=m,n=lengths(members),denominator=length(pop),proportion=if(length(pop))lengths(members)/length(pop)else NA_real_)]
    counts[[length(counts)+1L]]<<-c0
    for(tm in cfg$times){
      z<-m10_partition(x,pop,members,tm)
      for(j in seq_along(members))for(metric in c('risk','capture')){
        v<-if(metric=='risk')z$risks[[j]]$Q else z$capture[j]
        iid<-if(metric=='risk')z$risks[[j]]$iid else z$capture_if[,j]
        r<-m09_aj_row(v,iid,m,tm,meta$scenario[j],meta$scope[j],paste0(type,'_',metric))
        for(nm in setdiff(names(meta),c('scenario','scope')))r[,(nm):=meta[[nm]][j]]
        r[,`:=`(cell_n=length(members[[j]]),cell_status=if(!length(members[[j]]))'empty'else if(length(members[[j]])<2)'insufficient_n'else if(!any(x$event_status[members[[j]]]==1L & x$followup_years[members[[j]]]<=tm))'no_ASCVD_events'else 'estimable',
                 denominator_method='sum_stratum_proportion_times_AJ',event_burden=z$burden[j],total_event_burden=z$denominator)]
        risk[[length(risk)+1L]]<<-r
      }
      if('region'%in%names(meta))for(reg in unique(meta$region)){
        jj<-which(meta$region==reg)
        r<-m09_aj_row(sum(z$capture[jj]),rowSums(z$capture_if[,jj,drop=FALSE]),m,tm,meta$scenario[1],meta$scope[1],'joint_region_capture')
        r[,`:=`(guideline='joint',cell_id=paste0('region__',reg),region=reg,cell_n=sum(lengths(members)[jj]),cell_status='aggregate',
          denominator_method='sum_stratum_proportion_times_AJ',event_burden=sum(z$burden[jj]),total_event_burden=z$denominator)]
        risk[[length(risk)+1L]]<<-r
      }
      qa[[length(qa)+1L]]<<-data.table(imputation_id=m,type=type,scope=meta$scope[1],guideline=meta$guideline[1],scenario=meta$scenario[1],time=tm,
        n=sum(lengths(members)),denominator=length(pop),capture_sum=sum(z$capture),max_abs_IF_sum=z$if_sum_error,status=if(z$passed)'passed'else 'not_estimable')
    }
  }
  # H/T membership comes from the one authoritative per-MI partition index.
  for(group in ht){
    m07_log('10','MI ',m,' H/T partition ',group$meta$scope[1],
      ' / ',group$meta$guideline[1])
    add_partition(group$population,group$members,group$meta,'HT_stratum')
  }
  lev<-c('L','M_not_upgraded','M_upgraded_H','H','VH')
  for(k in c('SH','ST','SHT')){
    m07_log('10','MI ',m,' joint classification ',k,' (25 cells per horizon)')
    a<-m10_origin_category(x,'2019',k);b<-m10_origin_category(x,'2025',k)
    labels[[k]]<-data.table(imputation_id=m,participant_id=x$participant_id,scenario=k,origin_2019=a,origin_2025=b)
    meta<-CJ(origin_2019=lev,origin_2025=lev);meta[,`:=`(scope='overall',guideline='joint',scenario=k,cell_id=paste(origin_2019,origin_2025,sep='__'))]
    hv<-c('M_upgraded_H','H','VH')
    meta[,region:=ifelse(origin_2019%in%hv,ifelse(origin_2025%in%hv,'both_HVH','removed_HVH'),ifelse(origin_2025%in%hv,'new_HVH','neither_HVH'))]
    members<-lapply(seq_len(25),function(j)which(a==meta$origin_2019[j]&b==meta$origin_2025[j]))
    add_partition(seq_len(N),members,meta,'joint_category')
  }
  list(risks=rbindlist(risk,fill=TRUE),counts=rbindlist(counts,fill=TRUE),qa=rbindlist(qa,fill=TRUE),labels=rbindlist(labels))
}

m10_v5_export <- function(z,cfg){
  d<-rbindlist(lapply(z,`[[`,'risks'),fill=TRUE)
  keys<-c('scenario','scope','guideline','cell_id','origin_2019','origin_2025','region','h','t','subgroup_var','reference','comparison','estimand','pool_scale','time')
  pooled<-m09_pool(d,cfg,intersect(keys,names(d)))
  info<-d[,.(min_cell_n=min(cell_n),max_cell_n=max(cell_n),mean_cell_n=mean(cell_n),cell_status=paste(sort(unique(cell_status)),collapse=';'),
    event_burden=mean(event_burden),total_event_burden=mean(total_event_burden),denominator_method=denominator_method[1]),by=intersect(keys,names(d))]
  pooled<-merge(pooled,info,by=intersect(keys,names(d)))
  pooled[min_cell_n==0,`:=`(lcl=NA_real_,ucl=NA_real_,status='empty_in_some_MI_CI_unavailable')]
  pooled[grepl('no_ASCVD_events|insufficient_n',cell_status),`:=`(lcl=NA_real_,ucl=NA_real_,status='boundary_or_sparse_cell_CI_unavailable')]
  for(type in c('HT_stratum','joint_category')){
    name<-if(type=='HT_stratum')'HT_stratum_risks_and_normalized_event_capture'else 'joint_origin_category_risks_and_event_capture'
    m07_export(cfg,'observed_risk_and_event_capture',name,d[startsWith(estimand,type)],pooled[startsWith(estimand,type)])
  }
  ct<-rbindlist(lapply(z,`[[`,'counts'),fill=TRUE);kk<-intersect(c('scope','guideline','scenario','cell_id','origin_2019','origin_2025','region','h','t'),names(ct))
  m07_export(cfg,'risk_classification','origin_category_partition_counts',ct,m10_descriptive_pool(ct,kk,c('n','denominator','proportion')))
  lab<-rbindlist(lapply(z,`[[`,'labels'));m07_result(cfg,'strategy_assignments','risk_category_origins',lab,'per_imputation')
  m07_result(cfg,'strategy_assignments','risk_category_origin_dictionary',data.table(code=c('L','M_not_upgraded','M_upgraded_H','H','VH'),definition=c('Original low','Original moderate, not upgraded','Original moderate, upgraded to high','Original high','Original very high')))
  m07_result(cfg,'analysis_quality_control','normalized_capture_partition_checks',rbindlist(lapply(z,`[[`,'qa')))
  # Aggregate the same cell burdens, so region capture and cell capture agree exactly.
  m07_export(cfg,'observed_risk_and_event_capture','joint_HVH_region_capture',d[estimand=='joint_region_capture'],pooled[estimand=='joint_region_capture'])
}




# Verify the per-MI classification and AJ cell populations BEFORE any export.
m10_check_one <- function(z,cfg){
  ct<-z$partition$counts[scenario=='S0' & guideline%in%c('2019','2025') & !is.na(h) & !is.na(t)]
  check<-ct[,.(n=sum(n),denominator=unique(denominator)[1L],n_den=uniqueN(denominator),cells=.N,
    share=sum(proportion)),by=.(imputation_id,scope,guideline)]
  if(!nrow(check)||any(check$cells!=4L|check$n_den!=1L|
      abs(check$n-check$denominator)>1e-6|abs(check$share-1)>1e-9))
    stop('10: H/T counts do not exhaust the original moderate-risk population')
  r<-unique(z$partition$risks[estimand%in%c('HT_stratum_risk','HT_stratum_capture'),
      .(imputation_id,scope,guideline,scenario,h,t,time,estimand,cell_n)])
  want<-merge(r,ct[,.(imputation_id,scope,guideline,scenario,h,t,n)],
    by=c('imputation_id','scope','guideline','scenario','h','t'),all.x=TRUE)
  if(!nrow(want)||anyNA(want$n)||any(abs(want$cell_n-want$n)>1e-6))
    stop('10: H/T AJ population does not equal its authoritative classification count')
  invisible(TRUE)
}

m10_one <- function(m,cfg){
  x<-m07_read(cfg,m);m07_unique(x)
  scopes<-m07_scope_spec(x)
  ht<-m10_ht_groups(x,scopes)
  partition<-m10_v5_one(m,cfg,x,ht)
  z<-m10_core(m,cfg,x,scopes)
  hvh<-m10_hvh_outputs(m,cfg,x,scopes)
  z$observed<-rbindlist(list(z$observed,hvh$observed),fill=TRUE)
  ht_risks<-copy(partition$risks[estimand=='HT_stratum_risk'])
  ht_risks[,estimand:='observed_HT_stratum_risk']
  z$paired_observed<-rbindlist(list(hvh$paired,ht_risks),fill=TRUE)
  z$domains<-m10_domains(m,x,scopes,ht)
  z$overlap<-m10_overlap_from_counts(partition$counts)
  z$reclassification<-m10_reclassification_grid(z$reclassification,m,nrow(x))
  z$partition<-partition
  m10_check_one(z,cfg)
  z
}

# A new versioned name prevents accidental reading of either retired cache.
# Cache contents and signatures are saved atomically. Old files are never deleted.
m10_task <- function(m,cfg){
  cache_dir<-m07_dir(cfg,'10')
  p<-file.path(cache_dir,sprintf('classification_unified_v6_m%02d.rds',m))
  sigp<-paste0(p,'.signature.rds')
  code_files<-file.path(cfg$code_dir,c('Function_07_Guideline_strategies_M01.R',
    'Function_08_Counterfactual_risk_M01.R','Function_09_Main_analysis_M01.R',
    'Function_10_RM_classification_M01.R'))
  if(!all(file.exists(code_files)))stop('10: code signature source missing: ',
    paste(basename(code_files[!file.exists(code_files)]),collapse=', '))
  optional_dependencies<-file.path(cfg$code_dir,
    c('Function_RM_sources_M01.R','Function_M01_Analysis_common.R'))
  code_files<-c(code_files,optional_dependencies[file.exists(optional_dependencies)])
  input<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
  if(!file.exists(input))stop('10: missing 07 input: ',input)
  sig<-list(contract='classification_unified_v6',input=m07_md5(input),
    times=cfg$times,scientific_key=m07_scientific_key(cfg),
    code=vapply(code_files,m07_md5,character(1)))
  rebuild<-tolower(Sys.getenv('M01_CLASSIFICATION_REBUILD','0'))%in%c('1','true','yes')
  if(!rebuild&&file.exists(p)&&file.exists(sigp)&&identical(readRDS(sigp),sig)){
    m07_log('10','MI ',m,' unified cache hit')
    z<-readRDS(p)
    m10_check_one(z,cfg)
    return(z)
  }
  z<-m07_timed('10',paste('unified classification and AJ MI',m),m10_one(m,cfg))
  m07_atomic(z,p);m07_atomic(sig,sigp)
  z
}

m10_export_unified <- function(z,cfg){
  if(length(z)!=cfg$m)stop('10: expected exactly ',cfg$m,' MI results')
  for(o in z)m10_check_one(o,cfg)
  # Same per-MI objects are the sole source for ALL related CSVs.
  m10_v5_export(lapply(z,`[[`,'partition'),cfg)
  m09_emit(z,cfg,'observed_risk_and_event_capture',
    c(observed='high_risk_and_upgraded_risks_capture',
      paired_observed='high_risk_paired_changes_and_HT_risks'))
  for(part in c('description','overlap','reclassification','domains')){
    d<-rbindlist(lapply(z,`[[`,part),fill=TRUE)
    if(part=='description'){
      for(board in c('intervention_requirements','risk_classification')){
        dd<-if(board=='intervention_requirements')d[!metric%in%c('n_upgrade','n_hvh')]
          else d[metric%in%c('n_upgrade','n_hvh')]
        m07_export(cfg,board,'classification_and_requirements',dd,
          m10_mean(dd,c('scope','guideline','scenario','metric')))
      }
    }else{
      keys<-switch(part,overlap=c('scope','guideline','h','t'),
        reclassification=c('score_2019_category','score2_2025_category'),
        domains=c('scope','guideline','domain','level'))
      vals<-setdiff(names(d),c(keys,'imputation_id'))
      m07_export(cfg,'risk_classification',part,d,m10_descriptive_pool(d,keys,vals))
    }
  }
  d<-rbindlist(lapply(z,`[[`,'description'),fill=TRUE)
  changes<-m10_reduction_changes(d)
  for(nm in names(changes)){
    a<-changes[[nm]];keys<-intersect(c('scope','guideline','scenario'),names(a))
    vals<-setdiff(names(a),c(keys,'imputation_id'))
    m07_export(cfg,'intervention_requirements',paste0('required_reduction_',nm),
      a,m10_descriptive_pool(a,keys,vals))
  }
  r<-rbindlist(lapply(z,`[[`,'reclassification'))
  levels<-c('L','M','H','VH')
  r[,direction:=c('lower','unchanged','higher')[sign(match(score2_2025_category,levels)-
    match(score_2019_category,levels))+2L]]
  a<-r[,.(n=sum(n),denominator=unique(denominator)),by=.(imputation_id,direction)]
  a[,proportion:=n/denominator]
  m07_export(cfg,'risk_classification','reclassification_direction',a,
    m10_descriptive_pool(a,'direction',c('n','denominator','proportion')))
  # Commit marker is written LAST. A failed partial export leaves MD5 mismatch
  # (or no manifest), and 13 refuses to combine mixed-version output files.
  m10_export_manifest(cfg)
  invisible(TRUE)
}

m10_export_manifest <- function(cfg){
  rel<-c(
    ht_counts='risk_classification/pooled/origin_category_partition_counts.csv',
    ht_risk='observed_risk_and_event_capture/pooled/HT_stratum_risks_and_normalized_event_capture.csv',
    rm_domains='risk_classification/pooled/domains.csv',
    rm_overlap_compat='risk_classification/pooled/overlap.csv',
    risk_requirements='risk_classification/pooled/classification_and_requirements.csv',
    intervention_requirements='intervention_requirements/pooled/classification_and_requirements.csv',
    hvh_risk='observed_risk_and_event_capture/pooled/high_risk_and_upgraded_risks_capture.csv',
    hvh_paired='observed_risk_and_event_capture/pooled/high_risk_paired_changes_and_HT_risks.csv',
    reclassification='risk_classification/pooled/reclassification.csv')
  full<-file.path(cfg$out_dir,unname(rel))
  if(!all(file.exists(full)))stop('10: cannot finalize unified manifest; missing: ',
    paste(rel[!file.exists(full)],collapse=', '))
  data<-data.table(source=names(rel),relative_path=unname(rel),
    md5=vapply(full,m07_md5,character(1L)),bytes=as.numeric(file.info(full)$size),
    contract='classification_unified_v6',m=cfg$m)
  target<-file.path(cfg$out_dir,'analysis_quality_control','pooled',
    'classification_unified_export_manifest.csv')
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  temporary<-tempfile('classification_manifest_',tmpdir=dirname(target),fileext='.csv')
  on.exit(unlink(temporary),add=TRUE)
  fwrite(data,temporary,bom=TRUE)
  # Publish only after every other output was written and checked.
  if(file.exists(target)&&!file.remove(target))stop('10: unable to replace old manifest: ',target)
  if(!file.rename(temporary,target))stop('10: cannot commit unified export manifest: ',target)
  invisible(data)
}



# Publication summaries and exports (same-stage functions).
# Independent read-only supplement; no model bundle and no legacy cache mutation.
m10_publication_summaries <- function(cfg) {
  risks<-counts<-inputs<-list()
  for(m in seq_len(cfg$m)) {
    x<-as.data.table(m07_read(cfg,m));m07_unique(x)
    groups<-m10_ht_groups(x,list(overall=rep(TRUE,nrow(x))))
    for(year in c('2019','2025')) {
      partition<-groups[[paste('overall',year,sep='__')]]
      pop<-partition$population
      a0<-x[[paste0('A_',year,'_S0')]]
      positive<-unlist(partition$members[partition$meta$h==1L|partition$meta$t==1L],use.names=FALSE)
      joint_new<-which(x[[paste0('A_',year,'_SHT')]]==1L & a0==0L)
      emit<-function(ix,up,sc,cell){
        if(!all(ix%in%pop)||!all(ix%in%up)||length(ix)>length(up))stop('New eligibility subset QA failed')
        counts[[length(counts)+1L]]<<-data.table(imputation_id=m,guideline=year,
          scenario=sc,cell=cell,n_new=length(ix),n_upgrade=length(up),n_initial_M=length(pop),
          proportion_new=if(length(up))length(ix)/length(up)else NA_real_)
        z<-if(length(ix))m09_aj_embed(x,ix,15)else list(Q=NA_real_,iid=rep(NA_real_,nrow(x)))
        r<-m09_aj_row(z$Q,z$iid,m,15,sc,'overall','observed_newly_eligible_risk')
        r[,`:=`(guideline=year,cell=cell,n_new=length(ix),n_upgrade=length(up))]
        risks[[length(risks)+1L]]<<-r
      }
      for(sc in c('SH','ST','SHT')){
        a<-x[[paste0('A_',year,'_',sc)]]
        ix<-which(a==1L & a0==0L);up<-which(x[[paste0('upgrade_',year,'_',sc)]]==1L)
        if(any(a<a0)||length(ix)!=sum(a)-sum(a0))stop('Intervention monotonicity QA failed')
        emit(ix,up,sc,'all')
      }
      sum_n<-0L
      for(j in which(partition$meta$h==1L|partition$meta$t==1L)){
        up<-partition$members[[j]];ix<-intersect(up,joint_new);sum_n<-sum_n+length(ix)
        emit(ix,up,'SHT',partition$meta$cell_id[j])
      }
      if(sum_n!=length(joint_new)||!setequal(positive,which(x[[paste0('upgrade_',year,'_SHT')]]==1L)))
        stop('SNT union/partition QA failed')
    }
    p<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
    inputs[[m]]<-data.table(imputation_id=m,path=normalizePath(p,winslash='/'),md5=m07_md5(p))
    rm(x);gc(FALSE);m07_log('10','Publication summaries MI ',m,'/',cfg$m)
  }
  d<-rbindlist(risks);ct<-rbindlist(counts)
  keys<-c('scenario','scope','subgroup_var','reference','comparison','estimand','pool_scale','time','guideline','cell')
  pooled<-m09_pool(d,cfg,keys)
  empty<-ct[,.(all_empty=all(n_new==0)),by=.(scenario,guideline,cell)]
  pooled<-merge(pooled,empty,by=c('scenario','guideline','cell'))
  pooled[all_empty==TRUE,status:='not_applicable'];pooled[,all_empty:=NULL]
  m07_export(cfg,'observed_risk_and_event_capture','rm_newly_eligible_observed_risks',d,pooled)
  m07_export(cfg,'intervention_requirements','rm_newly_eligible_counts',ct,
    ct[,lapply(.SD,function(v)if(all(is.na(v)))NA_real_ else mean(v)),
      by=.(guideline,scenario,cell),.SDcols=c('n_new','n_upgrade','n_initial_M','proportion_new')])
  m07_result(cfg,'analysis_quality_control','publication_summary_input_manifest',rbindlist(inputs))
  invisible(pooled)
}
