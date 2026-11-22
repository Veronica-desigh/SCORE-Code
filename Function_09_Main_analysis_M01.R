# MODULE 09 | full population / subgroups / IF contrasts / Rubin / P
source(file.path(Sys.getenv("M01_CODE_DIR","."),"Function_08_Counterfactual_risk_M01.R"),local=FALSE)

m09_one_effect <- function(pred,weights,estimand,scope,scenario,m,dfcom,group="",reference="") {
  states<-names(weights);if(!all(states%in%names(pred)))stop("Missing predictions: ",paste(setdiff(states,names(pred)),collapse=","))
  nr<-length(pred[[1]]$estimate);n<-nrow(pred[[1]]$iid)
  is_ratio<-grepl("^rr_",estimand)
  if(is_ratio){
    a<-pred[[states[1]]];b<-pred[[states[2]]]
    if(any(a$estimate<=0|b$estimate<=0))return(data.table())
    q<-log(a$estimate/b$estimate)
    mat<-sweep(a$iid,2,a$estimate,"/")-sweep(b$iid,2,b$estimate,"/")
  } else {
    q<-Reduce(`+`,Map(function(state,weight)pred[[state]]$estimate*weight,states,unname(weights)))
    mat<-Reduce(`+`,Map(function(state,weight)pred[[state]]$iid*weight,states,unname(weights)))
  }
  if(!identical(dim(mat),c(n,as.integer(nr))))stop("IF length/order mismatch")
  if(any(!is.finite(mat)))stop("Nonfinite IF ",estimand)
  data.table(imputation_id=m,time=pred[[1]]$times,scenario=scenario,scope=scope,
    subgroup_var=group,reference=reference,estimand=estimand,pool_scale=if(is_ratio)"log" else "natural",
    Q=as.numeric(q),U=colSums(mat^2),dfcom=as.numeric(dfcom),n_standardized=if(is.null(pred[[1]]$n_standardized))NA_real_ else pred[[1]]$n_standardized)
}
m09_pred <- function(bundle,scope,scenario,transitions=NULL){
  key<-if(is.null(transitions))paste0("scope__",scope) else paste0("transition__",scenario,"__",scope,"__",transitions)
  if(!key%in%names(bundle$targets$index))return(NULL)
  required<-c("natural",paste0("2019_",scenario),paste0("2025_",scenario))
  preds<-setNames(lapply(required,function(a)m08_prediction(bundle,a,key)),required)
  preds
}
m09_contrasts <- function(pred,scope,scenario,m,dfcom,extra=TRUE){
  s19<-paste0("2019_",scenario);s25<-paste0("2025_",scenario)
  terms<-list(risk_natural=c(natural=1),risk_2019=setNames(1,s19),risk_2025=setNames(1,s25),
              rd_2019_vs_natural=setNames(c(1,-1),c(s19,"natural")),
              rd_2025_vs_natural=setNames(c(1,-1),c(s25,"natural")),
              rd_2025_vs_2019=setNames(c(1,-1),c(s25,s19)),
              arr_2019_vs_natural=setNames(c(1,-1),c("natural",s19)),
              arr_2025_vs_natural=setNames(c(1,-1),c("natural",s25)),arr_difference_2025_vs_2019=setNames(c(1,-1),c(s19,s25)),
              rr_2019_vs_natural=setNames(c(1,-1),c(s19,"natural")),
              rr_2025_vs_natural=setNames(c(1,-1),c(s25,"natural")),rr_2025_vs_2019=setNames(c(1,-1),c(s25,s19)))
  if(!extra)terms<-terms[!grepl("^arr_|^rr_",names(terms))]
  out<-lapply(names(terms),function(v)m09_one_effect(pred,terms[[v]],v,scope,scenario,m,dfcom))
  ans<-data.table::rbindlist(out,fill=TRUE)
  key<-paste0("scope__",scope)
  ans
}
m09_het <- function(bundle,scenario,m,dfcom){
  spec<-m07_subgroups();out<-list()
  for(i in seq_len(nrow(spec))){
    row<-spec[i];ref<-m09_pred(bundle,row$reference,scenario);cmp<-m09_pred(bundle,row$comparison,scenario)
    if(is.null(ref)||is.null(cmp))stop("Missing subgroup prediction")
    a<-paste0("2025_",scenario);b<-paste0("2019_",scenario)
    H<-length(bundle$engine$times)
    Q<-cmp[[a]]$estimate-cmp[[b]]$estimate-ref[[a]]$estimate+ref[[b]]$estimate
    IF<-cmp[[a]]$iid-cmp[[b]]$iid-ref[[a]]$iid+ref[[b]]$iid
    out[[i]]<-data.table(imputation_id=m,time=bundle$engine$times,scenario=scenario,scope=paste0("heterogeneity_",row$subgroup_var),
      subgroup_var=row$subgroup_var,reference=row$reference,comparison=row$comparison,
      estimand="delta_rd_2025_vs_2019",pool_scale="natural",Q=Q,U=colSums(IF^2),dfcom=dfcom,
      n_standardized=NA_real_)
  }
  rbindlist(out,fill=TRUE)
}
# Barnard-Rubin t df, same df used for confidence intervals and two-sided P values.
m09_rubin <- function(Q,U,dfcom=Inf,confidence=.95){
  M<-length(Q);if(M<2 || any(!is.finite(Q))||any(!is.finite(U)|U<0))stop("Invalid Rubin inputs")
  qbar<-mean(Q);ubar<-mean(U);B<-stats::var(Q);T<-ubar+(1+1/M)*B
  if(!is.finite(T)||T<0)stop("Rubin variance invalid")
  lambda<-if(T>0)(1+1/M)*B/T else 0
  df_old<-if(lambda <= .Machine$double.eps)Inf else (M-1)/(lambda^2)
  if(is.finite(dfcom)&&dfcom>0){dfobs<-((dfcom+1)/(dfcom+3))*dfcom*(1-lambda);df<-if(dfobs>0)1/(1/df_old+1/dfobs) else df_old}else df<-df_old
  se<-sqrt(T);if(se>0){stat<-qbar/se;p<-if(is.finite(df))2*stats::pt(-abs(stat),df) else 2*stats::pnorm(-abs(stat))}else {stat<-NA_real_;p<-NA_real_}
  crit<-if(is.finite(df))stats::qt((1+confidence)/2,df) else stats::qnorm((1+confidence)/2)
  list(estimate=qbar,se=se,lcl=qbar-crit*se,ucl=qbar+crit*se,df=df,t_statistic=stat,p_value=p,
       within_variance=ubar,between_variance=B,total_variance=T,lambda=lambda,M=M)
}
m09_pool <- function(data,cfg,keys=c("scenario","scope","subgroup_var","reference","comparison","estimand","pool_scale","time")) {
  d<-data.table::copy(data.table::as.data.table(data))
  if(!nrow(d))return(d)
  keys<-unique(c(keys,intersect(c("guideline","rm_comparison","domain","h","t"),names(d))))
  if(!"comparison"%in%names(d))d[,comparison:=NA_character_]
  if("transition"%in%names(d))keys<-c(keys,"transition")
  if(anyDuplicated(d,by=c(keys,"imputation_id")))stop("Duplicate estimand per MI")
  d[,{
    valid<-(.N==cfg$m && uniqueN(imputation_id)==cfg$m && all(is.finite(Q)&is.finite(U)&U>=0))
    if(!valid){list(estimate=NA_real_,se=NA_real_,lcl=NA_real_,ucl=NA_real_,df=NA_real_,t_statistic=NA_real_,p_value=NA_real_,
         within_variance=NA_real_,between_variance=NA_real_,total_variance=NA_real_,
         m=uniqueN(imputation_id),m_valid=sum(is.finite(Q)&is.finite(U)&U>=0),status="invalid_incomplete_MI")}
    else{
      r<-m09_rubin(Q,U,min(dfcom),cfg$confidence_level)
      islog<-pool_scale[1]=="log";list(estimate=if(islog)exp(r$estimate)else r$estimate,
        se=if(islog)exp(r$estimate)*r$se else r$se,
        lcl=if(islog)exp(r$lcl)else r$lcl,ucl=if(islog)exp(r$ucl)else r$ucl,
        df=r$df,t_statistic=r$t_statistic,p_value=if(grepl("heterogeneity|delta_rd",estimand[1]))r$p_value else NA_real_,
        within_variance=r$within_variance,between_variance=r$between_variance,total_variance=r$total_variance,
        m=r$M,m_valid=r$M,status="estimable")}
  },by=keys]
}
m09_full_one <- function(m,cfg,scenarios="S0"){
  bundle<-m08_bundle(cfg,m);dfcom<-bundle$engine$complete_df
  out<-list();hets<-list()
  scopes<-names(m07_scope_spec(m07_read(cfg,m)))
  for(k in scenarios) {
    for(scope in scopes){
      preds<-m09_pred(bundle,scope,k)
      if(is.null(preds))stop("Missing 08 scope prediction ",scope)
      out[[length(out)+1L]]<-m09_contrasts(preds,scope,k,m,dfcom)
    }
    hets[[k]]<-m09_het(bundle,k,m,dfcom)
  }
  list(effects=rbindlist(out,fill=TRUE),heterogeneity=rbindlist(hets,fill=TRUE))
}

# Aalen-Johansen selected-time CIF and empirical influence contributions.
# Exact AJ point; selected-time influence from the two-cause counting process.
m09_aj <- function(times,status,t,full_n=length(times)) {
  times<-as.numeric(times);status<-as.integer(status)
  n<-length(times)
  if(n<2L)return(list(Q=NA_real_,iid=rep(NA_real_,n)))
  if(anyNA(times)||anyNA(status)||any(times<=0)||any(!status%in%0:2))stop("Invalid AJ inputs")
  evtime<-sort(unique(times[status!=0 & times<=t]))
  if(!length(evtime))return(list(Q=0,iid=numeric(n)))
  etid<-match(times[status!=0&times<=t],evtime)
  d1<-tabulate(etid[status[status!=0&times<=t]==1],length(evtime))
  d2<-tabulate(etid[status[status!=0&times<=t]==2],length(evtime))
  # findInterval depends on sorted data; instead use reverse cumulative distribution.
  sorted<-sort(times)
  # At an event time, all individuals with T >= event time are at risk.
  Y<-n-findInterval(evtime-1e-9,sorted)
  if(any(Y<=0)||any(d1+d2>Y))stop("Invalid competing-risk risk set")
  h1<-d1/Y;h2<-d2/Y; htot<-h1+h2
  s_after<-cumprod(1-htot);sprev<-c(1,head(s_after,-1L))
  inc<-sprev*h1;ft<-sum(inc);f_by<-cumsum(inc)
  tail<-ft-f_by;den<-pmax(1-htot,.Machine$double.eps)
  g1<-sprev-tail/den;g2<--tail/den
  v<-(g1*h1+g2*h2)/Y
  pos<-findInterval(times,evtime);base<-c(0,cumsum(v))[pos+1L]
  event_own<-numeric(n);sel<-which(status!=0 & times<=t)
  j<-match(times[sel],evtime)
  event_own[sel]<-ifelse(status[sel]==1L,g1[j],g2[j])/Y[j]
  iid<-event_own-base
  # selected-time IF uses risk-set finite sample formulas; must center for summed variance
  iid<-iid-mean(iid)
  list(Q=as.numeric(ft),iid=iid,events=sum(d1),n=n)
}
m09_aj_embed <- function(x,ix,t){
  z<-m09_aj(x$followup_years[ix],x$event_status[ix],t)
  iid<-numeric(nrow(x));iid[ix]<-z$iid
  list(Q=z$Q,iid=iid,n=length(ix))
}
m09_aj_capture <- function(x,group_ix,t,scope_ix=seq_len(nrow(x))){
  n<-length(scope_ix);g<-intersect(scope_ix,group_ix)
  if(!length(g))return(list(risk=NA_real_,risk_if=rep(NA_real_,nrow(x)),capture=NA_real_,capture_if=rep(NA_real_,nrow(x))))
  f<-m09_aj_embed(x,g,t);tot<-m09_aj_embed(x,scope_ix,t)
  prop<-length(g)/n;ip<-numeric(nrow(x));ip[scope_ix]<- (as.integer(scope_ix%in%g)-prop)/n
  num<-prop*f$Q;inum<-prop*f$iid+f$Q*ip
  if(tot$Q<=0){capture<-NA_real_;icap<-rep(NA_real_,nrow(x))}
  else{capture<-num/tot$Q;icap<-inum/tot$Q-num/(tot$Q^2)*tot$iid}
  list(risk=f$Q,risk_if=f$iid,capture=capture,capture_if=icap,p=prop,denominator=tot$Q)
}
m09_aj_row <- function(value,iid,m,time,scenario,scope,estimand){
  data.table(imputation_id=m,time=time,scenario=scenario,scope=scope,
    subgroup_var="",reference="",comparison="",estimand=estimand,pool_scale="natural",
    Q=value,U=if(all(is.finite(iid)))sum(iid^2) else NA_real_,dfcom=Inf,n_standardized=NA_real_)
}

m09_task <- function(m,cfg){
  m07_cached(cfg,'09',m,c(file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m)),file.path(m07_dir(cfg,'08'),sprintf('bundle_m%02d.rds',m))),function(m){
    z<-m09_full_one(m,cfg);x<-m07_read(cfg,m);obs<-list()
    for(scope in names(m07_scope_spec(x))){
      ix<-which(m07_scope_spec(x)[[scope]])
      for(t in cfg$times){o<-m09_aj_embed(x,ix,t);obs[[length(obs)+1L]]<-m09_aj_row(o$Q,o$iid,m,t,'S0',scope,'observed_scope_risk')}
    };z$observed<-rbindlist(obs);z
  })
}
m09_emit <- function(z,cfg,board,parts){
  for(part in names(parts)){
    d<-rbindlist(lapply(z,`[[`,part),fill=TRUE);if(!nrow(d))next
    pooled<-m09_pool(d,cfg);m07_export(cfg,board,parts[[part]],d,pooled)
  }
}
m09_evalues <- function(d){
  d<-copy(d[grepl('^rr_',estimand)]);if(!nrow(d))return(d)
  ev<-function(r){r<-pmax(r,1/r);r+sqrt(r*(r-1))}
  d[,ci_nearest_null:=ifelse(lcl<=1&ucl>=1,1,ifelse(ucl<1,ucl,lcl))]
  d[,`:=`(evalue_point=ev(estimate),evalue_ci=ev(ci_nearest_null))];d
}



