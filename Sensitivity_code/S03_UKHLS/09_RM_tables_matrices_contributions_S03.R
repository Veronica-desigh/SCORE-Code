# 09 | Descriptive reports, two matrices and fixed-S0 six-group contribution decomposition.
# Standalone; requires completed 08 (no re-imputation, no risk estimation, no confidence intervals).
# ========================== USER CONFIG ==========================
DATASET_ID <- "S03_UKHLS"   # Fixed UKHLS only
THREADS <- 5L
TEST_N <- 0L
SMOKING_MODE <- "same_wave"
PROJECT_ROOT <- ""
# ================================================================
DATASET_ID <- "S03_UKHLS"
THREADS <- as.integer(Sys.getenv("S03_RM_THREADS",Sys.getenv("S03_THREADS",as.character(THREADS))))
TEST_N <- as.integer(Sys.getenv("S03_RM_TEST_N",Sys.getenv("S03_ANALYSIS_TEST_N",as.character(TEST_N))))
PROJECT_ROOT <- Sys.getenv("S03_RM_PROJECT_ROOT",Sys.getenv("S03_PROJECT_ROOT",PROJECT_ROOT))
SMOKING_MODE <- Sys.getenv("S03_SMOKING_MODE",SMOKING_MODE)
if(!is.finite(THREADS)||THREADS<1||!is.finite(TEST_N)||TEST_N<0)stop("Invalid settings")
args <- grep("^--file=",commandArgs(),value=TRUE)
frames <- Filter(Negate(is.null),lapply(sys.frames(),function(z)z$ofile))
start <- if(length(frames))dirname(tail(frames,1)[[1]]) else if(length(args))dirname(sub("^--file=","",args[1])) else getwd()
root <- if(nzchar(PROJECT_ROOT))PROJECT_ROOT else normalizePath(start,winslash="/",mustWork=TRUE)
if(!nzchar(PROJECT_ROOT))while(!dir.exists(file.path(root,"LJW_outcome"))){up<-dirname(root);if(up==root)stop("Project root missing");root<-up}
out<-file.path(root,"LJW_outcome","data_intervention",DATASET_ID)
if(TEST_N>0)out<-file.path(out,paste0("test_",TEST_N))
if(SMOKING_MODE=="wave2_prior")out<-file.path(out,"wave2_smoking_sensitivity")
out<-file.path(out,"external_RM_analyses")
for(p in c("02_per_imputation","03_pooled","quality_control","logs"))dir.create(file.path(out,p),recursive=TRUE,showWarnings=FALSE)
t0<-Sys.time();logfile<-file.path(out,"logs","09_RM_tables.log")
log<-function(...){x<-sprintf("[%s +%.1fs] %s",format(Sys.time(),"%H:%M:%S"),as.numeric(difftime(Sys.time(),t0,units="secs")),paste0(...,collapse=""));message(x);cat(x,"\n",file=logfile,append=TRUE)}
write_bom<-function(z,path){tmp<-tempfile(fileext=".csv");on.exit(unlink(tmp),add=TRUE);utils::write.csv(z,tmp,row.names=FALSE,na="",fileEncoding="UTF-8");raw<-readBin(tmp,"raw",n=file.info(tmp)$size);con<-file(path,"wb");on.exit(close(con),add=TRUE);writeBin(as.raw(c(239,187,191)),con);writeBin(raw,con);invisible(path)}
groups<-c("unchanged_untreated","newly_intervention","intervention_removed", "treatment_intensified","unchanged_intensity","treatment_deintensified")
orig_levels<-c("L_score","M_score","H_rm","H_score","VH_score")
scenarios<-c("S0","SN","ST","SNT")

table_one <- function(task,outdir) {
  obj<-readRDS(task$path);imp<-obj$imputation_id
  tabs<-list(level=list(),overall=list(),year_diff=list(),group=list(),riskmat=list(),groupmat=list(),fixed=list(),change=list(),qa=list())
  add<-function(topic,df){
    updated<-tabs
    updated[[topic]][[length(updated[[topic]])+1L]]<-df
    tabs<<-updated
  }
  for(setname in names(obj$sets)){
    states<-obj$sets[[setname]];s0<-states[["S0"]];n<-nrow(s0)
    if(anyNA(s0$baseline_S0_group)||any(!s0$baseline_S0_group%in%groups))stop("S0 group invalid")
    if(length(unique(s0$pidp))!=n)stop("ID not unique")
    for(sc in scenarios){
      d<-states[[sc]];stopifnot(nrow(d)==n,identical(d$pidp,s0$pidp),identical(d$baseline_S0_group,s0$baseline_S0_group))
      for(year in c("2019","2025")){
        if(year=="2019"&&sc=="ST")next # no 2019 ST; never mislabel reference S0 as ST
        y<-if(year=="2019")"19" else "25"
        cat<-d[[paste0("risk",y)]];A<-d[[paste0("A",y)]];D<-d[[paste0("D",y)]]
        base<-s0[[paste0("risk",y)]];a0<-s0[[paste0("A",y)]]
        for(lvl in c("L","M","H","VH","H+VH")) {
          nn<-sum(cat%in%if(lvl=="H+VH")c("H","VH") else lvl)
          add("level",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,guideline=year,category=lvl,N=n,n=nn,pct_total=100*nn/n))
        }
        hv<-sum(cat%in%c("H","VH"));eligible<-sum(A==1L)
        add("overall",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,guideline=year,N=n,
          n_hvh=hv,pct_hvh=100*hv/n,n_eligible=eligible,pct_eligible=100*eligible/n,
          n_upgraded=sum(d[[paste0("upgrade",y)]]),
          n_additional_hvh=hv-sum(base%in%c("H","VH")),
          pp_additional_hvh=100*(hv-sum(base%in%c("H","VH")))/n,
          n_additional_eligible=eligible-sum(a0==1L),
          pp_additional_eligible=100*(eligible-sum(a0==1L))/n,
          mean_minimum_D=mean(D)))
      }
      # 2019 comparator for ST is S0; all tables below explicitly identify the comparison pair.
      pair19<-if(sc=="ST")"S0" else if(sc=="SNT")"SN" else sc
      add("year_diff",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,
        comparison_2019=pair19,comparison_2025=sc,N=n,
        n_hvh_2019=sum(d$hv19),n_hvh_2025=sum(d$hv25),
        difference_hvh_n=sum(d$hv25)-sum(d$hv19),
        difference_hvh_pp=100*(sum(d$hv25)-sum(d$hv19))/n,
        n_eligible_2019=sum(d$A19),n_eligible_2025=sum(d$A25),
        difference_eligible_n=sum(d$A25)-sum(d$A19),
        difference_eligible_pp=100*(sum(d$A25)-sum(d$A19))/n,
        difference_mean_minimum_D=mean(d$D25-d$D19)))
      for(g in groups){
        sel<-d$transition==g;ng<-sum(sel)
        add("group",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,
          comparison_2019=pair19,comparison_2025=sc,transition=g,N=n,n=ng,pct_total=100*ng/n,
          n_eligible_2019=sum(d$A19[sel]),n_eligible_2025=sum(d$A25[sel]),
          eligible_pct_total_2019=100*sum(d$A19[sel])/n,
          eligible_pct_total_2025=100*sum(d$A25[sel])/n,
          eligible_pct_group_2019=if(ng)100*mean(d$A19[sel])else NA_real_,
          eligible_pct_group_2025=if(ng)100*mean(d$A25[sel])else NA_real_,
          n_hvh_2019=sum(d$hv19[sel]),n_hvh_2025=sum(d$hv25[sel]),
          mean_D_2019=if(ng)mean(d$D19[sel])else NA_real_,
          mean_D_2025=if(ng)mean(d$D25[sel])else NA_real_,
          mean_D_difference=if(ng)mean(d$D25[sel]-d$D19[sel])else NA_real_))
      }
      for(a in orig_levels)for(b in orig_levels){
        sel<-d$origin19==a;cnt<-sum(sel & d$origin25==b)
        add("riskmat",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,
          comparison_2019=pair19,comparison_2025=sc,origin_2019=a,origin_2025=b,N=n,n=cnt,
          pct_total=100*cnt/n,pct_row=if(sum(sel))100*cnt/sum(sel)else NA_real_))
      }
      for(a in groups)for(b in groups){
        sel<-s0$transition==a;cnt<-sum(sel&d$transition==b)
        add("groupmat",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,
          score_based_group=a,RM_based_group=b,N=n,n=cnt,pct_total=100*cnt/n,
          pct_row=if(sum(sel))100*cnt/sum(sel)else NA_real_))
      }
    }
    # FIXED contribution: ALWAYS group on original S0 six groups, even under SNT.
    for(sc in c("S0","SNT")){
      d<-states[[sc]]
      for(g in c(groups,"Overall")){
        sel<-if(g=="Overall")rep(TRUE,n) else s0$transition==g
        ng<-sum(sel)
        nh19<-sum(d$hv19[sel]);nh25<-sum(d$hv25[sel])
        na19<-sum(d$A19[sel]);na25<-sum(d$A25[sel])
        sumd19<-sum(d$D19[sel]);sumd25<-sum(d$D25[sel])
        add("fixed",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,fixed_S0_group=g,N=n,n_fixed=ng,
          n_hvh_2019=nh19,n_hvh_2025=nh25,
          pct_hvh_total_2019=100*nh19/n,pct_hvh_total_2025=100*nh25/n,
          pct_hvh_group_2019=if(ng)100*nh19/ng else NA_real_,
          pct_hvh_group_2025=if(ng)100*nh25/ng else NA_real_,
          n_eligible_2019=na19,n_eligible_2025=na25,
          pct_eligible_total_2019=100*na19/n,pct_eligible_total_2025=100*na25/n,
          pct_eligible_group_2019=if(ng)100*na19/ng else NA_real_,
          pct_eligible_group_2025=if(ng)100*na25/ng else NA_real_,
          mean_D_group_2019=if(ng)sumd19/ng else NA_real_,
          mean_D_group_2025=if(ng)sumd25/ng else NA_real_,
          contribution_hvh_n=nh25-nh19,contribution_hvh_pp=100*(nh25-nh19)/n,
          contribution_eligible_n=na25-na19,contribution_eligible_pp=100*(na25-na19)/n,
          contribution_mean_D=(sumd25-sumd19)/n))
      }
    }
    f<-do.call(rbind,tabs$fixed)
    f<-f[f$RM_set==setname,,drop=FALSE]
    s0f<-f[f$strategy=="S0",,drop=FALSE];sntf<-f[f$strategy=="SNT",,drop=FALSE]
    ix<-match(s0f$fixed_S0_group,sntf$fixed_S0_group)
    if(anyNA(ix)||nrow(s0f)!=nrow(sntf)||any(s0f$n_fixed!=sntf$n_fixed[ix]))stop("Fixed groups changed")
    for(j in seq_len(nrow(s0f))){
      x<-s0f[j,];y<-sntf[ix[j],]
      add("change",data.frame(imputation_id=imp,RM_set=setname,fixed_S0_group=x$fixed_S0_group,
        N=n,n_fixed=x$n_fixed,
        S0_hvh_contribution_n=x$contribution_hvh_n,SNT_hvh_contribution_n=y$contribution_hvh_n,
        change_hvh_contribution_n=y$contribution_hvh_n-x$contribution_hvh_n,
        change_hvh_contribution_pp=y$contribution_hvh_pp-x$contribution_hvh_pp,
        change_hvh_2019_n=y$n_hvh_2019-x$n_hvh_2019,change_hvh_2025_n=y$n_hvh_2025-x$n_hvh_2025,
        change_hvh_2019_pct_group=if(x$n_fixed)100*(y$n_hvh_2019-x$n_hvh_2019)/x$n_fixed else NA_real_,
        change_hvh_2025_pct_group=if(x$n_fixed)100*(y$n_hvh_2025-x$n_hvh_2025)/x$n_fixed else NA_real_,
        S0_eligible_contribution_n=x$contribution_eligible_n,SNT_eligible_contribution_n=y$contribution_eligible_n,
        change_eligible_contribution_n=y$contribution_eligible_n-x$contribution_eligible_n,
        change_eligible_contribution_pp=y$contribution_eligible_pp-x$contribution_eligible_pp,
        change_eligible_2019_n=y$n_eligible_2019-x$n_eligible_2019,
        change_eligible_2025_n=y$n_eligible_2025-x$n_eligible_2025,
        change_eligible_2019_pct_group=if(x$n_fixed)100*(y$n_eligible_2019-x$n_eligible_2019)/x$n_fixed else NA_real_,
        change_eligible_2025_pct_group=if(x$n_fixed)100*(y$n_eligible_2025-x$n_eligible_2025)/x$n_fixed else NA_real_,
        S0_mean_D_contribution=x$contribution_mean_D,SNT_mean_D_contribution=y$contribution_mean_D,
        change_mean_D_contribution=y$contribution_mean_D-x$contribution_mean_D))
    }
    # Additivity within each imputation: group contributions exactly reconstruct the overall contrast.
    for(sc in c("S0","SNT")) {
      u<-f[f$strategy==sc & f$fixed_S0_group!="Overall",,drop=FALSE]
      v<-f[f$strategy==sc & f$fixed_S0_group=="Overall",,drop=FALSE]
      err<-max(abs(c(sum(u$n_fixed)-n,
        sum(u$contribution_hvh_n)-v$contribution_hvh_n,
        sum(u$contribution_eligible_n)-v$contribution_eligible_n,
        sum(u$contribution_mean_D)-v$contribution_mean_D)))
      if(nrow(u)!=6L||nrow(v)!=1L||err>1e-8)stop("Fixed contribution additivity failed: ",setname," ",sc)
      add("qa",data.frame(imputation_id=imp,RM_set=setname,strategy=sc,N=n,
        n_risk_matrix=sum(vapply(tabs$riskmat,function(t)if(t$RM_set==setname&&t$strategy==sc)t$n else 0,numeric(1))),
        n_group_matrix=sum(vapply(tabs$groupmat,function(t)if(t$RM_set==setname&&t$strategy==sc)t$n else 0,numeric(1))),
        n_transition_groups=sum(vapply(tabs$group,function(t)if(t$RM_set==setname&&t$strategy==sc)t$n else 0,numeric(1))),
        contribution_max_abs_error=err))
    }
    for(sc in scenarios){
      sub<-do.call(rbind,tabs$group);sub<-sub[sub$RM_set==setname & sub$strategy==sc,,drop=FALSE]
      rr<-do.call(rbind,tabs$riskmat);rr<-rr[rr$RM_set==setname & rr$strategy==sc,,drop=FALSE]
      mm<-do.call(rbind,tabs$groupmat);mm<-mm[mm$RM_set==setname & mm$strategy==sc,,drop=FALSE]
      if(nrow(sub)!=6L||sum(sub$n)!=n||nrow(rr)!=25L||sum(rr$n)!=n||nrow(mm)!=36L||sum(mm$n)!=n)stop("Matrix/transition totals fail: ",setname," ",sc)
    }
  }
  result<-lapply(tabs,function(k)if(length(k))do.call(rbind,k) else data.frame())
  saveRDS(result,file.path(outdir,sprintf("MI_%03d_tables.rds",imp)),compress="xz")
  data.frame(imputation_id=imp,sets=length(obj$sets),N=nrow(obj$sets[[1]]$S0),pass=TRUE)
}

pool_table<-function(d,keys) {
  if(!nrow(d))return(d)
  values<-names(d)[vapply(d,is.numeric,logical(1))]
  values<-setdiff(values,c(keys,"imputation_id"))
  splits<-split(seq_len(nrow(d)),interaction(d[keys],drop=TRUE,lex.order=TRUE))
  do.call(rbind,lapply(splits,function(idx){
    x<-d[idx,,drop=FALSE];nums<-lapply(x[values],function(v)if(all(is.na(v)))NA_real_ else mean(v,na.rm=TRUE))
    data.frame(x[1L,keys,drop=FALSE],nums,m=nrow(x),m_valid_min=if(length(values))min(vapply(x[values],function(v)sum(!is.na(v)),integer(1))) else NA_integer_,row.names=NULL,check.names=FALSE)
  }))
}
manifest_path<-file.path(out,"quality_control","08_completed_chains.csv")
if(!file.exists(manifest_path))stop("Run step 08 first: ",manifest_path)
signature_path<-file.path(out,"quality_control","08_input_signatures.csv")
if(!file.exists(signature_path))stop("Missing step 08 input signatures; rerun step 08")
signatures<-read.csv(signature_path,fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
if(any(!file.exists(signatures$path)))stop("Step 08 original input no longer exists")
current_hash<-unname(tools::md5sum(signatures$path))
if(any(is.na(current_hash)|current_hash!=signatures$md5))
  stop("The step-05 or step-07 input has changed since 08. Rerun step 08 to avoid stale cross-step results.")
manifest<-read.csv(manifest_path,fileEncoding="UTF-8-BOM",check.names=FALSE)
if(!nrow(manifest)||anyDuplicated(manifest$imputation_id)||!all(manifest$completed))stop("Invalid step 08 manifest")
tasks<-lapply(manifest$imputation_id,function(i)list(id=i,path=file.path(out,"01_individual_states",sprintf("MI_%03d.rds",i))))
missing<-vapply(tasks,function(x)!file.exists(x$path),logical(1));if(any(missing))stop("Missing step 08 state(s): ",paste(manifest$imputation_id[missing],collapse=","))
log("09 started; ",length(tasks)," MI chains; ",THREADS," workers")
workers<-min(THREADS,length(tasks));res<-list()
if(workers>1L){cl<-parallel::makePSOCKcluster(workers);parallel::clusterExport(cl,c("table_one","groups","orig_levels","scenarios"),envir=.GlobalEnv)}
tryCatch({
for(ix in split(seq_along(tasks),ceiling(seq_along(tasks)/workers))){
  z<-if(workers>1L)parallel::parLapplyLB(cl,tasks[ix],table_one,outdir=file.path(out,"02_per_imputation")) else lapply(tasks[ix],table_one,outdir=file.path(out,"02_per_imputation"))
  res<-c(res,z);log("Per-MI tables: ",max(ix),"/",length(tasks)," (",round(100*max(ix)/length(tasks)),"%)")
}
}, finally=if(workers>1L)parallel::stopCluster(cl))
write_bom(do.call(rbind,res),file.path(out,"quality_control","09_computation_progress.csv"))
files<-file.path(out,"02_per_imputation",sprintf("MI_%03d_tables.rds",manifest$imputation_id))
topic_names<-c(level="01_risk_levels",overall="02_hvh_eligibility",year_diff="02b_2019_2025_differences",
               group="03_six_transitions",riskmat="04_risk_matrix",groupmat="05_six_group_matrix",
               fixed="06_fixed_contributions",change="07_contribution_change",qa="08_additivity_QA")
keys<-list(level=c("RM_set","strategy","guideline","category"),overall=c("RM_set","strategy","guideline"),
           year_diff=c("RM_set","strategy","comparison_2019","comparison_2025"),
           group=c("RM_set","strategy","comparison_2019","comparison_2025","transition"),
           riskmat=c("RM_set","strategy","comparison_2019","comparison_2025","origin_2019","origin_2025"),
           groupmat=c("RM_set","strategy","score_based_group","RM_based_group"),
           fixed=c("RM_set","strategy","fixed_S0_group"),change=c("RM_set","fixed_S0_group"),
           qa=c("RM_set","strategy"))
for(k in names(topic_names)){
  rows<-lapply(files,function(p)readRDS(p)[[k]])
  d<-do.call(rbind,rows)
  if(anyDuplicated(d[c("imputation_id",keys[[k]])]))stop("Duplicate keys in ",k)
  if(k=="qa"&&any(d$n_risk_matrix!=d$N|d$n_group_matrix!=d$N|d$n_transition_groups!=d$N|d$contribution_max_abs_error>1e-8))stop("QA failure")
  write_bom(d,file.path(out,"02_per_imputation",paste0(topic_names[[k]],".csv")))
  if(k!="qa"){
    pooled<-pool_table(d,keys[[k]])
    # S38-S45 reuse these results keyed by RM_set. Verify source coverage.
    expected <- c(groupmat = 144L, group = 24L, overall = 7L)
    if (k %in% names(expected)) {
      coverage <- table(pooled$RM_set)
      if (any(coverage != expected[[k]]))
        stop("Incomplete supplementary source grid for ", k)
      if (k == "groupmat") {
        totals <- aggregate(pct_total ~ RM_set + strategy, pooled, sum)
        if (any(abs(totals$pct_total - 100) > 1e-8))
          stop("Supplementary matrix percentages must sum to 100")
      }
    }
    # The MI mean of counts is an expected count, potentially fractional; do not label as observed integer N.
    write_bom(pooled,file.path(out,"03_pooled",paste0(topic_names[[k]],"_MI_mean.csv")))
    # Auditable across-imputation variability; not inferential confidence limits.
    fields <- setdiff(names(d)[vapply(d,is.numeric,logical(1))],c("imputation_id",keys[[k]]))
    idx <- split(seq_len(nrow(d)),interaction(d[keys[[k]]],drop=TRUE,lex.order=TRUE))
    spreads <- do.call(rbind,lapply(idx,function(ii){
      x <- d[ii,,drop=FALSE]
      v <- lapply(x[fields],function(a)c(
        sd=if(sum(!is.na(a))>1L)stats::sd(a,na.rm=TRUE) else NA_real_,
        min=if(all(is.na(a)))NA_real_ else min(a,na.rm=TRUE),
        max=if(all(is.na(a)))NA_real_ else max(a,na.rm=TRUE)))
      data.frame(x[1L,keys[[k]],drop=FALSE],as.list(unlist(v)),check.names=FALSE)
    }))
    write_bom(spreads,file.path(out,"quality_control",paste0("09_",topic_names[[k]],"_MI_spread.csv")))
  }
  log("Wrote ",k," | ",nrow(d)," per-MI rows")
}

# S32-S35 use the primary analytical cohort and its saved RM imputations.
if (TEST_N == 0L && SMOKING_MODE == "same_wave" && all(manifest$sets == 2L)) local({
  elsa <- DATASET_ID == "S04_ELSA"
  code <- if (elsa) "S04" else "S03"
  key <- if (elsa) "idauniq" else "pidp"
  clean <- file.path(root, "LJW_outcome", "data_clean", DATASET_ID)
  intervention <- dirname(out)
  mi_dir <- file.path(intervention, "rm_imputation")
  read_csv <- function(path) read.csv(path, fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
  states <- lapply(tasks, function(task) readRDS(task$path))
  base <- states[[1]]$sets$expanded$S0
  ids <- as.character(base[[key]])
  align <- function(x) {
    ix <- match(ids, as.character(x[[key]]))
    stopifnot(!anyNA(ix), !anyDuplicated(x[[key]]))
    x[ix, , drop = FALSE]
  }
  raw <- align(read_csv(file.path(clean, paste0("04_eligible_cohort_", code, ".csv"))))
  mi <- readRDS(file.path(mi_dir, if (elsa) "04_RM_imputed_cohorts_S04.rds" else "07_RM_imputed_cohorts_S03.rds"))
  stopifnot(setequal(seq_along(mi), manifest$imputation_id))
  mi <- lapply(mi, align)
  mapping <- read_csv(file.path(out, "quality_control",
                               if (elsa) "05_rm_mapping.csv" else "08_rm_mapping.csv"))
  summarize <- function(x, binary, applicable = rep(TRUE, length(x))) {
    observed <- x[applicable & !is.na(x)]
    missing <- sum(applicable & is.na(x))
    if (binary) stopifnot(all(as.character(observed) %in% c("0", "1")))
    q <- if (!binary && length(observed)) quantile(observed, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3)
    positive <- if (binary) sum(as.character(observed) == "1") else NA_real_
    data.frame(n_total = length(x), n_applicable = sum(applicable),
               n_observed = length(observed), n_missing = missing,
               missing_pct = 100 * missing / length(x), n_positive = positive,
               positive_pct = if (length(observed)) 100 * positive / length(observed) else NA_real_,
               q25 = q[1], median = q[2], q75 = q[3])
  }
  by_variable <- function(x, variables, binary) {
    stopifnot(all(variables %in% names(x)))
    do.call(rbind, lapply(variables, function(v) {
      applicable <- if (v == "rm_natural_menopause_lt40_proxy") x$sex == "female" else rep(TRUE, nrow(x))
      cbind(variable = v, summarize(x[[v]], v %in% binary, applicable))
    }))
  }
  if (anyNA(raw$sex) || any(!raw$sex %in% c(1L, 2L))) {
    stop("UKHLS analytical cohort: expected sex codes 1 (male) and 2 (female)")
  }
  raw$female <- as.integer(raw$sex == 2L)
  raw$moderate_2019 <- as.integer(base$risk19 == "M")
  raw$moderate_2025 <- as.integer(base$risk25 == "M")
  binary_base <- c("female", "current_smoking", "moderate_2019", "moderate_2025")
  baseline <- by_variable(raw, c("age", "female", "current_smoking", "sbp", "dbp",
                                "tc", "ldl", "hdl", if (elsa) "trig" else "tg",
                                "moderate_2019", "moderate_2025"), binary_base)
  variables <- c("bmi", if (elsa) "waist" else "waist_cm",
                 if (elsa) "hscrp_w6" else "hscrp", mapping$column)
  after <- do.call(rbind, lapply(seq_along(mi), function(i) {
    cbind(imputation_id = i, by_variable(mi[[i]], variables, mapping$column))
  }))
  average <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  pooled <- aggregate(after[setdiff(names(after), c("variable", "imputation_id"))],
                      after["variable"], average)
  if (elsa) {
    before <- read_csv(file.path(mi_dir, "04_variable_distribution_missingness_before_S04.csv"))
    before <- before[match(variables, before$variable), ]
    stopifnot(identical(before$variable, variables), all(before$n_total == length(ids)))
    before$missing_pct <- 100 * before$n_missing / before$n_total
    before$positive_pct <- before$positive_pct_observed
  } else {
    before <- by_variable(raw, variables, mapping$column)
  }
  distribution <- rbind(cbind(stage = "before", m = 0L, before[names(pooled)]),
                        cbind(stage = "after", m = length(mi), pooled))
  prevalence <- do.call(rbind, lapply(states, function(state) {
    i <- state$imputation_id
    do.call(rbind, lapply(c("strict", "expanded"), function(set) {
      cases <- state$sets[[set]]
      do.call(rbind, lapply(c("2019", "2025"), function(year) {
        suffix <- substr(year, 3L, 4L)
        moderate <- align(cases$S0)[[paste0("risk", suffix)]] == "M"
        stopifnot(sum(moderate) > 0)
        keep <- mapping$guideline %in% c("2019_and_2025", paste0(year, "_only"))
        if (set == "strict") keep <- keep & mapping$feasibility == "checkmark"
        md <- mapping[keep, ]
        values <- setNames(lapply(md$column, function(v) {
          x <- as.character(mi[[i]][[v]])
          if (v == "rm_natural_menopause_lt40_proxy") x[mi[[i]]$sex == "male"] <- "0"
          x <- x[moderate]
          stopifnot(all(x %in% c("0", "1")))
          x == "1"
        }), md$domain)
        strategies <- if (year == "2019") c("SN", "SNT") else c("SN", "ST", "SNT")
        for (strategy in strategies) {
          domain <- switch(strategy, SN = md$domain[md$domain != "hscrp"],
                           ST = "hscrp", SNT = md$domain)
          positive <- Reduce(function(a, b) a | b, values[domain])
          upgraded <- align(cases[[strategy]])[[paste0("upgrade", suffix)]][moderate]
          stopifnot(length(positive) == length(upgraded), all(positive == upgraded))
          if (strategy != "SNT" || set == "expanded") values[[strategy]] <- positive
        }
        n <- vapply(values, sum, integer(1))
        data.frame(imputation_id = i, RM_set = set, guideline = year, modifier = names(n),
                   n = unname(n), n_moderate = sum(moderate), pct = 100 * unname(n) / sum(moderate))
      }))
    }))
  }))
  pooled_prevalence <- aggregate(prevalence[c("n", "n_moderate", "pct")],
                                 prevalence[c("RM_set", "guideline", "modifier")], mean)
  pooled_prevalence$m <- length(mi)
  flow_path <- if (elsa) file.path(intervention, "04_ELSA_eligibility_flow_S03.csv") else
    Sys.glob(file.path(clean, "*", "04_eligibility_flow_S03.csv"))
  stopifnot(length(flow_path) == 1L)
  flow <- read_csv(flow_path)
  stopifnot(all(flow$n_before - flow$n_excluded == flow$n_remaining),
            identical(head(flow$n_remaining, -1L), tail(flow$n_before, -1L)),
            tail(flow$n_remaining, 1L) == length(ids))
  flow <- cbind(source_file = basename(flow_path), flow)
  tables <- list(S32_selection_flow = flow, S33_baseline = baseline,
                 S34_imputation = distribution, S34_per_imputation = after,
                 S35_RM_prevalence = pooled_prevalence, S35_per_imputation = prevalence)
  for (name in names(tables)) {
    layer <- if (endsWith(name, "_per_imputation")) "02_per_imputation" else "03_pooled"
    write_bom(cbind(cohort = DATASET_ID, tables[[name]]), file.path(out, layer, paste0(name, ".csv")))
  }
  log("Wrote supplementary sources S32-S35")
})

# Main Table 2: guideline-specific results and cross-guideline category changes.
if (TEST_N == 0L && SMOKING_MODE == "same_wave" && all(manifest$sets == 2L)) local({
  read_result <- function(stem) {
    read.csv(file.path(out,"03_pooled",paste0(stem,"_MI_mean.csv")),
             fileEncoding="UTF-8-BOM",check.names=FALSE)
  }
  overall <- read_result("02_hvh_eligibility")
  transitions <- read_result("03_six_transitions")
  annual <- list()
  paired <- list()
  for (set in c("strict","expanded")) {
    baseline <- subset(overall, RM_set==set & strategy=="S0")
    rm <- subset(overall, RM_set==set & strategy=="SNT")
    stopifnot(nrow(baseline)==2L,nrow(rm)==2L,!anyDuplicated(rm$guideline))
    baseline <- baseline[match(rm$guideline,baseline$guideline),]
    stopifnot(all(rm$N==baseline$N))
    rm$delta_mean_D_vs_S0 <- rm$mean_minimum_D-baseline$mean_minimum_D
    annual[[set]] <- rbind(transform(baseline,delta_mean_D_vs_S0=0),rm)

    base_groups <- subset(transitions,RM_set==set & strategy=="S0")
    rm_groups <- subset(transitions,RM_set==set & strategy=="SNT")
    stopifnot(nrow(base_groups)==6L,nrow(rm_groups)==6L,!anyDuplicated(rm_groups$transition))
    base_groups <- base_groups[match(rm_groups$transition,base_groups$transition),]
    stopifnot(all(rm_groups$N==base_groups$N))
    paired[[set]] <- data.frame(RM_set=set,transition=rm_groups$transition,N=rm_groups$N,
      n_S0=base_groups$n,pct_S0=base_groups$pct_total,
      delta_n=rm_groups$n-base_groups$n,delta_pp=rm_groups$pct_total-base_groups$pct_total)
    stopifnot(abs(sum(paired[[set]]$delta_n))<1e-7, abs(sum(paired[[set]]$delta_pp))<1e-7)
  }
  write_bom(do.call(rbind,annual),file.path(out,"03_pooled","Table_2_guideline_results.csv"))
  write_bom(do.call(rbind,paired),file.path(out,"03_pooled","Table_2_transition_net_changes.csv"))
})
log("09 completed | descriptives only; no risk models, CIs, or Rubin inferential intervals")






