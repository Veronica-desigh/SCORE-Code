# 08 | External-cohort RM upgrade strategies: strict (checkmark) and expanded (checkmark + triangle)
# Standalone script; reads EXISTING completed RM imputation; does NOT refit MI or rescore SCORE/SCORE2.
# Reuses and validates the frozen M01/S03 treatment rules (M->H only; LDL-C targets; six transitions).
# ============================= USER CONFIG =============================
DATASET_ID <- "S03_UKHLS"   # Fixed UKHLS only
THREADS <- 5L
TEST_N <- 0L
SMOKING_MODE <- "same_wave"  # S03 only: same_wave / wave2_prior
PROJECT_ROOT <- ""           # Empty: search parents of script directory
SCORE_PATH <- ""             # Optional: full path to 05 individual SCORE/SCORE2 results
MI_PATH <- ""                # Optional: full path to 07 completed RM imputation .rds
RM_MAP_CSV <- ""             # Optional CSV with domain,column; required to remap ELSA names
RUN_SETS <- c("strict", "expanded")
# =======================================================================
DATASET_ID <- "S03_UKHLS"
THREADS <- as.integer(Sys.getenv("S03_RM_THREADS", Sys.getenv("S03_THREADS", as.character(THREADS))))
TEST_N <- as.integer(Sys.getenv("S03_RM_TEST_N", Sys.getenv("S03_ANALYSIS_TEST_N", as.character(TEST_N))))
SMOKING_MODE <- Sys.getenv("S03_SMOKING_MODE", SMOKING_MODE)
PROJECT_ROOT <- Sys.getenv("S03_RM_PROJECT_ROOT", Sys.getenv("S03_PROJECT_ROOT", PROJECT_ROOT))
SCORE_PATH <- Sys.getenv("S03_RM_SCORE_PATH", SCORE_PATH)
MI_PATH <- Sys.getenv("S03_RM_MI_PATH", MI_PATH)
RM_MAP_CSV <- Sys.getenv("S03_RM_MAP", RM_MAP_CSV)
if (!identical(DATASET_ID,"S03_UKHLS")) stop("Invalid S03 cohort configuration")
if (!is.finite(THREADS) || THREADS < 1L || !is.finite(TEST_N) || TEST_N < 0L) stop("Invalid THREADS or TEST_N")
if (!all(RUN_SETS %in% c("strict", "expanded")) || anyDuplicated(RUN_SETS)) stop("Invalid RUN_SETS")
if (!SMOKING_MODE %in% c("same_wave", "wave2_prior")) stop("Invalid SMOKING_MODE")
args_file <- grep("^--file=", commandArgs(), value=TRUE)
frames <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
start <- if(length(frames)) dirname(tail(frames,1L)[[1L]]) else if(length(args_file)) dirname(sub("^--file=", "",args_file[1L])) else getwd()
root <- if(nzchar(PROJECT_ROOT)) normalizePath(PROJECT_ROOT, winslash="/", mustWork=TRUE) else normalizePath(start,winslash="/",mustWork=TRUE)
if (!nzchar(PROJECT_ROOT)) while(!dir.exists(file.path(root,"LJW_outcome"))) {
  up <- dirname(root); if(identical(up,root)) stop("Project root not found; specify PROJECT_ROOT")
  root <- up
}
base_out <- file.path(root,"LJW_outcome","data_intervention",DATASET_ID)
if(TEST_N>0L) base_out <- file.path(base_out,paste0("test_",TEST_N))
if(SMOKING_MODE=="wave2_prior") base_out <- file.path(base_out,"wave2_smoking_sensitivity")
prefix <- "S03"
if(!nzchar(SCORE_PATH)) SCORE_PATH <- file.path(base_out,"complete_case",paste0("05_score_strategy_participants_",prefix,".csv"))
if(!nzchar(MI_PATH)) MI_PATH <- file.path(base_out,"rm_imputation",paste0("07_RM_imputed_cohorts_",prefix,".rds"))
output <- file.path(base_out,"external_RM_analyses")
for(x in c("01_individual_states","quality_control","logs")) dir.create(file.path(output,x),recursive=TRUE,showWarnings=FALSE)
start_time <- Sys.time(); log_file <- file.path(output,"logs","08_RM_strategies.log")
log <- function(...) {msg <- sprintf("[%s | +%.1fs] %s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),as.numeric(difftime(Sys.time(),start_time,units="secs")),paste0(...,collapse="")); message(msg); cat(msg,"\n",file=log_file,append=TRUE)}
write_bom <- function(x,path) {tmp<-tempfile(fileext=".csv");on.exit(unlink(tmp),add=TRUE);utils::write.csv(x,tmp,row.names=FALSE,na="",fileEncoding="UTF-8");bytes<-readBin(tmp,"raw",n=file.info(tmp)$size);con<-file(path,"wb");on.exit(close(con),add=TRUE);writeBin(as.raw(c(239,187,191)),con);writeBin(bytes,con);invisible(path)}

# Semantic domains, not arbitrary approximations. ELSA mapping must be verified against its codebook.
maps <- c(
  obesity="rm_obesity_central", hscrp="rm_hscrp_gt2",
  deprivation="rm_social_deprivation", inactivity="rm_physical_inactivity_sport_proxy",
  arthritis="rm_chronic_inflammatory_arthritis_proxy", south_asian="rm_south_asian"
)
if(nzchar(RM_MAP_CSV)) {
  if(!file.exists(RM_MAP_CSV)) stop("RM_MAP_CSV does not exist: ",RM_MAP_CSV)
  mapping<-read.csv(RM_MAP_CSV, fileEncoding="UTF-8-BOM", stringsAsFactors=FALSE)
  if(!all(c("domain","column")%in%names(mapping))||anyDuplicated(mapping$domain))stop("RM map needs unique domain,column")
  if(!setequal(names(maps),mapping$domain))stop("RM map must define exactly: ",paste(names(maps),collapse=", "))
  maps<-setNames(mapping$column[match(names(maps),mapping$domain)],names(maps))
}
if(anyNA(maps)||any(!nzchar(maps))||anyDuplicated(unname(maps)))stop("Invalid RM column map")
nonlab19 <- c("obesity","deprivation","inactivity","arthritis")
nonlab25 <- c(nonlab19,"south_asian")
sets <- list(strict=list(`2019`=c("obesity"),`2025`=c("obesity")),
             expanded=list(`2019`=nonlab19,`2025`=nonlab25))
sets <- sets[RUN_SETS]

# Three-valued union: known positive dominates; fully observed negatives => 0; otherwise unresolved NA.
rm_binary <- function(z,var) {
  if(!var%in%names(z))stop("Missing RM column: ",var)
  a<-z[[var]]
  if(is.factor(a)) a<-as.character(a)
  if(is.logical(a)) a<-as.integer(a)
  a<-suppressWarnings(as.integer(a))
  if(any(!is.na(a)&!a%in%0:1))stop("RM column is not binary: ",var)
  a
}
rm_union <- function(z, domains, mapping) {
  if(!length(domains))return(rep(0L,nrow(z)))
  v<-lapply(domains,function(nm)rm_binary(z,unname(mapping[[nm]])))
  mat<-do.call(cbind,v)
  positive<-rowSums(mat==1L,na.rm=TRUE)>0L
  negative<-rowSums(mat==0L,na.rm=TRUE)==ncol(mat)
  result<-rep(NA_integer_,nrow(z));result[positive]<-1L;result[negative]<-0L
  result
}
# 2019 ST is not applicable. For cross-guideline ST matrices use 2019 S0 as EXPLICIT comparator.
year_domains <- function(setname,scenario,year,sets) {
  if(scenario=="S0")return(character())
  if(scenario=="ST")return(if(year=="2025")"hscrp" else character())
  if(scenario=="SN")return(sets[[setname]][[year]])
  if(scenario=="SNT")return(if(year=="2019") sets[[setname]][[year]] else unique(c(sets[[setname]][[year]],"hscrp")))
  stop("Unknown scenario")
}
treat <- function(ldl,level) {
  target <- unname(c(L=3.0,M=2.6,H=1.8,VH=1.4)[as.character(level)])
  if(anyNA(target)||anyNA(ldl)||any(!is.finite(ldl))||any(ldl<=0))stop("Invalid LDL/level")
  upper <- ifelse(level%in%c("H","VH"),pmin(target,ldl/2),target)
  eligible <- as.integer(ldl>target)
  reduction <- ifelse(eligible==1L,ldl-upper,0)
  list(target=target,A=eligible,D=reduction)
}
transition <- function(a19,a25,d19,d25) {
  x<-rep(NA_character_,length(a19));delta<-d25-d19;tol<-1e-10
  x[a19==0&a25==0]<-"unchanged_untreated"
  x[a19==0&a25==1]<-"newly_intervention"
  x[a19==1&a25==0]<-"intervention_removed"
  x[a19==1&a25==1&delta>tol]<-"treatment_intensified"
  x[a19==1&a25==1&abs(delta)<=tol]<-"unchanged_intensity"
  x[a19==1&a25==1&delta< -tol]<-"treatment_deintensified"
  if(anyNA(x))stop("Six-transition partition not exhaustive")
  x
}
origins <- function(b,u) {
  x<-paste0(as.character(b),"_score")
  x[u]<-"H_rm"
  if(!all(x%in%c("L_score","M_score","H_rm","H_score","VH_score")))stop("Invalid origin levels")
  x
}

# Small deterministic semantic self-tests; no study data or package dependencies.
unit <- data.frame(a=c(1L,0L,NA_integer_,0L),b=c(NA_integer_,0L,0L,NA_integer_))
if(!identical(rm_union(unit,c("a","b"),c(a="a",b="b")),c(1L,0L,NA_integer_,NA_integer_)))
  stop("RM union self-test failed: unknown must not become negative")
selftest_tr <- transition(c(0,0,1,1,1,1),c(0,1,0,1,1,1),
                          c(0,0,1.0,1.0,1.0,2.0),c(0,1.0,0,1.4,1.0,1.3))
if(!identical(selftest_tr,c("unchanged_untreated","newly_intervention","intervention_removed",
                            "treatment_intensified","unchanged_intensity","treatment_deintensified")))
  stop("Six-group self-test failed")
if(abs(treat(2,"H")$D-1.0)>1e-10||treat(1.7,"H")$A!=0L)stop("H-group LDL treatment self-test failed")

if(!file.exists(SCORE_PATH))stop("Missing 05 SCORE file: ",SCORE_PATH)
if(!file.exists(MI_PATH))stop("Missing existing 07 RM imputations: ",MI_PATH)
score<-read.csv(SCORE_PATH,fileEncoding="UTF-8-BOM",na.strings="",check.names=FALSE,stringsAsFactors=FALSE)
needed<-c("pidp","ldl","score_2019_category","score2_2025_category","trigger_2019","trigger_2025",
          "minimum_reduction_2019","minimum_reduction_2025","transition")
if(length(setdiff(needed,names(score))))stop("Missing 05 columns: ",paste(setdiff(needed,names(score)),collapse=", "))
if(!nrow(score)||anyNA(score[needed])||anyDuplicated(as.character(score$pidp)))stop("Invalid 05 score dataset")
score$pidp<-as.character(score$pidp)
if(!all(as.character(score$score_2019_category)%in%c("L","M","H","VH"))||
   !all(as.character(score$score2_2025_category)%in%c("L","M","H","VH"))) stop("05 risk categories invalid")
a0<-treat(score$ldl,as.character(score$score_2019_category)); b0<-treat(score$ldl,as.character(score$score2_2025_category))
if(any(a0$A!=as.integer(score$trigger_2019))||any(b0$A!=as.integer(score$trigger_2025))||
   any(abs(a0$D-score$minimum_reduction_2019)>1e-8)||any(abs(b0$D-score$minimum_reduction_2025)>1e-8)||
   any(transition(a0$A,b0$A,a0$D,b0$D)!=score$transition))stop("Frozen S0 treatment rules disagree with 05: halt to prevent drift")
mi<-readRDS(MI_PATH)
if(!is.list(mi)||!length(mi)||!all(vapply(mi,is.data.frame,logical(1)))) stop("MI input must be a list of completed data frames")
required_rm<-unname(maps[unique(c(unlist(lapply(sets,function(q)c(q[["2019"]],q[["2025"]]))),"hscrp"))])
if(anyNA(required_rm))stop("RM domain map unresolved")
# Explicit preflight across imputation chains. Unresolved M-class status fails; never reinterpret missing as RM-negative.
# Non-RM essential baseline covariates must be identical across MI chains (eligible variables were NOT imputed).
fixed_covariates <- c("age","sex","wave","current_smoking","tc","hdl","ldl","sbp","dbp","egfr")
fixed_covariates <- intersect(fixed_covariates,names(mi[[1]]))
ref_ids <- match(score$pidp,as.character(mi[[1]]$pidp))
if(anyNA(ref_ids))stop("MI1 cannot be aligned to the SCORE cohort")
reference_fixed <- mi[[1]][ref_ids,fixed_covariates,drop=FALSE]
coverage<-list(); tasks<-vector("list",length(mi));k<-1L
for(i in seq_along(mi)) {
  z<-mi[[i]]
  if(!"pidp"%in%names(z)||anyDuplicated(as.character(z$pidp)))stop("MI ",i," has duplicate/missing ID")
  if(anyNA(z$pidp)||length(setdiff(required_rm,names(z))))stop("MI ",i," RM fields missing: ",paste(setdiff(required_rm,names(z)),collapse=", "))
  loc<-match(score$pidp,as.character(z$pidp))
  if(length(loc)!=nrow(score)||anyNA(loc)||nrow(z)!=nrow(score))stop("MI ",i," does not match fixed 05 cohort exactly")
  z<-z[loc,,drop=FALSE]
  for(v in fixed_covariates){
    if(!v%in%names(z))stop("MI fixed covariate missing: ",v)
    if(!identical(as.character(z[[v]]),as.character(reference_fixed[[v]])))stop("MI changed eligible/score covariate: ",v," in chain ",i)
  }
  if("ldl"%in%names(z) && (anyNA(z$ldl)||any(abs(z$ldl-score$ldl)>1e-8))) stop("MI changed or lost observed LDL-C")
  for(domain in unique(names(maps))) rm_binary(z,unname(maps[[domain]]))
  for(setname in names(sets))for(sc in c("S0","SN","ST","SNT"))for(year in c("2019","2025")) {
    if(sc=="ST"&&year=="2019")next
    base<-as.character(score[[if(year=="2019")"score_2019_category" else "score2_2025_category"]])
    result<-rm_union(z,year_domains(setname,sc,year,sets),maps)
    n_undetermined<-sum(base=="M" & is.na(result))
    coverage[[k]]<-data.frame(imputation_id=i,set=setname,strategy=sc,guideline=year,
      n_initial_m=sum(base=="M"),n_unresolved_m=n_undetermined,n_unresolved_all=sum(is.na(result)),
      n_upgraded_known=sum(base=="M" & result==1L,na.rm=TRUE),stringsAsFactors=FALSE);k<-k+1L
  }
  tasks[[i]]<-list(id=i,data=z)
  if(i%%5L==0L||i==length(mi))log("Preflight: ",i,"/",length(mi)," imputations checked")
}
coverage<-do.call(rbind,coverage)
write_bom(coverage,file.path(output,"quality_control","08_rm_identifiability_by_imputation.csv"))
feasibility <- ifelse(names(maps)%in%c("obesity","hscrp"),"checkmark","triangle_proxy")
guideline <- ifelse(names(maps)%in%c("south_asian","hscrp","menopause"),"2025_only","2019_and_2025")
write_bom(data.frame(dataset=DATASET_ID,domain=names(maps),column=unname(maps),
  feasibility=feasibility,guideline=guideline),file.path(output,"quality_control","08_rm_mapping.csv"))
domain_manifest_rows<-list();manifest_k<-1L
for(setname in names(sets)) for(sc in c("S0","SN","ST","SNT")) for(y in c("2019","2025")) {
  applicable<-!(sc=="ST"&&y=="2019")
  domains<-if(applicable)paste(year_domains(setname,sc,y,sets),collapse="|")else ""
  domain_manifest_rows[[manifest_k]]<-data.frame(RM_set=setname,strategy=sc,guideline=y,
    applicable=applicable,RM_domains=domains,
    comparison_2019_for_ST=if(sc=="ST")"S0"else "",stringsAsFactors=FALSE)
  manifest_k<-manifest_k+1L
}
domain_manifest<-do.call(rbind,domain_manifest_rows)
write_bom(domain_manifest,file.path(output,"quality_control","08_domain_scenario_manifest.csv"))
write_bom(data.frame(input=c("score","RM_imputed"),
  path=normalizePath(c(SCORE_PATH,MI_PATH),winslash="/",mustWork=TRUE),
  md5=unname(tools::md5sum(c(SCORE_PATH,MI_PATH)))),
  file.path(output,"quality_control","08_input_signatures.csv"))
# Validate the strict and expanded domains independently. Incomplete expanded RM
# does not invalidate a strict calculation; never convert unknown RM to 0.
set_status <- do.call(rbind,lapply(names(sets), function(nm) {
  tab <- coverage[coverage$set==nm,,drop=FALSE]
  data.frame(RM_set=nm,imputations=length(mi),
    n_unresolved_scenario_cells=sum(tab$n_unresolved_m>0L),
    max_unresolved_m=max(tab$n_unresolved_m),runnable=!any(tab$n_unresolved_m>0L))
}))
write_bom(set_status,file.path(output,"quality_control","08_RM_set_run_status.csv"))
skipped <- set_status$RM_set[!set_status$runnable]
if(length(skipped))log("RM sets not identifiable: ",paste(skipped,collapse=", "),
 " (see 08_RM_set_run_status.csv; missingness remains unknown)")
sets <- sets[set_status$RM_set[set_status$runnable]]
if(!length(sets))stop("No identifiable RM sets; inspect 08_RM_set_run_status.csv")
write_bom(data.frame(RM_set=names(sets),runnable=TRUE),
  file.path(output,"quality_control","08_RM_executed_sets.csv"))
log("Preflight OK: ",length(tasks)," MI chains | N=",nrow(score),
 " | sets=",paste(names(sets),collapse=",")," | threads=",THREADS)

# Only core individual state is kept: it is sufficient for all requested tables, matrices and contributions.
process_one <- function(task,score,sets,maps,outdir) {
  z<-task$data;n<-nrow(z)
  b19<-as.character(score$score_2019_category);b25<-as.character(score$score2_2025_category)
  ans<-list()
  for(setname in names(sets)) {
    cases<-list()
    for(sc in c("S0","SN","ST","SNT")) {
      v19<-rm_union(z,year_domains(setname,sc,"2019",sets),maps)
      v25<-rm_union(z,year_domains(setname,sc,"2025",sets),maps)
      if(sc=="ST")v19<-rep(0L,n) # comparison is 2019 S0; NOT 2019 ST
      up19<-b19=="M"&v19==1L;up25<-b25=="M"&v25==1L
      if(anyNA(up19)||anyNA(up25))stop("Undetermined moderate risk RM status")
      r19<-b19;r25<-b25;r19[up19]<-"H";r25[up25]<-"H"
      t19<-treat(score$ldl,r19);t25<-treat(score$ldl,r25)
      cases[[sc]]<-data.frame(pidp=score$pidp,baseline_S0_group=score$transition,
        risk19=r19,risk25=r25,origin19=origins(b19,up19),origin25=origins(b25,up25),
        upgrade19=as.integer(up19),upgrade25=as.integer(up25),
        hv19=as.integer(r19%in%c("H","VH")),hv25=as.integer(r25%in%c("H","VH")),
        A19=t19$A,A25=t25$A,D19=t19$D,D25=t25$D,
        transition=transition(t19$A,t25$A,t19$D,t25$D),stringsAsFactors=FALSE)
    }
    if(!identical(cases[["S0"]]$transition,score$transition))stop("S0 transition drift")
    if(!identical(cases[["SN"]]$A19,cases[["SNT"]]$A19)||
       !identical(cases[["SN"]]$risk19,cases[["SNT"]]$risk19))stop("2019 SNT must equal SN")
    for(year in c("19","25")) {
      if(any(cases[["SNT"]][[paste0("A",year)]]<cases[["SN"]][[paste0("A",year)]])||
         any(cases[["SNT"]][[paste0("A",year)]]<cases[["ST"]][[paste0("A",year)]]))stop("SNT intervention monotonicity failed")
    }
    ans[[setname]]<-cases
  }
  saveRDS(list(imputation_id=task$id,N=length(score$pidp),sets=ans),
          file.path(outdir,sprintf("MI_%03d.rds",task$id)),compress="xz")
  data.frame(imputation_id=task$id,N=n,sets=length(ans),completed=TRUE)
}
workers<-min(THREADS,length(tasks));progress<-list(); idx_batches<-split(seq_along(tasks),ceiling(seq_along(tasks)/workers))
if(workers>1L) {
  cl<-parallel::makePSOCKcluster(workers)
  parallel::clusterExport(cl,c("rm_binary","rm_union","year_domains","treat","transition","origins","process_one"),envir=.GlobalEnv)
}
tryCatch({
for(j in seq_along(idx_batches)) {
  ix<-idx_batches[[j]]
  result<-if(workers>1L) parallel::parLapplyLB(cl,tasks[ix],process_one,score=score,sets=sets,maps=maps,outdir=file.path(output,"01_individual_states")) else
    lapply(tasks[ix],process_one,score=score,sets=sets,maps=maps,outdir=file.path(output,"01_individual_states"))
  progress<-c(progress,result)
  log("Core strategies: ",max(ix),"/",length(tasks)," chains (",round(100*max(ix)/length(tasks)),"%)")
}
}, finally=if(workers>1L)parallel::stopCluster(cl))
write_bom(do.call(rbind,progress),file.path(output,"quality_control","08_completed_chains.csv"))
write_bom(data.frame(dataset=DATASET_ID,RM_set=rep(names(sets),each=4),
  strategy=rep(c("S0","SN","ST","SNT"),length(sets)),
  strategy_2019=rep(c("S0","SN","S0 (ST unavailable)","SN (same as SNT)"),length(sets)),
  strategy_2025=rep(c("S0","SN","ST","SNT"),length(sets)),
  comparison_note=rep(c("S0 vs S0","SN vs SN","2019 S0 vs 2025 ST","2019 SN vs 2025 SNT"),length(sets))),
  file.path(output,"quality_control","08_strategy_pairing.csv"))
log("08 completed | ",length(tasks)," imputations, elapsed seconds=",round(as.numeric(difftime(Sys.time(),start_time,units="secs")),1))
