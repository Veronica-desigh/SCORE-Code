# S04 | 03 | S0 guideline risk distributions and six LDL-C treatment transitions
# Self-contained descriptive summary: no rescoring, no reclassification, no MI.
# Run after 02_Guideline_strategies_S04.R

THREADS <- 5L # configured globally; intentionally no parallelism for full-cohort sums
THREADS <- as.integer(Sys.getenv("S04_WORKERS", as.character(THREADS)))
if(is.na(THREADS) || THREADS < 1L) stop("Invalid S04_WORKERS")
start_time <- Sys.time()
log <- function(...) message(sprintf("[%s +%.1fs] %s",format(Sys.time(),"%H:%M:%S"),
  as.numeric(difftime(Sys.time(),start_time,units="secs")),paste0(...,collapse="")))
script_path <- function(){
  a <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
  f <- Filter(Negate(is.null),lapply(sys.frames(),function(z)z$ofile))
  if(length(f)) return(dirname(normalizePath(tail(f,1)[[1]],mustWork=FALSE)))
  if(length(a)) return(dirname(normalizePath(sub("^--file=","",a[1]),mustWork=FALSE)))
  normalizePath(getwd(),mustWork=TRUE)
}
script_dir <- script_path()
root <- Sys.getenv("S04_PROJECT_ROOT", normalizePath(file.path(script_dir,"../../.."),mustWork=FALSE))
analysis_dir <- Sys.getenv("S04_ANALYSIS_DIR",file.path(root,"LJW_outcome","data_intervention","S04_ELSA","complete_case"))
write_csv <- function(d,name){
  dest <- file.path(analysis_dir,paste0("03_",name,"_S04.csv"))
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  write.csv(d,dest,row.names=FALSE,na="",fileEncoding="UTF-8")
  log("Saved ",basename(dest),": ",nrow(d)," rows")
}
p <- file.path(analysis_dir,"02_score_strategy_participants_S04.csv")
if(!file.exists(p)) stop("Step 02 required: ",p)
d <- read.csv(p,colClasses=c(idauniq="character"),check.names=FALSE,stringsAsFactors=FALSE)
required <- c("idauniq","ldl","score_2019_category","score2_2025_category",
  "trigger_2019","trigger_2025","minimum_reduction_2019","minimum_reduction_2025","transition")
missing <- setdiff(required,names(d))
if(length(missing)) stop("Step 02 missing fields: ",paste(missing,collapse=", "))
if(!nrow(d)||anyNA(d[required])||anyDuplicated(d$idauniq)) stop("Invalid 02 scoring data")
levels <- c("L","M","H","VH")
groups <- c("unchanged_untreated","newly_intervention","intervention_removed",
            "treatment_intensified","unchanged_intensity","treatment_deintensified")
if(any(!d$transition %in% groups)||any(!as.character(d$score_2019_category)%in%levels)||
   any(!as.character(d$score2_2025_category)%in%levels)) stop("Invalid risk or transition labels")
if(any(!d$trigger_2019 %in% 0:1)||any(!d$trigger_2025 %in% 0:1)||
   any(!is.finite(d$minimum_reduction_2019))||any(!is.finite(d$minimum_reduction_2025))) stop("Invalid treatment fields")
n <- nrow(d)
log("Loaded fixed cohort n=",n,"; starting summary")
risk <- do.call(rbind,lapply(c("2019","2025"),function(yr){
  v <- if(yr=="2019")as.character(d$score_2019_category) else as.character(d$score2_2025_category)
  res <- table(factor(v,levels=levels))
  data.frame(guideline=yr,category=c(levels,"H+VH"),N=n,
    n=c(as.integer(res),sum(v%in%c("H","VH"))),
    pct_total=100*c(as.integer(res),sum(v%in%c("H","VH")))/n)
}))
write_csv(risk,"risk_level_H_VH")
log("Risk distribution complete (1/3)")
transitions <- do.call(rbind,lapply(groups,function(g){
 ix <- d$transition==g; ng<-sum(ix)
 data.frame(transition=g,N=n,n=ng,pct_total=100*ng/n,
   n_eligible_2019=sum(d$trigger_2019[ix]==1L),
   n_eligible_2025=sum(d$trigger_2025[ix]==1L),
   mean_minimum_reduction_2019=if(ng)mean(d$minimum_reduction_2019[ix])else NA_real_,
   mean_minimum_reduction_2025=if(ng)mean(d$minimum_reduction_2025[ix])else NA_real_,
   mean_minimum_reduction_difference=if(ng)mean(d$minimum_reduction_2025[ix]-d$minimum_reduction_2019[ix])else NA_real_)
}))
if(sum(transitions$n)!=n || length(unique(transitions$transition))!=6L)stop("Six-group exhaustiveness failed")
write_csv(transitions,"six_transitions_LDL_reduction")
log("Treatment transitions complete (2/3)")
qa <- data.frame(item=c("n_fixed","n_in_six_groups","n_eligible_2019","n_eligible_2025",
                           "n_HVH_2019","n_HVH_2025","score_csv_md5"),
 value=c(n,sum(transitions$n),sum(d$trigger_2019),sum(d$trigger_2025),
         sum(d$score_2019_category %in% c("H","VH")),
         sum(d$score2_2025_category %in% c("H","VH")),unname(tools::md5sum(p))))
write_csv(qa,"verification")
log("All S0 summaries complete (3/3); total elapsed seconds=",
 round(as.numeric(difftime(Sys.time(),start_time,units="secs")),1))
