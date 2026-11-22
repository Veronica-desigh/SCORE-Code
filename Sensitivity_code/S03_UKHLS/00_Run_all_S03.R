# S03 UKHLS: one launcher for all ten independent R scripts.
# Every script is self-contained; no old S03 Function_* scripts are sourced.
# ---------------------- EDIT HERE ----------------------
THREADS <- 5L
TEST_N <- 0L                # 0=full sample; >0=test (01-04 linked sample)
RUN_STEPS <- 1:10           # e.g. 7:10 after 01-06 were run
M <- 5L
MAXIT <- 5L
SMOKING_MODE <- "same_wave" # same_wave / wave2_prior
PROJECT_ROOT <- ""         # absolute project root (contains LJW_outcome and DATA)
SCRIPT_DIR <- ""           # absolute folder of 01-10 .R files
# -------------------------------------------------------
files <- c("01_Clean_interview_S03.R","02_Clean_nurse_S03.R",
           "03_Clean_blood_S03.R","04_Eligibility_criteria_S03.R",
           "05_Guideline_strategies_S03.R","06_Treatment_transitions_S03.R",
           "07_RM_imputation_distribution_S03.R","08_RM_strategies_S03.R",
           "09_RM_tables_matrices_contributions_S03.R","10_RM_vector_matrices_S03.R")
if(length(RUN_STEPS)<1L || anyNA(RUN_STEPS) || any(!RUN_STEPS %in% 1:10) ||
   !identical(as.integer(RUN_STEPS),sort(unique(as.integer(RUN_STEPS))))) stop("RUN_STEPS must be sorted unique 1:10")
if(!is.finite(THREADS) || THREADS<1 || !is.finite(TEST_N) || TEST_N<0 ||
   !is.finite(M) || M<1 || !is.finite(MAXIT) || MAXIT<1)stop("Invalid numeric settings")
if(!SMOKING_MODE %in% c("same_wave","wave2_prior"))stop("Unknown SMOKING_MODE")
find_scripts <- function(){
  f <- grep("^--file=",commandArgs(),value=TRUE)
  contexts <- Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
  bases <- c(if(nzchar(SCRIPT_DIR))SCRIPT_DIR,
             if(length(contexts))dirname(tail(contexts,1)[[1]]),
             if(length(f))dirname(sub("^--file=","",f[1])),getwd())
  if(nzchar(SCRIPT_DIR)){
    if(!dir.exists(SCRIPT_DIR))stop("SCRIPT_DIR missing: ",SCRIPT_DIR)
    return(normalizePath(SCRIPT_DIR,winslash="/",mustWork=TRUE))
  }
  candidates <- character()
  for(base in bases){
    cur<-normalizePath(base,winslash="/",mustWork=FALSE)
    for(j in 1:12){
      candidates<-c(candidates,cur,file.path(cur,"LJW_code","Sensitivity_code","S03_UKHLS"),
                    file.path(cur,"S03_UKHLS"))
      upper<-dirname(cur);if(upper==cur)break;cur<-upper
    }
  }
  candidates<-unique(candidates)
  found<-candidates[vapply(candidates,function(d)all(file.exists(file.path(d,files[RUN_STEPS]))),logical(1))]
  if(!length(found))stop("S03 scripts cannot be located. Set SCRIPT_DIR explicitly. Working directory: ",getwd())
  normalizePath(found[1],winslash="/",mustWork=TRUE)
}
code_dir<-find_scripts()
missing<-files[RUN_STEPS][!file.exists(file.path(code_dir,files[RUN_STEPS]))]
if(length(missing))stop("Scripts missing: ",paste(missing,collapse=", "))
root<-if(nzchar(PROJECT_ROOT))normalizePath(PROJECT_ROOT,winslash="/",mustWork=TRUE) else code_dir
if(!nzchar(PROJECT_ROOT)){
  while(!dir.exists(file.path(root,"LJW_outcome")) && !dir.exists(file.path(root,"DATA","UKHLS","UKHLS","csv"))){
    upper<-dirname(root);if(upper==root)stop("Cannot detect project root. Specify PROJECT_ROOT")
    root<-upper
  }
}
if(any(RUN_STEPS<=4) && !dir.exists(file.path(root,"DATA","UKHLS","UKHLS","csv")))stop("Missing UKHLS source directory")
if(5 %in% RUN_STEPS && !file.exists(file.path(root,"LJW_code","Main_code","M01_ASCVD","Function_07_Guideline_strategies_M01.R")))stop("Missing M01 SCORE engine")
for(i in RUN_STEPS){
  path<-file.path(code_dir,files[i]);tryCatch(parse(path,encoding="UTF-8"),
    error=function(e)stop("Parse error in ",files[i],": ",conditionMessage(e)))
}
base_out<-file.path(root,"LJW_outcome","data_intervention","S03_UKHLS")
if(TEST_N>0L)base_out<-file.path(base_out,paste0("test_",TEST_N))
if(SMOKING_MODE=="wave2_prior")base_out<-file.path(base_out,"wave2_smoking_sensitivity")
logdir<-file.path(base_out,"logs");dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
logfile<-file.path(logdir,"00_Run_all_S03_runtime.log")
t0<-Sys.time()
log<-function(...){line<-sprintf("[%s +%.1fs] %s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
  as.numeric(difftime(Sys.time(),t0,units="secs")),paste0(...,collapse=""))
  message(line);cat(line,"\n",file=logfile,append=TRUE)}
keys<-c("S03_THREADS","S03_SAMPLE_N","S03_ANALYSIS_TEST_N","S03_M","S03_MAXIT",
        "S03_PROJECT_ROOT","S03_RM_PROJECT_ROOT","S03_RM_THREADS","S03_RM_TEST_N",
        "S03_SMOKING_MODE","S03_RM_SCORE_PATH","S03_RM_MI_PATH","S03_RM_MAP")
original<-Sys.getenv(keys,unset=NA_character_)
restore<-function(){for(j in seq_along(keys))if(is.na(original[j]))Sys.unsetenv(keys[j])
  else do.call(Sys.setenv,setNames(list(original[j]),keys[j]))}
Sys.setenv(S03_THREADS=as.character(THREADS),S03_SAMPLE_N=as.character(TEST_N),
 S03_ANALYSIS_TEST_N=as.character(TEST_N),S03_M=as.character(M),
 S03_MAXIT=as.character(MAXIT),S03_PROJECT_ROOT=root,
 S03_RM_PROJECT_ROOT=root,S03_RM_THREADS=as.character(THREADS),
 S03_RM_TEST_N=as.character(TEST_N),S03_SMOKING_MODE=SMOKING_MODE,
 S03_RM_SCORE_PATH="",S03_RM_MI_PATH="",S03_RM_MAP="")
rs<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
if(!file.exists(rs))rs<-"Rscript"
summary<-list()
log("S03 START | steps=",paste(RUN_STEPS,collapse=",")," | workers=",THREADS,
  " | M=",M," | iterations=",MAXIT," | scripts=",code_dir," | root=",root)
tryCatch({
 for(k in seq_along(RUN_STEPS)){
   i<-RUN_STEPS[k];st<-Sys.time();path<-file.path(code_dir,files[i])
   log("STEP ",i,"/10 START ",basename(path))
   result<-tryCatch(system2(rs,args=shQuote(path),stdout="",stderr=""),error=function(e)e)
   code<-if(inherits(result,"error"))1L else if(length(result)==1L)as.integer(result)else 1L
   summary[[k]]<-data.frame(step=i,script=files[i],exit_status=code,
       elapsed_seconds=round(as.numeric(difftime(Sys.time(),st,units="secs")),2))
   utils::write.csv(do.call(rbind,summary),file.path(logdir,"00_Run_all_S03_summary.csv"),row.names=FALSE)
   log("STEP ",i,"/10 END | status=",code," | seconds=",summary[[k]]$elapsed_seconds)
   if(is.na(code)||code!=0L)stop("S03 step ",i," failed. Check step-specific log.")
 }
 log("S03 FINISHED | all requested steps succeeded")
},finally=restore())
