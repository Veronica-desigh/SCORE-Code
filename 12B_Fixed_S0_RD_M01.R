# 12B | Fixed S0 groups: guideline RD, RM-minus-S0 changes and weighted contributions.
# Run after stages 07 and 08. No refit; no plots.
options(stringsAsFactors=FALSE)
project_root <- Sys.getenv('M01_PROJECT_ROOT',
  'D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir <- Sys.getenv('M01_CODE_DIR',
  file.path(project_root,'LJW_code','Main_code','M01_ASCVD'))
Sys.setenv(M01_PROJECT_ROOT=project_root,M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_12B_Fixed_S0_RD_M01.R'),encoding='UTF-8')
cfg <- m07_config()
cfg$code_dir <- code_dir
cfg$out_dir <- m07_env()$out_dir
cfg$clean_dir <- Sys.getenv('M01_CLEAN_DIR',cfg$clean_dir)
cfg$root <- project_root
cfg$workers <- m12a_check_workers() # default 3, set M01_WORKERS to override
if(!requireNamespace('survival',quietly=TRUE))stop('Install survival')
if(!requireNamespace('Matrix',quietly=TRUE))stop('Install Matrix')
data.table::setDTthreads(1L)
t0 <- proc.time()[3]
m07_log('12B','START fixed S0 risk-difference analysis | workers=',cfg$workers,
        ' | times=',paste(cfg$times,collapse=','))
invisible(m07_check_manifest(cfg))
result <- m12x_parallel(seq_len(cfg$m),cfg,'12B',
  'Function_12B_Fixed_S0_RD_M01.R','m12b_task',workers=cfg$workers)
m07_timed('12B','Rubin pooling, plotting CSV and QA',m12b_export(result,cfg))
m07_log('12B','DONE | elapsed ',round(proc.time()[3]-t0,1),' sec')
