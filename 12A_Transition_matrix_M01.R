# 12A | Three S0 -> RM six-by-six matrices (SN, ST, SNT).
# Run after stage 07. No plots, outcome model fitting or MI re-imputation.
options(stringsAsFactors=FALSE)
project_root <- Sys.getenv('M01_PROJECT_ROOT',
  'D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir <- Sys.getenv('M01_CODE_DIR',
  file.path(project_root,'LJW_code','Main_code','M01_ASCVD'))
Sys.setenv(M01_PROJECT_ROOT=project_root,M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_12A_Transition_matrix_M01.R'),encoding='UTF-8')
cfg <- m07_config()
cfg$code_dir <- code_dir
cfg$out_dir <- m07_env()$out_dir
cfg$clean_dir <- Sys.getenv('M01_CLEAN_DIR',cfg$clean_dir)
cfg$root <- project_root
cfg$workers <- m12a_check_workers() # default 3, set M01_WORKERS to override
if(!requireNamespace('data.table',quietly=TRUE))stop('Install data.table')
data.table::setDTthreads(1L)
t0 <- proc.time()[3]
m07_log('12A','START | 3 S0->RM matrices | workers=',cfg$workers)
invisible(m07_check_manifest(cfg))
result <- m12x_parallel(seq_len(cfg$m),cfg,'12A',
  'Function_12A_Transition_matrix_M01.R','m12a_task',workers=cfg$workers)
m07_timed('12A','pooled summaries, plotting CSV and QA',m12a_export(result,cfg))
m07_log('12A','DONE | elapsed ',round(proc.time()[3]-t0,1),' sec')
