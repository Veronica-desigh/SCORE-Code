# 12 | Six-way treatment transitions / IF contribution decomposition.
# Only the optional summaries-only mode reads already saved 10/12 CSVs.
project_root<-Sys.getenv('M01_PROJECT_ROOT','D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir<-Sys.getenv('M01_CODE_DIR',file.path(project_root,'LJW_code/Main_code/M01_ASCVD'))
Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_12_Treatment_transitions_M01.R'),encoding='UTF-8')
cfg<-m07_config()
cfg$out_dir<-m07_env()$out_dir
cfg$code_dir<-code_dir
cfg$clean_dir<-Sys.getenv('M01_CLEAN_DIR',cfg$clean_dir)
cfg$workers<-as.integer(Sys.getenv('M01_WORKERS',as.character(cfg$workers)))
if(identical(Sys.getenv('M01_REDUCTION_SUMMARIES_ONLY','0'),'1')){
  m07_timed('12','saved mean reduction summaries only',m12_mean_reduction_exports(cfg))
}else{
  invisible(m07_check_manifest(cfg))
  z<-m07_worker(seq_len(cfg$m),function(m)m12_task(m,cfg),cfg$workers,
    '12','six-category transitions and contribution decomposition')
  m12_export(z,cfg)
}
