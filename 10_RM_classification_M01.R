# 10 | Unified RM classification and normalized observed event capture.
# Do NOT read legacy stage-10 caches. All classification/count/risk outputs
# are exported from the same per-MI result (new unified cache v6).
project_root<-Sys.getenv('M01_PROJECT_ROOT','D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir<-Sys.getenv('M01_CODE_DIR',file.path(project_root,'LJW_code/Main_code/M01_ASCVD'))
Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_10_RM_classification_M01.R'),encoding='UTF-8')
if(identical(Sys.getenv('M01_CLASSIFICATION_MODE','all'),'partitions'))
  stop('10: retired partitions-only mode. This version must export ALL outputs together; unset M01_CLASSIFICATION_MODE.')
cfg<-m07_config()
cfg$out_dir<-m07_env()$out_dir
cfg$code_dir<-code_dir
cfg$clean_dir<-Sys.getenv('M01_CLEAN_DIR',cfg$clean_dir)
cfg$workers<-as.integer(Sys.getenv('M01_WORKERS',as.character(cfg$workers)))
if(identical(Sys.getenv('M01_PUBLICATION_SUMMARIES_ONLY','0'),'1')){
  invisible(m07_check_manifest(cfg,check_code=FALSE))
  m10_publication_summaries(cfg);quit(save='no',status=0)
}
invisible(m07_check_manifest(cfg))
z<-m07_worker(seq_len(cfg$m),function(m)m10_task(m,cfg),cfg$workers,
  '10','unified RM classification, H/T partitions, normalized event capture')
m10_export_unified(z,cfg)
m10_publication_summaries(cfg)
rm(z);invisible(gc())
m07_log('10','Unified outputs completed; retired caches were neither read nor deleted')
