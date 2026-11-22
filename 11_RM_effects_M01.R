# 11 | RM counterfactual effects and paired heterogeneity; scientific formulas unchanged.
project_root<-Sys.getenv('M01_PROJECT_ROOT','D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir<-Sys.getenv('M01_CODE_DIR',file.path(project_root,'LJW_code/Main_code/M01_ASCVD'))
Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_11_RM_effects_M01.R'),encoding='UTF-8')
cfg<-m07_config()
cfg$out_dir<-m07_env()$out_dir
cfg$code_dir<-code_dir
cfg$clean_dir<-Sys.getenv('M01_CLEAN_DIR',cfg$clean_dir)
cfg$workers<-as.integer(Sys.getenv('M01_WORKERS',as.character(cfg$workers)))
invisible(m07_check_manifest(cfg))
z<-m07_worker(seq_len(cfg$m),function(m)m11_task(m,cfg),cfg$workers,'11','RM effects')
m11_validate_outputs(z,cfg)
m09_emit(z,cfg,'strategy_risks_and_effects',c(
  effects='all_scenario_risks_and_effects',
  heterogeneity='all_scenario_guideline_heterogeneity'))
m09_emit(z,cfg,'risk_modifier_incremental_effects',c(
  incremental='incremental_benefits_and_guideline_difference',
  incremental_heterogeneity='incremental_benefit_heterogeneity'))
