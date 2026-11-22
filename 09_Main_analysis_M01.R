# S0 standardized effects and observed competing-risk risks.
# v5 tables 11/12 and figure 2 reuse the saved scope NC and paired RD estimates.
project_root<-Sys.getenv('M01_PROJECT_ROOT','D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code')
code_dir<-Sys.getenv('M01_CODE_DIR',file.path(project_root,'LJW_code/Main_code/M01_ASCVD'));Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,'Function_09_Main_analysis_M01.R'))
cfg<-m07_config();m07_check_manifest(cfg)
z<-m07_worker(seq_len(cfg$m),function(m)m09_task(m,cfg),cfg$workers,'09','S0 effects and observed risks')
m09_emit(z,cfg,'strategy_risks_and_effects',c(effects='strategy_risks_and_risk_differences',heterogeneity='guideline_risk_difference_heterogeneity'))
m09_emit(z,cfg,'observed_risk_and_event_capture',c(observed='observed_population_risks'))
d<-m09_pool(rbindlist(lapply(z,`[[`,'effects')),cfg)
m07_result(cfg,'unmeasured_confounding_sensitivity','risk_ratio_evalues',m09_evalues(d))
