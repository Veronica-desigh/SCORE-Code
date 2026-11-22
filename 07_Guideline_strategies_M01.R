# 07 | Prepare guideline-RM strategies from completed MI; no outcome-model refit.
# One main script sources one corresponding function script.
options(stringsAsFactors=FALSE)
project_root<-Sys.getenv("M01_PROJECT_ROOT","D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
code_dir<-Sys.getenv("M01_CODE_DIR",file.path(project_root,"LJW_code","Main_code","M01_ASCVD"))
Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,"Function_07_Guideline_strategies_M01.R"),encoding="UTF-8")
if(identical(Sys.getenv('M01_MAIN_TABLE1_ONLY','0'),'1')){
  cfg<-m07_config();cfg$out_dir<-m07_env()$out_dir;cfg$code_dir<-code_dir
  cfg$clean_dir<-m07_env()$clean_dir
  m07_main_table1(cfg);quit(save='no',status=0)
}
if(identical(Sys.getenv('M01_REPORT_SUMMARIES_ONLY','0'),'1')){
  cfg<-m07_config();cfg$out_dir<-m07_env()$out_dir;cfg$code_dir<-code_dir
  cfg$clean_dir<-m07_env()$clean_dir
  invisible(m07_check_manifest(cfg,check_code=FALSE))
  m07_reporting_summaries(cfg);m07_main_table1(cfg);quit(save='no',status=0)
}
if(!requireNamespace("RiskScorescvd",quietly=TRUE))stop("Install RiskScorescvd")
cfg<-m07_env()
cfg$m<-as.integer(Sys.getenv("M01_MI_M","5"))
cfg$workers<-as.integer(Sys.getenv("M01_WORKERS","2"))
cfg$score_workers<-as.integer(Sys.getenv("M01_SCORE_WORKERS","2"))
cfg$score_region<-Sys.getenv("M01_SCORE_REGION","Low")
cfg$score_seed<-2026L
cfg$score_cutpoints<-c(1,5,10);cfg$score2_cutpoints<-c(2,10,20)
cfg$ldl_targets<-c(L=3,M=2.6,H=1.8,VH=1.4)
cfg$ldl_delta<-0.4;cfg$ldl_min<-0.6
cfg$times<-c(10,15);cfg$grid_days<-as.integer(Sys.getenv("M01_GRID_DAYS","14"))
cfg$prediction_chunk<-as.integer(Sys.getenv("M01_PREDICTION_CHUNK","2048"))
cfg$confidence_level<-.95;cfg$qa_tolerance<-1e-6
cfg$lipid_lookback_days<-365L
if(anyNA(c(cfg$m,cfg$workers,cfg$score_workers,cfg$grid_days,cfg$prediction_chunk)) || cfg$m<2 || cfg$workers<1 || cfg$grid_days<1)stop("Invalid configuration")
dir.create(cfg$out_dir,recursive=TRUE,showWarnings=FALSE)
data.table::setDTthreads(1)
mi<-m07_manifest(cfg)
base_path<-file.path(cfg$clean_dir,"data_analysis_base_incomplete_M01.rds")
if(!file.exists(base_path))base_path<-file.path(cfg$clean_dir,"data_analysis_base_M01.csv")
if(!file.exists(base_path))stop("Cannot find 05 analysis base: ",base_path)
base<-if(grepl("[.]rds$",base_path))readRDS(base_path) else data.table::fread(base_path)
base<-data.table::as.data.table(base);base[,participant_id:=as.character(participant_id)];m07_unique(base)
m07_assert(base,c("age","sex","smoking_status","systolic_blood_pressure","diastolic_blood_pressure","total_cholesterol","ldl_cholesterol","hdl_cholesterol"))
if(!"current_smoker"%in%names(base))base<-derive_current_smoker(base)
out<-m07_dir(cfg,"07")
master_path<-file.path(out,"observed_score_master.rds")
master_sig_path<-paste0(master_path,".signature.rds")
score_signature<-paste(m07_md5(base_path),m07_md5(file.path(code_dir,"Function_07_Guideline_strategies_M01.R")),
  cfg$score_region,cfg$score_seed,paste(cfg$score_cutpoints,collapse=":"),paste(cfg$score2_cutpoints,collapse=":"),sep="|")
if(file.exists(master_path)&&file.exists(master_sig_path)&&identical(readRDS(master_sig_path),score_signature)) {
  master<-readRDS(master_path)
  m07_log("07","SCORE/SCORE2 frozen master cache hit")
} else {
  scored<-m07_timed("07","SCORE + SCORE2 score calculation",calculate_score_values(base,cfg$score_region,cfg$score_workers))
  scored<-add_score_categories(scored,cfg$score_cutpoints,cfg$score2_cutpoints)
  scored<-add_stable_random_quantile(scored,seed=cfg$score_seed)
  master<-scored[,.(participant_id,age,sex,current_smoker,systolic_blood_pressure,diastolic_blood_pressure,total_cholesterol,ldl_cholesterol,hdl_cholesterol,
                    score_2019_percent,score2_2025_percent,score_2019_category,score2_2025_category,z_common)]
  master[,score_2019_category:=as.character(score_2019_category)]
  master[,score2_2025_category:=as.character(score2_2025_category)]
  m07_atomic(master,master_path);m07_atomic(score_signature,master_sig_path)
}
m07_atomic(mi,file.path(out,"input_manifest.rds"))
m07_atomic(list(path=base_path,md5=m07_md5(base_path)),file.path(out,"observed_base_signature.rds"))
m07_saveconfig(cfg)
fn<-function(m){
  p<-file.path(out,sprintf("m%02d.rds",m));sigfile<-paste0(p,".sig")
  sig<-paste(mi[imputation_id==m,md5],m07_md5(file.path(code_dir,"Function_07_Guideline_strategies_M01.R")),
    m07_md5(file.path(out,"observed_score_master.rds")),m07_scientific_key(cfg),sep="|")
  if(file.exists(p)&&file.exists(sigfile)&&identical(readRDS(sigfile),sig)){m07_log("07","MI",m," cache hit");return(data.table(imputation_id=m,reused=TRUE))}
  completed<-readRDS(mi[imputation_id==m,path]);z<-m07_timed("07",paste0("MI",m," scenarios + QA"),m07_strategy_one(completed,master,cfg,m))
  m07_atomic(z,p);m07_atomic(sig,sigfile)
  data.table(imputation_id=m,reused=FALSE,n=nrow(z))
}
# Sequential processing avoids duplicating large completed datasets in memory.
z<-m07_worker(seq_len(cfg$m),fn,workers=1L,stage="07",context="MI strategy")
m07_csv(rbindlist(z,fill=TRUE),file.path(out,"strategy_assignment_progress.csv"))
m07_log("07","DONE | source lookback 365 days | ",cfg$m," imputations")
m07_reporting_summaries(cfg)
m07_main_table1(cfg)
