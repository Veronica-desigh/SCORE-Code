# 08 | One CSC model per imputation; natural+8 strategy CIF/IF with resumable caches.
# v5 joint matrices use observed AJ in 10; no new counterfactual target or Cox refit.
project_root<-Sys.getenv("M01_PROJECT_ROOT","D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
code_dir<-Sys.getenv("M01_CODE_DIR",file.path(project_root,"LJW_code","Main_code","M01_ASCVD"))
Sys.setenv(M01_CODE_DIR=code_dir)
source(file.path(code_dir,"Function_08_Counterfactual_risk_M01.R"))
if(!requireNamespace("Matrix",quietly=TRUE))stop("Install Matrix")
cfg<-m07_config();m07_check_manifest(cfg)
data.table::setDTthreads(1)
# PSOCK on Windows needs closure loaded explicitly; by default use serial on Windows, MI fork on Unix.
workers<-cfg$workers
z<-m07_worker(seq_len(cfg$m),function(m)m08_predict_mi(m,cfg),workers=workers,stage="08",context="MI prediction")
m07_csv(rbindlist(z,fill=TRUE),file.path(m07_dir(cfg,"08"),"prediction_progress.csv"))
m07_log("08","DONE; competing death hazards frozen across hypothetical LDL states (legacy implementation)")

