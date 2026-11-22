# 仅原始资料确认阳性且符合原版定义的 RM 可触发。与相应原版比较升级损失、新增干预、降幅及配对风险差。
# 人数和比例逐份插补计算，风险和差值使用同一批个体的配对影响函数。
# ============================= 用户可修改参数 =============================
workers <- 3L  # 并行工作进程数，默认3；可设1为串行。总脚本可用S01_WORKERS_OVERRIDE覆盖。
s01_default_root <- "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code"
# ========================================================================
root <- Sys.getenv("S01_PROJECT_ROOT", unset="")
if(!nzchar(root)) root <- s01_default_root
root <- gsub("\\", "/", root, fixed=TRUE)
root <- sub("/+$", "", root)
s01_suffix <- "/LJW_code/Sensitivity_code/S01_different_RM_strategies"
if(endsWith(tolower(root),tolower(s01_suffix)))
  root <- substr(root,1L,nchar(root)-nchar(s01_suffix))
code <- Sys.getenv("S01_CODE_DIR",unset="")
if(!nzchar(code)) code <- file.path(root,"LJW_code","Sensitivity_code","S01_different_RM_strategies")
code <- gsub("\\", "/", code, fixed=TRUE)
code <- sub("/+$", "", code)
if(!dir.exists(root)) stop("S01项目根目录不存在: ",root)
if(!dir.exists(code)) stop("S01代码目录不存在: ",code)
Sys.setenv(S01_PROJECT_ROOT=root,S01_CODE_DIR=code,
           OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
worker_override <- Sys.getenv("S01_WORKERS_OVERRIDE",unset="")
if(nzchar(worker_override)) workers <- suppressWarnings(as.numeric(worker_override))
function_file <- file.path(code,"Function_03_Observed_confirmed_RM_S01.R")
if(!file.exists(function_file)) stop("缺少对应函数文件: ",function_file)
options(S01_CURRENT_FUNCTION=basename(function_file))
script_started <- proc.time()[["elapsed"]]
message(sprintf("[%s] 开始执行 03_Observed_confirmed_RM_S01 | workers=%s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),workers))
on.exit_msg <- function(){
  spent <- proc.time()[["elapsed"]]-script_started
  message(sprintf("[%s] 03_Observed_confirmed_RM_S01 已运行 %.1f秒（%.2f分钟）",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),spent,spent/60))
}
source(function_file,encoding="UTF-8")
definitions <- s01_definitions(code)
if(is.null(definitions$observed)) definitions$observed <- s01_observed_rule(root)
s01_run(root,"OBS",s01_obs_scenarios,definitions,workers=workers)
on.exit_msg()
