# S1：仅删除五类年度适用proxy RM；S2：仅删除家族史/心理社会应激；S3：两个独立阈值调整。
# STRICT现为固化的集合删除，不依赖任何旧版definitions$strict确认函数；缺失P90直接报错。
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
function_file <- file.path(code,"Function_02_Strict_RM_S01.R")
if(!file.exists(function_file)) stop("缺少对应函数文件: ",function_file)
options(S01_CURRENT_FUNCTION=basename(function_file))
script_started <- proc.time()[["elapsed"]]
message(sprintf("[%s] 开始执行 02_Strict_RM_S01 | workers=%s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),workers))
on.exit_msg <- function(){
  spent <- proc.time()[["elapsed"]]-script_started
  message(sprintf("[%s] 02_Strict_RM_S01 已运行 %.1f秒（%.2f分钟）",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),spent,spent/60))
}
source(function_file,encoding="UTF-8")
definitions <- s01_definitions(code)
definitions$strict <- NULL  # 不使用旧式严格确认函数
p90_frozen <- s01_townsend_p90(root)
if(!is.null(definitions$p90) && !isTRUE(all.equal(as.numeric(definitions$p90),p90_frozen,tolerance=1e-10)))
  stop("旧配置中的Townsend P90与M01固定基线队列不一致，请核对")
definitions$p90 <- p90_frozen
s01_run(root,"STRICT",s01_strict_scenarios,definitions,workers=workers)
on.exit_msg()
