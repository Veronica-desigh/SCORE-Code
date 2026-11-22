# FIX 固定 S0 干预强度减弱者，LOST 固定原版 SNT 升级但 strict 不升级者。
# 输出组内分类、最低降幅、观察风险、模型估计风险及配对差的总体加权贡献。
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
function_file <- file.path(code,"Function_08_Fixed_population_mechanisms_S01.R")
if(!file.exists(function_file)) stop("缺少对应函数文件: ",function_file)
options(S01_CURRENT_FUNCTION=basename(function_file))
script_started <- proc.time()[["elapsed"]]
message(sprintf("[%s] 开始执行 08_Fixed_population_mechanisms_S01 | workers=%s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),workers))
on.exit_msg <- function(){
  spent <- proc.time()[["elapsed"]]-script_started
  message(sprintf("[%s] 08_Fixed_population_mechanisms_S01 已运行 %.1f秒（%.2f分钟）",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),spent,spent/60))
}
source(function_file,encoding="UTF-8")
definitions <- s01_definitions(code)
# FIX 的额外情景允许由 definitions$fix_scenarios 指定；
# LOST 必须始终包含新版（纯删除）SNT_STRICT，不依赖旧版 definitions$strict。
fixed <- definitions$fix_scenarios
if (is.null(fixed)) fixed <- list()
if (!is.list(fixed)) stop("definitions$fix_scenarios 必须是情景列表")
fix_ids <- vapply(fixed, function(z) z$id, character(1))
if (anyDuplicated(fix_ids)) stop("FIX 额外情景 ID 重复")
if (!"SNT_STRICT" %in% fix_ids) {
  fixed <- c(fixed, list(s01_scenario("SNT_STRICT", "SNT", "strict")))
}
s01_run(root,"FIX_LOST",fixed,definitions,mechanisms=TRUE,workers=workers)
on.exit_msg()
