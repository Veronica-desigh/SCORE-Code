# 10_Run_all_S01.R | 按00 -> 09严格顺序逐一运行，每个子脚本单独R进程，不共享全局状态。
# 仅此一个总运行脚本；不新增辅助脚本。失败立即停止，禁止09在上游失败时继续汇总。
# ============================= 用户可修改参数 =============================
workers_override <- NA_integer_ # 默认NA:每个主脚本使用各自workers <- 3L；例如改为2L统一覆盖
start_step <- 0L          # 如需断点续跑可设为 1L...9L；前置步骤必须已有完整结果
end_step <- 9L            # 可设置提前结束；正式完整运行保持9L
s01_default_root <- "D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code"
# ========================================================================
root <- Sys.getenv("S01_PROJECT_ROOT", unset="")
if(!nzchar(root)) root <- s01_default_root
root <- gsub("\\", "/", root, fixed=TRUE)
root <- sub("/+$", "", root)
suffix <- "/LJW_code/Sensitivity_code/S01_different_RM_strategies"
if(endsWith(tolower(root),tolower(suffix))) root <- substr(root,1L,nchar(root)-nchar(suffix))
code <- Sys.getenv("S01_CODE_DIR",unset="")
if(!nzchar(code)) code <- file.path(root,"LJW_code","Sensitivity_code","S01_different_RM_strategies")
code <- normalizePath(code,winslash="/",mustWork=TRUE)
if(!dir.exists(root)) stop("项目根目录不存在: ",root)
if(!is.na(workers_override)){
  if(length(workers_override)!=1L || !is.finite(workers_override) ||
     workers_override<1L || workers_override!=as.integer(workers_override))
    stop("workers_override 必须为正整数或 NA_integer_")
  Sys.setenv(S01_WORKERS_OVERRIDE=as.character(workers_override))
} else Sys.unsetenv("S01_WORKERS_OVERRIDE")
Sys.setenv(S01_PROJECT_ROOT=root,S01_CODE_DIR=code,
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
script_names <- c(
  "00_Validate_inputs_S01.R",
  "01_Higher_moderate_risk_S01.R",
  "02_Strict_RM_S01.R",
  "03_Observed_confirmed_RM_S01.R",
  "04_Common_RM_S01.R",
  "05_Leave_one_out_RM_S01.R",
  "06_Multidomain_RM_S01.R",
  "07_One_at_a_time_RM_S01.R",
  "08_Fixed_population_mechanisms_S01.R",
  "09_Pool_and_export_S01.R")
if(length(start_step)!=1L || length(end_step)!=1L ||
   is.na(start_step) || is.na(end_step) ||
   start_step<0L || end_step>9L || start_step>end_step ||
   start_step!=as.integer(start_step) || end_step!=as.integer(end_step))
  stop("start_step/end_step 必须满足0 <= start_step <= end_step <= 9")
# 在Windows上Rscript.exe通常位于R.home('bin')；不同安装结构做兼容检查。
rscript <- Sys.which("Rscript")
if(!nzchar(rscript)){
  choices <- c(file.path(R.home("bin"),"Rscript.exe"),
               file.path(R.home("bin"),"Rscript"),
               file.path(R.home("bin"),"x64","Rscript.exe"))
  if(!any(file.exists(choices))) stop("找不到Rscript可执行文件；请先安装R或配置PATH")
  rscript <- choices[which(file.exists(choices))[1L]]
}
steps <- seq.int(start_step,end_step)
run_started <- Sys.time()
cat(sprintf("[%s] S01总运行开始 | %d个步骤 | 每脚本线程=%s\n",
  format(Sys.time(),"%Y-%m-%d %H:%M:%S"),length(steps),
  if(is.na(workers_override)) "各自参数 (默认3)" else as.character(workers_override)))
timings <- list()
for(step in steps){
  name <- script_names[[step+1L]]
  path <- file.path(code,name)
  if(!file.exists(path)) stop("找不到主脚本: ",path)
  t0 <- Sys.time()
  cat(sprintf("\n[%s] [%02d/%02d] 开始: %s\n",
      format(t0,"%H:%M:%S"),match(step,steps),length(steps),name))
  flush.console()
  status <- tryCatch(system2(rscript,args=c("--vanilla",shQuote(path)),
                            stdout="",stderr="",wait=TRUE),
    error=function(e) {message("执行异常：",conditionMessage(e)); 1L})
  took <- as.numeric(difftime(Sys.time(),t0,units="secs"))
  timings[[length(timings)+1L]] <- data.frame(
    step=step,script=name,exit_status=as.integer(status),seconds=took,
    finished_at=format(Sys.time(),"%Y-%m-%d %H:%M:%S"))
  cat(sprintf("[%s] [%02d/%02d] %s | 退出码=%d | 耗时 %.1fs | 总进度 %.1f%%\n",
    format(Sys.time(),"%H:%M:%S"),match(step,steps),length(steps),name,
    as.integer(status),took,100*match(step,steps)/length(steps)))
  flush.console()
  if(is.na(status) || status!=0L) stop("脚本失败：",name,"；停止后续运行，保护汇总完整性")
}
total_seconds <- as.numeric(difftime(Sys.time(),run_started,units="secs"))
cat(sprintf("\n[%s] S01全部完成 | %.1f秒 (%.2f分钟; %.2f小时)\n",
  format(Sys.time(),"%Y-%m-%d %H:%M:%S"),total_seconds,total_seconds/60,total_seconds/3600))
print(do.call(rbind,timings),row.names=FALSE)
