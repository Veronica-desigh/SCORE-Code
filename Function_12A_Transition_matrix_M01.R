# MODULE 12A | S0 -> SN/ST/SNT six-way treatment-decision transition matrices.
# Read-only reuse of stage-07 completed MI strategy assignments.
# Exactly six rows x six columns for each RM strategy; no inference/plots here.
source(file.path(Sys.getenv('M01_CODE_DIR', '.'),
                 'Function_07_Guideline_strategies_M01.R'), local=FALSE)
suppressPackageStartupMessages(library(data.table))

m12a_version <- function() '12A_fixed_S0_transition_matrix_v1.1'
m12a_scenarios <- function() c(SN='SH', ST='ST', SNT='SHT')
m12a_levels <- function() m07_transitions()
m12a_group_codes <- function() setNames(paste0('G',seq_along(m12a_levels())),m12a_levels())
m12a_board <- function(cfg) {
  p <- file.path(cfg$out_dir,'treatment_transition_matrix_S0_vs_RM')
  dir.create(p,recursive=TRUE,showWarnings=FALSE)
  p
}
m12a_check_workers <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv('M01_WORKERS','3')))
  if(length(n)!=1L || is.na(n) || n<1L)stop('M01_WORKERS must be a positive integer (default 3)')
  n
}
m12a_labels <- function() {
  data.table(order=seq_along(m12a_levels()),group=m12a_levels(),
             code=unname(m12a_group_codes()),
             label=c('Unchanged untreated','Newly intervention',
                     'Intervention removed','Treatment intensified',
                     'Unchanged intensity','Treatment deintensified'))
}
m12a_signature <- function(cfg,m) {
  files <- file.path(cfg$code_dir,c('Function_07_Guideline_strategies_M01.R',
                                     'Function_12A_Transition_matrix_M01.R'))
  input <- file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
  config <- file.path(cfg$out_dir,'study_config.rds')
  dep <- c(input,config,files)
  if(any(!file.exists(dep)))stop('12A missing dependencies: ',paste(dep[!file.exists(dep)],collapse=', '))
  list(version=m12a_version(),inputs=setNames(vapply(dep,m07_md5,character(1L)),basename(dep)))
}
m12a_assignment_check <- function(x,m) {
  x <- as.data.table(x)
  m07_assert(x,c('participant_id',paste0('transition_',m07_rms())),
             sprintf('12A MI %02d input',m))
  m07_unique(x)
  lev <- m12a_levels()
  for(k in m07_rms()) {
    s <- as.character(x[[paste0('transition_',k)]])
    if(anyNA(s)||any(!s%in%lev))stop('12A invalid/NA transition_',k,' in MI ',m)
  }
  if('imputation_id'%in%names(x) && any(x$imputation_id!=m))
    stop('12A MI indexing mismatch')
  invisible(TRUE)
}
m12a_one <- function(m,cfg) {
  x <- as.data.table(m07_read(cfg,m))
  m12a_assignment_check(x,m)
  n <- nrow(x)
  if(n<1L)stop('12A empty cohort')
  fixed <- as.character(x[['transition_S0']])
  lev <- m12a_levels(); lab <- m12a_labels()
  rows <- vector('list',length(m12a_scenarios()))
  checks <- vector('list',length(m12a_scenarios()))
  reductions <- vector('list',length(m12a_scenarios()))
  i <- 0L
  for(public in names(m12a_scenarios())) {
    i <- i+1L
    raw <- unname(m12a_scenarios()[[public]])
    current <- as.character(x[[paste0('transition_',raw)]])
    # CJ sorts (row order, column order); row-major flattening is exact.
    pos <- (match(fixed,lev)-1L)*6L+match(current,lev)
    if(anyNA(pos)||any(pos<1L|pos>36L))stop('12A invalid 6x6 index')
    mat <- CJ(row_order=seq_len(6L),column_order=seq_len(6L))
    mat[,n:=tabulate(pos,nbins=36L)]
    mat[,`:=`(imputation_id=m,scenario=public,scenario_internal=raw,
               n_total=nrow(x),row_group=lev[row_order],column_group=lev[column_order])]
    mat[,`:=`(row_n=sum(n),column_n=NA_integer_),by=row_order]
    mat[,column_n:=sum(n),by=column_order]
    mat[,`:=`(pct_total=100*n/n_total,
               pct_within_row=ifelse(row_n>0,100*n/row_n,NA_real_),
               pct_within_column=ifelse(column_n>0,100*n/column_n,NA_real_),
               diagonal=row_order==column_order)]
    mat[,`:=`(row_code=paste0('G',row_order),column_code=paste0('G',column_order))]
    rows[[i]] <- mat
    reductions[[i]] <- rbindlist(lapply(c("2019","2025"), function(year) {
      delta <- x[[paste0("D_",year,"_",raw)]] - x[[paste0("D_",year,"_S0")]]
      if(length(delta)!=n || any(!is.finite(delta))) stop("12A invalid paired LDL-C reduction")
      data.table(fixed_group=fixed, delta=delta)[,
        .(n=.N, mean_delta_mmol_l=mean(delta)), by=fixed_group][,
        `:=`(imputation_id=m, scenario=public, guideline=year)]
    }))
    row_actual <- tabulate(match(fixed,lev),nbins=6L)
    col_actual <- tabulate(match(current,lev),nbins=6L)
    row_found <- mat[,sum(n),by=row_order][order(row_order)]$V1
    col_found <- mat[,sum(n),by=column_order][order(column_order)]$V1
    pass <- identical(as.integer(row_actual),as.integer(row_found)) &&
      identical(as.integer(col_actual),as.integer(col_found)) &&
      sum(mat$n)==n && abs(sum(mat$pct_total)-100)<1e-8 && nrow(mat)==36L
    checks[[i]] <- data.table(imputation_id=m,scenario=public,n=n,
                              cells=nrow(mat),total_count=sum(mat$n),
                              total_pct=sum(mat$pct_total),
                              row_margins_equal=identical(as.integer(row_actual),as.integer(row_found)),
                              column_margins_equal=identical(as.integer(col_actual),as.integer(col_found)),
                              passed=pass)
    if(!pass)stop('12A matrix failed marginal QA: MI ',m,' scenario ',public)
    m07_log('12A','MI ',m,' ',public,' matrix 36/36 | diagonal ',
            sum(mat$n[mat$diagonal]),'/',n,' | ',round(100*sum(mat$n[mat$diagonal])/n,2),'%')
  }
  ans <- list(cells=rbindlist(rows),qa=rbindlist(checks),
              reductions=rbindlist(reductions),
              imputation_id=m,n=n,
              id_signature=m07_hash_object(as.character(x$participant_id)),
              fixed_signature=m07_hash_object(fixed),
              fixed_counts=tabulate(match(fixed,lev),nbins=6L))
  ans
}
m12a_task <- function(m,cfg) {
  root <- m12a_board(cfg); cp <- file.path(root,'checkpoints',sprintf('mi%02d.rds',m))
  signature_file <- paste0(cp,'.signature.rds')
  sig <- m12a_signature(cfg,m)
  if(file.exists(cp)&&file.exists(signature_file) &&
     identical(readRDS(signature_file),sig)) {
    z <- readRDS(cp)
    if(!isTRUE(all(z$qa$passed)))stop('12A cached QA failed')
    m07_log('12A','MI ',m,' cache hit')
    return(z)
  }
  z <- m07_timed('12A',sprintf('MI %02d three 6x6 matrices',m),m12a_one(m,cfg))
  m07_atomic(z,cp);m07_atomic(sig,signature_file)
  z
}
m12a_pooled_cells <- function(all) {
  d <- rbindlist(lapply(all,`[[`,'cells'))
  keys <- c('scenario','scenario_internal','row_order','column_order','row_group',
            'column_group','row_code','column_code','diagonal')
  z <- d[,.(mean_n=mean(n),sd_n=if(.N>1L)sd(n) else 0,
            min_n=min(n),max_n=max(n),
            mean_row_n=mean(row_n),mean_column_n=mean(column_n),
            mean_n_total=mean(n_total),
            pct_total=mean(pct_total),pct_within_row=if(all(is.na(pct_within_row)))NA_real_ else mean(pct_within_row,na.rm=TRUE),
            pct_within_column=if(all(is.na(pct_within_column)))NA_real_ else mean(pct_within_column,na.rm=TRUE),
            m=.N),by=keys]
  z[,scenario_order:=match(scenario,names(m12a_scenarios()))]
  setorderv(z,c('scenario_order','row_order','column_order'))
  z
}
m12a_plot_data <- function(pool,cfg) {
  p <- copy(pool)
  p[,`:=`(row_label=m12a_labels()$label[row_order],
           column_label=m12a_labels()$label[column_order],
           fill_value=pct_total,fill_unit='percent_of_total_cohort',
           display_n=as.integer(round(mean_n)),
           display_pct=pct_total,
           cell_label=sprintf('%s (%.1f%%)',formatC(round(mean_n),format='f',digits=0,big.mark=','),pct_total),
           axis_y='Score-based (S0) treatment-decision groups',
           axis_x='RM-based treatment-decision groups',
           figure_type='six_by_six_S0_to_RM')]
  if(any(p$m!=cfg$m)||nrow(p)!=108L)stop('12A: expected 108 pooled plot cells')
  p
}
m12a_row_summary <- function(d) {
  out <- d[,.(n_original=sum(n),n_unchanged=sum(n[diagonal]),
               n_moved=sum(n[!diagonal]),
               pct_stay=if(sum(n)>0)100*sum(n[diagonal])/sum(n) else NA_real_,
               pct_moved=if(sum(n)>0)100*sum(n[!diagonal])/sum(n) else NA_real_),
           by=.(imputation_id,scenario,row_order,row_code,row_group)]
  out
}
m12a_column_summary <- function(d) {
  d[,.(n_rm=sum(n),pct_of_cohort=100*sum(n)/unique(n_total)),
    by=.(imputation_id,scenario,column_order,column_code,column_group)]
}
m12a_export <- function(all,cfg) {
  if(length(all)!=cfg$m)stop('12A MI count mismatch')
  if(!all(vapply(all,function(z)all(z$qa$passed),logical(1))))stop('12A per-MI QA failed')
  if(length(unique(vapply(all,`[[`,character(1),'id_signature')))!=1L ||
     length(unique(vapply(all,`[[`,character(1),'fixed_signature')))!=1L)
    stop('12A participant IDs or S0 fixed group membership differ across imputations')
  root <- m12a_board(cfg)
  d <- rbindlist(lapply(all,`[[`,'cells'))
  if(anyDuplicated(d,by=c('imputation_id','scenario','row_order','column_order')))
    stop('12A duplicated matrix cells')
  if(nrow(d)!=108L*cfg$m)stop('12A MI cells missing')
  pool <- m12a_pooled_cells(all)
  plot <- m12a_plot_data(pool,cfg)
  rs <- m12a_row_summary(d)
  css <- m12a_column_summary(d)
  rs_pool <- rs[,.(mean_n_original=mean(n_original),mean_n_unchanged=mean(n_unchanged),
                    mean_n_moved=mean(n_moved),pct_stay=mean(pct_stay),
                    pct_moved=mean(pct_moved),m=.N),
                by=.(scenario,row_order,row_code,row_group)]
  cs_pool <- css[,.(mean_n_rm=mean(n_rm),pct_of_cohort=mean(pct_of_cohort),m=.N),
                 by=.(scenario,column_order,column_code,column_group)]
  if(any(abs(pool[,sum(mean_n),by=scenario]$V1-mean(vapply(all,`[[`,numeric(1),'n')))>1e-7))
    stop('12A pooled cell totals mismatch')
  for(s in names(m12a_scenarios())) {
    z <- plot[scenario==s]
    if(nrow(z)!=36L)stop('12A incomplete plot scenario: ',s)
    m07_csv(z,file.path(root,paste0('12A_plot_matrix_',s,'.csv')))
  }
  reductions <- rbindlist(lapply(all,`[[`,"reductions"))
  pooled_reductions <- reductions[, .(n=mean(n), mean_delta_mmol_l=mean(mean_delta_mmol_l),
                                      m=uniqueN(imputation_id)),
                                  by=.(fixed_group, scenario, guideline)]
  if (any(pooled_reductions$m != cfg$m) || nrow(pooled_reductions)!=36L)
    stop("12A incomplete fixed-group paired LDL-C results")
  m07_csv(reductions,file.path(root,"12A_fixed_LDL_change_per_imputation.csv"))
  m07_csv(pooled_reductions,file.path(root,"12A_fixed_LDL_change_MI_mean.csv"))
  m07_csv(d,file.path(root,'12A_matrix_cells_per_imputation.csv'))
  m07_csv(pool,file.path(root,'12A_matrix_cells_MI_mean.csv'))
  m07_csv(plot,file.path(root,'12A_plot_transition_matrix.csv'))
  m07_csv(rs,file.path(root,'12A_row_summary_per_imputation.csv'))
  m07_csv(rs_pool,file.path(root,'12A_plot_row_transition_summary.csv'))
  m07_csv(css,file.path(root,'12A_column_summary_per_imputation.csv'))
  m07_csv(cs_pool,file.path(root,'12A_plot_column_summary.csv'))
  m07_csv(rbindlist(lapply(all,`[[`,'qa')),file.path(root,'12A_matrix_quality_control.csv'))
  m07_csv(m12a_labels(),file.path(root,'12A_group_dictionary.csv'))
  manifests <- rbindlist(lapply(seq_along(all),function(i){
    m <- all[[i]]$imputation_id
    s <- m12a_signature(cfg,m)
    data.table(imputation_id=m,dependency=names(s$inputs),md5=unname(s$inputs),version=s$version)
  }))
  m07_csv(manifests,file.path(root,'12A_input_manifest.csv'))
  m07_log('12A','EXPORT DONE | 3 matrices x 36 cells; MI=',cfg$m,
          ' | output=',normalizePath(root,winslash='/',mustWork=FALSE))
  invisible(list(cells=pool,plot=plot,row=rs_pool))
}

# Independent PSOCK runner for stages 12A/12B: 3 processes by default,
# stage-specific elapsed logs, continuous parent-side progress, resumable MI caches.
# Do NOT call m07_worker(): it intentionally only dispatches numbered stages 08-12.
m12x_parallel <- function(indices,cfg,stage,module,worker_function,workers=3L) {
  if(!length(indices))return(list())
  workers <- min(as.integer(workers),length(indices))
  if(is.na(workers)||workers<1L)stop('Invalid workers')
  dispatch_serial <- function(m){
    m07_timed(stage,paste0('MI ',m,' ',stage),get(worker_function,mode='function')(m,cfg))
  }
  if(workers==1L)return(lapply(indices,dispatch_serial))
  logdir <- file.path(cfg$out_dir,'run_logs')
  dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
  taskdir <- tempfile(pattern=paste0('tasks_',stage,'_'),tmpdir=logdir)
  dir.create(taskdir,showWarnings=FALSE)
  on.exit(unlink(taskdir,recursive=TRUE,force=TRUE),add=TRUE)
  lib <- .libPaths()
  env <- c(M01_PROJECT_ROOT=Sys.getenv('M01_PROJECT_ROOT',cfg$root),
           M01_CODE_DIR=cfg$code_dir,M01_ANALYSIS_OUTPUT=cfg$out_dir,
           M01_CLEAN_DIR=cfg$clean_dir)
  cl <- parallel::makePSOCKcluster(workers,setup_timeout=60L)
  on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::clusterCall(cl,function(path,libs,envvars){
    .libPaths(libs)
    do.call(Sys.setenv,as.list(envvars))
    Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    source(path,local=FALSE)
    data.table::setDTthreads(1L)
    NULL
  },file.path(cfg$code_dir,module),lib,env)
  task <- function(m,cfg,stage,taskdir,worker_function){
    p <- file.path(taskdir,sprintf('mi%02d.rds',m))
    result <- tryCatch({
      value <- m07_timed(stage,paste0('MI ',m,' ',stage),
                         get(worker_function,mode='function')(m,cfg))
      list(ok=TRUE,value=value)
    },error=function(e)list(ok=FALSE,message=conditionMessage(e)))
    saveRDS(result,p,compress=FALSE)
    TRUE
  }
  environment(task) <- .GlobalEnv
  n <- length(indices); output <- vector('list',n)
  busy <- rep(FALSE,workers); slot <- integer(workers)
  next_j <- 1L; completed <- 0L
  submit <- function(w,j){
    parallel:::sendCall(cl[[w]],task,
      list(indices[j],cfg,stage,taskdir,worker_function),tag=j)
    busy[w] <<- TRUE; slot[w] <<- j
  }
  for(w in seq_len(workers)){
    if(next_j<=n){submit(w,next_j);next_j<-next_j+1L}
  }
  started <- Sys.time(); tick <- Sys.time()
  m07_log(stage,'Parallel MI workers=',workers,' | jobs=',n,' | stage=',stage)
  while(completed<n){
    for(w in which(busy)){
      j <- slot[w]; p <- file.path(taskdir,sprintf('mi%02d.rds',indices[j]))
      if(!file.exists(p))next
      z <- tryCatch(readRDS(p),error=function(e)NULL)
      if(is.null(z))next
      parallel:::recvResult(cl[[w]])
      if(!isTRUE(z$ok))stop(stage,' MI ',indices[j],': ',z$message)
      output[[j]] <- z$value; completed <- completed+1L
      busy[w] <- FALSE
      m07_log(stage,'Completed MI ',indices[j],' | ',completed,'/',n,
              ' | elapsed ',round(as.numeric(difftime(Sys.time(),started,units='secs')),1),' s')
      if(next_j<=n){submit(w,next_j);next_j<-next_j+1L}
    }
    if(as.numeric(difftime(Sys.time(),tick,units='secs'))>=5){
      m07_log(stage,'In progress MI=',paste(indices[slot[busy]],collapse=','),
              ' | ',completed,'/',n,' completed | elapsed ',
              round(as.numeric(difftime(Sys.time(),started,units='secs')),1),' s')
      tick <- Sys.time()
    }
    if(completed<n)Sys.sleep(.25)
  }
  output
}

