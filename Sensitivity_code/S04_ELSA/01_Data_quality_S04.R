# S04 | 01 | Quality check, variable distributions and missingness
# Run: Rscript 01_Data_quality_S04.R
# Read-only audit of already cleaned data; no exclusions, recoding or imputation.
# ============================= COMMON FUNCTIONS (LOCAL TO THIS SCRIPT) =============================
S04_WORKERS <- suppressWarnings(as.integer(Sys.getenv("S04_WORKERS", "5")))
if (is.na(S04_WORKERS) || S04_WORKERS < 1L) stop("S04_WORKERS must be a positive integer")
S04_BATCH_SIZE <- suppressWarnings(as.integer(Sys.getenv("S04_BATCH_SIZE", "1000")))
if (is.na(S04_BATCH_SIZE) || S04_BATCH_SIZE < 1L) stop("S04_BATCH_SIZE must be positive")

script_path <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) return(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE))
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(e) e$ofile))
  if (length(frames)) return(normalizePath(tail(frames, 1L)[[1L]], mustWork=FALSE))
  normalizePath(getwd(), mustWork=FALSE)
}
script_dir <- if (dir.exists(script_path())) script_path() else dirname(script_path())
root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
input_dir <- Sys.getenv("S04_INPUT_DIR", file.path(root, "DATA", "ELSA"))
analysis_dir <- Sys.getenv("S04_ANALYSIS_DIR", file.path(root, "LJW_outcome", "data_intervention", "S04_ELSA", "complete_case"))
qc_dir <- Sys.getenv("S04_QC_DIR", file.path(root, "LJW_outcome", "data_clean", "S04_ELSA", "quality_control"))

elapsed_text <- function(seconds) sprintf("%02d:%02d:%05.2f", floor(seconds/3600), floor(seconds%%3600/60), seconds%%60)
start_clock <- Sys.time()
log_message <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep="")
progress <- function(done,total,label="Completed") {
  pct <- if (total) 100*done/total else 100
  log_message(sprintf("%s: %d/%d (%.1f%%); elapsed %s",label,done,total,pct,elapsed_text(as.numeric(difftime(Sys.time(),start_clock,units="secs")))))
}
finish <- function(script_label) log_message(script_label," finished; elapsed ",elapsed_text(as.numeric(difftime(Sys.time(),start_clock,units="secs"))))
write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  write.csv(x,path,row.names=FALSE,na="")
  log_message("Wrote ",path," [",nrow(x)," rows]")
}
read_elsa <- function() {
  dictionary_file <- file.path(input_dir,"ELSA_variable.csv")
  data_file <- file.path(input_dir,"ELSA_data.rds")
  for (p in c(dictionary_file,data_file)) if (!file.exists(p)) stop("Missing ELSA input: ",p)
  dictionary <- read.csv(dictionary_file,fileEncoding="UTF-8-BOM",check.names=FALSE,stringsAsFactors=FALSE)
  needed <- c("output_column_name","source_field")
  if (!all(needed %in% names(dictionary))) stop("Dictionary missing: ",paste(setdiff(needed,names(dictionary)),collapse=", "))
  if (anyNA(dictionary[needed]) || anyDuplicated(dictionary$output_column_name) || anyDuplicated(dictionary$source_field)) stop("Dictionary has missing or duplicate mappings")
  d <- as.data.frame(readRDS(data_file))
  if (!all(names(d) %in% dictionary$output_column_name)) stop("Unmapped ELSA data fields: ",paste(setdiff(names(d),dictionary$output_column_name),collapse=", "))
  names(d) <- dictionary$source_field[match(names(d),dictionary$output_column_name)]
  names(d)[names(d)=="__derived_high_bp__"] <- "high_bp"
  if (!"idauniq" %in% names(d)) stop("Missing participant ID idauniq")
  d$idauniq <- as.character(d$idauniq)
  if (anyNA(d$idauniq) || any(!nzchar(d$idauniq)) || anyDuplicated(d$idauniq)) stop("Missing or duplicate participant identifiers")
  list(data=d, dictionary=dictionary)
}
write_qc <- function(d, name) write_csv(d,file.path(qc_dir,paste0("01_",name,"_S04.csv")))
log_message("Step 01 starting (configured workers=",S04_WORKERS,"). No need for parallelism in small QC tasks.")
loaded <- read_elsa()
d <- loaded$data
log_message("Loaded ",nrow(d)," participants and ",ncol(d)," fields; workers setting = ",S04_WORKERS)
# Analyze every column actually supplied; group by the original variable dictionary.
dictionary <- loaded$dictionary
names(dictionary)[names(dictionary) == 'source_field'] <- 'source_name'
dictionary$variable <- dictionary$source_name
dictionary$variable[dictionary$variable == '__derived_high_bp__'] <- 'high_bp'
# Wave 4 hs-CRP and persistence are outside the approved S04 analysis scope.
cols <- setdiff(names(d), c("hscrp_w4", "w4bldwt", "crp_persistent_w4w6"))
meta <- dictionary[match(cols, dictionary$variable), , drop = FALSE]
stopifnot(all(meta$variable == cols))

# Numeric variables with continuous interpretations: do not treat binary proxies as continuous.
continuous <- c('age','sbp','dbp','tc','hdl','ldl','trig','hba1c','bmi','waist',
                'hscrp_w6','w6nurwt','w6bldwt','w4bldwt')
summary_rows <- lapply(seq_along(cols), function(i) {
  if (i==1L || i%%5L==0L || i==length(cols)) progress(i,length(cols),"Variable summaries")
  v <- cols[i]
  x <- d[[v]]
  applicable <- if (identical(v,'rm_natural_menopause_lt40_proxy')) {
    if (!('sex' %in% names(d))) stop('sex is required for menopause applicability')
    d$sex %in% c('female', 'Female', 2L, '2')
  } else rep(TRUE,nrow(d))
  x <- x[applicable]
  nobs <- sum(!is.na(x))
  is_numeric <- v %in% continuous
  if (is_numeric && !is.numeric(x)) {
    stop('Expected numeric variable: ',v,'; got ',class(x)[1L])
  }
  obs <- x[!is.na(x)]
  q <- if (is_numeric && nobs) as.numeric(quantile(obs, c(.25,.5,.75))) else rep(NA_real_,3L)
  data.frame(group = meta$panel[i], variable = v, label = meta$output_column_name[i],
             n_total = nrow(d), n_applicable = length(x),
             n_not_applicable = nrow(d)-length(x), n_observed = nobs,
             n_missing = sum(is.na(x)),
             missing_pct = if (length(x)) round(mean(is.na(x))*100,3) else NA_real_,
             n_distinct_observed = length(unique(obs)),
             min = if (is_numeric && nobs) min(obs) else NA_real_,
             max = if (is_numeric && nobs) max(obs) else NA_real_,
             mean = if (is_numeric && nobs) mean(obs) else NA_real_,
             sd = if (is_numeric && nobs>1L) sd(obs) else NA_real_,
             q25=q[1L], median=q[2L], q75=q[3L],
             stringsAsFactors=FALSE)
})
summary <- do.call(rbind, summary_rows)
stopifnot(all(summary$n_observed + summary$n_missing + summary$n_not_applicable == nrow(d)))
write_qc(summary,'distribution_missingness')

# All observed category levels, without assuming that every factor uses 0/1.
cat_vars <- setdiff(cols, c(continuous, 'idauniq','baseline_date'))
category_list <- lapply(cat_vars, function(v) {
  applicable <- if (v == 'rm_natural_menopause_lt40_proxy') d$sex %in% c('female','Female',2L,'2') else rep(TRUE,nrow(d))
  x <- as.character(d[[v]][applicable])
  x <- x[!is.na(x)]
  tab <- table(x)
  if (!length(tab)) return(data.frame(variable=v, category=NA_character_,n_observed=0L,n=0L,pct_observed=NA_real_))
  data.frame(variable=v,category=names(tab),n_observed=length(x),n=as.integer(tab),
             pct_observed=100*as.integer(tab)/length(x),stringsAsFactors=FALSE)
})
write_qc(do.call(rbind,category_list),'categorical_distribution')

# Data validity findings are reports only: never overwrite the supplied data.
checks <- list(
  missing_identifier = sum(is.na(d$idauniq)|!nzchar(d$idauniq)),
  duplicate_identifier = anyDuplicated(d$idauniq),
  age_outside_50_69 = if ('age'%in%cols) sum(!is.na(d$age)&(d$age<50|d$age>69)) else NA_integer_,
  sex_missing = if ('sex'%in%cols) sum(is.na(d$sex)) else NA_integer_,
  smoking_missing = if ('current_smoking'%in%cols) sum(is.na(d$current_smoking)) else NA_integer_
)
write_qc(data.frame(check=names(checks), n=unlist(checks)),'input_checks')

# Wave 6 hs-CRP only. Wave 4 and persistence are outside S04 RM scope.
if ('hscrp_w6' %in% cols) {
  x <- d$hscrp_w6
  observed <- !is.na(x)
  write_qc(data.frame(variable='hscrp_w6', n_observed=sum(observed),
                      n_missing=sum(!observed),
                      n_above_2=sum(x[observed] > 2),
                      positive_pct_observed=if(sum(observed))100*mean(x[observed]>2)else NA_real_),
           'hscrp_wave6_only')
}

# Complete-case availability counts; no exclusions are applied to the QC dataset.
score_fields <- c('age','sex','current_smoking','sbp','dbp','tc','hdl','ldl')
missing_fields <- setdiff(score_fields,cols)
if (length(missing_fields)) stop('Missing essential SCORE/treatment fields: ',paste(missing_fields,collapse=', '))
cc <- complete.cases(d[,score_fields,drop=FALSE])
write_qc(data.frame(item=c('source_rows','complete_scoring_inputs','incomplete_scoring_inputs'),
                    n=c(nrow(d),sum(cc),sum(!cc))), 'analysis_readiness')
message('S04 QC complete: ',nrow(d),' records, ',sum(cc),' complete SCORE/treatment inputs.')

finish("01 data quality")
