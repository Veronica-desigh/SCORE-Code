# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_05_Eligibility_format_M01.R
# 任务：供05调用的变量字典、随访、QA
# 输入：04合格队列和变量字典
# 输出：分析基础数据和通用QA
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function_05_Eligibility_format_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   为 05_Eligibility_format_M01.R 提供分析基础数据整理函数。
#
# 设计原则：
#   1. 04 已输出一人一行的最终 analytic cohort；05 只核验，不去重、不改变样本量。
#   2. 05 统一变量名/取值，同时明确变量角色：required observed、to be imputed、
#      post-MI derived、outcome-model/subgroup/strategy use。
#   3. 05 不拟合最终 SES、不填补 physical activity，也不派生最终 hypertension；
#      这些变量在 06 的每套 completed dataset 中重新生成。
#   4. ASCVD/竞争死亡结局直接复用 02 的冻结定义，并构造 cause-specific
#      Nelson-Aalen cumulative hazards 作为 MI auxiliary predictors。
#   5. 普通分析协变量允许缺失；只有 required_observed 变量在 05 中必须完整。
# ============================================================================== #

library(data.table)

# ============================================================================== #
# 0. 通用检查、计时和格式化函数 ####
# ============================================================================== #

fmt_time <- function(sec) {
  if (is.null(sec) || is.na(sec)) {
    return(NA_character_)
  }
  if (sec < 60) {
    return(sprintf("%.2f sec", sec))
  }
  if (sec < 3600) {
    return(sprintf("%.2f min", sec / 60))
  }
  sprintf("%.2f hr", sec / 3600)
}

step_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

run_with_timer <- function(expr, label = NULL) {
  t0 <- Sys.time()
  if (!is.null(label)) {
    step_message("开始：", label)
  }

  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(label)) {
    step_message("完成：", label, " | 耗时 = ", fmt_time(elapsed))
  }

  if (exists("m01_record_timing", mode = "function", inherits = TRUE))
    m01_record_timing(if (is.null(label)) "unnamed_run_with_timer" else label, elapsed)
  list(value = value, elapsed_sec = elapsed)
}

assert_file_exists <- function(path, label = "文件") {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(label, "不存在：", path, call. = FALSE)
  }
  invisible(path)
}

assert_required_columns <- function(dt, cols, label = "数据表") {
  dt <- as.data.table(dt)
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  missing_cols <- setdiff(cols, names(dt))

  if (length(missing_cols) > 0L) {
    stop(
      label, "缺少关键列：", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_unique_id <- function(dt, id_col, label = "数据表") {
  dt <- as.data.table(dt)
  assert_required_columns(dt, id_col, label)

  n_missing <- dt[, sum(is.na(get(id_col)) | trimws(as.character(get(id_col))) == "")]
  n_duplicate_rows <- dt[, .N, by = id_col][N > 1L, sum(N)]
  if (length(n_duplicate_rows) == 0L || is.na(n_duplicate_rows)) {
    n_duplicate_rows <- 0L
  }

  if (n_missing > 0L || n_duplicate_rows > 0L) {
    stop(
      label, "的 ID 检查失败：缺失 ID = ", n_missing,
      "，重复记录行数 = ", n_duplicate_rows,
      call. = FALSE
    )
  }

  invisible(TRUE)
}

clean_character_missing <- function(x) {
  x <- trimws(as.character(x))
  missing_strings <- c(
    "", "NA", "N/A", "NaN", "NULL",
    "Do not know", "Prefer not to answer"
  )
  x[is.na(x) | x %in% missing_strings] <- NA_character_
  x
}

to_snake_value <- function(x) {
  x <- clean_character_missing(x)
  x <- tolower(x)
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- NA_character_
  x
}

typed_na <- function(type, n) {
  switch(
    type,
    date = as.IDate(rep(NA_character_, n)),
    integer = rep(NA_integer_, n),
    numeric = rep(NA_real_, n),
    character = rep(NA_character_, n),
    logical = rep(NA, n),
    rep(NA_character_, n)
  )
}

coerce_declared_type <- function(x, type) {
  switch(
    type,
    date = as.IDate(x),
    integer = suppressWarnings(as.integer(x)),
    numeric = suppressWarnings(as.numeric(x)),
    character = clean_character_missing(x),
    logical = as.logical(x),
    x
  )
}

# ============================================================================== #
# 1. 读取、纳排一致性审计和变量字典 ####
# ============================================================================== #

read_eligible_analysis_base <- function(path, id_col = "Participant_ID") {
  assert_file_exists(path, "04 最终合格人群数据")

  dt <- fread(path)
  setDT(dt)
  assert_unique_id(dt, id_col, "04 最终合格人群数据")
  dt[]
}

audit_expected_values <- function(dt, expected_values, stop_on_failure = TRUE) {
  dt <- as.data.table(dt)

  if (is.null(names(expected_values)) || any(names(expected_values) == "")) {
    stop("expected_values 必须是带变量名的向量或列表。", call. = FALSE)
  }

  assert_required_columns(dt, names(expected_values), "纳排一致性审计数据")

  out <- rbindlist(lapply(names(expected_values), function(v) {
    expected <- expected_values[[v]]
    x <- dt[[v]]
    x_chr <- as.character(x)
    expected_chr <- as.character(expected)
    observed <- sort(unique(x_chr[!is.na(x_chr)]))

    data.table(
      variable = v,
      expected_value = expected_chr,
      n_total = length(x),
      n_missing = sum(is.na(x)),
      n_unexpected = sum(!is.na(x_chr) & x_chr != expected_chr),
      observed_values = if (length(observed) == 0L) NA_character_ else paste(observed, collapse = "; "),
      passed = sum(is.na(x)) == 0L && sum(!is.na(x_chr) & x_chr != expected_chr) == 0L
    )
  }))

  if (stop_on_failure && any(!out$passed)) {
    failed <- out[passed == FALSE, variable]
    stop(
      "04 纳排结果一致性检查失败：", paste(failed, collapse = ", "),
      "。请先检查 04 输出，05 不会自动删除这些记录。",
      call. = FALSE
    )
  }

  out[]
}

validate_variable_dictionary <- function(dictionary) {
  dictionary <- as.data.table(dictionary)
  required_dictionary_cols <- c(
    "source_name", "standard_name", "label", "unit", "role", "type",
    "required_column", "required_observed", "impute", "imputation_method",
    "derive_after_mi", "used_for_eligibility", "used_for_strategy",
    "used_for_outcome_model", "used_for_subgroup"
  )
  assert_required_columns(dictionary, required_dictionary_cols, "05 变量字典")

  if (anyDuplicated(dictionary$source_name)) {
    stop("变量字典中 source_name 存在重复。", call. = FALSE)
  }
  if (anyDuplicated(dictionary$standard_name)) {
    stop("变量字典中 standard_name 存在重复。", call. = FALSE)
  }

  valid_name <- grepl("^[a-z][a-z0-9_]*$", dictionary$standard_name)
  if (any(!valid_name)) {
    stop(
      "以下标准变量名不符合 lower_snake_case：",
      paste(dictionary[!valid_name, standard_name], collapse = ", "),
      call. = FALSE
    )
  }

  valid_types <- c("date", "integer", "numeric", "character", "logical")
  invalid_types <- setdiff(unique(dictionary$type), valid_types)
  if (length(invalid_types) > 0L) {
    stop("变量字典包含未支持的类型：", paste(invalid_types, collapse = ", "), call. = FALSE)
  }

  if (any(dictionary$required_observed & dictionary$impute)) {
    bad <- dictionary[required_observed == TRUE & impute == TRUE, standard_name]
    stop(
      "变量不能同时 required_observed=TRUE 且 impute=TRUE：",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(dictionary$derive_after_mi & dictionary$impute)) {
    bad <- dictionary[derive_after_mi == TRUE & impute == TRUE, standard_name]
    stop(
      "派生变量不能直接插补；应插补其组成变量：",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

new_variable_spec <- function(source_name,
                              standard_name,
                              label,
                              unit = NA_character_,
                              role,
                              type,
                              required_column = TRUE,
                              required_observed = FALSE,
                              impute = FALSE,
                              imputation_method = NA_character_,
                              derive_after_mi = FALSE,
                              used_for_eligibility = FALSE,
                              used_for_strategy = FALSE,
                              used_for_outcome_model = FALSE,
                              used_for_subgroup = FALSE) {
  data.table(
    source_name = source_name,
    standard_name = standard_name,
    label = label,
    unit = unit,
    role = role,
    type = type,
    required_column = as.logical(required_column),
    required_observed = as.logical(required_observed),
    impute = as.logical(impute),
    imputation_method = as.character(imputation_method),
    derive_after_mi = as.logical(derive_after_mi),
    used_for_eligibility = as.logical(used_for_eligibility),
    used_for_strategy = as.logical(used_for_strategy),
    used_for_outcome_model = as.logical(used_for_outcome_model),
    used_for_subgroup = as.logical(used_for_subgroup)
  )
}

apply_variable_dictionary <- function(dt, dictionary) {
  dt <- as.data.table(dt)
  dictionary <- copy(as.data.table(dictionary))
  validate_variable_dictionary(dictionary)

  required_sources <- dictionary[required_column == TRUE, source_name]
  assert_required_columns(dt, required_sources, "04 最终合格人群数据")

  out <- data.table(row_index___ = seq_len(nrow(dt)))

  for (i in seq_len(nrow(dictionary))) {
    src <- dictionary$source_name[i]
    dst <- dictionary$standard_name[i]
    declared_type <- dictionary$type[i]

    if (src %in% names(dt)) {
      out[, (dst) := coerce_declared_type(dt[[src]], declared_type)]
    } else {
      out[, (dst) := typed_na(declared_type, nrow(dt))]
      warning("可选变量来源列不存在：", src, "；已生成缺失变量 ", dst, call. = FALSE)
    }
  }

  out[, row_index___ := NULL]

  required_vars <- dictionary[required_observed == TRUE, standard_name]
  missing_required <- vapply(required_vars, function(v) sum(is.na(out[[v]])), integer(1))
  if (any(missing_required > 0L)) {
    details <- paste0(
      names(missing_required)[missing_required > 0L], "=",
      missing_required[missing_required > 0L]
    )
    stop(
      "04 required-information gate 后仍存在 required_observed 缺失：",
      paste(details, collapse = ", "),
      call. = FALSE
    )
  }

  out[]
}

standardize_analysis_values <- function(dt) {
  dt <- copy(as.data.table(dt))

  if ("sex" %in% names(dt)) {
    dt[, sex := fcase(
      tolower(sex) == "female", "female",
      tolower(sex) == "male", "male",
      default = NA_character_
    )]
  }

  if ("ethnicity" %in% names(dt)) {
    dt[, ethnicity := fcase(
      tolower(ethnicity) == "white", "white",
      tolower(ethnicity) %in% c("non-white", "non_white"), "non_white",
      default = to_snake_value(ethnicity)
    )]
  }

  if ("education" %in% names(dt)) {
    dt[, education := fcase(
      tolower(education) == "low education", "low",
      tolower(education) == "medium education", "medium",
      tolower(education) == "high education", "high",
      default = NA_character_
    )]
  }

  if ("household_income" %in% names(dt)) {
    dt[, household_income := fcase(
      household_income == "Less than 18,000", "less_than_18000",
      household_income == "18,000 to 30,999", "18000_to_30999",
      household_income == "31,000 to 51,999", "31000_to_51999",
      household_income == "52,000 to 100,000", "52000_to_100000",
      household_income == "Greater than 100,000", "greater_than_100000",
      default = NA_character_
    )]
  }

  if ("employment_status" %in% names(dt)) {
    dt[, employment_status := fcase(
      employment_status == "Other not employed", "other_not_employed",
      employment_status == "Retired", "retired",
      employment_status == "Employed/self-employed", "employed_or_self_employed",
      default = NA_character_
    )]
  }

  if ("smoking_status" %in% names(dt)) {
    dt[, smoking_status := fcase(
      tolower(smoking_status) == "never", "never",
      tolower(smoking_status) %in% c("previous", "former"), "former",
      tolower(smoking_status) == "current", "current",
      default = NA_character_
    )]
  }

  if ("region" %in% names(dt)) {
    region_clean <- tolower(clean_character_missing(dt$region))
    dt[, region := fcase(
      region_clean == "england", "England",
      region_clean == "wales", "Wales",
      region_clean == "scotland", "Scotland",
      default = NA_character_
    )]
    invalid_region <- !is.na(region_clean) & is.na(dt$region)
    if (any(invalid_region)) {
      stop(
        "region 存在无法识别的非缺失水平；允许真实缺失进入 MI，但不允许未知编码。",
        call. = FALSE
      )
    }
  }

  if ("physical_activity" %in% names(dt)) {
    dt[, physical_activity := to_snake_value(physical_activity)]
  }
  if ("censor_reason" %in% names(dt)) {
    dt[, censor_reason := to_snake_value(censor_reason)]
  }

  dt[]
}

# ============================================================================== #
# 2. 二分类 LCA 个体层面 SES ####
# ============================================================================== #

derive_binary_ses_lca <- function(data,
                                  id_var = "participant_id",
                                  income_var = "household_income",
                                  education_var = "education",
                                  employment_var = "employment_status",
                                  nrep = 30L,
                                  maxiter = 5000L,
                                  seed = 123L) {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    stop(
      "二分类 SES LCA 需要 R 包 'poLCA'。请先安装后再运行 06 post-MI derivation：install.packages('poLCA')。",
      call. = FALSE
    )
  }

  dt <- copy(as.data.table(data))
  required_vars <- c(id_var, income_var, education_var, employment_var)
  assert_required_columns(dt, required_vars, "SES LCA 输入数据")
  assert_unique_id(dt, id_var, "SES LCA 输入数据")

  income_levels <- c(
    "less_than_18000",
    "18000_to_30999",
    "31000_to_51999",
    "52000_to_100000",
    "greater_than_100000"
  )
  education_levels <- c("low", "medium", "high")
  employment_levels <- c("other_not_employed", "retired", "employed_or_self_employed")

  lca_dt <- dt[, .(
    participant_id___ = get(id_var),
    income_code = match(get(income_var), income_levels),
    education_code = match(get(education_var), education_levels),
    employment_code = match(get(employment_var), employment_levels)
  )]

  if (any(!complete.cases(lca_dt[, .(income_code, education_code, employment_code)]))) {
    stop(
      "SES 组成变量存在缺失或未识别类别。该函数只应在 06 completed MI dataset 上调用；请检查插补和变量标准化。",
      call. = FALSE
    )
  }

  set.seed(seed)
  lca_formula <- cbind(income_code, education_code, employment_code) ~ 1
  lca_fit <- poLCA::poLCA(
    lca_formula,
    data = as.data.frame(lca_dt[, .(income_code, education_code, employment_code)]),
    nclass = 2L,
    nrep = as.integer(nrep),
    maxiter = as.integer(maxiter),
    verbose = FALSE,
    graphs = FALSE,
    na.rm = FALSE
  )

  ordinal_score <-
    (lca_dt$income_code - 1) / 4 +
    (lca_dt$education_code - 1) / 2 +
    (lca_dt$employment_code - 1) / 2

  posterior <- lca_fit$posterior
  class_mean_score <- colSums(posterior * ordinal_score) / colSums(posterior)
  ordered_classes <- order(class_mean_score)

  class_label_map <- rep(NA_character_, 2L)
  class_label_map[ordered_classes] <- c("low", "high")
  posterior_ordered <- posterior[, ordered_classes, drop = FALSE]
  max_posterior <- apply(posterior_ordered, 1L, max)

  ses_by_id <- data.table(
    participant_id___ = lca_dt$participant_id___,
    ses_class = class_label_map[lca_fit$predclass],
    ses_probability_low = posterior_ordered[, 1L],
    ses_probability_high = posterior_ordered[, 2L]
  )
  setnames(ses_by_id, "participant_id___", id_var)

  posterior_entropy <- -sum(posterior_ordered * log(pmax(posterior_ordered, .Machine$double.eps)))
  relative_entropy <- 1 - posterior_entropy / (nrow(posterior_ordered) * log(2))

  safe_scalar <- function(x) {
    if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x[1L])
  }

  n_iterations <- safe_scalar(lca_fit$numiter)
  model_summary <- data.table(
    metric = c(
      "n_participants", "n_classes", "nrep", "maxiter", "seed",
      "log_likelihood", "aic", "bic", "g_squared", "chi_squared",
      "n_iterations_best_model", "converged_before_maxiter",
      "relative_entropy", "mean_max_posterior_probability",
      "proportion_max_posterior_below_0_70",
      "low_class_mean_ordinal_score", "high_class_mean_ordinal_score"
    ),
    value = c(
      nrow(lca_dt), 2, nrep, maxiter, seed,
      safe_scalar(lca_fit$llik), safe_scalar(lca_fit$aic), safe_scalar(lca_fit$bic),
      safe_scalar(lca_fit$Gsq), safe_scalar(lca_fit$Chisq),
      n_iterations, as.numeric(!is.na(n_iterations) && n_iterations < maxiter),
      relative_entropy, mean(max_posterior), mean(max_posterior < 0.70),
      class_mean_score[ordered_classes[1L]], class_mean_score[ordered_classes[2L]]
    )
  )

  class_distribution <- ses_by_id[, .(
    n = .N,
    percent = 100 * .N / nrow(ses_by_id),
    mean_assigned_probability = if (ses_class[1L] == "low") {
      mean(ses_probability_low)
    } else {
      mean(ses_probability_high)
    }
  ), by = ses_class][order(match(ses_class, c("low", "high")))]

  list(
    ses_by_id = ses_by_id[],
    model_summary = model_summary[],
    class_distribution = class_distribution[],
    model = lca_fit
  )
}

# ============================================================================== #
# 3. ASCVD 结局与随访接口 ####
# ============================================================================== #

standardize_primary_ascvd_followup <- function(dt) {
  dt <- copy(as.data.table(dt))
  required <- c(
    "baseline_date", "followup_end_date", "event_status", "ascvd_event",
    "ascvd_source_date", "competing_death_source_date", "censor_date",
    "source_followup_years"
  )
  assert_required_columns(dt, required, "ASCVD 随访接口数据")

  date_vars <- intersect(
    c(
      "baseline_date", "followup_end_date", "ascvd_source_date",
      "competing_death_source_date", "censor_date", "loss_to_followup_date",
      "administrative_censor_date"
    ),
    names(dt)
  )
  dt[, (date_vars) := lapply(.SD, as.IDate), .SDcols = date_vars]

  dt[, event_status := as.integer(event_status)]
  invalid_status <- setdiff(unique(dt$event_status[!is.na(dt$event_status)]), 0:2)
  if (length(invalid_status) > 0L || anyNA(dt$event_status)) {
    stop("event_status 必须完整且只能取 0、1、2。", call. = FALSE)
  }

  source_event <- as.integer(dt$ascvd_event)
  expected_event <- as.integer(dt$event_status == 1L)
  if (anyNA(source_event) || any(source_event != expected_event)) {
    stop("primary_ascvd_event 与 primary_ascvd_status 不一致。", call. = FALSE)
  }

  dt[, `:=`(
    ascvd_event = as.integer(event_status == 1L),
    competing_death_event = as.integer(event_status == 2L),
    ascvd_event_date = as.IDate(NA),
    competing_death_date = as.IDate(NA)
  )]
  dt[event_status == 1L, ascvd_event_date := ascvd_source_date]
  dt[event_status == 2L, competing_death_date := competing_death_source_date]

  dt[, followup_days := as.numeric(followup_end_date - baseline_date)]
  dt[, followup_years := followup_days / 365.25]

  if (anyNA(dt$baseline_date) || anyNA(dt$followup_end_date)) {
    stop("baseline_date 或 followup_end_date 存在缺失。", call. = FALSE)
  }
  if (any(dt$followup_days < 0)) {
    stop("存在 followup_end_date 早于 baseline_date 的记录。", call. = FALSE)
  }
  if (any(dt$followup_end_date > dt$censor_date, na.rm = TRUE)) {
    stop("存在 followup_end_date 晚于 censor_date 的记录。", call. = FALSE)
  }
  if (any(dt[event_status == 0L, followup_end_date != censor_date], na.rm = TRUE)) {
    stop("删失者的 followup_end_date 与 censor_date 不一致。", call. = FALSE)
  }
  if (any(dt[event_status == 1L, followup_end_date != ascvd_source_date], na.rm = TRUE)) {
    stop("ASCVD 事件者的随访终止日期与 primary_ascvd_time 不一致。", call. = FALSE)
  }
  if (any(dt[event_status == 2L, followup_end_date != competing_death_source_date], na.rm = TRUE)) {
    stop("竞争死亡者的随访终止日期与 competing_death_time 不一致。", call. = FALSE)
  }
  if (any(dt$event_status == 1L & is.na(dt$ascvd_source_date))) {
    stop("ASCVD 事件者缺少 primary_ascvd_time。", call. = FALSE)
  }
  if (any(dt$event_status == 2L & is.na(dt$competing_death_source_date))) {
    stop("竞争死亡者缺少 competing_death_time。", call. = FALSE)
  }

  diff_years <- abs(dt$followup_years - dt$source_followup_years)
  if (any(diff_years > 1e-8, na.rm = TRUE)) {
    stop("重新计算的 followup_years 与 02 输出不一致。", call. = FALSE)
  }

  dt[, c("ascvd_source_date", "competing_death_source_date", "source_followup_years") := NULL]
  dt[]
}


# ============================================================================== #
# 4. MI auxiliary predictors：cause-specific Nelson-Aalen cumulative hazards ####
# ============================================================================== #

# 对每位参与者在其 observed follow-up time 处评估 marginal cause-specific
# Nelson-Aalen cumulative hazard。两条 hazard 仅作为 MICE auxiliary predictors，
# 不替代后续正式 cause-specific Cox / g-computation outcome model。
add_mi_nelson_aalen_predictors <- function(dt) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("05 需要 survival 包来计算 MI Nelson-Aalen auxiliary predictors。", call. = FALSE)
  }

  out <- copy(as.data.table(dt))
  assert_required_columns(
    out,
    c("followup_years", "ascvd_event", "competing_death_event"),
    "MI Nelson-Aalen 输入"
  )

  if (anyNA(out$followup_years) || any(!is.finite(out$followup_years)) ||
      any(out$followup_years < 0)) {
    stop("followup_years 必须完整、有限且 >=0，才能构造 Nelson-Aalen predictors。", call. = FALSE)
  }

  cumulative_hazard_at_time <- function(time, event) {
    event <- as.integer(event)
    if (anyNA(event) || any(!event %in% 0:1)) {
      stop("cause-specific event indicator 必须完整且为 0/1。", call. = FALSE)
    }
    if (sum(event) == 0L) return(rep(0, length(time)))

    fit <- survival::coxph(
      survival::Surv(time, event) ~ 1,
      ties = "breslow",
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    bh <- survival::basehaz(fit, centered = FALSE)
    if (nrow(bh) == 0L) return(rep(0, length(time)))

    idx <- findInterval(time, bh$time)
    ans <- numeric(length(time))
    positive <- idx > 0L
    ans[positive] <- bh$hazard[idx[positive]]
    ans
  }

  out[, na_hazard_ascvd := cumulative_hazard_at_time(followup_years, ascvd_event)]
  out[, na_hazard_death :=
        cumulative_hazard_at_time(followup_years, competing_death_event)]

  if (anyNA(out$na_hazard_ascvd) || anyNA(out$na_hazard_death) ||
      any(out$na_hazard_ascvd < 0) || any(out$na_hazard_death < 0)) {
    stop("Nelson-Aalen auxiliary predictor 生成失败。", call. = FALSE)
  }

  out[]
}


# ============================================================================== #
# 5. 随访、结局和缺失情况描述 ####
# ============================================================================== #

estimate_reverse_km_followup <- function(dt) {
  dt <- as.data.table(dt)
  assert_required_columns(dt, c("followup_years", "event_status"), "反向 Kaplan-Meier 输入数据")

  if (!requireNamespace("survival", quietly = TRUE)) {
    warning("未安装 survival 包，反向 Kaplan-Meier 随访时间未计算。", call. = FALSE)
    return(data.table(
      metric = c("reverse_km_median_years", "reverse_km_lcl_years", "reverse_km_ucl_years"),
      value = NA_real_
    ))
  }

  reverse_event <- as.integer(dt$event_status == 0L)
  fit <- survival::survfit(survival::Surv(dt$followup_years, reverse_event) ~ 1)
  fit_table <- summary(fit)$table

  get_fit_value <- function(name) {
    if (name %in% names(fit_table)) as.numeric(fit_table[[name]]) else NA_real_
  }

  data.table(
    metric = c("reverse_km_median_years", "reverse_km_lcl_years", "reverse_km_ucl_years"),
    value = c(
      get_fit_value("median"),
      get_fit_value("0.95LCL"),
      get_fit_value("0.95UCL")
    )
  )
}

build_followup_outcome_summary <- function(dt) {
  dt <- as.data.table(dt)
  required <- c(
    "participant_id", "baseline_date", "followup_end_date", "followup_years",
    "event_status", "ascvd_event", "competing_death_event", "loss_to_followup_event"
  )
  assert_required_columns(dt, required, "随访和结局描述数据")

  n_total <- nrow(dt)
  total_person_years <- sum(dt$followup_years)
  n_ascvd <- sum(dt$ascvd_event == 1L)
  n_competing <- sum(dt$competing_death_event == 1L)
  n_censored <- sum(dt$event_status == 0L)
  n_lost <- sum(dt$loss_to_followup_event == 1L, na.rm = TRUE)
  followup_quantiles <- quantile(
    dt$followup_years,
    probs = c(0, 0.25, 0.50, 0.75, 1),
    na.rm = TRUE,
    names = FALSE
  )

  incidence <- if (total_person_years > 0) {
    stats::poisson.test(n_ascvd, T = total_person_years)
  } else {
    NULL
  }

  incidence_rate <- if (is.null(incidence)) NA_real_ else as.numeric(incidence$estimate) * 1000
  incidence_lcl <- if (is.null(incidence)) NA_real_ else incidence$conf.int[1L] * 1000
  incidence_ucl <- if (is.null(incidence)) NA_real_ else incidence$conf.int[2L] * 1000

  summary_dt <- data.table(
    section = c(
      rep("cohort", 4),
      rep("followup", 7),
      rep("outcome", 11)
    ),
    metric = c(
      "n_participants",
      "baseline_date_min",
      "baseline_date_max",
      "followup_end_date_max",
      "total_person_years",
      "minimum_followup_years",
      "q1_followup_years",
      "median_followup_years",
      "q3_followup_years",
      "maximum_followup_years",
      "mean_followup_years",
      "n_ascvd_events",
      "percent_ascvd_events",
      "n_competing_deaths",
      "percent_competing_deaths",
      "n_censored",
      "percent_censored",
      "n_lost_to_followup",
      "percent_lost_to_followup",
      "ascvd_incidence_rate_per_1000_person_years",
      "ascvd_incidence_rate_95ci_lower",
      "ascvd_incidence_rate_95ci_upper"
    ),
    value = as.character(c(
      n_total,
      format(min(dt$baseline_date), "%Y-%m-%d"),
      format(max(dt$baseline_date), "%Y-%m-%d"),
      format(max(dt$followup_end_date), "%Y-%m-%d"),
      total_person_years,
      followup_quantiles[1L],
      followup_quantiles[2L],
      followup_quantiles[3L],
      followup_quantiles[4L],
      followup_quantiles[5L],
      mean(dt$followup_years),
      n_ascvd,
      100 * n_ascvd / n_total,
      n_competing,
      100 * n_competing / n_total,
      n_censored,
      100 * n_censored / n_total,
      n_lost,
      100 * n_lost / n_total,
      incidence_rate,
      incidence_lcl,
      incidence_ucl
    )),
    unit = c(
      "participants", "date", "date", "date",
      "person_years", rep("years", 6),
      "events", "percent", "events", "percent", "participants", "percent",
      "participants", "percent", rep("per_1000_person_years", 3)
    )
  )

  reverse_km <- estimate_reverse_km_followup(dt)
  reverse_km[, `:=`(section = "followup", unit = "years")]
  reverse_km[, value := as.character(value)]
  setcolorder(reverse_km, names(summary_dt))

  rbindlist(list(summary_dt, reverse_km), use.names = TRUE)[]
}

build_missingness_summary <- function(dt, dictionary = NULL) {
  dt <- as.data.table(dt)

  out <- data.table(
    variable = names(dt),
    class = vapply(dt, function(x) class(x)[1L], character(1)),
    n_missing = vapply(dt, function(x) sum(is.na(x)), integer(1)),
    percent_missing = vapply(dt, function(x) 100 * mean(is.na(x)), numeric(1)),
    n_unique = vapply(dt, function(x) uniqueN(x, na.rm = TRUE), integer(1))
  )

  if (!is.null(dictionary)) {
    dictionary <- as.data.table(dictionary)
    meta <- dictionary[, .(
      variable = standard_name,
      label,
      unit,
      role,
      required_observed,
      impute,
      imputation_method,
      derive_after_mi,
      used_for_eligibility,
      used_for_strategy,
      used_for_outcome_model,
      used_for_subgroup
    )]
    out <- merge(meta, out, by = "variable", all.y = TRUE, sort = FALSE)
    out <- out[match(names(dt), variable)]
  }

  out[]
}

validate_final_analysis_base <- function(dt, input_n, input_ids = NULL) {
  dt <- as.data.table(dt)
  assert_unique_id(dt, "participant_id", "05 最终分析基础数据")

  if (nrow(dt) != input_n) {
    stop(
      "05 改变了 04 的样本量：输入 = ", input_n, "，输出 = ", nrow(dt),
      call. = FALSE
    )
  }

  if (!is.null(input_ids)) {
    input_ids <- as.character(input_ids)
    output_ids <- as.character(dt$participant_id)
    if (!setequal(input_ids, output_ids)) {
      stop("05 输入与输出 participant_id 集合不一致。", call. = FALSE)
    }
    if (!identical(input_ids, output_ids)) {
      stop("05 改变了 participant_id 的行顺序。", call. = FALSE)
    }
  }

  invalid_names <- names(dt)[!grepl("^[a-z][a-z0-9_]*$", names(dt))]
  if (length(invalid_names) > 0L) {
    stop("最终分析数据存在非 lower_snake_case 变量：", paste(invalid_names, collapse = ", "), call. = FALSE)
  }

  if (any(dt$ascvd_event + dt$competing_death_event > 1L)) {
    stop("存在同时记为 ASCVD 事件和竞争死亡的记录。", call. = FALSE)
  }

  invisible(TRUE)
}


# ============================================================================== #
# 5. 最终数据契约审计 ####
# ============================================================================== #

build_row_identity_audit <- function(input_ids, output_ids) {
  input_ids <- as.character(input_ids)
  output_ids <- as.character(output_ids)
  out <- data.table(
    check = c(
      "same_row_count", "input_id_complete", "output_id_complete",
      "input_id_unique", "output_id_unique", "same_id_set", "same_id_order"
    ),
    n_violations = c(
      abs(length(input_ids) - length(output_ids)),
      sum(is.na(input_ids) | trimws(input_ids) == ""),
      sum(is.na(output_ids) | trimws(output_ids) == ""),
      length(input_ids) - uniqueN(input_ids),
      length(output_ids) - uniqueN(output_ids),
      length(setdiff(input_ids, output_ids)) + length(setdiff(output_ids, input_ids)),
      if (length(input_ids) == length(output_ids)) sum(input_ids != output_ids) else max(length(input_ids), length(output_ids))
    )
  )
  out[, passed := n_violations == 0L]
  if (any(!out$passed)) stop("05 输入输出行或 ID 恒等性检查失败。", call. = FALSE)
  out[]
}

build_type_contract_audit <- function(dt, dictionary) {
  dt <- as.data.table(dt)
  dictionary <- as.data.table(dictionary)
  assert_required_columns(dictionary, c("standard_name", "type"), "05 类型契约字典")
  assert_required_columns(dt, dictionary$standard_name, "05 类型契约数据")

  out <- dictionary[, .(variable = standard_name, expected_type = type)]
  out[, observed_class := vapply(variable, function(v) class(dt[[v]])[1L], character(1))]
  out[, passed := mapply(function(expected, observed) {
    switch(
      expected,
      date = observed %in% c("IDate", "Date"),
      integer = observed == "integer",
      numeric = observed == "numeric",
      character = observed == "character",
      logical = observed == "logical",
      FALSE
    )
  }, expected_type, observed_class)]
  if (any(!out$passed)) {
    bad <- out[passed == FALSE, paste0(variable, "=", observed_class, " (expected ", expected_type, ")")]
    stop("05 最终变量类型不符合字典：", paste(bad, collapse = "; "), call. = FALSE)
  }
  out[]
}

build_required_observed_audit <- function(dt, dictionary) {
  dt <- as.data.table(dt)
  dictionary <- as.data.table(dictionary)
  required_vars <- dictionary[required_observed == TRUE, standard_name]
  assert_required_columns(dt, required_vars, "05 required_observed 检查")
  out <- data.table(
    variable = required_vars,
    n_total = nrow(dt),
    n_missing = vapply(required_vars, function(v) sum(is.na(dt[[v]])), integer(1))
  )
  out[, passed := n_missing == 0L]
  if (any(!out$passed)) {
    bad <- out[passed == FALSE, paste0(variable, "=", n_missing)]
    stop("05 required_observed 变量存在缺失：", paste(bad, collapse = "; "), call. = FALSE)
  }
  out[]
}

build_value_domain_audit <- function(dt) {
  dt <- as.data.table(dt)
  binary_vars <- c(
    "family_history_cvd", "antihypertensive_medication_baseline",
    "alcohol_score", "diet_score",
    "loss_to_followup_event", "ascvd_event", "competing_death_event",
    "chronic_respiratory_disease", "chronic_liver_disease"
  )
  assert_required_columns(dt, c(binary_vars, "event_status", "region", "followup_years"), "05 取值域检查")

  binary_audit <- rbindlist(lapply(binary_vars, function(v) {
    x <- dt[[v]]
    data.table(
      variable = v,
      allowed_values = "0;1",
      observed_values = paste(sort(unique(x[!is.na(x)])), collapse = ";"),
      n_missing = sum(is.na(x)),
      n_invalid = sum(!is.na(x) & !x %in% 0:1)
    )
  }))
  status_audit <- data.table(
    variable = "event_status",
    allowed_values = "0;1;2",
    observed_values = paste(sort(unique(dt$event_status[!is.na(dt$event_status)])), collapse = ";"),
    n_missing = sum(is.na(dt$event_status)),
    n_invalid = sum(!is.na(dt$event_status) & !dt$event_status %in% 0:2)
  )
  region_levels <- c("England", "Wales", "Scotland")
  region_audit <- data.table(
    variable = "region",
    allowed_values = paste(region_levels, collapse = ";"),
    observed_values = paste(sort(unique(dt$region[!is.na(dt$region)])), collapse = ";"),
    n_missing = sum(is.na(dt$region)),
    n_invalid = sum(!is.na(dt$region) & !dt$region %in% region_levels)
  )
  followup_audit <- data.table(
    variable = "followup_years",
    allowed_values = ">=0 and finite",
    observed_values = paste0(
      "min=", min(dt$followup_years, na.rm = TRUE),
      ";max=", max(dt$followup_years, na.rm = TRUE)
    ),
    n_missing = sum(is.na(dt$followup_years)),
    n_invalid = sum(!is.na(dt$followup_years) &
                      (!is.finite(dt$followup_years) | dt$followup_years < 0))
  )
  out <- rbindlist(list(binary_audit, status_audit, region_audit, followup_audit), use.names = TRUE)
  out[, passed := n_invalid == 0L]
  if (any(!out$passed)) {
    bad <- out[passed == FALSE, paste0(variable, ": missing=", n_missing, ", invalid=", n_invalid)]
    stop("05 变量取值域检查失败（真实缺失允许进入 MI；非法非缺失值不允许）：", paste(bad, collapse = "; "), call. = FALSE)
  }
  out[]
}

build_followup_event_consistency_audit <- function(dt) {
  dt <- as.data.table(dt)
  checks <- data.table(
    check = c(
      "status_matches_ascvd_and_competing_death",
      "ascvd_and_competing_death_mutually_exclusive",
      "ascvd_event_date_present",
      "competing_death_date_present",
      "followup_end_not_after_censor",
      "censored_followup_end_matches_censor"
    ),
    n_violations = c(
      sum(dt$event_status != fifelse(dt$ascvd_event == 1L, 1L, fifelse(dt$competing_death_event == 1L, 2L, 0L))),
      sum(dt$ascvd_event == 1L & dt$competing_death_event == 1L),
      sum(dt$ascvd_event == 1L & is.na(dt$ascvd_event_date)),
      sum(dt$competing_death_event == 1L & is.na(dt$competing_death_date)),
      sum(dt$followup_end_date > dt$censor_date, na.rm = TRUE),
      sum(dt$event_status == 0L & dt$followup_end_date != dt$censor_date, na.rm = TRUE)
    )
  )
  checks[, passed := n_violations == 0L]
  if (any(!checks$passed)) stop("05 随访与事件一致性检查失败。", call. = FALSE)
  checks[]
}

run_analysis_base_unit_tests <- function() {
  region_test <- data.table(region = c("England", "wales", "SCOTLAND"))
  region_result <- standardize_analysis_values(region_test)
  stopifnot(identical(region_result$region, c("England", "Wales", "Scotland")))

  row_audit <- build_row_identity_audit(c("1", "2"), c("1", "2"))
  stopifnot(all(row_audit$passed))

  event_test <- data.table(
    event_status = 0:2,
    ascvd_event = c(0L, 1L, 0L),
    competing_death_event = c(0L, 0L, 1L),
    ascvd_event_date = as.IDate(c(NA, "2020-01-02", NA)),
    competing_death_date = as.IDate(c(NA, NA, "2020-01-03")),
    followup_end_date = as.IDate(c("2020-01-04", "2020-01-02", "2020-01-03")),
    censor_date = as.IDate(rep("2020-01-04", 3L))
  )
  stopifnot(all(build_followup_event_consistency_audit(event_test)$passed))

  na_test <- data.table(
    followup_years = c(1, 2, 3),
    ascvd_event = c(0L, 1L, 0L),
    competing_death_event = c(0L, 0L, 1L)
  )
  na_result <- add_mi_nelson_aalen_predictors(na_test)
  stopifnot(
    all(is.finite(na_result$na_hazard_ascvd)),
    all(is.finite(na_result$na_hazard_death)),
    all(na_result$na_hazard_ascvd >= 0),
    all(na_result$na_hazard_death >= 0)
  )

  step_message(
    "05 合成边界测试通过：地区允许真实缺失、ID 恒等、随访状态一致、",
    "Nelson-Aalen MI auxiliary predictors 可生成。"
  )
  invisible(TRUE)
}

# Explicit three-valued OR: a known positive wins, two negatives yield zero, otherwise NA.
rm01_tri_or <- function(a,b) {
  if (length(a) != length(b)) stop("RM binary OR length mismatch")
  ifelse((!is.na(a) & a == 1L) | (!is.na(b) & b == 1L), 1L,
         ifelse(!is.na(a) & !is.na(b) & a == 0L & b == 0L,0L,NA_integer_))
}
rm01_version <- "M01_RM_v1.8_20261101"
rm01_manifest <- function() {
  path <- file.path(code_dir,"config","RM_domains_v1.7.csv")
  if (!file.exists(path)) stop("Missing RM domain manifest: ",path)
  x <- data.table::fread(path, encoding="UTF-8")
  required <- c("rm_id","information_type","include_2019","include_2025","rule")
  if (!all(required %in% names(x)) || nrow(x)!=16L || anyDuplicated(x$rm_id))
    stop("Invalid or incomplete RM domains manifest; expected 16 measurable domains")
  for (v in c("include_2019","include_2025")) {
    x[, (v) := as.logical(get(v))]
    if (anyNA(x[[v]])) stop("Invalid RM flag ",v)
  }
  x[, version := rm01_version]
  x
}
rm01_derive <- function(d, cutoff, preimputed_premature_menopause=FALSE) {
  if(length(cutoff)!=1L || !is.finite(cutoff)) stop("Invalid frozen Townsend Q80")
  out <- data.table::copy(d); female <- out$sex=="female"
  manifest <- rm01_manifest()
  for(i in seq_len(nrow(manifest))) {
    value <- eval(parse(text=manifest$rule[i]),envir=out,enclos=environment())
    out[, (paste0("rm_",manifest$rm_id[i])) := as.integer(value)]
  }
  # Reportability can establish Lp(a) status without inventing a concentration.
  flag <- as.character(out$rm_source_lpa_flag)
  low <- !is.na(flag)&(flag %in% c("2","4")|grepl("too low",flag,ignore.case=TRUE))
  high <- !is.na(flag)&(flag %in% c("3","5")|grepl("too high",flag,ignore.case=TRUE))
  out[is.na(rm_lpa)&low,rm_lpa:=0L];out[is.na(rm_lpa)&high,rm_lpa:=1L]
  if("rm_lpa_status" %in% names(out)) out[!is.na(rm_lpa_status),rm_lpa:=as.integer(as.character(rm_lpa_status))]
  early <- (out$rm_source_menopause==1L & out$rm_source_menopause_age<40) |
    (out$rm_source_oophorectomy==1L & out$rm_source_oophorectomy_age<40)
  out[, rm_menopause_uncertain := as.integer(female & is.na(early))]
  out[, rm_premature_menopause := as.integer(female & !is.na(early) & early)]
  out[!female,rm_hdp:=0L]
  for(year in c("2019","2025")) for(domain in c("H","T","any")) {
    use <- manifest[[paste0("include_",year)]] & (domain=="any" | manifest$information_type==domain)
    columns <- paste0("rm_",manifest$rm_id[use]);a<-as.matrix(out[,..columns])
    value<-ifelse(rowSums(a==1L,na.rm=TRUE)>0L,1L,ifelse(rowSums(is.na(a))==0L,0L,NA_integer_))
    out[, (paste0("rm_",tolower(domain),"_",year)) := as.integer(value)]
  }
  out
}
rm06_reproductive_components <- function() paste0("rm_source_", c(
  "menopause", "menopause_age", "hysterectomy", "hysterectomy_age",
  "oophorectomy", "oophorectomy_age"))
rm06_raw_binary_sources <- function() paste0("rm_source_",c(
  "south_asian","inactive","family_proxy","stressor_proxy","inflammatory",
  "psychiatric","af","hdp","hiv","sleep_apnoea","fatty_liver",
  "menopause","hysterectomy","oophorectomy"))
rm06_binary_sources <- function() setdiff(rm06_raw_binary_sources(),
  paste0("rm_source_",c("south_asian","menopause","hysterectomy","oophorectomy")))
rm06_binary_targets <- function() c(rm06_binary_sources(),"rm_lpa_status")
rm06_continuous_sources <- function() c(
  "waist_circumference","rm_source_crp")
rm06_imputation_targets <- function() c(rm06_continuous_sources(),rm06_binary_targets())
# 05/06共享描述与缺失审计。输入数据及变量列表；只返回报告，不插补、不写文件。
# 从M01母版的调用契约恢复该已缺失依赖；统计估计仍在原05/06函数中执行。
m01qa_assert_vars <- function(x,vars,label="data") {
  absent <- setdiff(vars,names(x))
  if(length(absent))stop(label," missing: ",paste(absent,collapse=", "))
  invisible(TRUE)
}
m01qa_rm_vars <- function() c(paste0("rm_",rm01_manifest()$rm_id),"rm_any_2019","rm_any_2025","rm_h_2019","rm_t_2019","rm_h_2025","rm_t_2025")
m01qa_missingness <- function(x,vars,stage,dataset_id=NA_integer_) {
  m01qa_assert_vars(x,vars)
  data.table::rbindlist(lapply(vars,function(v)data.table::data.table(
    stage=stage,dataset_id=dataset_id,variable=v,n=nrow(x),n_missing=sum(is.na(x[[v]])),
    n_observed=sum(!is.na(x[[v]])),pct_missing=100*mean(is.na(x[[v]])))))
}
m01qa_distribution <- function(x,vars,stage,dataset_id=NA_integer_) {
  m01qa_assert_vars(x,vars)
  data.table::rbindlist(lapply(vars,function(v) {
    z<-x[[v]];known<-z[!is.na(z)]
    if(is.numeric(z)) {
      q<-if(length(known))as.numeric(stats::quantile(known,c(0,.25,.5,.75,1)))else rep(NA_real_,5)
      data.table::data.table(stage=stage,dataset_id=dataset_id,variable=v,type="numeric",n=nrow(x),
        n_missing=sum(is.na(z)),mean=if(length(known))mean(known)else NA_real_,sd=stats::sd(known),
        min=q[1],q1=q[2],median=q[3],q3=q[4],max=q[5])
    } else {
      tab<-table(as.character(known));if(!length(tab))tab<-setNames(0L,"no observed values")
      data.table::data.table(stage=stage,dataset_id=dataset_id,variable=v,type="categorical",n=nrow(x),
        n_missing=sum(is.na(z)),level=names(tab),count=as.integer(tab),pct_total=100*as.integer(tab)/nrow(x))
    }
  }),fill=TRUE)
}
m01qa_rm_binary <- function(x,vars,stage,dataset_id=NA_integer_) {
  m01qa_assert_vars(x,vars)
  data.table::rbindlist(lapply(vars,function(v) {
    z<-x[[v]];if(any(!is.na(z)&!z%in%0:1))stop("RM not binary: ",v)
    k<-sum(z==1,na.rm=TRUE);known<-sum(!is.na(z))
    data.table::data.table(stage=stage,dataset_id=dataset_id,variable=v,n=nrow(x),n_missing=sum(is.na(z)),
      n_positive=k,n_negative=sum(z==0,na.rm=TRUE),pct_positive_total=100*k/nrow(x),
      pct_positive_observed=if(known)100*k/known else NA_real_)
  }))
}
m01qa_pre_post <- function(before,after,vars,dataset_id=NA_integer_,pre_suffix="") {
  m01qa_assert_vars(before,paste0(vars,pre_suffix));m01qa_assert_vars(after,vars)
  stopifnot(nrow(before)==nrow(after))
  data.table::rbindlist(lapply(vars,function(v) {
    a<-before[[paste0(v,pre_suffix)]];b<-after[[v]];known<-!is.na(a)
    changed<-known&(is.na(b)|as.character(a)!=as.character(b))
    data.table::data.table(stage="pre_post",dataset_id=dataset_id,variable=v,n=nrow(after),
      n_missing_before=sum(is.na(a)),n_missing_after=sum(is.na(b)),
      n_imputed=sum(is.na(a)&!is.na(b)),n_observed_changed=sum(changed),
      n_changed_observed=sum(known & !is.na(b) & as.character(a)!=as.character(b)),
      n_new_missing=sum(known & is.na(b)),
      pct_missing_before=100*mean(is.na(a)),pct_missing_after=100*mean(is.na(b)))
  }))
}
m01qa_cross_imputation_missingness <- function(x) {
  x[,.(n_imputations=data.table::uniqueN(dataset_id),
    n_missing_before=mean(n_missing_before),n_missing_after_min=min(n_missing_after),
    n_missing_after_max=max(n_missing_after),n_imputed_mean=mean(n_imputed),
    n_observed_changed_max=max(n_observed_changed)),by=variable]
}

# Final eligible-cohort audit, calculated before MI because follow-up is fixed.
build_followup_support_audit <- function(d) {
 stopifnot(!anyNA(d$followup_years),all(d$followup_years>=0))
 groups<-c("Overall",sort(unique(as.character(d$region))))
 summary<-risk<-annual<-list()
 for(g in groups) {
  z<-if(g=="Overall")d else d[region==g];t<-z$followup_years
  q<-quantile(t,c(.25,.5,.75),names=FALSE)
  summary[[g]]<-data.table(region=g,n=nrow(z),q1_years=q[1],median_years=q[2],q3_years=q[3],iqr_years=q[3]-q[1],
   ascvd_events=sum(z$ascvd_event),competing_deaths=sum(z$competing_death_event),max_years=max(t))
  risk[[g]]<-data.table(region=g,year=c(10,15),n_risk=c(sum(t>=10),sum(t>=15)),denominator=nrow(z))
  annual[[g]]<-rbindlist(lapply(10:14,function(a)data.table(region=g,start_year=a,end_year=a+1,
   n_risk_start=sum(t>=a),n_ascvd=sum(t>=a&t<a+1&z$ascvd_event==1),
   n_competing_death=sum(t>=a&t<a+1&z$competing_death_event==1),
   n_censored=sum(t>=a&t<a+1&z$ascvd_event==0&z$competing_death_event==0),
   person_years=sum(pmax(0,pmin(t,a+1)-a)))))
 }
 list(observed_followup=rbindlist(summary),risk_set=rbindlist(risk)[,pct_baseline:=100*n_risk/denominator],annual_tail_events=rbindlist(annual))
}

