# M01 v1.8 | 2026-11-01 | one main script / one function module
# ============================================================================== #
# 文件模块说明 | Function_06_Multiple_imputation_M01.R
# 任务：供06调用的MICE类型/变量准备、SES LCA和MI QA
# 输入：05 incomplete数据、MICE完成数据
# 输出：MICE配置、SES/活动/血压派生与诊断
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# ============================================================================== #
# Function_06_Multiple_imputation_M01_v1.1.R
# ------------------------------------------------------------------------------ #
# 用途：
#   为 06_Multiple_imputation_M01.R 提供 MICE/FCS、多重插补后变量派生、
#   completed dataset QA、SES 2-class LCA 和 MI 稳定性汇总函数。
#
# 核心原则：
#   1. 联合插补 outcome-model / subgroup / RM源变量；评分和纳排输入保持观测值。
#   2. 不直接插补 physical_activity、ses_class、hypertension。
#   3. 每套 completed dataset 中重新派生上述变量。
#   4. SCORE/SCORE2、风险等级、LDL target/trigger/counterfactual LDL 不进入 MI。
# ============================================================================== #

library(data.table)

# ------------------------------------------------------------------------------ #
# 0. 通用工具
# ------------------------------------------------------------------------------ #
mi06_step_message <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

mi06_assert_file <- function(path, label = "文件") {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(label, "不存在：", path, call. = FALSE)
  }
  invisible(path)
}

mi06_assert_columns <- function(dt, cols, label = "数据") {
  dt <- data.table::as.data.table(dt)
  missing_cols <- setdiff(unique(cols), names(dt))
  if (length(missing_cols) > 0L) {
    stop(label, "缺少变量：", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

mi06_assert_unique_id <- function(dt, id_col = "participant_id", label = "数据") {
  dt <- data.table::as.data.table(dt)
  mi06_assert_columns(dt, id_col, label)
  if (anyNA(dt[[id_col]]) || any(trimws(as.character(dt[[id_col]])) == "")) {
    stop(label, "存在缺失 participant_id。", call. = FALSE)
  }
  if (data.table::uniqueN(dt[[id_col]]) != nrow(dt)) {
    stop(label, "存在重复 participant_id。", call. = FALSE)
  }
  invisible(TRUE)
}

mi06_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(as.character(tools::md5sum(path)))
}

# ------------------------------------------------------------------------------ #
# 1. 变量类型和 MICE 输入准备
# ------------------------------------------------------------------------------ #
mi06_factor_levels <- function() {
  list(
    sex = c("female", "male"),
    smoking_status = c("never", "former", "current"),
    region = c("England", "Wales", "Scotland"),
    ethnicity = c("0", "1"),
    education = c("low", "medium", "high"),
    household_income = c(
      "less_than_18000", "18000_to_30999", "31000_to_51999",
      "52000_to_100000", "greater_than_100000"
    ),
    employment_status = c(
      "other_not_employed", "retired", "employed_or_self_employed"
    )
  )
}

mi06_binary_variables <- function() {
  c(
    rm06_binary_targets(),
    "alcohol_score", "diet_score", "family_history_cvd",
    "antihypertensive_medication_baseline",
    "chronic_liver_disease", "chronic_respiratory_disease",
    "ascvd_event", "competing_death_event"
  )
}

mi06_ordered_variables <- function() {
  c("education", "household_income")
}

mi06_nominal_variables <- function() {
  c("sex", "smoking_status", "region", "ethnicity", "employment_status")
}

mi06_prepare_mice_types <- function(dt) {
  out <- as.data.frame(data.table::copy(data.table::as.data.table(dt)))
  lvls <- mi06_factor_levels()

  for (v in intersect(names(lvls), names(out))) {
    if (v %in% mi06_ordered_variables()) {
      out[[v]] <- ordered(as.character(out[[v]]), levels = lvls[[v]])
    } else {
      out[[v]] <- factor(as.character(out[[v]]), levels = lvls[[v]])
    }
  }

  for (v in intersect(mi06_binary_variables(), names(out))) {
    x <- out[[v]]
    if (is.factor(x)) x <- as.character(x)
    x <- suppressWarnings(as.integer(x))
    bad <- !is.na(x) & !x %in% 0:1
    if (any(bad)) stop(v, " 存在非 0/1 取值。", call. = FALSE)
    out[[v]] <- factor(x, levels = c(0, 1))
  }

  out
}

mi06_restore_completed_types <- function(dt) {
  out <- data.table::copy(data.table::as.data.table(dt))
  char_vars <- intersect(
    c("sex", "smoking_status", "region", "ethnicity", "education",
      "household_income", "employment_status"),
    names(out)
  )
  for (v in char_vars) data.table::set(out, j = v, value = as.character(out[[v]]))

  for (v in intersect(mi06_binary_variables(), names(out))) {
    x <- out[[v]]
    if (is.factor(x)) x <- as.character(x)
    data.table::set(out, j = v, value = suppressWarnings(as.integer(x)))
  }
  out[]
}

# 构建“丰富但不把派生变量/日期塞入 MICE”的 predictor pool。
mi06_default_predictor_pool <- function() {
  c(
    rm06_imputation_targets(),
    # demographic / fixed strategy covariates
    "age", "sex", "smoking_status",
    "systolic_blood_pressure", "diastolic_blood_pressure",
    "total_cholesterol", "ldl_cholesterol", "hdl_cholesterol",
    "hba1c", "egfr_2021",

    # variables to be imputed / SES components
    "region", "ethnicity", "townsend_deprivation_index",
    "education", "household_income", "employment_status",
    "bmi", "alcohol_score", "diet_score", "triglycerides",
    "family_history_cvd", "antihypertensive_medication_baseline",
    "chronic_liver_disease", "chronic_respiratory_disease",

    # outcome information for survival-data MI
    "ascvd_event", "competing_death_event",
    "na_hazard_ascvd", "na_hazard_death", "followup_years"
  )
}

mi06_build_mice_spec <- function(data_analysis_base, variable_dictionary) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("06 需要 R 包 mice：install.packages('mice')。", call. = FALSE)
  }

  dt <- data.table::copy(data.table::as.data.table(data_analysis_base))
  dict <- data.table::copy(data.table::as.data.table(variable_dictionary))
  mi06_assert_unique_id(dt)

  # UKB raw MET 可能保留负值特殊编码。负值不代表真实活动量，不能作为 PMM donor，
  # 因此在 MICE 输入层统一转为 NA；05 原始 incomplete base 本身不被修改。
  if ("physical_activity_met_minutes_week" %in% names(dt)) {
    dt[
      !is.na(physical_activity_met_minutes_week) &
        physical_activity_met_minutes_week < 0,
      physical_activity_met_minutes_week := NA_real_
    ]
  }
  mi06_assert_columns(
    dict,
    c("standard_name", "impute", "imputation_method"),
    "05 variable dictionary"
  )

  impute_flag <- tolower(as.character(dict$impute)) %in% c("true", "1")
  impute_vars <- dict[impute_flag, as.character(standard_name)]
  impute_vars <- intersect(impute_vars, names(dt))
  if (length(impute_vars) == 0L) stop("变量字典中没有可执行的 MI 变量。", call. = FALSE)

  predictor_pool <- intersect(mi06_default_predictor_pool(), names(dt))
  input_vars <- unique(c(impute_vars, predictor_pool))
  mice_input <- dt[, ..input_vars]
  mice_input <- mi06_prepare_mice_types(mice_input)

  method <- rep("", length(input_vars))
  names(method) <- input_vars
  method_map <- dict[standard_name %in% impute_vars, .(
    standard_name = as.character(standard_name),
    imputation_method = as.character(imputation_method)
  )]
  method[method_map$standard_name] <- method_map$imputation_method

  # 只有真实存在缺失的变量才执行插补；完整变量保留为 predictor，但 method=""。
  n_missing <- vapply(mice_input, function(x) sum(is.na(x)), integer(1))
  method[n_missing == 0L] <- ""

  allowed_methods <- c("", "pmm", "logreg", "polr", "polyreg")
  invalid_methods <- setdiff(unique(method), allowed_methods)
  if (length(invalid_methods) > 0L) {
    stop("发现未支持的 MI method：", paste(invalid_methods, collapse = ", "), call. = FALSE)
  }

  # 方法与数据类型的一致性 QA。
  for (v in names(method)[method != ""]) {
    meth <- method[[v]]
    x <- mice_input[[v]]
    if (meth == "logreg" && (!is.factor(x) || nlevels(x) != 2L)) {
      stop(v, " 使用 logreg 时必须是二分类 factor。", call. = FALSE)
    }
    if (meth == "polr" && !is.ordered(x)) {
      stop(v, " 使用 polr 时必须是 ordered factor。", call. = FALSE)
    }
    if (meth == "polyreg" && (!is.factor(x) || is.ordered(x))) {
      stop(v, " 使用 polyreg 时必须是 nominal factor。", call. = FALSE)
    }
    if (meth == "pmm" && !is.numeric(x)) {
      stop(v, " 使用 pmm 时必须是 numeric。", call. = FALSE)
    }

    observed <- x[!is.na(x)]
    if (length(observed) == 0L) {
      stop(v, " 全部缺失，无法执行 ", meth, " 插补。", call. = FALSE)
    }
    if (is.factor(x) && data.table::uniqueN(as.character(observed)) < 2L) {
      stop(v, " 有缺失但已观测值不足两个类别，无法稳定执行 ", meth, "。", call. = FALSE)
    }
    if (meth == "pmm" && data.table::uniqueN(as.numeric(observed)) < 2L) {
      stop(v, " 有缺失但已观测连续值不足两个不同取值，无法稳定执行 PMM。", call. = FALSE)
    }
  }

  predictor_matrix <- mice::make.predictorMatrix(mice_input)
  predictor_matrix[,] <- 0L
  predictor_cols <- intersect(predictor_pool, colnames(predictor_matrix))
  target_rows <- names(method)[method != ""]
  predictor_matrix[target_rows, predictor_cols] <- 1L
  diag(predictor_matrix) <- 0L

  # 不用两个 event indicator 去互相“预测”彼此的确定性组合 event_status；
  # event_status 本身不进入 MICE 输入，减少不必要的确定性共线。
  list(
    mice_input = mice_input,
    impute_vars = impute_vars,
    method = method,
    predictor_matrix = predictor_matrix,
    missingness = data.table(
      variable = names(n_missing),
      n_missing = as.integer(n_missing),
      percent_missing = 100 * as.numeric(n_missing) / nrow(mice_input),
      method = unname(method[names(n_missing)])
    )
  )
}

# ------------------------------------------------------------------------------ #
# 2. Post-MI deterministic / model-based derivation
# ------------------------------------------------------------------------------ #
mi06_classify_physical_activity <- function(met_value) {
  data.table::fcase(
    is.na(met_value) | met_value < 0, NA_character_,
    met_value <= 600, "level1",
    met_value <= 3000, "level2",
    met_value <= 6000, "level3",
    met_value > 6000, "level4",
    default = NA_character_
  )
}

mi06_derive_current_smoker <- function(smoking_status) {
  data.table::fcase(
    smoking_status == "current", 1L,
    smoking_status %chin% c("never", "former"), 0L,
    default = NA_integer_
  )
}

mi06_derive_hypertension <- function(sbp, dbp, medication) {
  medication <- suppressWarnings(as.integer(medication))
  out <- rep(NA_integer_, length(sbp))

  # 冻结规则：只要 observed SBP/DBP 已达到阈值，就直接判定 HTN=1，
  # 不受降压药原始缺失影响。06 completed data 中 medication 应已完整。
  bp_positive <- !is.na(sbp) & !is.na(dbp) & (sbp >= 140 | dbp >= 90)
  out[bp_positive] <- 1L

  bp_below <- !is.na(sbp) & !is.na(dbp) & sbp < 140 & dbp < 90
  out[bp_below & !is.na(medication) & medication == 1L] <- 1L
  out[bp_below & !is.na(medication) & medication == 0L] <- 0L
  out
}

mi06_derive_binary_ses_lca <- function(data,
                                        id_var = "participant_id",
                                        income_var = "household_income",
                                        education_var = "education",
                                        employment_var = "employment_status",
                                        nrep = 30L,
                                        maxiter = 5000L,
                                        seed = 123L) {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    stop("06 2-class SES LCA 需要 R 包 poLCA：install.packages('poLCA')。", call. = FALSE)
  }

  dt <- data.table::copy(data.table::as.data.table(data))
  required <- c(id_var, income_var, education_var, employment_var)
  mi06_assert_columns(dt, required, "SES LCA 输入")
  mi06_assert_unique_id(dt, id_var, "SES LCA 输入")

  income_levels <- c(
    "less_than_18000", "18000_to_30999", "31000_to_51999",
    "52000_to_100000", "greater_than_100000"
  )
  education_levels <- c("low", "medium", "high")
  employment_levels <- c(
    "other_not_employed", "retired", "employed_or_self_employed"
  )

  lca_dt <- dt[, .(
    participant_id___ = get(id_var),
    income_code = match(get(income_var), income_levels),
    education_code = match(get(education_var), education_levels),
    employment_code = match(get(employment_var), employment_levels)
  )]
  if (any(!stats::complete.cases(lca_dt[, .(income_code, education_code, employment_code)]))) {
    stop("06 SES 组成变量在 completed dataset 中仍有缺失/非法水平。", call. = FALSE)
  }

  set.seed(seed)
  fit <- poLCA::poLCA(
    cbind(income_code, education_code, employment_code) ~ 1,
    data = as.data.frame(lca_dt[, .(income_code, education_code, employment_code)]),
    nclass = 2L,
    nrep = as.integer(nrep),
    maxiter = as.integer(maxiter),
    verbose = FALSE,
    graphs = FALSE,
    calc.se = FALSE,
    na.rm = FALSE
  )

  # 对 poLCA 的任意 class number 做跨 MI label alignment。
  # ordinal_score 越高代表总体社会经济状况越高。
  ordinal_score <-
    (lca_dt$income_code - 1) / 4 +
    (lca_dt$education_code - 1) / 2 +
    (lca_dt$employment_code - 1) / 2
  posterior <- fit$posterior
  class_mean_score <- colSums(posterior * ordinal_score) / colSums(posterior)
  ordered_classes <- order(class_mean_score)
  label_map <- rep(NA_character_, 2L)
  label_map[ordered_classes] <- c("low", "high")
  posterior_ordered <- posterior[, ordered_classes, drop = FALSE]

  ses <- data.table(
    participant_id___ = lca_dt$participant_id___,
    ses_class = label_map[fit$predclass],
    ses_probability_low = posterior_ordered[, 1L],
    ses_probability_high = posterior_ordered[, 2L]
  )
  data.table::setnames(ses, "participant_id___", id_var)

  max_posterior <- apply(posterior_ordered, 1L, max)
  posterior_entropy <- -sum(posterior_ordered * log(pmax(posterior_ordered, .Machine$double.eps)))
  relative_entropy <- 1 - posterior_entropy / (nrow(posterior_ordered) * log(2))

  safe_scalar <- function(x) {
    if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x[1L])
  }

  model_summary <- data.table(
    metric = c(
      "n_participants", "log_likelihood", "aic", "bic",
      "relative_entropy", "mean_max_posterior_probability",
      "proportion_max_posterior_below_0_70",
      "low_class_mean_ordinal_score", "high_class_mean_ordinal_score"
    ),
    value = c(
      nrow(lca_dt), safe_scalar(fit$llik), safe_scalar(fit$aic), safe_scalar(fit$bic),
      relative_entropy, mean(max_posterior), mean(max_posterior < 0.70),
      class_mean_score[ordered_classes[1L]], class_mean_score[ordered_classes[2L]]
    )
  )

  class_distribution <- ses[, .(
    n = .N,
    percent = 100 * .N / nrow(ses),
    mean_assigned_probability = if (ses_class[1L] == "low") {
      mean(ses_probability_low)
    } else {
      mean(ses_probability_high)
    }
  ), by = ses_class][order(match(ses_class, c("low", "high")))]

  component_profiles <- merge(
    ses[, c(id_var, "ses_class"), with = FALSE],
    dt[, c(id_var, income_var, education_var, employment_var), with = FALSE],
    by = id_var,
    all.x = TRUE,
    sort = FALSE
  )
  profile_rows <- list()
  for (v in c(income_var, education_var, employment_var)) {
    tmp <- component_profiles[, .N, by = .(ses_class, value = as.character(get(v)))]
    tmp[, percent := 100 * N / sum(N), by = ses_class]
    tmp[, variable := v]
    profile_rows[[v]] <- tmp
  }
  component_profile <- rbindlist(profile_rows, use.names = TRUE, fill = TRUE)

  list(
    ses_by_id = ses[],
    model_summary = model_summary[],
    class_distribution = class_distribution[],
    component_profile = component_profile[],
    model = fit
  )
}

mi06_post_imputation_derivation <- function(completed_dt,
                                             imputation_id,
                                             lca_nrep = 30L,
                                             lca_maxiter = 5000L,
                                             lca_seed_base = 205000L) {
  out <- data.table::copy(data.table::as.data.table(completed_dt))
  mi06_assert_unique_id(out)
  mi06_assert_columns(
    out,
    c(
      "rm_source_inactive", "smoking_status",
      "education", "household_income", "employment_status",
      "systolic_blood_pressure", "diastolic_blood_pressure",
      "antihypertensive_medication_baseline"
    ),
    "06 completed dataset"
  )

  # Updated study plan: direct Field 22035 binary inactivity; raw MET is QA only.
  mi06_assert_columns(out,"rm_source_inactive","completed RM source")
  out[, physical_activity := as.integer(rm_source_inactive)]
  out[, physical_inactivity := as.integer(rm_source_inactive)]
  out[, current_smoker := mi06_derive_current_smoker(smoking_status)]
  out[, hypertension := mi06_derive_hypertension(
    systolic_blood_pressure,
    diastolic_blood_pressure,
    antihypertensive_medication_baseline
  )]

  lca <- mi06_derive_binary_ses_lca(
    out,
    nrep = lca_nrep,
    maxiter = lca_maxiter,
    seed = as.integer(lca_seed_base + imputation_id)
  )
  out[lca$ses_by_id, on = "participant_id", `:=`(
    ses_class = i.ses_class,
    ses_probability_low = i.ses_probability_low,
    ses_probability_high = i.ses_probability_high
  )]

  if (anyNA(out$physical_activity)) stop("06 派生 physical_activity 后仍有缺失。", call. = FALSE)
  if (anyNA(out$current_smoker)) stop("06 派生 current_smoker 后仍有缺失。", call. = FALSE)
  if (anyNA(out$hypertension)) stop("06 派生 hypertension 后仍有缺失。", call. = FALSE)
  if (anyNA(out$ses_class)) stop("06 派生 ses_class 后仍有缺失。", call. = FALSE)

  list(data = out[], lca = lca)
}

# ------------------------------------------------------------------------------ #
# 3. Completed dataset QA
# ------------------------------------------------------------------------------ #
mi06_outcome_model_variables <- function() {
  c(
    "age", "sex", "ethnicity", "region", "ses_class", "rm_deprivation",
    "family_history_cvd", "current_smoker", "rm_inactive",
    "alcohol_score", "diet_score", "rm_obesity",
    "systolic_blood_pressure", "diastolic_blood_pressure",
    "antihypertensive_medication_baseline",
    "ldl_cholesterol", "hdl_cholesterol", "triglycerides", "egfr_2021",
    "chronic_liver_disease", "chronic_respiratory_disease",
    "rm_inflammatory", "rm_psychiatric", "rm_af", "rm_lpa", "rm_crp",
    "followup_years", "ascvd_event", "competing_death_event", "event_status"
  )
}

mi06_validate_completed_dataset <- function(dt, reference_ids, imputation_id) {
  dt <- data.table::as.data.table(dt)
  mi06_assert_unique_id(dt)
  required <- mi06_outcome_model_variables()
  mi06_assert_columns(dt, required, paste0("MI ", imputation_id))

  ids <- as.character(dt$participant_id)
  reference_ids <- as.character(reference_ids)
  n_missing_required <- sum(vapply(required, function(v) sum(is.na(dt[[v]])), integer(1)))
  invalid_hypertension <- if ("hypertension" %in% names(dt)) {
    sum(!is.na(dt$hypertension) & !dt$hypertension %in% 0:1)
  } else Inf

  audit <- data.table(
    imputation_id = imputation_id,
    check = c(
      "same_row_count", "same_id_set", "same_id_order",
      "outcome_model_variables_complete", "ses_class_valid",
      "physical_activity_valid", "hypertension_valid"
    ),
    n_violations = c(
      abs(length(reference_ids) - nrow(dt)),
      length(setdiff(reference_ids, ids)) + length(setdiff(ids, reference_ids)),
      if (length(reference_ids) == length(ids)) sum(reference_ids != ids) else max(length(reference_ids), length(ids)),
      n_missing_required,
      sum(is.na(dt$ses_class) | !dt$ses_class %in% c("low", "high")),
      sum(is.na(dt$physical_activity) | !dt$physical_activity %in% 0:1),
      invalid_hypertension
    )
  )
  audit[, passed := n_violations == 0L]
  if (any(!audit$passed)) {
    bad <- audit[passed == FALSE, paste0(check, "=", n_violations)]
    stop("MI ", imputation_id, " completed dataset QA 失败：", paste(bad, collapse = "; "), call. = FALSE)
  }
  audit[]
}

# ------------------------------------------------------------------------------ #
# 4. 分布与稳定性 QA
# ------------------------------------------------------------------------------ #
mi06_summarize_completed_distribution <- function(dt, imputation_id, variables) {
  dt <- data.table::as.data.table(dt)
  rows <- lapply(intersect(variables, names(dt)), function(v) {
    x <- dt[[v]]
    if (is.numeric(x) || is.integer(x)) {
      x_num <- as.numeric(x)
      data.table(
        imputation_id = imputation_id,
        variable = v,
        statistic = c("n", "mean", "sd", "q1", "median", "q3", "min", "max"),
        level = NA_character_,
        value = c(
          sum(!is.na(x_num)), mean(x_num, na.rm = TRUE), stats::sd(x_num, na.rm = TRUE),
          stats::quantile(x_num, 0.25, na.rm = TRUE, names = FALSE),
          stats::median(x_num, na.rm = TRUE),
          stats::quantile(x_num, 0.75, na.rm = TRUE, names = FALSE),
          min(x_num, na.rm = TRUE), max(x_num, na.rm = TRUE)
        )
      )
    } else {
      tab <- dt[, .N, by = .(level = as.character(get(v)))]
      tab[, value := 100 * N / sum(N)]
      tab[, `:=`(
        imputation_id = imputation_id,
        variable = v,
        statistic = "percent"
      )]
      tab[, .(imputation_id, variable, statistic, level, value)]
    }
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}


# 比较原始 observed values 与每套 MI 中“原本缺失位置”的 imputed values。
# 该输出仅用于插补分布 QA，不参与正式估计。
mi06_compare_observed_imputed <- function(original_dt,
                                           completed_list,
                                           imputed_vars) {
  original <- data.table::as.data.table(original_dt)
  if (length(completed_list) < 1L) stop("completed_list 为空。", call. = FALSE)
  mi06_assert_unique_id(original)

  rows <- list()
  k <- 0L
  for (v in intersect(imputed_vars, names(original))) {
    x0 <- original[[v]]
    # raw MET 的负值特殊编码在 MICE 输入层按缺失处理，这里 QA 使用同一规则。
    if (identical(v, "physical_activity_met_minutes_week")) {
      x0 <- suppressWarnings(as.numeric(x0))
      x0[!is.na(x0) & x0 < 0] <- NA_real_
    }
    miss_idx <- is.na(x0)
    obs_idx <- !miss_idx

    summarize_numeric <- function(x, source, m) {
      x <- suppressWarnings(as.numeric(x))
      x <- x[is.finite(x)]
      data.table::data.table(
        variable = v,
        source = source,
        imputation_id = m,
        statistic = c("n", "mean", "sd", "q1", "median", "q3", "min", "max"),
        level = NA_character_,
        value = c(
          length(x),
          if (length(x)) mean(x) else NA_real_,
          if (length(x) > 1L) stats::sd(x) else NA_real_,
          if (length(x)) stats::quantile(x, 0.25, names = FALSE) else NA_real_,
          if (length(x)) stats::median(x) else NA_real_,
          if (length(x)) stats::quantile(x, 0.75, names = FALSE) else NA_real_,
          if (length(x)) min(x) else NA_real_,
          if (length(x)) max(x) else NA_real_
        )
      )
    }

    summarize_categorical <- function(x, source, m) {
      x <- as.character(x)
      x <- x[!is.na(x)]
      if (length(x) == 0L) {
        return(data.table::data.table(
          variable = v, source = source, imputation_id = m,
          statistic = "percent", level = NA_character_, value = NA_real_
        ))
      }
      tab <- data.table::data.table(level = x)[, .N, by = level]
      tab[, value := 100 * N / sum(N)]
      tab[, `:=`(
        variable = v, source = source, imputation_id = m,
        statistic = "percent"
      )]
      tab[, .(variable, source, imputation_id, statistic, level, value)]
    }

    is_num <- is.numeric(x0) || is.integer(x0)
    k <- k + 1L
    rows[[k]] <- if (is_num) {
      summarize_numeric(x0[obs_idx], "observed", 0L)
    } else {
      summarize_categorical(x0[obs_idx], "observed", 0L)
    }

    if (any(miss_idx)) {
      for (m in seq_along(completed_list)) {
        comp <- data.table::as.data.table(completed_list[[m]])
        if (!identical(as.character(comp$participant_id), as.character(original$participant_id))) {
          stop("observed-vs-imputed QA 时 participant_id 顺序不一致。", call. = FALSE)
        }
        k <- k + 1L
        rows[[k]] <- if (is_num) {
          summarize_numeric(comp[[v]][miss_idx], "imputed", as.integer(m))
        } else {
          summarize_categorical(comp[[v]][miss_idx], "imputed", as.integer(m))
        }
      }
    }
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

mi06_build_subgroup_stability <- function(completed_list) {
  if (length(completed_list) < 1L) stop("completed_list 为空。", call. = FALSE)
  ref_ids <- as.character(completed_list[[1L]]$participant_id)
  m <- length(completed_list)

  ses_mat <- sapply(completed_list, function(dt) {
    dt <- data.table::as.data.table(dt)
    if (!identical(as.character(dt$participant_id), ref_ids)) {
      stop("构造 subgroup stability 时 MI ID 顺序不一致。", call. = FALSE)
    }
    as.integer(dt$ses_class == "low")
  }, simplify = "matrix")

  htn_mat <- sapply(completed_list, function(dt) {
    dt <- data.table::as.data.table(dt)
    as.integer(dt$hypertension == 1L)
  }, simplify = "matrix")

  participant <- data.table(
    participant_id = ref_ids,
    n_low_ses = rowSums(ses_mat),
    n_hypertension = rowSums(htn_mat),
    m = m
  )

  summary <- rbindlist(list(
    participant[, .N, by = .(count = n_low_ses)][, variable := "ses_low_count"],
    participant[, .N, by = .(count = n_hypertension)][, variable := "hypertension_count"]
  ))
  summary[, percent := 100 * N / sum(N), by = variable]
  setorder(summary, variable, count)
  list(participant = participant[], summary = summary[])
}

# ------------------------------------------------------------------------------ #
# 5. Branch A/B 汇合前的不变量 QA helper
# ------------------------------------------------------------------------------ #
mi06_strategy_required_observed_variables <- function() {
  c(
    "participant_id", "age", "sex", "smoking_status",
    "systolic_blood_pressure", "diastolic_blood_pressure",
    "total_cholesterol", "ldl_cholesterol", "hdl_cholesterol"
  )
}

mi06_validate_strategy_inputs_unchanged <- function(original_dt,
                                                     completed_dt,
                                                     imputation_id,
                                                     tolerance = 1e-12) {
  original <- data.table::as.data.table(original_dt)
  completed <- data.table::as.data.table(completed_dt)
  vars <- mi06_strategy_required_observed_variables()
  mi06_assert_columns(original, vars, "05 original incomplete base")
  mi06_assert_columns(completed, vars, paste0("MI ", imputation_id))
  if (!identical(
    as.character(original$participant_id),
    as.character(completed$participant_id)
  )) {
    stop("MI ", imputation_id, " participant_id 顺序与 05 base 不一致。", call. = FALSE)
  }

  rows <- lapply(setdiff(vars, "participant_id"), function(v) {
    x <- original[[v]]
    y <- completed[[v]]
    if (is.numeric(x) || is.integer(x)) {
      bad <- xor(is.na(x), is.na(y)) |
        (!is.na(x) & !is.na(y) & abs(as.numeric(x) - as.numeric(y)) > tolerance)
    } else {
      bad <- xor(is.na(x), is.na(y)) |
        (!is.na(x) & !is.na(y) & as.character(x) != as.character(y))
    }
    data.table::data.table(
      imputation_id = as.integer(imputation_id),
      variable = v,
      n_mismatch = sum(bad),
      passed = sum(bad) == 0L
    )
  })
  out <- data.table::rbindlist(rows)
  if (any(!out$passed)) {
    bad <- out[passed == FALSE, paste0(variable, "=", n_mismatch)]
    stop(
      "MI 改变了 Branch A strategy-required observed inputs：",
      paste(bad, collapse = "; "),
      call. = FALSE
    )
  }
  out[]
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
# ============================================================================== #
# 文件模块说明 | Function_06_RM_imputation_M01.R
# 任务：供06调用的RM特殊插补及二分类完成派生
# 输入：05观测RM和连续指标、MICE状态
# 输出：限制条件下的完成RM、原始观测保持不变
# 规则：函数默认不写文件；具体输出路径由01/02/05/06主脚本统一指定。
# ============================================================================== #

# M01 v1.7 FCS: continuous sources, unresolved Lp(a) binary status; menopause fixed.
# Only continuous non-ratio RM sources and already-classified binary RMs
# are imputed. Reproductive status, ages and urine ACR sources are NEVER imputed.
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
rm06_flag_limits <- function(flag,low,high) {
  code <- as.character(flag)
  lo <- code %in% c("2","4") | grepl("too low",code,ignore.case=TRUE)
  hi <- code %in% c("3","5") | grepl("too high",code,ignore.case=TRUE)
  lo[is.na(lo)]<-FALSE;hi[is.na(hi)]<-FALSE
  list(lower=ifelse(hi,high,0),upper=ifelse(lo,low,Inf))
}
rm06_prepare <- function(d) {
  x <- data.table::copy(data.table::as.data.table(d))
  x[, rm_lpa_status := as.integer(rm_lpa_observed)]
  x[, family_history_cvd := as.integer(rm_source_family_proxy)]
  needed <- c(rm06_imputation_targets(),
    paste0("rm_source_",c("lpa_flag","crp_flag","tg_flag")))
  absent <- setdiff(needed,names(x))
  if(length(absent)) stop("Re-run revised 05: missing ",paste(absent,collapse=", "))
  female <- x$sex == "female"
  if(anyNA(female)) stop("Missing sex in MI input")
  # HDP is undefined in men. Standardize that structural non-applicability
  # before MICE, as in the historical 06 implementation.
  x[!female, rm_source_hdp := 0L]
  pm <- x$rm_premature_menopause
  if("rm_premature_menopause_observed" %in% names(x)) {
    recorded <- x$rm_premature_menopause_observed
    same <- (is.na(recorded)&is.na(pm)) |
      (!is.na(recorded)&!is.na(pm)&recorded==pm)
    if(any(!same)) stop("Direct binary RM differs from observed 05 classification")
  }
  if(any(!is.na(pm) & !pm %in% 0:1)) stop("Invalid observed premature menopause 0/1")
  if(any(!female & (is.na(pm)|pm!=0L))) stop("Male premature menopause must be 0")
  for(v in c("rm_source_hdp","rm_premature_menopause")) {
    bad <- !female & (is.na(x[[v]])|x[[v]] != 0L)
    if(any(bad)) stop("Male female-specific variable must be 0: ",v)
  }
  bounds <- list(
    lipoprotein_a=rm06_flag_limits(x$rm_source_lpa_flag,3.8,189),
    rm_source_crp=rm06_flag_limits(x$rm_source_crp_flag,0.08,80),
    triglycerides=rm06_flag_limits(x$rm_source_tg_flag,0.1,11.3))
  labels <- c("Reportable at assay and after aliquot correction, if attempted",
    "Reportable at assay but not reportable after any corrections (too low)",
    "Reportable at assay but not reportable after any corrections (too high)",
    "Not reportable at assay (too low)","Not reportable at assay (too high)")
  reportability_audit <- data.table::rbindlist(lapply(c("lpa","crp","tg"),function(v) {
    flag <- as.character(x[[paste0("rm_source_",v,"_flag")]])
    unknown <- !is.na(flag) & !flag %in% c(as.character(1:5),labels)
    if(any(unknown)) stop("Unknown reportability flag for ",v)
    data.table::data.table(assay=v,reportability=flag)[,.N,by=.(assay,reportability)]
  }))
  lpa_value<-as.numeric(x$lipoprotein_a)
  lpa_known<-is.finite(lpa_value)&lpa_value>0
  if(any(lpa_known & (lpa_value<bounds$lipoprotein_a$lower | lpa_value>bounds$lipoprotein_a$upper)))
    stop("Observed Lp(a) concentration conflicts with verified assay bounds")
  bounds$lipoprotein_a <- NULL
  bounds$apolipoprotein_b <- NULL
  for(v in names(bounds)) {
    val <- as.numeric(x[[v]])
    conflict <- !is.na(val) &
      (val<bounds[[v]]$lower | val>bounds[[v]]$upper)
    if(any(conflict,na.rm=TRUE)) stop("Observed value conflicts with assay flag: ",v)
    val[!is.na(val)&val<=0] <- NA_real_
    data.table::set(x,j=v,value=val)
  }
  list(data=x,bounds=bounds,female=female,
       reportability_audit=reportability_audit)
}

rm06_truncated_normal <- function(mu,sd,lower,upper) {
  a <- (lower-mu)/sd;b <- (upper-mu)/sd;u <- runif(length(mu));z<-numeric(length(mu))
  right <- a>0
  if(any(right)) {
    la<-pnorm(a[right],lower.tail=FALSE,log.p=TRUE);lb<-pnorm(b[right],lower.tail=FALSE,log.p=TRUE)
    lp<-la+log(u[right]+(1-u[right])*exp(lb-la))
    z[right]<-qnorm(lp,lower.tail=FALSE,log.p=TRUE)
  }
  left<-!right
  if(any(left)) {
    la<-pnorm(a[left],log.p=TRUE);lb<-pnorm(b[left],log.p=TRUE)
    lp<-lb+log(u[left]+(1-u[left])*exp(la-lb))
    z[left]<-qnorm(lp,log.p=TRUE)
  }
  ans<-mu+sd*z
  if(any(!is.finite(ans)))stop("Nonfinite truncated draw; check bounds")
  pmin(upper,pmax(lower,ans))
}
rm06_continuous_draw <- function(y,ry,x,wy,lower,upper,fit_eligible=rep(TRUE,length(y))) {
  donors <- ry & fit_eligible & is.finite(y)&y>0
  if(sum(donors)<5L)stop("Continuous RM has fewer than 5 eligible donors")
  X<-as.matrix(x);if(anyNA(X))stop("Incomplete FCS predictors")
  if(ncol(X)) {
    keep<-apply(X[donors,,drop=FALSE],2,sd)>1e-10;X<-X[,keep,drop=FALSE]
    if(ncol(X))X<-scale(X,center=colMeans(X[donors,,drop=FALSE]),scale=apply(X[donors,,drop=FALSE],2,sd))
  }
  X<-cbind(1,X)
  censored <- !ry & fit_eligible & (lower>0 | is.finite(upper))
  if(any(censored)) {
    # Use the interval-censored likelihood, not a regression fitted solely to
    # detectable values (which would bias urine albumin and Lp(a) distributions).
    use<-donors|censored;lo<-hi<-log(y[use]);ci<-censored[use]
    lo[ci]<-log(lower[use][ci]);hi[ci]<-log(upper[use][ci])
    lo[!is.finite(lo)]<-NA_real_;hi[!is.finite(hi)]<-NA_real_
    q<-qr(X[use,,drop=FALSE]);cols<-q$pivot[seq_len(q$rank)]
    design<-X[,cols,drop=FALSE]
    frame<-as.data.frame(design[use,,drop=FALSE])
    names(frame)<-paste0("p",seq_len(ncol(design)))
    frame$.lo<-lo;frame$.hi<-hi
    form<-stats::as.formula(paste0("survival::Surv(.lo,.hi,type='interval2') ~ 0 + ",paste(names(frame)[seq_len(ncol(design))],collapse=" + ")))
    fit<-survival::survreg(form,data=frame,dist="gaussian",control=survival::survreg.control(maxiter=100))
    if(fit$iter>=100L)stop("Censored RM model did not converge")
    theta<-c(stats::coef(fit),log(fit$scale));vc<-fit$var
    if(any(!is.finite(theta))||any(!is.finite(vc)))stop("Censored RM model fit failed")
    eig<-eigen((vc+t(vc))/2,symmetric=TRUE)
    if(min(eig$values)< -1e-8)stop("Invalid censored RM covariance")
    theta<-theta+as.vector(eig$vectors%*%(sqrt(pmax(eig$values,0))*rnorm(length(theta))))
    sigma<-exp(tail(theta,1));mu<-as.vector(design[wy,,drop=FALSE]%*%head(theta,-1))
    return(exp(rm06_truncated_normal(mu,sigma,log(pmax(lower[wy],.Machine$double.xmin)),log(upper[wy]))))
  }
  A<-crossprod(X[donors,,drop=FALSE])+diag(1e-5,ncol(X))
  V<-solve(A);beta<-as.vector(V%*%crossprod(X[donors,,drop=FALSE],log(y[donors])))
  resid<-log(y[donors])-as.vector(X[donors,,drop=FALSE]%*%beta)
  df<-max(3L,sum(donors)-qr(X[donors,,drop=FALSE])$rank)
  sigma<-sqrt(max(sum(resid^2),1e-6)/rchisq(1,df))
  beta<-beta+as.vector(t(chol(V))%*%rnorm(length(beta)))*sigma
  mu<-as.vector(X[wy,,drop=FALSE]%*%beta)
  exp(rm06_truncated_normal(mu,sigma,log(pmax(lower[wy],.Machine$double.xmin)),log(upper[wy])))
}
rm06_impute_source <- function(y,ry,x,wy=NULL,target) {
  if(is.null(wy))wy<-!ry
  b<-rm_mice_context$bounds[[target]]
  rm06_continuous_draw(y,ry,x,wy,b$lower,b$upper)
}
mice.impute.rm_crp <- function(y,ry,x,wy=NULL,...) rm06_impute_source(y,ry,x,wy,"rm_source_crp")
mice.impute.rm_tg <- function(y,ry,x,wy=NULL,...) rm06_impute_source(y,ry,x,wy,"triglycerides")

mice.impute.rm_south_asian <- function(y,ry,x,wy=NULL,...) {
  if(is.null(wy))wy<-!ry
  z<-mice::mice.impute.logreg(y,ry,x,wy=wy,...)
  ethnic_column<-grep("^ethnicitynon_white$",colnames(x),value=TRUE)
  if(length(ethnic_column)!=1L)stop("South Asian imputation requires broad ethnicity predictor")
  z[x[wy,ethnic_column]<0.5]<-0L
  z
}


# Female-specific risk modifiers: males contribute structural zeros,
# NEVER donors for female-condition logistic regression.
rm06_impute_female_binary <- function(y,ry,x,wy=NULL,target,...) {
  if(is.null(wy)) wy <- !ry
  ctx <- get("rm_mice_context",envir=.GlobalEnv,inherits=FALSE)
  if(any(wy & !ctx$female))
    stop("Attempted female-specific RM imputation in male: ",target)
  female_ry <- ry & ctx$female
  yy <- as.character(y[female_ry])
  if(sum(female_ry)<20L || length(unique(yy))!=2L ||
     sum(yy=="1")<5L || sum(yy=="0")<5L)
    stop("Too few female binary donors/positives for ",target)
  z <- mice::mice.impute.logreg(y=y,ry=female_ry,x=x,wy=wy,...)
  if(length(z)!=sum(wy) || anyNA(z))stop("Invalid binary draws for ",target)
  z
}
mice.impute.rm_female_hdp <- function(y,ry,x,wy=NULL,...) {
  rm06_impute_female_binary(y,ry,x,wy,target="rm_source_hdp",...)
}
rm06_extend_spec <- function(spec,context) {
  forbidden <- c(rm06_reproductive_components(),
    "rm_source_urine_albumin", "rm_source_urine_creatinine",
    "rm_albuminuria", "rm_albuminuria_observed")
  included <- intersect(forbidden, colnames(spec$predictor_matrix))
  if(length(included)) stop("Forbidden RM components in MICE: ",
    paste(included,collapse=", "),". Re-run revised 05")
  spec$method["family_history_cvd"] <- ""
  drop <- intersect(c("family_history_cvd","lipoprotein_a","apolipoprotein_b"),names(spec$mice_input))
  spec$mice_input[drop] <- NULL
  spec$method <- spec$method[!names(spec$method)%in%drop]
  spec$predictor_matrix <- spec$predictor_matrix[!rownames(spec$predictor_matrix)%in%drop,!colnames(spec$predictor_matrix)%in%drop,drop=FALSE]
  spec$where <- is.na(spec$mice_input)
  for(v in intersect(rm06_binary_targets(),names(spec$mice_input))) {
    x <- spec$mice_input[[v]]
    spec$mice_input[[v]] <- factor(as.integer(as.character(x)),levels=0:1)
  }
  meth <- c(rm_source_crp="rm_crp",triglycerides="rm_tg")
  for(v in names(meth))if(v %in% names(spec$method) && any(spec$where[,v]))
    spec$method[[v]] <- meth[[v]]
  for(v in c("rm_source_hdp")) {
    if(any(spec$where[,v])) {
      spec$method[[v]] <- "rm_female_hdp"
    }
    if(any(spec$where[!context$female,v]))
      stop("Male structural zeros selected for MI in ",v)
  }
  spec$visitSequence <- names(spec$method)[spec$method!=""]
  spec
}
# ------------------------------------------------------------------------------ #
# 06-Progress | Real-time MICE progress without changing imputation draws
# Task: write a per-chain progress event after each completed FCS iteration.
# Input: chain number, iteration limit and progress log directory.
# Output: MI_progress_chain_XX.csv plus console progress in the master process.
# Note: mice::post is an observation-free callback; it must not update imp/data.
# ------------------------------------------------------------------------------ #
rm06_progress_event <- function(event, iteration = 0L, detail = "") {
  if (!exists("rm_mice_progress_file", envir = .GlobalEnv, inherits = FALSE))
    return(invisible(NULL))
  path <- get("rm_mice_progress_file", envir = .GlobalEnv, inherits = FALSE)
  chain <- get("rm_mice_progress_chain", envir = .GlobalEnv, inherits = FALSE)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  detail <- gsub("[\\r\\n,]", " ", as.character(detail))
  cat(sprintf("%s,%d,%d,%s,%s\n", stamp, as.integer(chain),
              as.integer(iteration), as.character(event), detail),
      file = path, append = TRUE)
  if (isTRUE(get("rm_mice_progress_echo", envir = .GlobalEnv,
                 inherits = FALSE))) {
    limit <- get("rm_mice_progress_maxit", envir = .GlobalEnv,
                 inherits = FALSE)
    if (event == "iteration") {
      cat(sprintf("[%s] MICE | MI %02d | iteration %d/%d (%.0f%% of this chain)\n",
                  format(Sys.time(), "%H:%M:%S"), chain, iteration,
                  limit, 100 * iteration / limit))
    } else if (event %in% c("completed", "error")) {
      cat(sprintf("[%s] MICE | MI %02d | %s %s\n",
                  format(Sys.time(), "%H:%M:%S"), chain, event, detail))
    }
    flush.console()
  }
  invisible(NULL)
}

rm06_progress_tick <- function() {
  step <- get("rm_mice_progress_iteration", envir = .GlobalEnv,
              inherits = FALSE) + 1L
  assign("rm_mice_progress_iteration", step, envir = .GlobalEnv)
  rm06_progress_event("iteration", step)
  invisible(NULL)
}

rm06_progress_monitor <- function(progress_files, iterations, M, maxit,
                                  silent = FALSE) {
  for (m in seq_len(M)) {
    path <- progress_files[[m]]
    if (!file.exists(path)) next
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    if (!length(lines)) next
    events <- strsplit(lines, ",", fixed = TRUE)
    good <- Filter(function(v) length(v) >= 4L &&
                     v[[4L]] == "iteration" &&
                     !is.na(suppressWarnings(as.integer(v[[3L]]))), events)
    if (!length(good)) next
    newest <- max(vapply(good, function(v) as.integer(v[[3L]]), integer(1)))
    newest <- min(newest, maxit)
    if (newest <= iterations[[m]]) next
    for (step in seq.int(iterations[[m]] + 1L, newest)) {
      iterations[[m]] <- step
      if (!silent) {
        done <- sum(iterations)
        cat(sprintf("[%s] MICE progress | MI %02d/%02d | iteration %d/%d | total %d/%d (%.1f%%)\n",
                    format(Sys.time(), "%H:%M:%S"), m, M, step, maxit,
                    done, M * maxit, 100 * done / (M * maxit)))
        flush.console()
      }
    }
  }
  iterations
}

# Parallel PSOCK controller: use the parallel package's own dispatch/receive
# primitives, with socketSelect() for nonblocking polling of progress files.
# This avoids splitting MICE iterations or restarting FCS, preserving the
# original MICE sampler and the per-chain L'Ecuyer random streams.
rm06_progress_parallel_jobs <- function(cl, run, M, maxit, progress_files,
                                        poll_interval = 1) {
  ns <- asNamespace("parallel")
  if (!exists("sendCall", envir = ns, inherits = FALSE) ||
      !exists("recvOneResult", envir = ns, inherits = FALSE)) {
    stop("Installed parallel package does not support live PSOCK dispatch", call. = FALSE)
  }
  send <- get("sendCall", envir = ns, inherits = FALSE)
  receive <- get("recvOneResult", envir = ns, inherits = FALSE)
  nworkers <- length(cl)
  results <- vector("list", M)
  progress <- integer(M)
  next_job <- 1L
  received <- 0L
  for (w in seq_len(min(nworkers, M))) {
    send(cl[[w]], run, list(next_job), tag = next_job)
    next_job <- next_job + 1L
  }
  sockets <- lapply(cl, function(worker) worker$con)
  while (received < M) {
    # socketSelect() is supported for PSOCK socket connections on Windows.
    ready <- socketSelect(sockets, timeout = poll_interval)
    progress <- rm06_progress_monitor(progress_files, progress, M, maxit)
    if (!any(ready)) next
    packet <- receive(cl)
    received <- received + 1L
    m <- as.integer(packet$tag)
    if (length(m) != 1L || is.na(m) || m < 1L || m > M)
      stop("Unexpected MICE worker task tag", call. = FALSE)
    if (inherits(packet$value, "try-error") ||
        inherits(packet$value, "error")) {
      stop(sprintf("MICE chain %02d failed: %s", m,
                   as.character(packet$value)), call. = FALSE)
    }
    results[[m]] <- packet$value
    progress <- rm06_progress_monitor(progress_files, progress, M, maxit)
    cat(sprintf("[%s] MICE | chain %02d/%02d finished (%d/%d chains)\n",
                format(Sys.time(), "%H:%M:%S"), m, M, received, M))
    flush.console()
    if (next_job <= M) {
      send(cl[[packet$node]], run, list(next_job), tag = next_job)
      next_job <- next_job + 1L
    }
  }
  progress <- rm06_progress_monitor(progress_files, progress, M, maxit)
  results
}

rm06_parallel_mice <- function(spec, context, M, maxit, seed, workers, code_dir,
                               progress_dir = NULL, show_progress = TRUE,
                               poll_interval = 1) {
  workers <- min(as.integer(workers), M)
  if (!workers %in% 1:3) stop("MI workers must be 1, 2 or 3")
  if (!is.finite(poll_interval) || poll_interval < 0.2)
    stop("MI progress poll_interval must be >= 0.2 seconds")
  if (is.null(progress_dir)) progress_dir <- file.path(tempdir(), "MI_progress")
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  progress_files <- file.path(progress_dir,
    sprintf("MI_progress_chain_%02d.csv", seq_len(M)))
  # Overwrite only the designated progress logs, not other QC files.
  for (path in progress_files)
    writeLines("timestamp,chain,iteration,event,detail", path, useBytes = TRUE)

  set.seed(seed, kind = "L'Ecuyer-CMRG")
  streams <- vector("list", M)
  streams[[1L]] <- .Random.seed
  if (M > 1L) for (m in 2:M)
    streams[[m]] <- parallel::nextRNGStream(streams[[m - 1L]])

  visits <- spec$visitSequence
  if (!length(visits)) stop("No active MICE imputation targets")
  # The post-processing expression executes ONCE per FCS iteration, after
  # the final visited univariate imputation. No data values are changed.
  post <- mice::make.post(spec$mice_input)
  post[[tail(visits, 1L)]] <-
    'get("rm06_progress_tick", envir = .GlobalEnv, inherits = FALSE)()'

  run <- function(m) {
    data.table::setDTthreads(1L)
    assign("rm_mice_context", context, envir = .GlobalEnv)
    assign(".Random.seed", streams[[m]], envir = .GlobalEnv)
    assign("rm_mice_progress_file", progress_files[[m]], envir = .GlobalEnv)
    assign("rm_mice_progress_chain", as.integer(m), envir = .GlobalEnv)
    assign("rm_mice_progress_iteration", 0L, envir = .GlobalEnv)
    assign("rm_mice_progress_maxit", as.integer(maxit), envir = .GlobalEnv)
    chain_started <- Sys.time()
    assign("rm_mice_progress_echo", isTRUE(show_progress) && workers == 1L,
           envir = .GlobalEnv)
    rm06_progress_event("started", 0L)
    tryCatch({
      result <- mice::mice(
        spec$mice_input, m = 1L, maxit = maxit, method = spec$method,
        predictorMatrix = spec$predictor_matrix, where = spec$where,
        visitSequence = spec$visitSequence, post = post,
        seed = NA, printFlag = FALSE,
        remove.constant = FALSE, remove.collinear = FALSE)
      if (get("rm_mice_progress_iteration", envir = .GlobalEnv) != maxit)
        stop("MICE progress callback count differs from maxit; check post/visitSequence")
      rm06_progress_event("completed", maxit,
        sprintf("elapsed_seconds=%.2f", as.numeric(difftime(Sys.time(), chain_started, units = "secs"))))
      result
    }, error = function(e) {
      rm06_progress_event("error", get("rm_mice_progress_iteration",
                                      envir = .GlobalEnv), conditionMessage(e))
      stop(e)
    })
  }
  if (show_progress) {
    cat(sprintf("[%s] MICE started: M=%d, iterations=%d, workers=%d; total=%d rounds\n",
                format(Sys.time(), "%H:%M:%S"), M, maxit, workers, M * maxit))
    cat("MICE progress logs: ", progress_dir, "\n", sep = "")
    flush.console()
  }
  if (workers == 1L) {
    chains <- lapply(seq_len(M), run)
  } else {
    cl <- parallel::makePSOCKcluster(workers, setup_timeout = 60,
                                    timeout = 43200)
    chains <- tryCatch({
      parallel::clusterCall(cl, .libPaths, .libPaths())
      parallel::clusterExport(cl,
        c("spec", "context", "streams", "code_dir", "maxit", "workers",
          "progress_files", "post", "run", "show_progress"),
        envir = environment())
      parallel::clusterEvalQ(cl, {
        Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                   MKL_NUM_THREADS = "1")
        source(file.path(code_dir, "Function_06_Multiple_imputation_M01.R"),
               encoding = "UTF-8")
        NULL
      })
      if (show_progress) {
        rm06_progress_parallel_jobs(cl, run, M, maxit, progress_files,
                                    poll_interval = poll_interval)
      } else {
        parallel::parLapply(cl, seq_len(M), run)
      }
    }, finally = parallel::stopCluster(cl))
  }
  if (show_progress)
    cat(sprintf("[%s] MICE all %d chains completed; combining mids objects.\n",
                format(Sys.time(), "%H:%M:%S"), M))
  Reduce(function(a, b) {
    joined <- mice::ibind(a, b)
    joined$imp <- setNames(lapply(names(a$data), function(v) {
      av <- a$imp[[v]]; bv <- b$imp[[v]]
      if (is.null(av) && is.null(bv)) return(NULL)
      z <- cbind(av, bv)
      names(z) <- as.character(seq_len(ncol(z)))
      z
    }), names(a$data))
    joined
  }, chains)
}
rm06_finalize <- function(d,original,cutoff) {
  female<-d$sex=="female"
  completed_ethnicity <- suppressWarnings(as.integer(as.character(d$ethnicity)))
  if(anyNA(completed_ethnicity) || any(!completed_ethnicity %in% 0:1))
    stop("Ethnicity high-risk binary not completed")
  observed_ethnicity <- as.integer(original$rm_source_south_asian)
  if(any(!is.na(observed_ethnicity) &
         completed_ethnicity[!is.na(observed_ethnicity)] != observed_ethnicity[!is.na(observed_ethnicity)]))
    stop("Observed South Asian ethnicity changed during MI")
  data.table::set(d,j="rm_source_south_asian",value=completed_ethnicity)
  for(v in rm06_binary_targets())
    data.table::set(d,j=v,value=as.integer(as.character(d[[v]])))
  d[, family_history_cvd := as.integer(rm_source_family_proxy)]
  out <- rm01_derive(d,cutoff)
  columns<-m01qa_rm_vars()
  if(any(vapply(out[,..columns],function(x)anyNA(x)||any(!x %in% 0:1),logical(1))))
    stop("Incomplete/invalid completed binary RMs")
  for(year in c("2019","2025")) {
    h <- out[[paste0("rm_h_",year)]]
    t <- out[[paste0("rm_t_",year)]]
    any_rm <- out[[paste0("rm_any_",year)]]
    if(any(any_rm != as.integer(h == 1L | t == 1L)))
      stop("Inconsistent H/T/any composite RM for ",year)
  }
  observed<-rm01_derive(original,cutoff)
  for(v in columns) {
    known<-!is.na(observed[[v]])
    if(any(out[[v]][known]!=observed[[v]][known]))
      stop("Observed RM classification changed: ",v)
  }
  if (any(out$rm_inactive != out$physical_inactivity))
    stop("Physical-inactivity source and model indicator diverged")
  if(any(out$rm_premature_menopause[!female]!=0L))
    stop("Male premature menopause nonzero after MI")
  out
}
rm06_binary_rm_audit <- function(completed, original, imputation_id) {
  female<-completed$sex=="female"
  initial<-original$rm_premature_menopause
  final<-completed$rm_premature_menopause
  if(anyNA(final))stop("Premature menopause remains missing")
  if(any(final[!is.na(initial)]!=initial[!is.na(initial)]))
    stop("Known premature menopause status changed")
  data.table::data.table(
    imputation_id=imputation_id,
    n_female=sum(female),
    n_observed_positive=sum(female & initial==1L,na.rm=TRUE),
    n_observed_negative=sum(female & initial==0L,na.rm=TRUE),
    n_original_unknown=sum(female & is.na(initial)),
    n_imputed_positive=sum(female & is.na(initial) & final==1L),
    n_completed_positive=sum(female & final==1L),
    n_male_nonzero=sum(!female & final!=0L)
  )
}
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
# Appended to the sole 06 function module by the implementation script.
mi06_chain_long <- function(a,statistic) {
 if(is.null(a)||length(a)==0L)return(data.table())
 z<-as.data.table(as.table(a));setnames(z,c("variable","iteration","imputation_id","value"))
 z[,`:=`(iteration=as.integer(gsub("[^0-9]","",as.character(iteration))),imputation_id=as.integer(gsub("[^0-9]","",as.character(imputation_id))),statistic=statistic)]
 z
}
mi06_model_hash <- function(env=.GlobalEnv) {
 nms<-ls(env,all.names=TRUE)
 nms<-nms[grepl("^(mi06_|rm06_|rm01_|mice[.]impute[.]|m01qa_)",nms)]
 nms<-setdiff(nms,c("mi06_model_hash","mi06_check_model_snapshot","mi06_export_workbooks","mi06_write_xlsx","mi06_chain_long"))
 nms<-sort(nms[vapply(nms,function(n)is.function(get(n,envir=env)),logical(1))])
 code<-setNames(lapply(nms,function(n)paste(deparse(get(n,envir=env),width.cutoff=500L),collapse="\n")),nms)
 p<-tempfile(fileext=".rds");on.exit(unlink(p));saveRDS(code,p,version=2)
 mi06_md5(p)
}
mi06_check_model_snapshot <- function(clean_dir,code_dir) {
 snapshot<-file.path(clean_dir,"MI_M01/imputation_function_snapshot_M01.R")
 if(!file.exists(snapshot))stop("Missing saved model-function snapshot; cannot verify export against fitted models")
 e<-new.env(parent=.GlobalEnv);sys.source(snapshot,envir=e)
 if(!identical(mi06_model_hash(e),mi06_model_hash(.GlobalEnv)))
   stop("Statistical model functions changed; saved results cannot be relabelled using new models")
 invisible(snapshot)
}
mi06_write_xlsx <- function(tables,path,metadata) {
 if(!requireNamespace("openxlsx",quietly=TRUE))stop("Install openxlsx")
 wb<-openxlsx::createWorkbook();index<-data.table(sheet=names(tables),rows=vapply(tables,nrow,integer(1)))
 all<-c(list(result_index=index,run_metadata=metadata),tables)
 for(nm in names(all)) {
  z<-as.data.frame(all[[nm]]);if(nrow(z)>1048575L)stop("Excel limit exceeded: ",nm)
  openxlsx::addWorksheet(wb,nm)
  openxlsx::writeData(wb,nm,z,withFilter=nrow(z)>0)
  openxlsx::freezePane(wb,nm,firstRow=TRUE)
  openxlsx::setColWidths(wb,nm,cols=seq_len(max(1,ncol(z))),widths=22)
  openxlsx::addStyle(wb,nm,openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#203864",textDecoration="bold",wrapText=TRUE),rows=1,cols=seq_len(max(1,ncol(z))),gridExpand=TRUE)
 }
 tmp<-paste0(path,".tmp.xlsx");openxlsx::saveWorkbook(wb,tmp,overwrite=TRUE)
 if(!file.copy(tmp,path,overwrite=TRUE))stop("Cannot replace workbook: ",path)
 unlink(tmp)
}
mi06_export_workbooks <- function(clean_dir,code_dir) {
 out<-file.path(clean_dir,"analysis_excel");dir.create(out,recursive=TRUE,showWarnings=FALSE)
 # Saved statistical results are the source of truth, including export-only mode.
 qc<-file.path(clean_dir,"MI_M01/quality_control")
 required<-file.path(qc,c("MI_parameter_manifest_M01.csv","completed_dataset_contract_audit_M01.csv","MI_code_input_provenance_M01.csv"))
 if(any(!file.exists(required)))stop("Cannot export incomplete MI run: ",paste(required[!file.exists(required)],collapse=", "))
 audit<-fread(required[2]);if(any(!audit$passed))stop("MI contract audit failed")
 prov<-fread(required[3]);inputs<-prov[!component%in%c("main_06","function_06")]
 current<-vapply(inputs$path,mi06_md5,character(1))
 if(!identical(unname(current),as.character(inputs$md5)))stop("Saved MI results do not match current inputs/definitions; rerun analysis")
 snapshot<-mi06_check_model_snapshot(clean_dir,code_dir)
 if(!mi06_md5(snapshot)%in%prov[component=="function_06",md5])stop("Model snapshot does not match run provenance")
 pars<-fread(required[1]);meta<-rbindlist(list(data.table(parameter=c("RM_version","export_time","clean_directory","aggregation"),value=c("1.7",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),normalizePath(clean_dir,winslash="/"),"Each completed dataset summarized separately; no stacked denominator")),pars),fill=TRUE)
 paths<-list.files(clean_dir,pattern="[.]csv$",recursive=TRUE,full.names=TRUE)
 paths<-paths[!grepl("MI_completed/|progress/|analysis_excel/",gsub("\\\\","/",paths)) & !dirname(paths)==clean_dir]
 grouped<-list(cohort=list(),definitions=list(),rm=list(),mi=list(),distributions=list(),ses=list())
 links<-list()
 for(p in paths) {
  # Participant-level audit and prescription lists remain external, never truncated.
  if(file.info(p)$size>2e6) {links[[length(links)+1L]]<-data.table(path=normalizePath(p,winslash="/"),reason="Participant-level/large audit retained as CSV");next}
  z<-fread(p);base<-tools::file_path_sans_ext(basename(p))
  group<-if(grepl("SES|subgroup",base))"ses" else if(grepl("observed_vs_imputed|distribution|pre_post|variation",base))"distributions" else if(grepl("RM_",base))"rm" else if(grepl("MI_M01",p,fixed=TRUE))"mi" else "cohort"
  key<-substr(gsub("_M01$","",base),1,31);original<-key;j<-1L
  while(key%in%names(grouped[[group]])){j<-j+1L;key<-paste0(substr(original,1,27),"_",j)}
  grouped[[group]][[key]]<-z
 }
 for(nm in c("eligibility","outcome","variable"))grouped$definitions[[paste0(nm,"_definitions")]]<-fread(file.path(code_dir,"config",paste0(nm,"_definitions_v1.7.csv")))
 grouped$definitions$references<-as.data.table(read.csv(file.path(code_dir,"config/references_v1.7.csv"),fileEncoding="UTF-8-BOM",check.names=FALSE))
 grouped$rm$RM_domains<-fread(file.path(code_dir,"config/RM_domains_v1.7.csv"))
 grouped$rm$RM_unavailable<-data.table(year=c(2019,2019,2019,2019,2025,2025,2025),domain=c("CAC","Carotid/femoral plaque","LVH","CKD","CAC","Carotid/femoral plaque","Coronary plaque"),status=c(rep("No validated baseline phenotype; not imputed",3),"Eligibility exclusion; not upgrading RM",rep("No validated baseline phenotype; not imputed",3)))
 if(length(links))grouped$mi$external_data_files<-rbindlist(links)
 chain<-rbindlist(lapply(c("mean","variance"),function(nm){p<-file.path(qc,paste0("mice_chain_",nm,"_M01.rds"));if(file.exists(p))mi06_chain_long(readRDS(p),nm)}),fill=TRUE)
 grouped$mi$chain_diagnostics<-chain
 # Prespecified mean and min/max of per-imputation summaries, not pooled people.
 p<-file.path(qc,"observed_vs_imputed_distribution_M01.csv")
 d<-fread(p)
 grouped$distributions$observed_imputed_summary<-d[,.(n_imputations=uniqueN(imputation_id),mean_value=mean(value),min_value=min(value),max_value=max(value)),by=.(variable,source,statistic,level)]
 methods<-fread(file.path(qc,"imputation_method_M01.csv"))
 methods[,n:=as.integer(pars[parameter=="analysis_n",value])]
 methods[,pct_missing:=100*n_missing/n]
 methods[,model_description:=fcase(method=="pmm","Predictive mean matching",method=="logreg","Logistic regression",method=="polr","Proportional odds regression",method=="polyreg","Multinomial logistic regression",method=="rm_crp","Interval-censored log-normal regression with posterior draws",method=="rm_tg","Interval-censored log-normal regression with posterior draws",method=="rm_south_asian","Logistic regression restricted to non-White completed ethnicity",method=="rm_female_hdp","Female-only logistic regression; male structural zero",method=="","Not imputed; predictor or fixed value",default=method)]
 methods[,post_MI_derivation:=fcase(variable=="townsend_deprivation_index","Social deprivation using frozen eligible-cohort Q80",variable%in%c("bmi","waist_circumference"),"One combined obesity RM (BMI OR sex-specific waist)",variable=="rm_source_crp","hs-CRP >2 mg/L in 2025; >10 sensitivity flag",variable%in%c("education","household_income","employment_status"),"Two-class SES LCA in each completed dataset",variable=="rm_lpa_status","Shared Lp(a) RM in both years",variable=="rm_source_inactive","Physical inactivity (Field 22035, binary)",default="Retained or mapped to RM domain")]
 grouped$mi$variable_imputation_models<-methods
 lca<-fread(file.path(qc,"SES_LCA_model_summary_M01.csv"))
 grouped$ses$SES_model_across_imputations<-lca[,.(mean_value=mean(value),min_value=min(value),max_value=max(value)),by=metric]
 classes<-fread(file.path(qc,"SES_LCA_class_distribution_M01.csv"))
 grouped$ses$SES_class_across_imputations<-classes[,.(percent_mean=mean(percent),percent_min=min(percent),percent_max=max(percent)),by=ses_class]
 profiles<-fread(file.path(qc,"SES_LCA_component_profile_M01.csv"))
 grouped$ses$SES_profile_across_imputations<-profiles[,.(percent_mean=mean(percent),percent_min=min(percent),percent_max=max(percent)),by=.(ses_class,variable,value)]
 membership<-fread(file.path(qc,"subgroup_membership_stability_summary_M01.csv"))
 total_m<-as.integer(pars[parameter=="M",value])
 membership[,pattern:=fifelse(count==0,"Stable negative/high SES",fifelse(count==total_m,"Stable positive/low SES","Varied across imputations"))]
 grouped$ses$subgroup_stability_categories<-membership[,.(n=sum(N),percent=sum(percent)),by=.(variable,pattern)]
 filenames<-c(cohort="M01_人群纳排与随访审计_v1.7.xlsx",definitions="M01_疾病排除与ASCVD结局定义_v1.7.xlsx",rm="M01_变量来源与RM定义_v1.7.xlsx",mi="M01_多重插补模型与诊断_v1.7.xlsx",distributions="M01_插补前后分布与RM结果_v1.7.xlsx",ses="M01_SES潜类模型与亚组稳定性_v1.7.xlsx")
 manifest<-rbindlist(lapply(names(filenames),function(g){
  mi06_write_xlsx(grouped[[g]],file.path(out,filenames[[g]]),meta)
  data.table(group=g,file=filenames[[g]],sheet=names(grouped[[g]]))
 }))
 fwrite(manifest,file.path(out,"workbook_manifest_M01.csv"),bom=TRUE)
 invisible(manifest)
}

