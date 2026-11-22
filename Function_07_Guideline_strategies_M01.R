# MODULE 07 | shared utilities + preserved SCORE/SCORE2 + scenario engine
suppressPackageStartupMessages(library(data.table))

.resolve_score2019_input <- function(data, x, arg_name) {
  if (is.null(data)) {
    return(x)
  }

  if (!is.character(x) || length(x) != 1L) {
    stop(sprintf(
      "When `data` is provided, `%s` must be a single column name string.",
      arg_name
    ))
  }

  if (!(x %in% names(data))) {
    stop(sprintf(
      "Column `%s` specified in `%s` was not found in `data`.",
      x, arg_name
    ))
  }

  data[[x]]
}

# ---- .recycle_score2019_inputs: retained shared primitive ----
.recycle_score2019_inputs <- function(...) {
  xs <- list(...)
  lens <- vapply(xs, length, integer(1))
  n <- max(lens)

  if (any(lens != 1L & lens != n)) {
    stop("All resolved inputs must have the same length or length 1.")
  }

  lapply(xs, rep, length.out = n)
}

# ---- normalize_score2019_sex: retained shared primitive ----
normalize_score2019_sex <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  out <- rep(NA_character_, length(x))

  if (is.character(x)) {
    x2 <- trimws(tolower(x))
    out[x2 %in% c("male", "man", "men", "m", "1")] <- "men"
    out[x2 %in% c("female", "woman", "women", "f", "0")] <- "women"
  }

  idx_na <- is.na(out) & !is.na(x)
  if (any(idx_na)) {
    x_num <- suppressWarnings(as.numeric(x))
    out[idx_na & !is.na(x_num) & x_num == 1] <- "men"
    out[idx_na & !is.na(x_num) & x_num == 0] <- "women"
  }

  out
}

# ---- normalize_score2019_smoking: retained shared primitive ----
normalize_score2019_smoking <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  out <- rep(NA_character_, length(x))

  if (is.character(x)) {
    x2 <- trimws(tolower(x))
    out[x2 %in% c("1", "yes", "y", "true", "current", "smoker", "smoking")] <- "smoker"
    out[x2 %in% c(
      "0", "no", "n", "false", "never", "former", "previous",
      "non-smoker", "nonsmoker", "non smoker"
    )] <- "non-smoker"
  }

  idx_na <- is.na(out) & !is.na(x)
  if (any(idx_na)) {
    x_num <- suppressWarnings(as.numeric(x))
    out[idx_na & !is.na(x_num) & x_num == 1] <- "smoker"
    out[idx_na & !is.na(x_num) & x_num == 0] <- "non-smoker"
  }

  out
}

# ---- normalize_score2019_region: retained shared primitive ----
normalize_score2019_region <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  out <- rep(NA_character_, length(x))

  if (is.character(x)) {
    x2 <- trimws(tolower(x))
    out[x2 %in% c("low", "low risk", "low-risk")] <- "low"
    out[x2 %in% c("high", "high risk", "high-risk")] <- "high"
  }

  idx_na <- is.na(out) & !is.na(x)
  if (any(idx_na)) {
    x_num <- suppressWarnings(as.numeric(x))
    out[idx_na & !is.na(x_num) & x_num == 0] <- "low"
    out[idx_na & !is.na(x_num) & x_num == 1] <- "high"
  }

  out
}

# ---- .score2019_map_nearest: retained shared primitive ----
.score2019_map_nearest <- function(x, breaks, tie = c("upper", "lower"), cap = TRUE) {
  tie <- match.arg(tie)

  if (!is.numeric(x)) stop("`x` must be numeric.")
  if (!is.numeric(breaks) || length(breaks) < 2L) {
    stop("`breaks` must be a numeric vector with length >= 2.")
  }
  if (is.unsorted(breaks, strictly = TRUE)) {
    stop("`breaks` must be strictly increasing.")
  }
  if (!is.logical(cap) || length(cap) != 1L || is.na(cap)) {
    stop("`cap` must be TRUE or FALSE.")
  }

  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)

  if (!any(ok)) {
    return(out)
  }

  x_work <- x

  if (cap) {
    x_work[ok] <- pmin(pmax(x_work[ok], min(breaks)), max(breaks))
  } else {
    out_of_range <- ok & (x_work < min(breaks) | x_work > max(breaks))
    ok[out_of_range] <- FALSE
  }

  mids <- (breaks[-length(breaks)] + breaks[-1]) / 2

  if (tie == "upper") {
    idx <- findInterval(x_work[ok], mids) + 1L
  } else {
    # 等于中点时落入下侧
    idx <- findInterval(x_work[ok], mids, left.open = TRUE) + 1L
  }

  idx <- pmax(1L, pmin(idx, length(breaks)))
  out[ok] <- breaks[idx]
  out
}

# ---- score2019_map_age: retained shared primitive ----
score2019_map_age <- function(age,
                              age_breaks = c(40, 50, 55, 60, 65),
                              tie = c("upper", "lower"),
                              cap = TRUE) {
  .score2019_map_nearest(x = age, breaks = age_breaks, tie = tie, cap = cap)
}

# ---- score2019_map_sbp: retained shared primitive ----
score2019_map_sbp <- function(sbp,
                              sbp_breaks = c(120, 140, 160, 180),
                              tie = c("upper", "lower"),
                              cap = TRUE) {
  .score2019_map_nearest(x = sbp, breaks = sbp_breaks, tie = tie, cap = cap)
}

# ---- score2019_map_tc: retained shared primitive ----
score2019_map_tc <- function(tc,
                             tc_breaks = c(4, 5, 6, 7),
                             tie = c("upper", "lower"),
                             cap = TRUE) {
  .score2019_map_nearest(x = tc, breaks = tc_breaks, tie = tie, cap = cap)
}

# ---- score2019_chart_db: retained shared primitive ----
score2019_chart_db <- function() {
  if (exists("chart_dt", envir = .score2019_cache_env, inherits = FALSE)) {
    return(get("chart_dt", envir = .score2019_cache_env, inherits = FALSE))
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required for SCORE_2019_calculator.R.")
  }

  data.table <- asNamespace("data.table")

  make_dt <- function(values_by_age_desc, region_type, sex, smoking) {
    age_levels_desc <- c("70", "65", "60", "55", "50", "40")
    sbp_levels_desc <- c(180, 160, 140, 120)
    tc_levels <- c(4, 5, 6, 7)

    out_list <- vector("list", length(age_levels_desc))

    for (i in seq_along(age_levels_desc)) {
      age_chr <- age_levels_desc[i]
      mat <- values_by_age_desc[[age_chr]]

      dt <- data.table$as.data.table(mat)
      data.table$setnames(dt, paste0("V", seq_len(ncol(dt))), paste0("tc_", tc_levels))
      dt[, sbp_group := sbp_levels_desc]
      dt[, age_group := as.numeric(age_chr)]

      dt <- data.table$melt(
        dt,
        id.vars = c("age_group", "sbp_group"),
        variable.name = "tc_group",
        value.name = "risk"
      )

      dt[, tc_group := as.numeric(sub("^tc_", "", tc_group))]
      dt[, `:=`(
        region_type = region_type,
        sex = sex,
        smoking = smoking
      )]

      out_list[[i]] <- dt
    }

    data.table$rbindlist(out_list, use.names = TRUE)
  }

  # ---------------------------------------------------------------------------.
  # Low-risk chart  ####
  # ---------------------------------------------------------------------------.
  low_women_nonsmoker <- list(
    "70" = rbind(c(7, 8, 8, 9), c(6, 6, 7, 7), c(5, 5, 6, 6), c(4, 4, 5, 5)),
    "65" = rbind(c(4, 4, 5, 5), c(3, 3, 4, 4), c(2, 3, 3, 3), c(2, 2, 2, 2)),
    "60" = rbind(c(2, 3, 3, 3), c(2, 2, 2, 2), c(1, 1, 1, 2), c(1, 1, 1, 1)),
    "55" = rbind(c(1, 1, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1), c(0, 0, 0, 1)),
    "50" = rbind(c(1, 1, 1, 1), c(0, 0, 1, 1), c(0, 0, 0, 0), c(0, 0, 0, 0)),
    "40" = rbind(c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0))
  )

  low_women_smoker <- list(
    "70" = rbind(c(11, 11, 12, 13), c(9, 9, 10, 11), c(7, 8, 8, 9), c(6, 6, 7, 7)),
    "65" = rbind(c(7, 7, 8, 9), c(5, 6, 6, 7), c(4, 4, 5, 5), c(3, 3, 3, 4)),
    "60" = rbind(c(4, 5, 5, 6), c(3, 3, 4, 4), c(2, 2, 3, 3), c(2, 2, 2, 2)),
    "55" = rbind(c(3, 3, 3, 4), c(2, 2, 2, 3), c(1, 1, 1, 2), c(1, 1, 1, 1)),
    "50" = rbind(c(2, 2, 2, 3), c(1, 1, 1, 2), c(1, 1, 1, 1), c(0, 0, 0, 1)),
    "40" = rbind(c(1, 1, 1, 1), c(0, 0, 0, 1), c(0, 0, 0, 0), c(0, 0, 0, 0))
  )

  low_men_nonsmoker <- list(
    "70" = rbind(c(12, 14, 15, 17), c(10, 11, 13, 14), c(8, 9, 10, 12), c(7, 8, 9, 10)),
    "65" = rbind(c(8, 9, 10, 12), c(6, 7, 8, 9), c(5, 5, 6, 7), c(3, 4, 5, 5)),
    "60" = rbind(c(5, 6, 7, 8), c(4, 4, 5, 5), c(3, 3, 3, 4), c(2, 2, 2, 3)),
    "55" = rbind(c(3, 4, 4, 5), c(2, 2, 3, 3), c(1, 2, 2, 2), c(1, 1, 1, 2)),
    "50" = rbind(c(2, 2, 3, 3), c(1, 1, 2, 2), c(1, 1, 1, 1), c(0, 1, 1, 1)),
    "40" = rbind(c(1, 1, 1, 1), c(0, 0, 1, 1), c(0, 0, 0, 0), c(0, 0, 0, 0))
  )

  low_men_smoker <- list(
    "70" = rbind(c(18, 20, 22, 24), c(15, 16, 18, 20), c(12, 13, 15, 17), c(10, 11, 12, 14)),
    "65" = rbind(c(12, 14, 16, 18), c(9, 11, 12, 14), c(7, 8, 9, 11), c(5, 6, 7, 8)),
    "60" = rbind(c(8, 10, 11, 13), c(6, 7, 8, 9), c(4, 5, 6, 7), c(3, 4, 4, 5)),
    "55" = rbind(c(6, 7, 8, 9), c(4, 4, 5, 6), c(3, 3, 3, 4), c(2, 2, 2, 3)),
    "50" = rbind(c(4, 5, 5, 6), c(2, 3, 3, 4), c(1, 2, 2, 3), c(1, 1, 1, 2)),
    "40" = rbind(c(2, 2, 3, 3), c(1, 1, 1, 2), c(1, 1, 1, 1), c(0, 0, 0, 1))
  )

  # ---------------------------------------------------------------------------.
  # High-risk chart ####
  # ---------------------------------------------------------------------------.
  high_women_nonsmoker <- list(
    "70" = rbind(c(12, 13, 14, 15), c(10, 11, 12, 13), c(8, 9, 10, 10), c(7, 7, 8, 9)),
    "65" = rbind(c(7, 8, 8, 9), c(5, 6, 6, 7), c(4, 4, 5, 5), c(3, 3, 4, 4)),
    "60" = rbind(c(4, 4, 5, 5), c(3, 3, 3, 4), c(2, 2, 2, 3), c(1, 1, 2, 2)),
    "55" = rbind(c(2, 2, 3, 3), c(1, 2, 2, 2), c(1, 1, 1, 1), c(1, 1, 1, 1)),
    "50" = rbind(c(1, 1, 2, 2), c(1, 1, 1, 1), c(0, 0, 1, 1), c(0, 0, 0, 0)),
    "40" = rbind(c(0, 0, 1, 1), c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0))
  )

  high_women_smoker <- list(
    "70" = rbind(c(17, 19, 20, 21), c(14, 15, 16, 18), c(12, 13, 14, 15), c(10, 10, 11, 12)),
    "65" = rbind(c(11, 12, 13, 15), c(9, 9, 10, 11), c(7, 7, 8, 9), c(5, 5, 6, 7)),
    "60" = rbind(c(7, 8, 9, 10), c(5, 6, 6, 7), c(4, 4, 4, 5), c(3, 3, 3, 3)),
    "55" = rbind(c(5, 5, 6, 7), c(3, 3, 4, 4), c(2, 2, 2, 3), c(1, 1, 2, 2)),
    "50" = rbind(c(3, 3, 4, 4), c(2, 2, 2, 3), c(1, 1, 1, 2), c(1, 1, 1, 1)),
    "40" = rbind(c(1, 1, 2, 2), c(1, 1, 1, 1), c(0, 0, 0, 1), c(0, 0, 0, 0))
  )

  high_men_nonsmoker <- list(
    "70" = rbind(c(24, 26, 30, 33), c(20, 22, 25, 28), c(16, 18, 21, 24), c(13, 15, 17, 20)),
    "65" = rbind(c(15, 17, 20, 23), c(12, 14, 16, 18), c(9, 11, 12, 14), c(7, 8, 10, 11)),
    "60" = rbind(c(10, 11, 13, 15), c(7, 8, 10, 11), c(5, 6, 7, 8), c(4, 4, 5, 6)),
    "55" = rbind(c(6, 7, 9, 10), c(4, 5, 6, 7), c(3, 3, 4, 5), c(2, 2, 3, 3)),
    "50" = rbind(c(4, 5, 6, 7), c(2, 3, 3, 4), c(2, 2, 2, 3), c(1, 1, 1, 2)),
    "40" = rbind(c(2, 2, 2, 3), c(1, 1, 1, 2), c(0, 1, 1, 1), c(0, 0, 0, 1))
  )

  high_men_smoker <- list(
    "70" = rbind(c(33, 36, 40, 45), c(27, 31, 34, 39), c(23, 26, 29, 33), c(19, 22, 25, 28)),
    "65" = rbind(c(23, 26, 30, 34), c(18, 21, 24, 27), c(14, 16, 19, 22), c(11, 13, 15, 17)),
    "60" = rbind(c(16, 19, 22, 25), c(12, 14, 16, 19), c(9, 10, 12, 14), c(6, 7, 9, 10)),
    "55" = rbind(c(11, 13, 16, 18), c(8, 9, 11, 13), c(5, 6, 7, 9), c(4, 4, 5, 6)),
    "50" = rbind(c(8, 9, 11, 13), c(5, 6, 7, 9), c(3, 4, 5, 6), c(2, 2, 3, 4)),
    "40" = rbind(c(4, 4, 5, 7), c(2, 2, 3, 4), c(1, 1, 2, 2), c(1, 1, 1, 1))
  )

  chart_dt <- data.table$rbindlist(list(
    make_dt(low_women_nonsmoker, "low", "women", "non-smoker"),
    make_dt(low_women_smoker, "low", "women", "smoker"),
    make_dt(low_men_nonsmoker, "low", "men", "non-smoker"),
    make_dt(low_men_smoker, "low", "men", "smoker"),
    make_dt(high_women_nonsmoker, "high", "women", "non-smoker"),
    make_dt(high_women_smoker, "high", "women", "smoker"),
    make_dt(high_men_nonsmoker, "high", "men", "non-smoker"),
    make_dt(high_men_smoker, "high", "men", "smoker")
  ), use.names = TRUE)

  data.table$setcolorder(
    chart_dt,
    c("region_type", "sex", "smoking", "age_group", "sbp_group", "tc_group", "risk")
  )

  chart_dt[, lookup_key := paste(region_type, sex, smoking, age_group, sbp_group, tc_group, sep = "|")]
  data.table$setkey(chart_dt, lookup_key)

  assign("chart_dt", chart_dt, envir = .score2019_cache_env)
  chart_dt
}

# ---- score2019_lookup_risk_discrete: retained shared primitive ----
score2019_lookup_risk_discrete <- function(region_type,
                                           sex,
                                           smoking,
                                           age_group,
                                           sbp_group,
                                           tc_group) {
  xs <- .recycle_score2019_inputs(
    region_type = region_type,
    sex = sex,
    smoking = smoking,
    age_group = age_group,
    sbp_group = sbp_group,
    tc_group = tc_group
  )

  region_type <- xs$region_type
  sex <- xs$sex
  smoking <- xs$smoking
  age_group <- xs$age_group
  sbp_group <- xs$sbp_group
  tc_group <- xs$tc_group

  valid_region <- c("low", "high")
  valid_sex <- c("women", "men")
  valid_smoke <- c("non-smoker", "smoker")
  # 官方原始图表数据库保留 70 岁行以便核对，但本研究查表入口不允许访问该行。
  valid_age <- c(40, 50, 55, 60, 65)
  valid_sbp <- c(120, 140, 160, 180)
  valid_tc <- c(4, 5, 6, 7)

  ok <- !(is.na(region_type) | is.na(sex) | is.na(smoking) |
    is.na(age_group) | is.na(sbp_group) | is.na(tc_group))

  if (any(ok & !(region_type %in% valid_region))) {
    stop("`region_type` must be 'low' or 'high'.")
  }
  if (any(ok & !(sex %in% valid_sex))) {
    stop("`sex` must be 'women' or 'men'.")
  }
  if (any(ok & !(smoking %in% valid_smoke))) {
    stop("`smoking` must be 'non-smoker' or 'smoker'.")
  }
  if (any(ok & !(age_group %in% valid_age))) {
    stop("For this study, `age_group` must be one of 40, 50, 55, 60, 65; age 70 is disabled.")
  }
  if (any(ok & !(sbp_group %in% valid_sbp))) {
    stop("`sbp_group` must be one of 120, 140, 160, 180.")
  }
  if (any(ok & !(tc_group %in% valid_tc))) {
    stop("`tc_group` must be one of 4, 5, 6, 7.")
  }

  out <- rep(NA_real_, length(region_type))
  if (!any(ok)) {
    return(out)
  }

  chart_dt <- score2019_chart_db()
  key_vec <- paste(region_type[ok], sex[ok], smoking[ok], age_group[ok], sbp_group[ok], tc_group[ok], sep = "|")
  idx <- match(key_vec, chart_dt$lookup_key)
  out[ok] <- chart_dt$risk[idx]

  out
}

# ---- SCORE_2019: retained shared primitive ----
SCORE_2019 <- function(data = NULL,
                       sex,
                       age,
                       smoking,
                       sbp,
                       tc,
                       risk_region,
                       age_breaks = c(40, 50, 55, 60, 65),
                       sbp_breaks = c(120, 140, 160, 180),
                       tc_breaks = c(4, 5, 6, 7),
                       tie = c("upper", "lower"),
                       cap = TRUE,
                       return_details = FALSE) {
  tie <- match.arg(tie)

  if (!is.logical(return_details) || length(return_details) != 1L || is.na(return_details)) {
    stop("`return_details` must be TRUE or FALSE.")
  }

  sex_vec <- .resolve_score2019_input(data, sex, "sex")
  age_vec <- .resolve_score2019_input(data, age, "age")
  smoke_vec <- .resolve_score2019_input(data, smoking, "smoking")
  sbp_vec <- .resolve_score2019_input(data, sbp, "sbp")
  tc_vec <- .resolve_score2019_input(data, tc, "tc")
  region_vec <- .resolve_score2019_input(data, risk_region, "risk_region")

  xs <- .recycle_score2019_inputs(
    sex = sex_vec,
    age = age_vec,
    smoking = smoke_vec,
    sbp = sbp_vec,
    tc = tc_vec,
    risk_region = region_vec
  )

  sex_vec <- xs$sex
  age_vec <- xs$age
  smoke_vec <- xs$smoking
  sbp_vec <- xs$sbp
  tc_vec <- xs$tc
  region_vec <- xs$risk_region

  sex_std <- normalize_score2019_sex(sex_vec)
  smoke_std <- normalize_score2019_smoking(smoke_vec)
  region_std <- normalize_score2019_region(region_vec)

  age_group <- score2019_map_age(
    age = suppressWarnings(as.numeric(age_vec)),
    age_breaks = age_breaks,
    tie = tie,
    cap = cap
  )

  if (any(age_group > 65, na.rm = TRUE)) {
    stop("本研究 SCORE chart 年龄格点不得超过 65 岁。", call. = FALSE)
  }

  sbp_group <- score2019_map_sbp(
    sbp = suppressWarnings(as.numeric(sbp_vec)),
    sbp_breaks = sbp_breaks,
    tie = tie,
    cap = cap
  )

  tc_group <- score2019_map_tc(
    tc = suppressWarnings(as.numeric(tc_vec)),
    tc_breaks = tc_breaks,
    tie = tie,
    cap = cap
  )

  risk <- score2019_lookup_risk_discrete(
    region_type = region_std,
    sex = sex_std,
    smoking = smoke_std,
    age_group = age_group,
    sbp_group = sbp_group,
    tc_group = tc_group
  )

  if (!return_details) {
    return(risk)
  }

  data.frame(
    sex_input = sex_vec,
    age = age_vec,
    smoking_input = smoke_vec,
    sbp = sbp_vec,
    tc = tc_vec,
    risk_region_input = region_vec,
    sex = sex_std,
    smoking = smoke_std,
    region_type = region_std,
    age_score2019 = age_group,
    sbp_score2019 = sbp_group,
    tc_score2019 = tc_group,
    score2019_risk = risk,
    stringsAsFactors = FALSE
  )
}

# ---- assert_file_exists: retained shared primitive ----
assert_file_exists <- function(path, label = "文件") {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop(label, "不存在：", path, call. = FALSE)
  }
  invisible(path)
}

# ---- assert_required_columns: retained shared primitive ----
assert_required_columns <- function(dt, cols, label = "数据表") {
  dt <- as.data.table(dt)
  missing_cols <- setdiff(unique(cols), names(dt))
  if (length(missing_cols) > 0L) {
    stop(label, "缺少关键列：", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- assert_unique_id: retained shared primitive ----
assert_unique_id <- function(dt, id_col = "participant_id", label = "数据表") {
  assert_required_columns(dt, id_col, label)
  if (anyNA(dt[[id_col]])) stop(label, "的 ID 存在缺失。", call. = FALSE)
  n_dup <- nrow(dt) - uniqueN(dt[[id_col]])
  if (n_dup > 0L) stop(label, "的 ID 存在重复，重复行数：", n_dup, call. = FALSE)
  invisible(TRUE)
}

# ---- read_score_analysis_inputs: retained shared primitive ----
read_score_analysis_inputs <- function(data_path, dictionary_path, eligibility_audit_path) {
  assert_file_exists(data_path, "05 分析基础数据")
  assert_file_exists(dictionary_path, "05 变量字典")
  assert_file_exists(eligibility_audit_path, "05 纳排常数审计")

  data <- fread(data_path)
  dictionary <- fread(dictionary_path)
  eligibility_audit <- fread(eligibility_audit_path)

  setDT(data)
  setDT(dictionary)
  setDT(eligibility_audit)
  assert_unique_id(data, "participant_id", "05 分析基础数据")

  list(
    data = data,
    dictionary = dictionary,
    eligibility_audit = eligibility_audit
  )
}

# ---- validate_eligibility_audit: retained shared primitive ----
validate_eligibility_audit <- function(audit) {
  audit <- as.data.table(audit)
  assert_required_columns(audit, c("variable", "passed"), "05 纳排常数审计")
  passed <- tolower(as.character(audit$passed)) %in% c("true", "1")
  if (any(!passed | is.na(passed))) {
    bad <- audit[!passed | is.na(passed), variable]
    stop("05 纳排常数审计未通过：", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- validate_score_unit_contract: retained shared primitive ----
validate_score_unit_contract <- function(dictionary) {
  dictionary <- copy(as.data.table(dictionary))
  assert_required_columns(dictionary, c("standard_name", "unit"), "05 变量字典")

  expected <- data.table(
    standard_name = c(
      "age", "systolic_blood_pressure", "total_cholesterol", "hdl_cholesterol"
    ),
    expected_unit = c("years", "mmHg", "mmol/L", "mmol/L")
  )

  observed <- dictionary[expected, on = "standard_name", .(
    standard_name = i.standard_name,
    expected_unit = i.expected_unit,
    observed_unit = unit
  )]
  observed[, passed := !is.na(observed_unit) & observed_unit == expected_unit]

  if (any(!observed$passed)) {
    bad <- observed[passed == FALSE, paste0(standard_name, ": ", observed_unit, " != ", expected_unit)]
    stop("评分输入单位契约不一致：", paste(bad, collapse = "; "), call. = FALSE)
  }

  observed[]
}

# ---- validate_score_inputs: retained shared primitive ----
validate_score_inputs <- function(dt) {
  dt <- as.data.table(dt)
  required <- c(
    "participant_id", "age", "sex", "smoking_status",
    "systolic_blood_pressure", "total_cholesterol", "hdl_cholesterol"
  )
  assert_required_columns(dt, required, "SCORE/SCORE2 输入数据")
  assert_unique_id(dt, "participant_id", "SCORE/SCORE2 输入数据")

  allowed_sex <- c("female", "male")
  allowed_smoking <- c("never", "former", "current")
  bad_sex <- setdiff(unique(na.omit(as.character(dt$sex))), allowed_sex)
  bad_smoking <- setdiff(unique(na.omit(as.character(dt$smoking_status))), allowed_smoking)
  if (length(bad_sex) > 0L) stop("sex 存在未知分类：", paste(bad_sex, collapse = ", "), call. = FALSE)
  if (length(bad_smoking) > 0L) stop("smoking_status 存在未知分类：", paste(bad_smoking, collapse = ", "), call. = FALSE)

  missing_n <- vapply(dt[, ..required], function(x) sum(is.na(x)), integer(1))
  if (any(missing_n > 0L)) {
    stop(
      "评分输入存在缺失：",
      paste0(names(missing_n)[missing_n > 0L], "=", missing_n[missing_n > 0L], collapse = ", "),
      call. = FALSE
    )
  }

  range_rules <- data.table(
    variable = c("age", "systolic_blood_pressure", "total_cholesterol", "hdl_cholesterol"),
    lower = c(40, 50, 1, 0.1),
    upper = c(69, 250, 15, 5)
  )
  range_audit <- rbindlist(lapply(seq_len(nrow(range_rules)), function(i) {
    v <- range_rules$variable[i]
    x <- suppressWarnings(as.numeric(dt[[v]]))
    data.table(
      variable = v,
      n = length(x),
      n_missing = sum(is.na(x)),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      expected_lower = range_rules$lower[i],
      expected_upper = range_rules$upper[i],
      n_outside_expected = sum(x < range_rules$lower[i] | x > range_rules$upper[i], na.rm = TRUE)
    )
  }))

  if (any(range_audit$n_outside_expected > 0L)) {
    bad <- range_audit[n_outside_expected > 0L, variable]
    stop("评分输入超出预设合理范围：", paste(bad, collapse = ", "), call. = FALSE)
  }

  range_audit[]
}

# ---- validate_score_eligibility_boundaries: retained shared primitive ----
validate_score_eligibility_boundaries <- function(dt) {
  dt <- as.data.table(dt)
  required <- c(
    "age", "systolic_blood_pressure", "diastolic_blood_pressure",
    "total_cholesterol", "ldl_cholesterol"
  )
  assert_required_columns(dt, required, "SCORE/SCORE2 适用人群边界检查")

  audit <- rbindlist(list(
    data.table(
      rule = "Age must be 40-69 years",
      variable = "age",
      threshold = "40 <= age <= 69",
      observed_min = min(dt$age, na.rm = TRUE),
      observed_max = max(dt$age, na.rm = TRUE),
      n_violation = sum(dt$age < 40 | dt$age > 69, na.rm = TRUE)
    ),
    data.table(
      rule = "Markedly elevated SBP excluded",
      variable = "systolic_blood_pressure",
      threshold = "SBP < 180 mmHg",
      observed_min = min(dt$systolic_blood_pressure, na.rm = TRUE),
      observed_max = max(dt$systolic_blood_pressure, na.rm = TRUE),
      n_violation = sum(dt$systolic_blood_pressure >= 180, na.rm = TRUE)
    ),
    data.table(
      rule = "Markedly elevated DBP excluded",
      variable = "diastolic_blood_pressure",
      threshold = "DBP < 110 mmHg",
      observed_min = min(dt$diastolic_blood_pressure, na.rm = TRUE),
      observed_max = max(dt$diastolic_blood_pressure, na.rm = TRUE),
      n_violation = sum(dt$diastolic_blood_pressure >= 110, na.rm = TRUE)
    ),
    data.table(
      rule = "Extreme total cholesterol excluded",
      variable = "total_cholesterol",
      threshold = "TC <= 8 mmol/L",
      observed_min = min(dt$total_cholesterol, na.rm = TRUE),
      observed_max = max(dt$total_cholesterol, na.rm = TRUE),
      n_violation = sum(dt$total_cholesterol > 8, na.rm = TRUE)
    ),
    data.table(
      rule = "Markedly elevated LDL-C excluded",
      variable = "ldl_cholesterol",
      threshold = "LDL-C <= 4.9 mmol/L",
      observed_min = min(dt$ldl_cholesterol, na.rm = TRUE),
      observed_max = max(dt$ldl_cholesterol, na.rm = TRUE),
      n_violation = sum(dt$ldl_cholesterol > 4.9, na.rm = TRUE)
    )
  ), use.names = TRUE)
  audit[, passed := n_violation == 0L]

  if (any(!audit$passed)) {
    bad <- audit[passed == FALSE, paste0(rule, " (n=", n_violation, ")")]
    stop(
      "存在不符合 SCORE/SCORE2 共同适用边界的参与者：",
      paste(bad, collapse = "; "),
      call. = FALSE
    )
  }
  audit[]
}

# ---- derive_current_smoker: retained shared primitive ----
derive_current_smoker <- function(dt) {
  dt <- copy(as.data.table(dt))
  dt[, current_smoker := fcase(
    smoking_status == "current", 1L,
    smoking_status %chin% c("never", "former"), 0L,
    default = NA_integer_
  )]
  if (anyNA(dt$current_smoker)) stop("current_smoker 重分类后仍存在缺失。", call. = FALSE)
  dt[]
}

# ---- calculate_score2_percent_parallel: retained shared primitive ----
calculate_score2_percent_parallel <- function(dt, risk_region = "Low", n_workers = 2L) {
  dt <- as.data.table(dt)
  n <- nrow(dt)
  if (n == 0L) return(numeric())

  score_input <- list(
    age = as.numeric(dt$age),
    sex = as.character(dt$sex),
    smoker = as.numeric(dt$current_smoker),
    sbp = as.numeric(dt$systolic_blood_pressure),
    diabetes = rep(0, n),
    tc = as.numeric(dt$total_cholesterol),
    hdl = as.numeric(dt$hdl_cholesterol),
    region = rep(risk_region, n)
  )

  worker_fun <- function(index, input) {
    vapply(index, function(i) {
      as.numeric(RiskScorescvd::SCORE2(
        Risk.region = input$region[i],
        Age = input$age[i],
        Gender = input$sex[i],
        smoker = input$smoker[i],
        systolic.bp = input$sbp[i],
        diabetes = input$diabetes[i],
        total.chol = input$tc[i],
        total.hdl = input$hdl[i],
        classify = FALSE
      ))[1L]
    }, numeric(1))
  }

  n_workers <- max(1L, min(as.integer(n_workers), 14L, n))
  n_chunks <- min(n, max(n_workers, n_workers * 40L))
  chunks <- split(seq_len(n), cut(seq_len(n), breaks = n_chunks, labels = FALSE))

  if (n_workers == 1L) {
    parts <- vector("list",length(chunks))
    for(i in seq_along(chunks)){
      parts[[i]]<-worker_fun(chunks[[i]],score_input)
      cat(sprintf("[SCORE2] completed chunks %d/%d | %s\n",i,length(chunks),format(Sys.time(),"%H:%M:%S")))
    }
  } else {
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, .libPaths, .libPaths())
    parallel::clusterEvalQ(cl, library(RiskScorescvd))
    parallel::clusterExport(
      cl,
      varlist = c("score_input", "worker_fun"),
      envir = environment()
    )
    parts <- vector("list",length(chunks))
    waves <- split(seq_along(chunks),ceiling(seq_along(chunks)/n_workers))
    for(w in seq_along(waves)){
      ii <- waves[[w]]
      parts[ii] <- parallel::parLapply(cl,chunks[ii],function(index)worker_fun(index,score_input))
      cat(sprintf("[SCORE2] completed chunks %d/%d | %s\n",max(ii),length(chunks),format(Sys.time(),"%H:%M:%S")))
    }
  }

  out <- numeric(n)
  for (i in seq_along(chunks)) out[chunks[[i]]] <- parts[[i]]
  out
}

# ---- calculate_score_values: retained shared primitive ----
calculate_score_values <- function(dt,
                                   risk_region = "Low",
                                   n_workers = 2L,
                                   score_age_grid = c(40, 50, 55, 60, 65),
                                   score_sbp_grid = c(120, 140, 160, 180),
                                   score_tc_grid = c(4, 5, 6, 7),
                                   score_tie = "upper") {
  dt <- copy(as.data.table(dt))
  if (length(risk_region) != 1L || is.na(risk_region) || !nzchar(as.character(risk_region))) {
    stop("risk_region 必须是单个非缺失评分地区。", call. = FALSE)
  }
  if (!exists("SCORE_2019", mode = "function")) {
    stop("未加载 SCORE_2019() 评分引擎。", call. = FALSE)
  }

  score_2019_details <- SCORE_2019(
    data = NULL,
    sex = dt$sex,
    age = dt$age,
    smoking = dt$current_smoker,
    sbp = dt$systolic_blood_pressure,
    tc = dt$total_cholesterol,
    risk_region = rep(risk_region, nrow(dt)),
    age_breaks = score_age_grid,
    sbp_breaks = score_sbp_grid,
    tc_breaks = score_tc_grid,
    tie = score_tie,
    cap = TRUE,
    return_details = TRUE
  )

  dt[, `:=`(
    score_2019_percent = as.numeric(score_2019_details$score2019_risk),
    score_2019_age_chart = as.numeric(score_2019_details$age_score2019),
    score_2019_sbp_chart = as.numeric(score_2019_details$sbp_score2019),
    score_2019_tc_chart = as.numeric(score_2019_details$tc_score2019)
  )]

  if (any(dt$score_2019_age_chart > 65, na.rm = TRUE)) {
    stop("SCORE 2019 计算错误：年龄不得读取 70 岁图表行。", call. = FALSE)
  }

  dt[, score2_2025_percent := calculate_score2_percent_parallel(
    dt,
    risk_region = risk_region,
    n_workers = n_workers
  )]

  # risk_region 是评分算法地区；region 是 assessment-centre 地区，必须分别保留。
  dt[, risk_region := as.character(risk_region)]

  for (v in c("score_2019_percent", "score2_2025_percent")) {
    if (anyNA(dt[[v]])) stop(v, " 计算后存在缺失。", call. = FALSE)
    if (any(dt[[v]] < 0 | dt[[v]] > 100, na.rm = TRUE)) {
      stop(v, " 不在预期百分比范围 [0,100]。", call. = FALSE)
    }
    if (max(dt[[v]], na.rm = TRUE) <= 1) {
      stop(v, " 疑似返回比例尺度；本脚本要求百分比尺度，请核对评分函数。", call. = FALSE)
    }
  }

  dt[]
}

# ---- classify_score_percent: retained shared primitive ----
classify_score_percent <- function(x, cutpoints, levels = c("L", "M", "H", "VH")) {
  if (length(cutpoints) != 3L) stop("风险分类切点必须包含 3 个值。", call. = FALSE)
  out <- fcase(
    is.na(x), NA_character_,
    x < cutpoints[1], levels[1],
    x < cutpoints[2], levels[2],
    x < cutpoints[3], levels[3],
    default = levels[4]
  )
  factor(out, levels = levels, ordered = TRUE)
}

# ---- add_score_categories: retained shared primitive ----
add_score_categories <- function(dt,
                                 score_cutpoints = c(1, 5, 10),
                                 score2_cutpoints = c(2, 10, 20)) {
  dt <- copy(as.data.table(dt))
  risk_levels <- c("L", "M", "H", "VH")
  dt[, score_2019_category := classify_score_percent(score_2019_percent, score_cutpoints, risk_levels)]
  dt[, score2_2025_category := classify_score_percent(score2_2025_percent, score2_cutpoints, risk_levels)]

  old_num <- match(as.character(dt$score_2019_category), risk_levels)
  new_num <- match(as.character(dt$score2_2025_category), risk_levels)
  dt[, score_reclassification_direction := factor(
    fcase(
      is.na(old_num) | is.na(new_num), NA_character_,
      new_num < old_num, "Lower category under SCORE2",
      new_num == old_num, "Same category",
      new_num > old_num, "Higher category under SCORE2"
    ),
    levels = c("Lower category under SCORE2", "Same category", "Higher category under SCORE2"),
    ordered = TRUE
  )]
  dt[, score_reclassified := as.integer(score_2019_category != score2_2025_category)]
  dt[]
}

# ---- stable_u01_from_id: retained shared primitive ----
stable_u01_from_id <- function(id, seed = 2026L) {
  id_chr <- enc2utf8(as.character(id))
  vapply(id_chr, function(one_id) {
    if (is.na(one_id)) return(NA_real_)

    # UKB participant_id 为不超过 2^31-1 的整数。对这类 ID 使用可逆线性同余映射，
    # 因而在 ID 唯一时不会因哈希截断而产生重复随机数。
    id_num <- suppressWarnings(as.numeric(one_id))
    modulus <- 2147483647
    if (!is.na(id_num) && id_num >= 0 && id_num < modulus && id_num == floor(id_num)) {
      h <- (id_num * 48271 + as.double(seed)) %% modulus
      return((h + 0.5) / modulus)
    }

    # 非整数 ID 的后备路径。
    ints <- utf8ToInt(one_id)
    h <- as.double(seed %% modulus)
    if (length(ints) > 0L) {
      for (x in ints) h <- (h * 131 + x + 17) %% modulus
    }
    (h + 0.5) / modulus
  }, numeric(1))
}

# ---- add_stable_random_quantile: retained shared primitive ----
add_stable_random_quantile <- function(dt, id_col = "participant_id", seed = 2026L) {
  dt <- copy(as.data.table(dt))
  assert_unique_id(dt, id_col, "生成 z_common 的数据")
  dt[, z_common := stable_u01_from_id(get(id_col), seed = seed)]
  if (anyNA(dt$z_common) || any(dt$z_common <= 0 | dt$z_common >= 1)) {
    stop("z_common 生成失败或超出 (0,1)。", call. = FALSE)
  }
  dt[]
}

# ---- strategy07_ldl_absolute_target: retained shared primitive ----
strategy07_ldl_absolute_target <- function(risk_level,
  targets = c(L = 3.0, M = 2.6, H = 1.8, VH = 1.4)) {
  if (!is.numeric(targets) || !identical(names(targets),c("L","M","H","VH")) ||
      any(!is.finite(targets)) || any(targets<=0) || any(diff(targets)>=0))
    stop("LDL-C 风险等级目标必须为 L/M/H/VH 的递减正数")
  unname(targets[as.character(risk_level)])
}

# ---- strategy07_apply_ldl_treat_to_target: retained shared primitive ----
strategy07_apply_ldl_treat_to_target <- function(ldl,
                                                 risk_level,
                                                 z_common,
                                                 delta = 0.4,
                                                 ldl_min = 0.6,
                                                 ldl_targets = c(L=3.0,M=2.6,H=1.8,VH=1.4)) {
  n <- length(ldl)
  if (length(risk_level) != n || length(z_common) != n) {
    stop("LDL-C、风险等级和 z_common 长度不一致。", call. = FALSE)
  }
  if (any(!is.finite(ldl)) || any(ldl <= 0)) {
    stop("Branch A 的 LDL-C 必须为正的有限观测值。", call. = FALSE)
  }
  if (any(!is.finite(z_common)) || any(z_common < 0 | z_common > 1)) {
    stop("z_common 必须位于 [0,1]。", call. = FALSE)
  }
  if (length(delta) != 1L || !is.finite(delta) || delta <= 0) {
    stop("delta 必须为正的有限标量。", call. = FALSE)
  }
  if (length(ldl_min) != 1L || !is.finite(ldl_min) || ldl_min <= 0) {
    stop("ldl_min 必须为正的有限标量。", call. = FALSE)
  }

  absolute_target <- strategy07_ldl_absolute_target(risk_level,ldl_targets)
  if (anyNA(absolute_target)) stop("风险等级必须为 L/M/H/VH。", call. = FALSE)
  high_or_very_high <- as.character(risk_level) %in% c("H", "VH")

  # 与既有 07/08 冻结规则完全一致：
  # H/VH 的 >=50% 降幅要求只约束干预 upper bound，不参与“是否触发”的判断。
  upper <- ifelse(
    high_or_very_high,
    pmin(absolute_target, 0.5 * ldl),
    absolute_target
  )
  lower <- pmax(ldl_min, upper - delta)
  trigger <- ldl > absolute_target
  ldl_new <- ifelse(trigger, lower + z_common * (upper - lower), ldl)

  tol <- 1e-12
  if (any(trigger & lower > upper + tol)) stop("LDL intervention lower > upper。", call. = FALSE)
  if (any(trigger & (ldl_new < lower - tol | ldl_new > upper + tol))) {
    stop("触发者的 counterfactual LDL 不在冻结窗口内。", call. = FALSE)
  }
  if (any(!trigger & abs(ldl_new - ldl) > tol)) {
    stop("未触发者的 LDL 被意外改变。", call. = FALSE)
  }
  if (any(ldl_new > ldl + tol)) stop("干预后 LDL 高于基线 LDL。", call. = FALSE)

  data.table(
    absolute_target = absolute_target,
    lower = lower,
    upper = upper,
    intervention_trigger = as.integer(trigger),
    ldl_intervened = ldl_new,
    absolute_target_achieved = as.integer(ldl_new <= absolute_target + tol),
    relative_50pct_achieved = as.integer(ldl_new <= 0.5 * ldl + tol)
  )
}


# ==== Refactored orchestration utilities ==== 
m07_log <- function(stage,...){
  line<-paste0(sprintf("[%s][%s][PID:%s] ",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),stage,Sys.getpid()),paste0(...,collapse=""))
  cat(line,"\n",sep="");flush.console()
  # PSOCK workers do not always forward stdout; persist real-time per-PID progress.
  root<-Sys.getenv("M01_PROJECT_ROOT","D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
  out<-Sys.getenv("M01_ANALYSIS_OUTPUT",file.path(root,"LJW_outcome","data_intervention","M01_ASCVD","analysis_results"))
  logdir<-file.path(out,"run_logs");dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
  cat(line,"\n",file=file.path(logdir,sprintf("%s_PID%s.log",stage,Sys.getpid())),append=TRUE,sep="")
}
m07_assert <- function(x,cols,label="dataset") {miss<-setdiff(cols,names(x));if(length(miss))stop(label," missing: ",paste(miss,collapse=", ")); invisible(TRUE)}
m07_unique <- function(x){if(anyNA(x$participant_id)||anyDuplicated(x$participant_id))stop("participant_id not unique or NA")}
m07_md5 <- function(path){if(!file.exists(path))stop("Missing ",path);unname(tools::md5sum(path))}
m07_atomic <- function(object,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);tmp<-paste0(path,".part_",Sys.getpid());on.exit(if(file.exists(tmp))unlink(tmp),add=TRUE);saveRDS(object,tmp,compress=FALSE);if(file.exists(path))unlink(path);if(!file.rename(tmp,path))stop("Cannot commit ",path)}
m07_csv <- function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);data.table::fwrite(x,path,bom=TRUE,na="NA")}
m07_env <- function() {
  root<-Sys.getenv("M01_PROJECT_ROOT","D:/01_博士生/01_小论文/01_Baseline_intervention/04_Code")
  code_dir<-Sys.getenv("M01_CODE_DIR",file.path(root,"LJW_code","Main_code","M01_ASCVD"))
  clean_dir<-Sys.getenv("M01_CLEAN_DIR",file.path(root,"LJW_outcome","data_clean","M01_ASCVD"))
  out_dir<-Sys.getenv("M01_ANALYSIS_OUTPUT",file.path(root,"LJW_outcome","data_intervention","M01_ASCVD","analysis_results"))
  list(root=root,code_dir=code_dir,clean_dir=clean_dir,out_dir=out_dir)
}
m07_saveconfig <- function(cfg){m07_atomic(cfg,file.path(cfg$out_dir,"study_config.rds"))}
m07_config <- function(){env<-m07_env();p<-file.path(env$out_dir,"study_config.rds");if(!file.exists(p))stop("Run 07 first: ",p);cfg<-readRDS(p); cfg}
m07_manifest <- function(cfg){paths<-file.path(cfg$clean_dir,"MI_M01","MI_completed",sprintf("MI_%02d_analysis_M01.rds",seq_len(cfg$m)));if(!all(file.exists(paths)))stop("Missing completed MI: ",paste(paths[!file.exists(paths)],collapse="; "));data.table(imputation_id=seq_len(cfg$m),path=normalizePath(paths,winslash="/"),md5=vapply(paths,m07_md5,character(1)))}
m07_read <- function(cfg,m){p<-file.path(m07_dir(cfg,"07"),sprintf("m%02d.rds",m));if(!file.exists(p))stop("Missing 07 MI ",m);readRDS(p)}
m07_scope_spec <- function(x){
  m07_assert(x,c("age","sex","ses_class","hypertension"))
  list(overall=rep(TRUE,nrow(x)),age_lt_60=x$age<60,age_ge_60=x$age>=60,
  female=as.character(x$sex)=="female",male=as.character(x$sex)=="male",
  high_ses=as.character(x$ses_class)=="high",low_ses=as.character(x$ses_class)=="low",
  no_hypertension=as.integer(x$hypertension)==0L,hypertension=as.integer(x$hypertension)==1L)
}
m07_subgroups <- function(){data.table(subgroup_var=c("age","sex","ses","hypertension"),reference=c("age_lt_60","female","high_ses","no_hypertension"),comparison=c("age_ge_60","male","low_ses","hypertension"))}
m07_rms <- function()c("S0","SH","ST","SHT")
m07_states <- function()c("natural",as.vector(outer(c("2019","2025"),m07_rms(),paste,sep="_")))
m07_transitions <- function()c("unchanged_untreated","newly_intervention","intervention_removed","treatment_intensified","unchanged_intensity","treatment_deintensified")
m07_one <- function(ldl,level,z,cfg){
  a<-strategy07_apply_ldl_treat_to_target(ldl,level,z,cfg$ldl_delta,cfg$ldl_min,cfg$ldl_targets)
  G<-a$absolute_target;A<-a$intervention_trigger;U<-a$upper
  D<-ifelse(A==1,ldl-U,0)
  data.table(G=G,A=A,U=U,B=a$lower,D=D,Lcf=a$ldl_intervened,D_achieved=ldl-a$ldl_intervened)
}
m07_derive_transition <- function(A19,A25,D19,D25,epsilon=1e-10){
  delta<-D25-D19
  v<-data.table::fcase(A19==0&A25==0,"unchanged_untreated",A19==0&A25==1,"newly_intervention",A19==1&A25==0,"intervention_removed",A19==1&A25==1&delta>epsilon,"treatment_intensified",A19==1&A25==1&abs(delta)<=epsilon,"unchanged_intensity",A19==1&A25==1&delta< -epsilon,"treatment_deintensified",default=NA_character_)
  if(anyNA(v))stop("Treatment-transition assignment not exhaustive");v
}
m07_rm_definitions <- function() {
  d<-Sys.getenv("M01_CODE_DIR", ".")
  paths<-c(file.path(d,"config","RM_domains_v1.8.csv"),file.path(dirname(d),"config","RM_domains_v1.8.csv"))
  p<-paths[file.exists(paths)][1L]
  if(is.na(p))stop("Missing RM_domains_v1.8.csv; copy the config directory with the 07–12 scripts")
  tab<-data.table::fread(p)
  if(nrow(tab)!=16L||anyDuplicated(tab$domain)||!all(c("domain","type","include_2019","include_2025")%in%names(tab)))stop("Invalid RM domain configuration")
  expected<-c(`2019_H`=11,`2019_T`=1,`2025_H`=12,`2025_T`=2)
  for(year in c("2019","2025"))for(type in c("H","T")) {
    target_type<-type
    n<-nrow(tab[get("type")==target_type & get(paste0("include_",year))==1L])
    if(n!=expected[[paste0(year,"_",type)]])stop("RM definition size mismatch ",year,"/",type)
  }
  tab
}
m07_rm_binary <- function(x,year,type){
  d<-m07_rm_definitions()
  target_type<-type
  cols<-paste0("rm_",d[get("type")==target_type & get(paste0("include_",year))==1L,domain])
  m07_assert(x,cols,paste0(year," ",type," RM source"))
  mat<-vapply(cols,function(v){val<-as.integer(as.character(x[[v]]));if(anyNA(val)||any(!val%in%0:1))stop("Completed MI RM must be 0/1: ",v);val},integer(nrow(x)))
  z<-as.integer(rowSums(mat)>0)
  composite<-paste0("rm_",tolower(type),"_",year)
  if(composite%in%names(x) && any(as.integer(as.character(x[[composite]]))!=z))
    stop("06 composite RM and frozen v1.8 definition disagree: ",composite)
  z
}

m07_strategy_one <- function(completed,master,cfg,m){
  x<-m07_standardize_rm(data.table::copy(data.table::as.data.table(completed)));m07_unique(x);x[,participant_id:=as.character(participant_id)]
  for(v in c("age","sex","systolic_blood_pressure","diastolic_blood_pressure","total_cholesterol","ldl_cholesterol","hdl_cholesterol","current_smoker")){
    if(!v%in%names(x) || !v%in%names(master))stop("Observed SCORE input not shared ",v)
    ix<-match(x$participant_id,master$participant_id)
    if(anyNA(ix) || any(as.character(x[[v]])!=as.character(master[[v]][ix])))stop("Observed SCORE input changed after MI: ",v)
  }
  ix<-match(x$participant_id,master$participant_id)
  new<-master[ix,.(score_2019_percent,score2_2025_percent,score_2019_category,score2_2025_category,z_common)]
  stopifnot(nrow(new)==nrow(x))
  overlap<-intersect(names(new),names(x));if(length(overlap))stop("Strategy columns already in MI: ",paste(overlap,collapse=","))
  x<-data.table::as.data.table(cbind(x,new))
  for(year in c("2019","2025")){
    base<-as.character(x[[if(year=="2019")"score_2019_category" else "score2_2025_category"]])
    H<-m07_rm_binary(x,year,"H");T<-m07_rm_binary(x,year,"T")
    for(k in m07_rms()){
      hit<-switch(k,S0=rep(FALSE,nrow(x)),SH=H==1,ST=T==1,SHT=H==1|T==1)
      upgrade<-base=="M"&hit;level<-ifelse(upgrade,"H",base)
      nm<-paste(year,k,sep="_")
      x[[paste0("risk_",nm)]]<-level;x[[paste0("upgrade_",nm)]]<-as.integer(upgrade)
      y<-m07_one(x$ldl_cholesterol,level,x$z_common,cfg)
      for(v in names(y))x[[paste0(v,"_",nm)]]<-y[[v]]
      if(any(x[[paste0("A_",nm)]]==0 & (x[[paste0("D_",nm)]]!=0 | x[[paste0("Lcf_",nm)]]!=x$ldl_cholesterol)))stop("Non-intervened LDL invariant failed")
    }
  }
  for(k in m07_rms()){
    a19<-paste0("2019_",k);a25<-paste0("2025_",k)
    x[[paste0("transition_",k)]]<-m07_derive_transition(x[[paste0("A_",a19)]],x[[paste0("A_",a25)]],x[[paste0("D_",a19)]],x[[paste0("D_",a25)]])
  }
  for(year in c("2019","2025")){
    q<-function(k)x[[paste0("A_",year,"_",k)]]
    if(any(q("SHT")<q("SH")|q("SHT")<q("ST")|q("SH")<q("S0")|q("ST")<q("S0")))stop("RM strategy monotonicity violated")
    for(k in c("SH","ST","SHT")) {
      newly<-sum(q(k)==1&q("S0")==0);if(newly!=sum(q(k))-sum(q("S0")))stop("Additional-intervention QA violated")
    }
  }
  x<-data.table::as.data.table(as.list(x));data.table::set(x,j="imputation_id",value=rep(as.integer(m),nrow(x)))
  x[]
}

# Content- and configuration-addressed stage cache; never reuse stale estimates.
m07_cached <- function(cfg,stage,m,dependencies,fun){
  out<-m07_dir(cfg,stage);target<-file.path(out,sprintf("mi%02d.rds",m));sig_path<-paste0(target,".signature.rds")
  main_functions<-c(`09`="Function_09_Main_analysis_M01.R",`10`="Function_10_RM_classification_M01.R",
    `11`="Function_11_RM_effects_M01.R",`12`="Function_12_Treatment_transitions_M01.R")
  # The stage may source earlier function modules: hash the complete dependency chain.
  fns<-c("Function_07_Guideline_strategies_M01.R","Function_08_Counterfactual_risk_M01.R",
         "Function_09_Main_analysis_M01.R","Function_10_RM_classification_M01.R",
         "Function_11_RM_effects_M01.R","Function_12_Treatment_transitions_M01.R")
  use_n<-match(stage,c("07","08","09","10","11","12"))
  deps<-c(dependencies,file.path(cfg$out_dir,"study_config.rds"),
          file.path(Sys.getenv("M01_CODE_DIR"),fns[seq_len(use_n)]))
  if(any(!file.exists(deps)))stop("Cannot construct cache key; missing: ",paste(deps[!file.exists(deps)],collapse="; "))
  sig<-data.table(path=normalizePath(deps,winslash="/"),md5=vapply(deps,m07_md5,character(1)))
  if(file.exists(target)&&file.exists(sig_path)&&identical(as.data.frame(readRDS(sig_path)),as.data.frame(sig))){
    m07_log(stage,"MI",m," verified checkpoint reuse");return(readRDS(target))
  }
  result<-fun(m)
  m07_atomic(result,target);m07_atomic(sig,sig_path)
  result
}

# Result-oriented paths; numbered stages are internal dispatch keys only.
m07_dir <- function(cfg,stage){
  names<-c(`07`='strategy_assignments',`08`='risk_prediction',`09`='strategy_risks_and_effects',`10`='risk_classification',`11`='risk_modifier_incremental_effects',`12`='treatment_transitions')
  p<-file.path(cfg$out_dir,names[[stage]]);dir.create(p,recursive=TRUE,showWarnings=FALSE);p
}
m07_hash_object <- function(x){p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=2);m07_md5(p)}
m07_scientific_key <- function(cfg) m07_hash_object(list(schema=3L,targets=cfg$ldl_targets,delta=cfg$ldl_delta,min=cfg$ldl_min,rm=m07_rm_definitions(),epsilon=1e-10))
m07_timed <- function(stage,label,expr){
  tm<-proc.time()[3];status<-'ERROR';m07_log(stage,'START ',label)
  on.exit({elapsed<-unname(proc.time()[3]-tm);m07_log(stage,status,' ',label,' | elapsed ',round(elapsed,2),' sec')
    e<-m07_env();p<-file.path(e$out_dir,'run_logs',paste0('timing_PID',Sys.getpid(),'.csv'));m07_csv_append(data.table(timestamp=format(Sys.time()),stage=stage,task=label,status=status,elapsed_seconds=elapsed,pid=Sys.getpid()),p)},add=TRUE)
  z<-force(expr);status<-'OK';z
}
m07_csv_append <- function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);fwrite(x,path,append=file.exists(path),col.names=!file.exists(path))}
m07_result <- function(cfg,board,name,data,layer='pooled'){
  p<-file.path(cfg$out_dir,board,layer,paste0(name,'.csv'));m07_csv(data,p)
}
m07_export <- function(cfg,board,name,d,pooled=NULL){
  m07_result(cfg,board,name,d,'per_imputation');if(!is.null(pooled))m07_result(cfg,board,name,pooled)
}
m07_worker <- function(indices,fun,workers=1L,stage='',context=''){
  workers<-max(1L,min(as.integer(workers),length(indices)))
  if(workers==1L)return(lapply(indices,function(i)m07_timed(stage,paste(context,'MI',i),fun(i))))
  # Workers receive only their own MI and a small configuration; parent polls status files.
  cfg0<-get('cfg',envir=.GlobalEnv);d<-cfg0$code_dir
  taskdir<-file.path(cfg0$out_dir,'run_logs',paste0('tasks_',stage,'_',Sys.getpid()));dir.create(taskdir,recursive=TRUE,showWarnings=FALSE)
  cl<-parallel::makePSOCKcluster(workers,setup_timeout=30);on.exit(parallel::stopCluster(cl),add=TRUE)
  paths<-setNames(paste0('Function_',c('07_Guideline_strategies','08_Counterfactual_risk','09_Main_analysis','10_RM_classification','11_RM_effects','12_Treatment_transitions'),'_M01.R'),sprintf('%02d',7:12))
  parallel::clusterCall(cl,function(d,f,lib,cfg0){.libPaths(lib);Sys.setenv(M01_CODE_DIR=d,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');source(file.path(d,f),local=FALSE);data.table::setDTthreads(1);assign('cfg',cfg0,.GlobalEnv)},d,paths[[stage]],.libPaths(),cfg0)
  # Strip global bindings from the dispatched wrapper; all stage tasks are defined in paired modules.
  task<-function(i,stage,taskdir){
    p<-file.path(taskdir,paste0(i,'.rds'));tryCatch({z<-m07_timed(stage,paste('MI',i,'analysis'),m07_stage_task(stage,i,cfg));saveRDS(list(ok=TRUE,value=z),p)},error=function(e)saveRDS(list(ok=FALSE,message=conditionMessage(e)),p));TRUE
  }
  environment(task)<-.GlobalEnv
  # Asynchronous sends let the parent report elapsed time during long Cox fits.
  next_i<-1L;done<-0L;out<-vector('list',length(indices));busy<-rep(FALSE,workers);which_i<-integer(workers)
  send<-function(w,j){parallel:::sendCall(cl[[w]],task,list(indices[j],stage,taskdir),tag=j);busy[w]<<-TRUE;which_i[w]<<-j}
  for(w in seq_len(workers)){send(w,next_i);next_i<-next_i+1L}
  started<-Sys.time();last<-Sys.time()
  while(done<length(indices)){
    for(w in which(busy)){
      j<-which_i[w];p<-file.path(taskdir,paste0(indices[j],'.rds'))
      if(file.exists(p)){
        z<-tryCatch(readRDS(p),error=function(e)NULL);if(is.null(z))next
        parallel:::recvResult(cl[[w]])
        if(!z$ok)stop('MI ',indices[j],': ',z$message)
        out[[j]]<-z$value;done<-done+1L;busy[w]<-FALSE;m07_log(stage,'completed ',done,'/',length(indices),' ',context)
        if(next_i<=length(indices)){send(w,next_i);next_i<-next_i+1L}
      }
    }
    if(as.numeric(difftime(Sys.time(),last,units='secs'))>=5){m07_log(stage,'running MI ',paste(indices[which_i[busy]],collapse=','),' | completed ',done,'/',length(indices),' | elapsed ',round(as.numeric(difftime(Sys.time(),started,units='secs'))),' sec');last<-Sys.time()}
    if(done<length(indices))Sys.sleep(.2)
  };out
}
m07_stage_task <- function(stage,m,cfg){
  switch(stage,`08`=m08_predict_mi(m,cfg),`09`=m09_task(m,cfg),`10`=m10_task(m,cfg),`11`=m11_task(m,cfg),`12`=m12_task(m,cfg),stop('Unsupported worker stage'))
}
.score2019_cache_env <- new.env(parent=emptyenv())



# Check input signatures even when only a downstream stage is requested.
m07_check_manifest <- function(cfg, check_code=TRUE){
 saved<-readRDS(file.path(m07_dir(cfg,'07'),'input_manifest.rds'));now<-m07_manifest(cfg)
 if(!identical(as.data.frame(saved),as.data.frame(now)))stop('MI data changed: rerun strategy assignment')
 base<-readRDS(file.path(m07_dir(cfg,'07'),'observed_base_signature.rds'))
 if(!identical(base$md5,m07_md5(base$path)))stop('Observed score input changed: rerun strategy assignment')
 for(m in seq_len(cfg$m)){
  p<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m));sig<-paste(now[imputation_id==m,md5],m07_md5(file.path(cfg$code_dir,'Function_07_Guideline_strategies_M01.R')),m07_md5(file.path(m07_dir(cfg,'07'),'observed_score_master.rds')),m07_scientific_key(cfg),sep='|')
  if(!file.exists(paste0(p,'.sig')))stop('Missing strategy signature')
  stored<-readRDS(paste0(p,'.sig'))
  if(!check_code){
   # Explicit reporting-only reuse: still check MI, observed master and scientific parameters.
   a<-strsplit(stored,'|',fixed=TRUE)[[1]];b<-strsplit(sig,'|',fixed=TRUE)[[1]]
   if(length(a)!=length(b)||length(a)<4L||!identical(a[-2L],b[-2L]))stop('Frozen reporting inputs changed')
  }else if(!identical(stored,sig))stop('Stale strategy assignment: rerun strategy stage')
 };now
}

# Explicit compatibility with 06 source-domain names; never infer missing negatives.
m07_standardize_rm <- function(x){
 aliases<-c(rm_stress='rm_stressor_proxy',rm_family='rm_family_proxy')
 for(v in names(aliases))if(aliases[[v]]%in%names(x)){
  src<-as.integer(as.character(x[[aliases[[v]]]]))
  if(v%in%names(x)&&!identical(as.integer(as.character(x[[v]])),src))stop('Conflicting RM aliases: ',v)
  if(!v%in%names(x))set(x,j=v,value=src)
 };x
}


# Descriptive reporting only: reuse frozen pre-MI states and completed datasets.
m07_reporting_summaries <- function(cfg){
 templates<-file.path(cfg$code_dir,'config/report_static_templates')
 map<-c('Age, years'='age','Sex'='sex','Region'='region','Educational attainment'='education','Household income'='household_income','Employment status'='employment_status','Socioeconomic status'='ses_class','Social deprivation'='rm_deprivation','High-risk ethnicity'='ethnicity','Current smoking'='current_smoker','Physical inactivity'='rm_inactive','Alcohol consumption ≤10 g/day'='alcohol_score','Alcohol consumption'='alcohol_score','Healthy diet'='diet_score','Obesity and central obesity'='rm_obesity','Systolic blood pressure, mmHg'='systolic_blood_pressure','Systolic blood pressure'='systolic_blood_pressure','Diastolic blood pressure, mmHg'='diastolic_blood_pressure','Diastolic blood pressure'='diastolic_blood_pressure','Antihypertensive medication'='antihypertensive_medication_baseline','LDL cholesterol, mmol/L'='ldl_cholesterol','LDL-C'='ldl_cholesterol','HDL cholesterol, mmol/L'='hdl_cholesterol','HDL-C'='hdl_cholesterol','Triglycerides, mmol/L'='triglycerides','Triglycerides'='triglycerides','eGFR, mL/min/1.73 m²'='egfr_2021','eGFR (from creatinine)'='egfr_2021','Family history of CVD'='rm_family_proxy','Chronic liver disease'='chronic_liver_disease','Chronic respiratory disease'='chronic_respiratory_disease','Chronic inflammatory disorders'='rm_inflammatory','Major psychiatric disorders'='rm_psychiatric','Atrial fibrillation'='rm_af','Psychosocial stress'='rm_stressor_proxy','HIV-related condition'='rm_hiv','Obstructive sleep apnoea syndrome'='rm_sleep_apnoea','Obstructive sleep apnoea'='rm_sleep_apnoea','Non-alcoholic fatty liver disease'='rm_fatty_liver','NAFLD'='rm_fatty_liver','History of premature menopause'='rm_premature_menopause','Premature menopause'='rm_premature_menopause','Hypertensive disorders of pregnancy'='rm_hdp','Elevated lipoprotein(a)'='rm_lpa','Elevated Lp(a)'='rm_lpa','Elevated hs-CRP'='rm_crp','Townsend deprivation index'='townsend_deprivation_index','BMI, kg/m²'='bmi','BMI'='bmi','Waist circumference, cm'='waist_circumference','Waist circumference'='waist_circumference','hs-CRP, mg/L'='rm_source_crp','hs-CRP concentration'='rm_source_crp','Total cholesterol'='total_cholesterol','Physical inactivity, positive'='rm_inactive','High-risk ethnicity, positive'='ethnicity','Elevated Lp(a), positive'='rm_lpa')
 levels<-list(sex=c('female','male'),region=c('England','Wales','Scotland'),education=c('low','medium','high'),household_income=c('less_than_18000','18000_to_30999','31000_to_51999','52000_to_100000','greater_than_100000'),employment_status=c('employed_or_self_employed','retired','other_not_employed'),ses_class=c('high','low'))
 specs<-list()
 for(tab in c(3,5,6)){
  t<-fread(file.path(templates,sprintf('Supplementary_Table_%02d_reference.csv',tab)),header=FALSE,fill=TRUE,sep=",",blank.lines.skip=FALSE)
  v<-NULL;idx<-0L
  for(i in 2:nrow(t)){
   label<-t[[1]][i];if(is.na(label)||grepl('^Panel',label))next
   if(label%in%names(map)){v<-unname(map[label]);idx<-0L
    continuous<-v%in%c('age','systolic_blood_pressure','diastolic_blood_pressure','ldl_cholesterol','hdl_cholesterol','total_cholesterol','triglycerides','egfr_2021','townsend_deprivation_index','bmi','waist_circumference','rm_source_crp')
    kind<-if(tab==5)'missing' else if(continuous)'continuous' else if(tab==6&&grepl('positive|^Social deprivation$|^Obesity|^Elevated hs',label))'category' else 'missing'
    level<-if(kind=='category')'1'else NA_character_
   }else{
    if(is.null(v))stop('No variable for template row ',label)
    idx<-idx+1L;kind<-'category';level<-if(v%in%names(levels))levels[[v]][idx] else as.character(idx-1L)
   }
   specs[[length(specs)+1L]]<-data.table(table=tab,row=i,label=label,variable=v,kind=kind,level=level)
  }
 }
 spec<-rbindlist(specs);base<-as.data.table(readRDS(file.path(cfg$clean_dir,'data_analysis_base_incomplete_M01.rds')))
 prepare<-function(d,pre=FALSE){
  if(!'current_smoker'%in%names(d))d[,current_smoker:=as.integer(smoking_status=='current')]
  for(v in unique(spec$variable))if(pre&&paste0(v,'_observed')%in%names(d))set(d,j=v,value=d[[paste0(v,'_observed')]])
  if(pre&&! 'ses_class'%in%names(d))d[,ses_class:=NA_character_]
  miss<-setdiff(unique(spec$variable),names(d));if(length(miss))stop('Reporting variables missing: ',paste(miss,collapse=','));d
 }
 base<-prepare(base,TRUE);rows<-list()
 summarize<-function(d,m,source){
  out<-vector('list',nrow(spec))
  for(i in seq_len(nrow(spec))){r<-spec[i];v<-r$variable;x<-d[[v]];den<-length(x);missing<-sum(is.na(x));mask<-rep(TRUE,den)
   if(source=='imputed_only')mask<-is.na(base[[v]])
   if(r$table==5&&v%in%c('rm_hdp','rm_premature_menopause')){mask<-mask&as.character(d$sex)=='female'}
   x<-x[mask];N<-length(x);ok<-x[!is.na(x)];q<-rep(NA_real_,3)
   if(r$kind=='continuous'&&length(ok))q<-as.numeric(quantile(as.numeric(as.character(ok)),c(.25,.5,.75)))
   count<-if(r$kind=='category')sum(as.character(ok)==r$level)else NA_real_
   out[[i]]<-cbind(r,data.table(imputation_id=m,source=source,n_total=N,n_observed=length(ok),n_missing=if(r$table==5&&v=="rm_premature_menopause")sum(d$rm_menopause_uncertain[mask]==1L,na.rm=TRUE)else sum(is.na(x)),count=count,proportion=if(length(ok))count/length(ok)else NA_real_,q1=q[1],median=q[2],q3=q[3]))
  };rbindlist(out)
 }
 rows[[1]]<-summarize(base,0,'observed');manifest<-m07_manifest(cfg)
 for(m in seq_len(cfg$m)){
  d<-prepare(as.data.table(readRDS(manifest[imputation_id==m,path])))
  if(!identical(as.character(d$participant_id),as.character(base$participant_id)))stop('Reporting MI ID order mismatch')
  rows[[length(rows)+1L]]<-m07_timed('reporting',paste('MI',m,'baseline summaries'),summarize(d,m,'completed'))
  rows[[length(rows)+1L]]<-summarize(d,m,'imputed_only');rm(d);gc(FALSE)
 }
 all<-rbindlist(rows);m07_result(cfg,'baseline_and_imputation','baseline_and_imputation_descriptive_summary',all,'per_imputation')
 m07_result(cfg,'baseline_and_imputation','reporting_variable_manifest',spec)
 method<-file.path(cfg$clean_dir,'MI_M01/quality_control/imputation_method_M01.csv');if(file.exists(method))m07_result(cfg,'baseline_and_imputation','actual_imputation_methods',fread(method))
 invisible(all)
}




# Publication summaries and exports (same-stage functions).
# Reporting-only additions within the single stage-07 function script.
# Explicit reporting-only mode reuses frozen assignments after input checks.
m07_main_table1 <- function(cfg) {
  invisible(m07_check_manifest(cfg,check_code=FALSE))
  variables<-c('age','sex','rm_south_asian','ses_class','rm_deprivation',
    'current_smoker','rm_inactive','rm_obesity','systolic_blood_pressure',
    'hypertension','antihypertensive_medication_baseline','egfr_2021',
    'total_cholesterol','ldl_cholesterol','hdl_cholesterol','triglycerides')
  positive<-c(sex='female',ses_class='low')
  continuous<-c(1L,9L,12:16)
  result<-inputs<-list();identity<-event<-NULL
  for(m in seq_len(cfg$m)) {
    path<-file.path(m07_dir(cfg,'07'),sprintf('m%02d.rds',m))
    d<-as.data.table(readRDS(path));m07_unique(d)
    m07_assert(d,c(variables,'ascvd_event','event_status'))
    if(anyNA(d$ascvd_event)||any(!d$ascvd_event%in%0:1)||
       any(d$ascvd_event!=as.integer(d$event_status==1L)))stop('Table 1 event contract failed')
    if(m==1L){identity<-as.character(d$participant_id);event<-d$ascvd_event}
    if(!identical(identity,as.character(d$participant_id))||!identical(event,d$ascvd_event))
      stop('Table 1 ID/event status changed across MI')
    groups<-list(overall=rep(TRUE,nrow(d)),no_observed_ascvd=d$ascvd_event==0L,
      observed_ascvd=d$ascvd_event==1L)
    for(g in names(groups))for(i in seq_along(variables)){
      v<-variables[i];x<-d[[v]][groups[[g]]];n<-length(x)
      if(!n||anyNA(x))stop('Table 1 missing completed values or empty column: ',v,' ',g)
      iscont<-i%in%continuous
      q<-if(iscont)as.numeric(quantile(as.numeric(as.character(x)),c(.25,.5,.75)))else rep(NA_real_,3)
      count<-if(iscont)NA_real_ else sum(as.character(x)==if(v%in%names(positive))positive[[v]]else '1')
      result[[length(result)+1L]]<-data.table(imputation_id=m,row=i,scope=g,variable=v,
        kind=if(iscont)'continuous'else 'category',n_total=n,count=count,
        proportion=count/n,q1=q[1],median=q[2],q3=q[3])
    }
    inputs[[m]]<-data.table(imputation_id=m,path=normalizePath(path,winslash='/'),md5=m07_md5(path))
    rm(d);gc(FALSE)
  }
  a<-rbindlist(result)
  pooled<-a[,lapply(.SD,mean),by=.(row,scope,variable,kind),
    .SDcols=c('n_total','count','proportion','q1','median','q3')]
  pooled[,m:=cfg$m]
  m07_export(cfg,'baseline_and_imputation','main_table1_by_observed_ASCVD',a,pooled)
  m07_result(cfg,'baseline_and_imputation','main_table1_input_manifest',rbindlist(inputs))
  invisible(pooled)
}
