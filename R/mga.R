#' Bootstrap PLS multigroup analysis (PLS-MGA)
#'
#' A thin wrapper around \code{\link[seminr]{estimate_pls_mga}} that compares
#' structural path coefficients across two groups and applies the two-tailed
#' PLS-MGA decision rule (significant when \eqn{p < \alpha} or
#' \eqn{p > 1 - \alpha}). Estimation is performed entirely by \pkg{seminr}; this
#' function selects the groups, standardises the output, and optionally attaches
#' a MICOM caution when invariance was not established.
#'
#' @param model A fitted \pkg{seminr} model.
#' @param data A data frame with the model indicators and the grouping variable.
#' @param group_var Name of the grouping column in \code{data}.
#' @param group_a,group_b The two group labels to compare.
#' @param nboot Number of bootstrap resamples passed to \pkg{seminr}.
#' @param alpha Significance level.
#' @param micom_result Optional \code{metis_micom} result; when its invariance is
#'   not established a caution is added to the output.
#' @return An object of class \code{metis_mga} with a \code{paths} data frame.
#' @seealso \code{\link{metis_micom}}, \code{\link{metis_perm_mga}}
#' @export
metis_pls_mga <- function(
  model,
  data,
  group_var,
  group_a = NULL,
  group_b = NULL,
  nboot = 5000,
  alpha = 0.05,
  micom_result = NULL
) {
  .metis_mga_require_seminr()
  .metis_mga_require_model(model)
  nboot <- .metis_mga_validate_positive_integer(nboot, "nboot")
  alpha <- .metis_mga_validate_alpha(alpha)

  selection <- .metis_mga_select_two_groups(data, group_var, group_a, group_b)
  selected_model <- seminr::rerun(model, data = selection$data)
  cores <- .metis_mga_cores()

  mga <- seminr::estimate_pls_mga(
    pls_model = selected_model,
    condition = selection$condition,
    nboot = nboot,
    cores = cores
  )
  mga_df <- as.data.frame(mga, stringsAsFactors = FALSE)

  source <- if ("source" %in% names(mga_df)) mga_df$source else rep(NA_character_, nrow(mga_df))
  target <- if ("target" %in% names(mga_df)) mga_df$target else rep(NA_character_, nrow(mga_df))
  path <- ifelse(!is.na(source) & !is.na(target), paste(source, "->", target), rownames(mga_df))
  p_value <- .metis_mga_column(mga_df, "pls_mga_p")

  paths <- data.frame(
    path = path,
    group_a_beta = .metis_mga_column(mga_df, "group1_beta"),
    group_b_beta = .metis_mga_column(mga_df, "group2_beta"),
    difference = .metis_mga_column(mga_df, "diff"),
    p_value = p_value,
    decision = .metis_mga_decision(p_value, alpha),
    stringsAsFactors = FALSE
  )

  warning_message <- .metis_mga_micom_warning(micom_result)
  note <- "For PLS-MGA, p < alpha or p > 1 - alpha indicates significant group difference."
  if (nzchar(warning_message)) note <- paste(note, warning_message)

  out <- list(
    method = "PLS-MGA",
    engine = "seminr::estimate_pls_mga",
    groups = list(group_var = selection$group_var, group_a = selection$group_a, group_b = selection$group_b),
    settings = list(nboot = nboot, alpha = alpha, cores = cores),
    paths = paths,
    note = note
  )
  if (nzchar(warning_message)) out$warning <- warning_message
  class(out) <- c("metis_mga", class(out))
  out
}

#' Permutation multigroup analysis
#'
#' Optional permutation test of structural path differences across two groups.
#' Group labels are reshuffled, the model is re-estimated for each permutation,
#' and observed path differences are compared against the permutation
#' distribution. This is distinct from bootstrap PLS-MGA
#' (\code{\link{metis_pls_mga}}).
#'
#' @param model_spec A fitted \pkg{seminr} model to be re-estimated per permutation.
#' @param data A data frame with the model indicators and the grouping variable.
#' @param group_var Name of the grouping column in \code{data}.
#' @param group_a,group_b The two group labels to compare.
#' @param permutations Number of permutations.
#' @param alpha Significance level.
#' @param seed Random seed for reproducibility.
#' @return An object of class \code{metis_mga} with a \code{paths} data frame.
#' @seealso \code{\link{metis_pls_mga}}
#' @export
metis_perm_mga <- function(
  model_spec,
  data,
  group_var,
  group_a = NULL,
  group_b = NULL,
  permutations = 5000,
  alpha = 0.05,
  seed = 123
) {
  .metis_mga_require_seminr()
  .metis_mga_require_model(model_spec)
  permutations <- .metis_mga_validate_positive_integer(permutations, "permutations")
  alpha <- .metis_mga_validate_alpha(alpha)
  seed <- .metis_mga_validate_seed(seed)

  selection <- .metis_mga_select_two_groups(data, group_var, group_a, group_b)
  selected_model <- seminr::rerun(model_spec, data = selection$data)
  path_index <- .metis_mga_path_index(selected_model)

  group_model_a <- seminr::rerun(selected_model, data = selection$data[selection$condition, , drop = FALSE])
  group_model_b <- seminr::rerun(selected_model, data = selection$data[!selection$condition, , drop = FALSE])
  observed <- .metis_mga_path_diffs(group_model_a, group_model_b, path_index)

  set.seed(seed)
  permutation_values <- matrix(NA_real_, nrow = permutations, ncol = length(observed))
  colnames(permutation_values) <- names(observed)

  for (i in seq_len(permutations)) {
    perm_condition <- .metis_mga_permute_condition(length(selection$condition), selection$n_a)
    permutation_values[i, ] <- tryCatch({
      perm_model_a <- seminr::rerun(selected_model, data = selection$data[perm_condition, , drop = FALSE])
      perm_model_b <- seminr::rerun(selected_model, data = selection$data[!perm_condition, , drop = FALSE])
      .metis_mga_path_diffs(perm_model_a, perm_model_b, path_index)
    }, error = function(err) {
      rep(NA_real_, length(observed))
    })
  }

  paths <- lapply(names(observed), function(path) {
    permuted <- permutation_values[, path]
    ci <- .metis_mga_ci(permuted, c(alpha / 2, 1 - alpha / 2))
    p_value <- .metis_mga_p_value_two_tailed(observed[[path]], permuted)
    parts <- strsplit(path, " -> ", fixed = TRUE)[[1]]
    beta_a <- group_model_a$path_coef[parts[[1]], parts[[2]]]
    beta_b <- group_model_b$path_coef[parts[[1]], parts[[2]]]
    data.frame(
      path = path,
      group_a_beta = beta_a,
      group_b_beta = beta_b,
      difference = observed[[path]],
      ci_lower = ci[[1]],
      ci_upper = ci[[2]],
      p_value = p_value,
      decision = if (!is.na(p_value) && p_value < alpha) "significant" else "not significant",
      stringsAsFactors = FALSE
    )
  })

  out <- list(
    method = "Permutation MGA",
    engine = "seminr::rerun with random group-label reassignment",
    groups = list(group_var = selection$group_var, group_a = selection$group_a, group_b = selection$group_b),
    settings = list(permutations = permutations, alpha = alpha, seed = seed),
    paths = do.call(rbind, paths),
    note = "Permutation MGA re-estimates group models after each group-label reassignment. This is not PLS-MGA."
  )
  rownames(out$paths) <- NULL
  attr(out$paths, "permutation_values") <- as.data.frame(permutation_values, stringsAsFactors = FALSE)
  class(out) <- c("metis_mga", class(out))
  out
}

#' Print a multigroup analysis result
#'
#' @param x An object of class \code{metis_mga}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.metis_mga <- function(x, ...) {
  cat(x$method, "\n")
  cat("Engine:", x$engine, "\n")
  cat("Groups:", x$groups$group_var, "=", x$groups$group_a, "vs", x$groups$group_b, "\n")
  if (!is.null(x$settings$nboot)) {
    cat("Bootstraps:", x$settings$nboot, " Alpha:", x$settings$alpha, "\n\n")
  } else {
    cat("Permutations:", x$settings$permutations, " Alpha:", x$settings$alpha, " Seed:", x$settings$seed, "\n\n")
  }
  if (nrow(x$paths)) print(x$paths, row.names = FALSE)
  cat("\nNote\n")
  cat(x$note, "\n")
  invisible(x)
}

.metis_mga_require_seminr <- function() {
  if (!requireNamespace("seminr", quietly = TRUE)) {
    stop("Package 'seminr' is required for METIS MGA.")
  }
}

.metis_mga_require_model <- function(model) {
  if (!is.list(model)) stop("model must be a fitted seminr model object.")
  required <- c("rawdata", "path_coef", "smMatrix")
  missing <- required[vapply(required, function(name) is.null(model[[name]]), logical(1))]
  if (length(missing)) {
    stop(sprintf("model is missing required seminr fields: %s.", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

.metis_mga_validate_positive_integer <- function(value, name) {
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) stop(sprintf("%s must be a positive integer.", name))
  value
}

.metis_mga_validate_alpha <- function(alpha) {
  alpha <- suppressWarnings(as.numeric(alpha))
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be a number greater than 0 and less than 1.")
  alpha
}

.metis_mga_validate_seed <- function(seed) {
  seed <- suppressWarnings(as.integer(seed))
  if (is.na(seed)) stop("seed must be an integer.")
  seed
}

.metis_mga_cores <- function() {
  cores <- suppressWarnings(as.integer(Sys.getenv("METIS_MICOM_MGA_CORES", "1")))
  if (is.na(cores) || cores < 1L) 1L else cores
}

.metis_mga_select_two_groups <- function(data, group_var, group_a = NULL, group_b = NULL) {
  data <- as.data.frame(data)
  group_var <- as.character(group_var)
  if (length(group_var) != 1L || !nzchar(group_var)) stop("group_var must be a single column name.")
  if (!group_var %in% names(data)) stop(sprintf("group_var '%s' was not found in data.", group_var))

  labels <- as.character(data[[group_var]])
  available <- unique(labels[!is.na(labels) & nzchar(labels)])
  if (length(available) < 2L) stop("group_var must contain at least two non-missing groups.")

  if (is.null(group_a) && is.null(group_b)) {
    if (length(available) > 2L) stop("group_a and group_b are required when group_var has more than two groups.")
    group_a <- available[[1]]
    group_b <- available[[2]]
  } else if (is.null(group_a) || is.null(group_b)) {
    stop("group_a and group_b must be supplied together.")
  }

  group_a <- as.character(group_a)
  group_b <- as.character(group_b)
  if (identical(group_a, group_b)) stop("group_a and group_b must be different.")
  if (!group_a %in% available) stop(sprintf("group_a '%s' was not found in group_var.", group_a))
  if (!group_b %in% available) stop(sprintf("group_b '%s' was not found in group_var.", group_b))

  keep <- labels %in% c(group_a, group_b)
  selected <- data[keep, , drop = FALSE]
  selected_labels <- labels[keep]
  condition <- selected_labels == group_a
  n_a <- sum(condition)
  n_b <- sum(!condition)
  if (n_a < 2L || n_b < 2L) stop("Both selected groups must contain at least two observations.")

  list(data = selected, condition = condition, group_var = group_var, group_a = group_a, group_b = group_b, n_a = n_a, n_b = n_b)
}

.metis_mga_column <- function(data, name) {
  if (!name %in% names(data)) return(rep(NA_real_, nrow(data)))
  suppressWarnings(as.numeric(data[[name]]))
}

.metis_mga_decision <- function(p_value, alpha) {
  ifelse(is.na(p_value), "undetermined", ifelse(p_value < alpha | p_value > 1 - alpha, "significant", "not significant"))
}

.metis_mga_micom_warning <- function(micom_result) {
  if (is.null(micom_result)) return("")
  partial <- tryCatch(isTRUE(micom_result$invariance$partial), error = function(err) FALSE)
  if (partial) return("")
  "MICOM partial measurement invariance was not established; interpret group path comparisons cautiously."
}

.metis_mga_path_index <- function(model) {
  path_matrix <- model$smMatrix
  if (is.null(path_matrix) || is.null(rownames(path_matrix)) || is.null(colnames(path_matrix))) {
    path_matrix <- model$path_coef != 0
  }
  positions <- which(path_matrix != 0, arr.ind = TRUE)
  if (!nrow(positions)) stop("No structural paths were found in the seminr model.")
  data.frame(
    source = rownames(path_matrix)[positions[, "row"]],
    target = colnames(path_matrix)[positions[, "col"]],
    stringsAsFactors = FALSE
  )
}

.metis_mga_path_diffs <- function(model_a, model_b, path_index) {
  out <- vapply(seq_len(nrow(path_index)), function(i) {
    source <- path_index$source[[i]]
    target <- path_index$target[[i]]
    model_a$path_coef[source, target] - model_b$path_coef[source, target]
  }, numeric(1))
  stats::setNames(out, paste(path_index$source, path_index$target, sep = " -> "))
}

.metis_mga_permute_condition <- function(n, n_a) {
  condition <- rep(FALSE, n)
  condition[sample.int(n, size = n_a, replace = FALSE)] <- TRUE
  condition
}

.metis_mga_ci <- function(values, probs) {
  values <- values[!is.na(values)]
  if (!length(values)) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(values, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

.metis_mga_p_value_two_tailed <- function(observed, values) {
  values <- values[!is.na(values)]
  if (is.na(observed) || !length(values)) return(NA_real_)
  mean(abs(values) >= abs(observed))
}
