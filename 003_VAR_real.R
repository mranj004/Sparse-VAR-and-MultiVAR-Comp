# Lasso VAR(1) for the Nestler daily measurements.
# Dependency: install.packages("glmnet")
#
# In an R Markdown chunk (source() defines functions without running a fit):
#   source("003_VAR_real.R")
#   var_result <- run_var_lasso()
#   tibble::as_tibble(var_result$metrics)
#   tibble::as_tibble(var_result$metrics_by_variable)
#   var_result$models[["9"]]$A
#   shared_result <- run_var_lasso(mode = "pooled")
#   tibble::as_tibble(rbind(var_result$metrics, shared_result$metrics))
#
# Defaults: fit on days 1:50, tune lambda on targets 51:66, refit on
# days 1:66, and evaluate once on targets 67:82. Every forecast uses the
# OBSERVED previous day, including within the test block (one-step forecasts).
# Coefficients stay fixed throughout each evaluation block.
#
# For each participant, solve exactly the entrywise lasso objective
#   ||Y - X %*% t(A)||_F^2 / (2*n) + lambda * sum(abs(A)), n = T - 1.
# Rows of A are tomorrow's responses; columns are today's predictors.
# This separates into six Gaussian lasso regressions with the SAME lambda.
# No intercept, centering, or scaling is added. In pooled mode, n is the total
# number of within-participant training pairs, and one A is shared by all IDs.
# See https://glmnet.stanford.edu/articles/glmnet.html for the Gaussian loss.
# The first five measurements have approximately zero means over all 82 days
# within each selected participant, suggesting full-series centering upstream.
# The split below prevents leakage in fitting/tuning, but cannot undo that
# preprocessing. Original uncentered data would be needed to rule it out.

.var_real_dir <- local({
  paths <- lapply(sys.frames(), function(frame) frame$ofile)
  paths <- Filter(function(path) is.character(path) && length(path) == 1L, paths)
  if (length(paths)) dirname(normalizePath(tail(paths, 1L)[[1L]])) else getwd()
})

# Fit the entire requested lambda path. Each array slice is a 6 x 6 A matrix.
.var_lasso_path <- function(X, Y, lambdas) {
  p <- ncol(Y)
  result <- array(0, c(p, p, length(lambdas)),
                  dimnames = list(colnames(Y), colnames(X), as.character(lambdas)))
  for (j in seq_len(p)) {
    # An identically zero response has the exact zero solution at every lambda.
    if (all(Y[, j] == 0)) next
    fit_args <- list(
      x = X, y = Y[, j], family = "gaussian", alpha = 1,
      lambda = lambdas, intercept = FALSE, standardize = FALSE
    )
    # glmnet >= 5 uses control; retain compatibility with older installations.
    controls <- list(thresh = 1e-12, maxit = 1000000L)
    if ("control" %in% names(formals(glmnet::glmnet))) {
      fit_args$control <- controls
    } else fit_args <- c(fit_args, controls)
    fit <- do.call(glmnet::glmnet, fit_args)
    if (fit$jerr != 0L || length(fit$lambda) != length(lambdas)) {
      stop("glmnet did not fit the complete lambda path.")
    }
    result[j, , ] <- as.matrix(stats::coef(fit))[-1L, , drop = FALSE]
  }
  result
}

# MSFE = average squared Euclidean error per six-variable forecast.
# MSE = average squared error per scalar measurement; MSFE = 6 * MSE.
.var_forecast_metrics <- function(predictions, by = character()) {
  groups <- c("mode", "model", by)
  key <- do.call(interaction, c(predictions[groups], list(drop = TRUE, lex.order = TRUE)))
  rows <- lapply(split(predictions, key), function(z) {
    error <- z$actual - z$predicted
    n_forecasts <- nrow(unique(z[c("id", "day")]))
    cbind(z[1L, groups, drop = FALSE], data.frame(
      n_forecasts = n_forecasts, n_values = nrow(z),
      MSFE = sum(error^2) / n_forecasts,
      MSE = mean(error^2), RMSE = sqrt(mean(error^2)), MAE = mean(abs(error))
    ))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

run_var_lasso <- function(
    data_path = file.path(.var_real_dir, "dataset", "flip", "0002_nestler_ts.tsv"),
    mode = c("individual", "pooled"),
    train_end = 50L, validation_end = 66L, total_days = 82L,
    lambda_grid = NULL, n_lambda = 60L, lambda_min_ratio = 1e-4,
    expected_n = 21L) {
  mode <- match.arg(mode)
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop('Package "glmnet" is required. Run install.packages("glmnet") first.')
  }
  is_integer <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x == floor(x)
  if (!all(vapply(list(train_end, validation_end, total_days), is_integer, logical(1))) ||
      train_end < 3L || train_end >= validation_end || validation_end >= total_days) {
    stop("Require integer days: 3 <= train_end < validation_end < total_days.")
  }
  if (!is.null(expected_n) && (!is_integer(expected_n) || expected_n < 1L)) {
    stop("expected_n must be a positive integer or NULL.")
  }
  if (is.null(lambda_grid)) {
    if (!is_integer(n_lambda) || n_lambda < 2L ||
        !is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1L ||
        !is.finite(lambda_min_ratio) || lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
      stop("Require n_lambda >= 2 and 0 < lambda_min_ratio < 1.")
    }
  } else if (!is.numeric(lambda_grid) || !length(lambda_grid) ||
             any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
    stop("lambda_grid must contain finite nonnegative values.")
  }

  variables <- c("sociable", "creative", "friendly", "organised", "self_esteem", "nervous")
  data <- utils::read.delim(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("id", "day", variables)
  if (!all(required %in% names(data))) {
    stop("Missing columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (anyNA(data$id) || anyNA(data$day) || !is.numeric(data$day) ||
      any(!is.finite(data$day)) || any(data$day != floor(data$day))) {
    stop("IDs must be present and days must be finite integers.")
  }
  if (!all(vapply(data[variables], is.numeric, logical(1)))) {
    stop("All six measurement columns must be numeric.")
  }
  data$id <- as.character(data$id)
  ids <- unique(data$id)
  audit <- do.call(rbind, lapply(ids, function(id) {
    z <- data[data$id == id & data$day %in% seq_len(total_days), required]
    data.frame(
      id = id, n_rows = nrow(z),
      complete_days = nrow(z) == total_days &&
        !anyDuplicated(z$day) && setequal(z$day, seq_len(total_days)),
      finite_measurements = all(is.finite(as.matrix(z[variables])))
    )
  }))
  audit$selected <- audit$complete_days & audit$finite_measurements
  selected_ids <- audit$id[audit$selected]
  if (!length(selected_ids)) stop("No participants have complete daily records.")
  if (!is.null(expected_n) && length(selected_ids) != expected_n) {
    stop("Expected ", expected_n, " complete participants; found ", length(selected_ids),
         ". Set expected_n = NULL to allow a different count.")
  }
  selected <- data[data$id %in% selected_ids & data$day %in% seq_len(total_days), required]
  selected <- selected[order(match(selected$id, selected_ids), selected$day), ]
  rownames(selected) <- NULL

  # Build pairs separately for every ID: never connect two participants.
  series <- lapply(selected_ids, function(id) {
    z <- selected[selected$id == id, ]
    values <- as.matrix(z[variables])
    list(X = values[-nrow(values), , drop = FALSE],
         Y = values[-1L, , drop = FALSE],
         index = data.frame(id = id, day = z$day[-1L]))
  })
  names(series) <- selected_ids
  X <- do.call(rbind, lapply(series, `[[`, "X"))
  Y <- do.call(rbind, lapply(series, `[[`, "Y"))
  index <- do.call(rbind, lapply(series, `[[`, "index"))
  rownames(index) <- NULL
  index$split <- ifelse(index$day <= train_end, "train",
                        ifelse(index$day <= validation_end, "validation", "test"))
  groups <- if (mode == "individual") {
    setNames(lapply(selected_ids, function(id) which(index$id == id)), selected_ids)
  } else list(pooled = seq_len(nrow(index)))

  models <- tuning <- forecasts <- list()
  for (group in names(groups)) {
    rows <- groups[[group]]
    training <- rows[index$split[rows] == "train"]
    validation <- rows[index$split[rows] == "validation"]
    development <- rows[index$split[rows] != "test"]
    test <- rows[index$split[rows] == "test"]
    if (is.null(lambda_grid)) {
      # At lambda_max, every coefficient is zero. Use training data ONLY.
      lambda_max <- max(abs(crossprod(X[training, , drop = FALSE],
                                     Y[training, , drop = FALSE]))) / length(training)
      lambdas <- if (lambda_max == 0) 1 else
        exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = n_lambda))
    } else lambdas <- sort(unique(lambda_grid), decreasing = TRUE)

    path <- .var_lasso_path(X[training, , drop = FALSE], Y[training, , drop = FALSE], lambdas)
    validation_msfe <- vapply(seq_along(lambdas), function(k) {
      error <- Y[validation, , drop = FALSE] - X[validation, , drop = FALSE] %*% t(path[, , k])
      mean(rowSums(error^2))
    }, numeric(1))
    best <- which.min(validation_msfe) # descending grid: ties favor more sparsity
    chosen_lambda <- lambdas[best]
    tuning[[group]] <- data.frame(
      group = group, lambda = lambdas, validation_MSFE = validation_msfe,
      selected = seq_along(lambdas) == best
    )
    A <- .var_lasso_path(X[development, , drop = FALSE],
                        Y[development, , drop = FALSE], chosen_lambda)[, , 1L]
    residuals <- Y[development, , drop = FALSE] - X[development, , drop = FALSE] %*% t(A)
    models[[group]] <- list(
      A = A, lambda = chosen_lambda, n_training_pairs = length(development),
      validation_MSFE = validation_msfe[best], nonzero = sum(A != 0),
      spectral_radius = max(Mod(eigen(A, only.values = TRUE)$values)),
      objective = sum(residuals^2) / (2 * length(development)) + chosen_lambda * sum(abs(A))
    )

    # Two reference forecasts, using the same test rows as the VAR.
    # The mean baseline always uses each participant's own days 1:validation_end.
    means <- t(vapply(index$id[test], function(id) {
      colMeans(selected[selected$id == id & selected$day <= validation_end, variables])
    }, numeric(length(variables))))
    predictions <- list(lasso_VAR1 = X[test, , drop = FALSE] %*% t(A),
                        persistence = X[test, , drop = FALSE], training_mean = means)
    forecasts[[group]] <- do.call(rbind, lapply(names(predictions), function(model) {
      data.frame(
        mode = mode, model = model, id = rep(index$id[test], each = length(variables)),
        day = rep(index$day[test], each = length(variables)),
        predictor_day = rep(index$day[test] - 1L, each = length(variables)),
        variable = rep(variables, times = length(test)),
        actual = as.vector(t(Y[test, , drop = FALSE])),
        predicted = as.vector(t(predictions[[model]]))
      )
    }))
  }
  forecasts <- do.call(rbind, forecasts)
  rownames(forecasts) <- NULL
  forecasts$error <- forecasts$actual - forecasts$predicted
  forecasts$squared_error <- forecasts$error^2
  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  list(
    metrics = .var_forecast_metrics(forecasts),
    metrics_by_participant = .var_forecast_metrics(forecasts, "id"),
    metrics_by_variable = .var_forecast_metrics(forecasts, "variable"),
    predictions = forecasts, models = models, tuning = tuning,
    participant_audit = audit, selected_ids = selected_ids,
    data = selected, pair_index = index,
    settings = list(
      mode = mode, variables = variables, data_path = normalizePath(data_path),
      train_end = train_end, validation_end = validation_end, total_days = total_days,
      n_participants = length(selected_ids), n_participant_days = nrow(selected),
      n_test_forecasts = sum(index$split == "test"),
      intercept = FALSE, standardize = FALSE,
      forecast_type = "Fixed coefficients; one step ahead using observed previous day",
      primary_metric = "MSFE: mean squared Euclidean error per six-variable forecast",
      glmnet_version = as.character(utils::packageVersion("glmnet")),
      R_version = R.version.string
    )
  )
}
