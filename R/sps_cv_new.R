#' Site-level Cross-Validation for Synthetic Purposive Sampling
#' @param out Output from function \code{sps()}
#' @param estimates_selected data.frame with two columns: the first column represents estimates of the site-specific ATEs for the selected sites and the second column represents its corresponding standard error. The number of rows is equal to the number of the selected sites and \code{rownames(estimates_selected)} should be the names of the selected sites.
#' @param K (Default = \code{2}) Fold of the cross-validation.
#' @param max_iter (Default = \code{100}) How many times the function repeats K-fold cross-validation.
#' @param seed Numeric. \code{seed} used internally. Default = \code{1234}.
#' @return \code{sps_cv} returns two objects.
#'  \itemize{
#'    \item \code{p_value}: P-value for the null hypothesis that the the average-site ATE among non-selected sites is equal to the weighted average estimator of the average-site ATE based on site-specific ATE estimates in selected sites.
#'    \item \code{internal}: Objects useful for internal use of the function.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
sps_cv_new <- function(out = NULL, estimates_selected = NULL, K = 2, max_iter = 100, seed = 1234){

  ## Housekeeping
  if(("sps" %in% class(out)) == FALSE){
    stop(" `out` should be an output from function `sps()` ")
  }
  selected_sites <- out$selected_sites

  ## estimates_selected
  if(setequal(rownames(estimates_selected), selected_sites) == FALSE){
    stop(" `rownames(estimates_selected)` should match to `selected_sites` ")
  }

  ## X_selected
  X_selected <- out$internal$X[out$internal$ss == 1, , drop = FALSE]

  # reordering
  estimates_selected <- estimates_selected[match(rownames(X_selected), rownames(estimates_selected)), ]

  estimate <- estimates_selected[, 1]
  se <- estimates_selected[, 2]



  s <- nrow(X_selected)
  ind <- seq(1:s)
  half <- ceiling(s/K)

  n_check <- min(choose(s, half), max_iter)
  all_comb <- combn(s, half)

  if(choose(s, half) >= max_iter){
    set.seed(seed)
    use_index <- sample(x = seq(1:ncol(all_comb)), size = max_iter, replace = FALSE)
  }else{
    use_index <- seq(1:ncol(all_comb))
  }

  cv_tab <- matrix(NA, nrow = n_check, ncol = 5)
  for(i in 1:n_check){
    use_non_selected_index  <- all_comb[, use_index[i]]
    s_use <- setdiff(ind, use_non_selected_index)
    ss_use <- rep(0, s)
    ss_use[s_use] <- 1

    # estimate
    out_cv_i <- sps_estimator_non_selected(X = X_selected,
                                           ss = ss_use,
                                           estimate = estimate[s_use], se = se[s_use])

    # truth
    # estimate_cv <- mean(estimate[s_use])
    est_cv <- rma(yi = estimate[use_non_selected_index],
                  vi = se[use_non_selected_index]^2, method = "REML")

    cv_tab[i, 1:5] <- c(est_cv$beta, est_cv$se,
                        out_cv_i$non_selected[1], out_cv_i$non_selected[2],
                        out_cv_i$mean_RMSE)
  }
  colnames(cv_tab) <- c("Obs:Estimate", "Obs:Std. Error",
                        "Pred:Estimate", "Pred:Std. Error",
                        "RMSE")

  diff    <- cv_tab[,1] - cv_tab[,3]
  diff_se <- sqrt(cv_tab[,2]^2 + cv_tab[,4]^2)
  p_value <- 2*(1 - pnorm(abs(diff/diff_se)))

  # Holms correction
  sorted_p_value <- sort(p_value)
  m <- length(p_value)
  holm_p_value <- (m + 1 - seq(from = 1, to = m))*sorted_p_value
  final_p_value <- min(holm_p_value)

  out <- list("p_value" = final_p_value, "internal" = cv_tab)

  return(out)
}

sps_estimator_non_selected <- function(X, ss, estimate, se){

  out_W <- sps_weights(X = X, site = ss, site_name = NULL)
  W <- out_W$W
  mean_RMSE <- mean(out_W$RMSE)

  # Estimate between-site variance
  bet_var <- between_var_sps(estimate = estimate, se = se, X_selected = X[ss == 1, , drop = FALSE])

  # ###############
  # Estimate for Average of non-selected Sites
  # ###############
  std_w <- apply(W, 1, mean)
  estimate_non <- sum(std_w * estimate)
  se_non <- sqrt(sum(std_w^2 *(se^2 + bet_var)))
  out_non <- c(estimate_non, se_non); names(out_non) <- c("Estimate", "Std. Error")

  out <- list("non_selected" = out_non,
              "mean_RMSE" = mean_RMSE,
              "bet_se" = sqrt(bet_var))

  return(out)
}
