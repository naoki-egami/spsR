#' @export
sps_cv <- function(estimate, se, X_selected, site_name = NULL,
                   K = 2, max_iter = 100, seed){

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
    out_cv_i <- sps_estimator_non_selected(X = X_selected, ss = ss_use,
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

  out <- list("p_value" = final_p_value, "table" = cv_tab)

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
