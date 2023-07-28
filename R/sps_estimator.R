sps_estimator <- function(out, estimate, se, W = NULL){

  if(is.null(W) == TRUE){
    # first estimate weights
    out_W <- sps_weights(X = out$X, site = out$ss, site_name = NULL)

    W <- out_W$W
  }
  N <- out$N

  # Estimate between-site variance
  bet_var <- between_var_sps(estimate = estimate, se = se,
                             X_selected = out$X[out$ss == 1, , drop = FALSE])

  N_S <- length(estimate)
  N_R <- N - N_S

  # Each Site
  estimate_each <- t(W) %*% estimate
  se_each <- sqrt(t(W^2) %*% c(se^2) + bet_var)
  out_each <- cbind(estimate_each, se_each); colnames(out_each) <- c("Estimate", "Std. Error")

  # ###############
  # Overall
  # ###############
  std_w <- (1 + apply(W, 1, sum))/sum(1 + apply(W, 1, sum))
  estimate_overall <- sum(std_w * estimate)
  se_overall <- sqrt(sum(std_w^2 *(se^2 + bet_var)))
  out_overall <- c(estimate_overall, se_overall); names(out_overall) <- c("Estimate", "Std. Error")

  out <- list("overall" = out_overall, "each" = out_each, "bet_se" = sqrt(bet_var))

  return(out)
}
