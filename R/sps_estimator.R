#' @export
sps_estimator <- function(out = NULL, estimates_selected = NULL, X_selected = NULL, X = NULL, ss = NULL, W = NULL){

  if(is.null(out) == FALSE){
    sps_used <- TRUE
    data_name <- out$selected_sites

    if(all(rownames(estimates_selected) %in% data_name) == FALSE){
      stop(" `rownames(estimates_selected)` should match to names of the selected sites (`out$selected`) ")
    }

    X_use  <- out$X
    ss_use <- out$ss
    N <- out$N
    X_selected <- X_use[ss_use == 1, , drop = FALSE]

  }else if(is.null(X_selected) == FALSE){
    X_selected <- X_selected
    data_name <- rownames(X_selected)

    if(all(rownames(estimates_selected) %in% data_name) == FALSE){
      stop(" `rownames(estimates_selected)` should match to `rownames(X_selected)` ")
    }

    X_use  <- X
    ss_use <- ss
    N <- nrow(X_use)

  }else{
    stop(" `out` or `X_selected` should be specified.  ")
  }

  # reordering
  estimates_selected <- estimates_selected[match(data_name, rownames(estimates_selected)),]

  estimate <- estimates_selected[, 1]
  se <- estimates_selected[, 2]

  if(is.null(W) == TRUE){
    # first estimate weights
    out_W <- sps_weights(X = X_use, site = ss_use, site_name = NULL)
    W <- out_W$W
  }

  # Estimate between-site variance
  bet_var <- between_var_sps(estimate = estimate, se = se, X_selected = X_selected)

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

  class(out) <- c(class(out), "sps_estimator")

  return(out)
}
