#' Synthetic Purposive Sampling Estimator for the Average-Site ATE
#' @param out Output from function \code{sps()}
#' @param estimates_selected data.frame with two columns: the first column represents estimates of the site-specific ATEs for the selected sites and the second column represents its corresponding standard error. The number of rows is equal to the number of the selected sites and \code{rownames(estimates_selected)} should be names of the selected sites.
#' @param X (Optional. Use this only when sites are selected based on \code{sps()}. Default = \code{NULL}) Site-level variables for the target population of sites. Row names should be names of sites.
#' @param X_selected (Optional. Use this only when sites are selected based on \code{sps()}. Default = \code{NULL}) Site-level variables for the selected sites. Row names should be names of sites.
#' @param site_selected (Optional. Use this only when sites are selected based on \code{sps()}. Default = \code{NULL}) Names of sites users selected.
#' @param W (Optional. Default = \code{NULL}) Estimated weights to be used for estimation. When \code{NULL}, the function will estimate weights given the selected sites.
#' @return \code{sps_estimator} returns an object of \code{sps_estimator} class.
#'  \itemize{
#'    \item \code{average_site_ATE}: An estimate of the average-site ATE and its corresponding standard error.
#'    \item \code{site_specific_ATE}: Estimates and standard errors for the site-specific ATEs in non-selected sites.
#'    \item \code{bet_se}: Estimated between-site standard errors.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
sps_estimator <- function(out = NULL, estimates_selected = NULL, X = NULL, X_selected = NULL, site_selected = NULL, W = NULL){

  if(is.null(out) == FALSE){
    sps_used <- TRUE
    data_name <- out$selected_sites

    if(all(rownames(estimates_selected) %in% data_name) == FALSE){
      stop(" `rownames(estimates_selected)` should match to names of the selected sites (`out$selected`) ")
    }

    X_use  <- out$internal$X
    ss_use <- out$internal$ss
    N <- out$internal$N
    X_selected <- X_use[ss_use == 1, , drop = FALSE]

  }else if(is.null(X_selected) == FALSE){
    X_selected <- X_selected
    data_name <- rownames(X_selected)

    if(all(rownames(estimates_selected) %in% data_name) == FALSE){
      stop(" `rownames(estimates_selected)` should match to `rownames(X_selected)` ")
    }
    if(all(is.null(site_selected)) == TRUE){
      stop(" `When users do not supply `out` from `sps()`, they should supply `site_selected` ")
    }
    X_use  <- X
    ss_use <- rep(0, nrow(X))
    ss_use[site_selected] <- 1
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

  out <- list("average_site_ATE" = out_overall, "site_specific_ATE" = out_each, "bet_se" = sqrt(bet_var), "estimates_selected" = estimates_selected)

  class(out) <- c(class(out), "sps_estimator")

  return(out)
}
