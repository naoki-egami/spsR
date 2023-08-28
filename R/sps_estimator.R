#' Synthetic Purposive Sampling Estimator for the Average-Site ATE
#' @param out Output from function \code{sps()}
#' @param estimates_selected data.frame with two columns: the first column represents estimates of the site-specific ATEs for the selected sites and the second column represents its corresponding standard error. The number of rows is equal to the number of the selected sites and \code{rownames(estimates_selected)} should be names of the selected sites.
#' @param X (Optional. Use this only when sites are not selected based on \code{sps()}. Default = \code{NULL}) Site-level variables for the target population of sites. Row names should be names of sites.
#' @param selected_sites (Optional. Use this only when sites are not selected based on \code{sps()}. Default = \code{NULL}) Names of sites users selected. This should be a subset of rownames(X).
#' @return \code{sps_estimator} returns an object of \code{sps_estimator} class.
#'  \itemize{
#'    \item \code{average_site_ATE}: An estimate of the average-site ATE and its corresponding standard error.
#'    \item \code{site_specific_ATE}: Estimates and standard errors for the site-specific ATEs in non-selected sites.
#'    \item \code{bet_se}: Estimated between-site standard errors.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export

sps_estimator <- function(out = NULL, estimates_selected = NULL, X = NULL, selected_sites = NULL){

  ## Housekeeping
  ### We order names in the order of X

  if(any(is.na(estimates_selected))){
    stop(" `estimates_selected` contains missing data. Please supply `estimates_selected` without missing values. ")
  }
  # cat(" `estimates_selected` should have point estimates and their standard errors ")

  if(is.null(out) == FALSE){  ## site selection is done with sps()

    ## out
    if(("sps" %in% class(out)) == FALSE){
      stop(" `out` should be an output from function `sps()` ")
    }
    selected_sites <- out$selected_sites

    ## estimates_selected
    if(setequal(rownames(estimates_selected), selected_sites) == FALSE){
      stop(" `rownames(estimates_selected)` should match to `selected_sites` ")
    }

    X_use  <- out$internal$X
    ss_use <- out$internal$ss
    N <- out$internal$N
    X_selected <- X_use[ss_use == 1, , drop = FALSE]

    estimates_selected <- estimates_selected[match(rownames(X_selected), rownames(estimates_selected)), ]
    # B[match(A, B)] will match the order of B to the order of A
    # estimates_selected <- estimates_selected[match(data_name, rownames(estimates_selected)),]

  }else if(is.null(X) == FALSE){  ## site selection is not done with sps()

    if(is.null(X) == TRUE){
      stop(" When sites are selected not using `sps()`, please supply `X`  ")
    }
    if(is.null(selected_sites) == TRUE){
      stop(" When sites are selected not using `sps()`, please supply `selected_sites`  ")
    }
    if(setequal(rownames(estimates_selected), selected_sites) == FALSE){
      stop(" `rownames(estimates_selected)` should match to `selected_sites` ")
    }
    if(all(selected_sites %in% rownames(X)) == FALSE){
      stop(" `selected_sites` should be a subset of `rownames(X)` ")
    }

    selected_sites_ind <- which(rownames(X) %in% selected_sites)

    X_use  <- X
    ss_use <- rep(0, nrow(X))
    ss_use[selected_sites_ind] <- 1
    N <- nrow(X_use)
    X_selected <- X_use[ss_use == 1, , drop = FALSE]

    estimates_selected <- estimates_selected[match(rownames(X_selected), rownames(estimates_selected)), ]

  }else{
    stop(" `out` or `X_selected` should be specified.  ")
  }

  ###
  estimate <- estimates_selected[, 1]
  se <- estimates_selected[, 2]

  # if(is.null(W) == TRUE){
  # first estimate weights
  out_W <- sps_weights(X = X_use, site = ss_use, site_name = NULL)
  W <- out_W$W
  RMSE_X  <- out_W$RMSE
  #}

  # Estimate between-site variance
  bet_var <- between_var_sps(estimate = estimate, se = se, X_selected = X_selected)

  N_S <- length(estimate)
  N_R <- N - N_S

  # Each Site
  estimate_each <- t(W) %*% estimate
  se_each <- sqrt(t(W^2) %*% c(se^2) + bet_var)
  out_each <- cbind(estimate_each, se_each); colnames(out_each) <- c("Estimate", "Std. Error")
  rownames(out_each) <- setdiff(rownames(X_use), rownames(X_selected))

  # ###############
  # Overall
  # ###############
  std_w <- (1 + apply(W, 1, sum))/sum(1 + apply(W, 1, sum))
  estimate_overall <- sum(std_w * estimate)
  se_overall <- sqrt(sum(std_w^2 *(se^2 + bet_var)))
  out_overall <- c(estimate_overall, se_overall); names(out_overall) <- c("Estimate", "Std. Error")

  out <- list("average_site_ATE" = out_overall, "site_specific_ATE" = out_each,
              "bet_se" = sqrt(bet_var), "estimates_selected" = estimates_selected,
              "RMSE_X" = RMSE_X)

  class(out) <- c(class(out), "sps_estimator")

  return(out)
}
