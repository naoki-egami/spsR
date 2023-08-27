#' Synthetic Purposive Sampling: Site Selection for External Validity
#' @param X Site-level variables for the target population of sites. Row names should be names of sites. X cannot contain missing data.
#' @param N_s Number of study sites to be selected.
#' @param stratify (Optional. Default = \code{NULL}) Output from function \code{stratify_sps()}. This argument helps users incorporate practical and logistical constraints. See examples on \url{http://naokiegami.com/spsR/articles/stratify_sps.html}
#' @param site_include (Optional. Default = \code{NULL}) Names of sites users want to always include (or have already selected).
#' @param site_exclude (Optional. Default = \code{NULL}) Names of sites users want to always exclude.
#' @param lambda Values of the tuning parameters. If users want to change how to balance three parts of the objective function, they can change \code{lambda}. Default values are \code{c(1, 1, 0)}. Users who want to fine-tune the tuning parameters, please see the methodological details in Egami and Lee (2023+) (https://naokiegami.com/paper/sps.pdf).
#' @param seed Numeric. \code{seed} used internally. Default = \code{1234}.
#' @import CVXR
#' @import ggplot2
#' @import GGally
#' @import spsRdata
#' @importFrom metafor rma
#' @importFrom grDevices adjustcolor
#' @importFrom stats pnorm qnorm sd as.formula model.matrix var
#' @importFrom utils combn
#' @importFrom dplyr case_when
#' @return \code{sps} returns an object of \code{sps} class.
#'  \itemize{
#'    \item \code{selected_sites}: Names of sites the SPS algorithm selected.
#'    \item \code{W}: Estimated weights to approximate non-selected sites using selected sites. \code{W} will be used in the subsequent estimation of the average-site ATE.
#'    \item \code{obj}: Estimated values of the objective function, separately for three parts.
#'    \item \code{internal}: Objects useful for internal use of the function.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
#'

sps <- function(X, N_s,
                stratify = NULL,
                site_include = NULL,
                site_exclude = NULL,
                lambda = c(1, 1, 0),
                seed = 1234){

  # ###############
  # Housekeeping
  # ###############
  ## X
  ## X contains NA
  if(any(is.na(X))){
    stop(" X contains missing data. Please supply X without missing values. Consider using `impute_var()` in R package `spsRdata`. ")
  }

  ## factor or character
  if(all(sapply(X, class) == "numeric") == FALSE){
    stop(" X contains `factor` or `character` variables. Before using sps(), please convert them into numeric or binary variables. ")
  }

  ## sd of each variable is so differnet
  X_sd <- apply(X, 2, sd)
  if(max(X_sd)/min(X_sd) >= 100){
    warning(" Some variables in X have standard deviation more than 100 times larger than other variables. This might cause estimation problem. Please consider using `scale()` to make standard deviations of variables comparable. ")
  }

  ## Transform to matrix
  X <- as.matrix(X)
  class(X) <- "matrix"

  ## If there is only one variable
  if(ncol(X) == 1){
    X_orig <- X
    X <- cbind(1, X)
  }else{
    X_orig <- X
  }

  ## row names
  if(any(is.null(rownames(X)))){
    warning(" rownames(X) is NULL. We use 1, 2, ... to represent names of sites. ")
    rownames(X) <- seq(from = 1, to = nrow(X))
  }

  ## N_s
  if(N_s > nrow(X)){
    stop(" N_s cannot be larger than the number of sites in X. ")
  }

  if(N_s <= length(site_include)){
    stop(" N_s should be larger than length(site_include) ")
    # N_s <- N_s + length(site_include)
    # N_s is defined including site_include
  }

  ## stratify
  if(is.null(stratify) == TRUE){
    C   <- NULL
    c0  <- NULL
  }else{
    if("stratify_sps" %in% class(stratify)){
      C  <- stratify$C
      c0 <- stratify$c0
    }else{
      C_l <- list()
      c0 <- c()
      for(z in 1:length(stratify)){
        C_l[[z]] <- stratify[[z]]$C
        c0 <- c(c0, stratify[[z]]$c0)
      }
      C  <- do.call("rbind", C_l)
    }
  }

  # site_include
  if(all(is.null(site_include)) == FALSE){
    if(all(site_include %in% rownames(X)) == FALSE){
      stop(" site_include should be a subset of rownames(X) ")
    }else{
      site_include_ind <- which(rownames(X) %in% site_include)
    }
  }else{
    site_include_ind <- NULL
  }

  # site_exclude
  if(all(is.null(site_exclude)) == FALSE){
    if(all(site_exclude %in% rownames(X)) == FALSE){
      stop(" site_exclude should be a subset of rownames(X) ")
    }else{
      site_exclude_ind <- which(rownames(X) %in% site_exclude)
    }
  }else{
    site_exclude_ind <- NULL
  }

  ###### (End of Housekeeping)

  N <- nrow(X)
  L <- ncol(X)

  S <- Variable(N, integer = TRUE)
  W <- Variable(N, N)
  Q <- Variable(N, N)
  Z <- Variable(N, L)

  # constraints
  co_1 <- W >= 0
  co_2 <- Q >= 0
  co_3 <- S >= 0
  co_4 <- S <= 1
  co_5 <- sum(S) == N_s # the number of selected sites
  co_6 <- sum_entries(Q, axis = 2) == (1 - S)
  co_7 <- list()
  for(j in 1:N){
    co_7[[j]] <- Q[j,] <= S[j]
  }
  co_8 <- Q <= W
  co_9 <- list()
  for(j in 1:N){
    co_9[[j]] <- Q[j,] >= W[j,] - (1 - S[j])
  }
  co_10 <- list()
  for(k in 1:N){
    co_10[[k]] <- W[,k] <= (1 - S[k])
  }
  co_11 <- list()
  for(l in 1:L){
    co_11[[l]] <- Z[,l] == (1-S)*X[,l] - t(Q) %*% X[,l]
  }
  if(all(is.null(site_include_ind)) == FALSE){
    co_12 <- S[site_include_ind] == 1
  }else{
    co_12 <- NULL
  }
  if(all(is.null(site_exclude_ind)) == FALSE){
    co_13 <- S[site_exclude_ind] == 0
  }else{
    co_13 <- NULL
  }

  # Note: W_{jk} = 0 when k is not selected.
  # Note: W_{jk} can be arbitrary when j is not selected

  constraints <- c(list(co_1, co_2, co_3, co_4, co_5, co_6, co_8),
                   co_7, co_9, co_10, co_11, co_12, co_13)

  # Stratification
  if(is.null(C) == FALSE){
    co_14 <- C%*% S >= c0
    constraints <- c(constraints, co_14)
  }

  # interpolation bias (Distance from each non-selected sites)
  Xs <- matrix(NA, nrow = N, ncol = N)
  for(j in 1:N){
    for(k in 1:N){
      Xs[j,k] <- mean((X[j, ] - X[k,])^2)
    }
  }

  # Objective Function
  # obj <- Minimize( sum_squares(Z)/(L*(N-s)))

  # Check Infeasibility
  check_inf <- Problem(Minimize(0), constraints)
  res_inf <- solve(check_inf)
  if(res_inf$status != "optimal"){
    stop(" Constraints are infeasible. Please check whether conditions in `stratify` are feasible. ")
  }

  cat("Selecting Study Sites...")
  if(lambda[3] == 0){
    obj <- Minimize( lambda[1]*sum_squares(Z)/(L*(N-N_s)) + lambda[2]*(sum_entries(Q*Xs)/(N-N_s)))
  }else{
    obj <- Minimize( lambda[1]*sum_squares(Z)/(L*(N-N_s)) + lambda[2]*(sum_entries(Q*Xs)/(N-N_s)) + lambda[3]*sum_squares(Q)/(N-N_s) )
  }
  p <- Problem(obj, constraints)

  set.seed(seed)
  res <- solve(p)
  ss0 <- res$getValue(S)
  Z_out <- res$getValue(Z)
  W_out <- res$getValue(W)
  Q_out <- res$getValue(Q)

  obj_1 <- sum(Z_out^2)/(L*(N-N_s))
  obj_2  <- sum(Q_out*Xs)/(N-N_s)
  obj_3  <- sum(Q_out^2)/(N-N_s)
  ss_out   <- round(ss0)

  if(any(is.null(rownames(X_orig)) == TRUE)){
    rownames(X_orig) <- seq(1:nrow(X_orig))
  }
  selected_sites <- rownames(X_orig)[ss_out == 1]

  internal_use <- list("ss" = ss_out,
                       "X" = X_orig,
                       "Q" = Q_out,
                       "N_s" = N_s,
                       "N" = N)

  if(any(is.na(selected_sites))){
    message("\n The optimization fails. \n")
    if(res$status == "infeasible" | res$status == "infeasible_inaccurate"){
      message("\n Constraints are infeasible. Please check whether conditions in`stratify` are feasible. \n")
    }else if(res$status == "unbounded" | res$status == "unbounded_inaccurate"){
      message("\n The objective function is unbounded. Some problems in X. \n")
    }else if(res$status == "solver_error"){
      message("\n Errors in the underlying package `CVXR`. \n")
    }
  }

  out <- list("selected_sites" = selected_sites,
              "W" = W_out,
              "obj" = c(obj_1, obj_2, obj_3),
              "internal" = internal_use)
  class(out) <- c(class(out), "sps")
  return(out)
}
