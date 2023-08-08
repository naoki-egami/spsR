#' Synthetic Purposive Sampling: Site Selection for External Validity
#' @param X Site-level variables for the target population of sites. Row names should be names of sites.
#' @param N_s Number of study sites to be selected.
#' @param stratify (Optional. Default = \code{NULL}) Output from function \code{stratify_sps()}. This argument helps users incorporate practical and logistical constraints. See examples on \url{http://naokiegami.com/spsR/articles/stratify_sps.html}
#' @param site_selected (Optional. Default = \code{NULL}) Names of sites users always want to select (or have already selected).
#' @param site_unavailable (Optional. Default = \code{NULL}) Names of sites users cannot select.
#' @param lambda Values of the tuning parameters. If users want to change how to balance three parts of the objective function, they can change \code{lambda}. Default values are \code{c(1, 1, 0.1)}. Users who want to fine-tune the tuning parameters, please see examples on \url{http://naokiegami.com/spsR/articles/methods_guides.html}.
#' @param seed Numeric. \code{seed} used internally. Default = \code{1234}.
#' @import CVXR
#' @import ggplot2
#' @import GGally
#' @import spsRdata
#' @importFrom metafor rma
#' @importFrom grDevices adjustcolor
#' @importFrom stats pnorm qnorm
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
                site_selected = NULL,
                site_unavailable = NULL,
                lambda = c(1, 1, 0.05),
                seed = 1234){

  # Housekeeping

  ## X
  ### Need to add Missing data
  X <- as.matrix(X)
  class(X) <- "matrix"
  if(ncol(X) == 1){
    X_orig <- X
    X <- cbind(1, X)
  }else{
    X_orig <- X
  }

  ## N_s
  if(N_s <= length(site_selected)){
    stop(" N_s should be larger than length(site_selected) ")
    # N_s <- N_s + length(site_selected)
    # N_s is defined including site_selected
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
  # site_selected
  if(all(is.null(site_selected)) == FALSE){
    if(all(site_selected %in% rownames(X)) == FALSE){
      stop(" site_selected should be a subset of rownames(X) ")
    }else{
      site_selected_ind <- which(site_selected %in% rownames(X))
    }
  }else{
    site_selected_ind <- NULL
  }

  # site_unavailable
  if(all(is.null(site_unavailable)) == FALSE){
    if(all(site_unavailable %in% rownames(X)) == FALSE){
      stop(" site_unavailable should be a subset of rownames(X) ")
    }else{
      site_unavailable_ind <- which(site_unavailable %in% rownames(X))
    }
  }else{
    site_unavailable_ind <- NULL
  }


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
  if(all(is.null(site_selected_ind)) == FALSE){
    co_12 <- S[site_selected_ind] == 1
  }else{
    co_12 <- NULL
  }
  if(all(is.null(site_unavailable_ind)) == FALSE){
    co_13 <- S[site_unavailable_ind] == 0
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

  obj <- Minimize( lambda[1]*sum_squares(Z)/(L*(N-N_s)) +lambda[2]*(sum_entries(Q*Xs)/(N-N_s)) + lambda[3]*sum_squares(Q)/(N-N_s) )

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

  out <- list("selected_sites" = selected_sites,
              "W" = W_out,
              "obj" = c(obj_1, obj_2, obj_3),
              "internal" = internal_use)
  class(out) <- c(class(out), "sps")
  return(out)
}
