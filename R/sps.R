#' Synthetic Purposive Sampling: Site Selection for External Validity
#' @param X Site-level covariates for the target population of sites
#' @param N_s Number of study sites to be selected
#' @param C Linear constraints (left-hand side)
#' @param c0  Linear constraints (right-hand side)
#' @param site_selected A
#' @param site_unavailable A
#' @param lambda A
#' @param seed Numeric. `seed` used internally. Default = 1234.
#' @import CVXR
#' @import ggplot2
#' @import GGally
#' @importFrom metafor rma
#' @importFrom grDevices adjustcolor
#' @importFrom stats pnorm
#' @importFrom utils combn
#' @return \code{sps} returns an object of \code{sps} class.
#'  \itemize{
#'    \item \code{ss}: Estimated external robustness.
#'    \item \code{W}: Estimates of the pAMCE for all factors in each bootstrap sample.
#'    \item \code{obj}: Estimated values of the objective function
#'    \item \code{X}: X
#'    \item \code{Q}: Q
#'    \item \code{N}: N
#'    \item \code{N_s}: N_s
#'    \item \code{...}: Values for internal use.
#'  }
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
#'

sps <- function(X, N_s,
                C = NULL, c0 = NULL,
                site_selected = NULL,
                site_unavailable = NULL,
                lambda = c(1, 1, 0.1),
                seed = 1234){

  # Housekeeping
  X <- as.matrix(X)
  class(X) <- "matrix"

  if(ncol(X) == 1){
    X_orig <- X
    X <- cbind(1, X)
  }else{
    X_orig <- X
  }

  N_s <- N_s + length(site_selected)

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
  if(is.null(site_selected) == FALSE){
    co_12 <- S[site_selected] == 1
  }else{
    co_12 <- NULL
  }
  if(is.null(site_unavailable) == FALSE){
    co_13 <- S[site_unavailable] == 0
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

  obj <- Minimize( lambda[1]*sum_squares(Z)/(L*(N-N_s)) +
                     lambda[2]*(sum_entries(Q*Xs)/(N-N_s)) +
                     lambda[3]*sum_squares(Q)/(N-N_s) )

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

  out <- list("ss" = ss_out,
              "W" = W_out,
              "obj" = c(obj_1, obj_2, obj_3),
              "X" = X_orig,
              "Q" = Q_out,
              "N_s" = N_s,
              "N" = N)
  class(out) <- c(class(out), "sps")
  return(out)
}
