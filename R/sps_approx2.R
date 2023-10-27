#' Synthetic Purposive Sampling: Site Selection for External Validity
#' @param X Site-level variables for the target population of sites. Row names should be names of sites. X cannot contain missing data.
#' @param N_s Number of study sites to be selected.
#' @param stratify (Optional. Default = \code{NULL}) Output from function \code{stratify_sps()}. This argument helps users incorporate practical and logistical constraints. See examples on \url{http://naokiegami.com/spsR/articles/stratify_sps.html}
#' @param site_include (Optional. Default = \code{NULL}) Names of sites users want to always include (or have already selected).
#' @param site_exclude (Optional. Default = \code{NULL}) Names of sites users want to always exclude.
#' @param lambda Values of the tuning parameters. If users want to change how to balance three parts of the objective function, they can change \code{lambda}. Default values are \code{c(1, 1, 0)}. Users who want to fine-tune the tuning parameters, please see the methodological details in Egami and Lee (2023+) (https://naokiegami.com/paper/sps.pdf).
#' @param seed Numeric. \code{seed} used internally. Default = \code{1234}.
#' @param max_iter Numeric. The number of iterations used in the optimization. Default = \code{10}.
#' @param solver Solver we use in the internal CVXR optimization. Default = \code{ECOS_BB}. See the CVXR website for information on other solvers.
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

sps_approx2 <- function(X, N_s,
                        stratify = NULL,
                        site_include = NULL,
                        site_exclude = NULL,
                        lambda = c(1, 1, 0),
                        seed = 1234,
                        max_iter = 10,
                        solver = "ECOS_BB"){

  # ###############
  # Housekeeping
  # ###############
  ## X
  ## X contains NA
  if(any(is.na(X))){
    stop(" X contains missing data. Please supply X without missing values. Consider using `impute_var()` in R package `spsRdata`. ")
  }

  ## factor or character
  if(any(sapply(X, class) %in% c("factor", "character")) == TRUE){
    stop(" X contains `factor` or `character` variables. Before using sps(), please convert them into numeric or binary variables. ")
  }

  ## sd of each variable is so differnet
  X_sd <- apply(X, 2, sd)
  if(max(X_sd)/min(X_sd) >= 100){
    warning(" \n Some variables in X have standard deviation more than 100 times larger than other variables. This might cause estimation problem. Please consider using `scale()` to make standard deviations of variables comparable. ")
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
    warning(" \n rownames(X) is NULL. We use 1, 2, ... to represent names of sites. ")
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
      if(is.list(stratify) == FALSE){
        stop(" If there is more than one stratification, `stratify` should be a list where each element is an output from `stratify_sps()`. ")
      }else{
        for(z in 1:length(stratify)){
          if(("stratify_sps" %in% class(stratify[[z]])) == FALSE){
            stop(" If there is more than one stratification, `stratify` should be a list where each element is an output from `stratify_sps()`. ")
          }
        }
      }

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
  # Q <- Variable(N, N)
  # Z <- Variable(N, L)
  # H <- Variable(N, N)

  # constraints
  co_1 <- W >= 0
  # co_2 <- Q >= 0
  co_3 <- S >= 0
  co_4 <- S <= 1
  co_5 <- sum(S) == N_s # the number of selected sites
  # co_6 <- sum_entries(Q, axis = 2) == (1 - S)
  co_6 <- sum_entries(W, axis = 2) == (1 - S)
  # co_7 <- list()
  # for(j in 1:N){
  #   co_7[[j]] <- Q[j,] <= S[j]
  # }
  # co_8 <- Q <= W
  # co_9 <- list()
  # for(j in 1:N){
  #   co_9[[j]] <- Q[j,] >= W[j,] - (1 - S[j])
  # }
  co_10 <- list()
  for(k in 1:N){
    co_10[[k]] <- W[,k] <= (1 - S[k])
  }
  # co_11 <- list()
  # for(l in 1:L){
  #   co_11[[l]] <- Z[,l] == (1-S)*X[,l] - t(Q) %*% X[,l]
  # }
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

  constraints <- c(list(co_1, co_3, co_4, co_5, co_6), co_10, co_12, co_13)

  # co_2,
  # co_6,
  # co_8),
  # co_7, co_9, co_10, # co_11,

  # Stratification
  if(is.null(C) == FALSE){
    co_14 <- C%*% S >= c0
    constraints <- c(constraints, co_14)
  }

  ## New Constraints
  # W_{jk} = 0 when S_j = 0
  co_15 <- list()
  for(j in 1:N){
    co_15[[j]] <- W[j,] <= S[j]
  }

  ## H_{ij}
  # co_16 <- co_17 <- co_18 <- co_19 <- list()
  # for(i in 1:N){
  #   co_16[[i]] <- H[i, ] <= S[i]
  # }
  # for(j in 1:N){
  #   co_17[[j]] <- H[,j] <= S[j]
  # }
  # for(i in 1:N){
  #   co_18[[i]] <- list()
  #   for(j in 1:N){
  #     co_18[[i]][[j]] <- H[i,j] >= S[i] + S[j] - 1
  #   }
  # }
  # for(i in 1:N){
  #   co_19[[i]] <- H[i,i] == 0
  # }

  constraints <- c(constraints, co_15)

  # interpolation bias (Distance from each non-selected sites)
  Xs <- matrix(NA, nrow = N, ncol = N)
  for(j in 1:N){
    for(k in 1:N){
      Xs[j,k] <- mean((X[j, ] - X[k,])^2)
    }
  }
  X_multi <- matrix(NA, nrow = N, ncol = N)
  for(j in 1:N){
    for(k in 1:N){
      X_multi[j,k] <-  sum(X[j, ]*X[k, ])
    }
  }

  # Objective Function
  # obj <- Minimize( sum_squares(Z)/(L*(N-s)))

  # Problem is big!
  if(choose(N, N_s) > 1000000){
    choose_p <- paste0("choose(", N, ", ", N_s, ") = ", choose(N, N_s))
    cat("Note 1: The optimization problem is big. There are more than 1 million unique combinations of sites (roughly ", choose_p, "combinations). Thus, it will take some computational time to solve the problem. \n")
    cat("Note 2: If users want to reduce the computational time, they can divide the population of sites into sub-populations and run sps() within each sub-populations. For example, if users want to select 9 sites out of 50 sites, they can divide 50 sites into three sub-populations and select 3 sites from each sub-population of sites by running sps() separately. \n\n")
  }

  # Check Infeasibility
  cat("Checking whether constraints specified in `stratify` are feasible...\n")
  if(is.null(stratify) == FALSE){
    check_inf <- Problem(Minimize(0), constraints = constraints)
    res_inf   <- solve(check_inf,  solver = solver)
    if(res_inf$status != "optimal"){
      stop(" Constraints are infeasible. Please check whether conditions in `stratify` are feasible. ")
    }
  }

  # New Objective Function
  X_j_2_mean <- apply(X^2, 1, sum)

  Main_1 <- sum_entries((1-S) * X_j_2_mean) - 2*sum_entries(W * X_multi) + N*sum_entries(t(W) %*% X_j_2_mean)
  Main_2 <- sum_entries(W*Xs)
  Main_3 <- sum_squares(W)

  cat("Selecting Study Sites...\n")
  if(lambda[3] == 0){
    obj <- Minimize( lambda[1]*Main_1/(L*(N-N_s)) + lambda[2]*(Main_2/(N-N_s)))
  }else{
    obj <- Minimize( lambda[1]*Main_1/(L*(N-N_s)) + lambda[2]*(Main_2/(N-N_s)) + lambda[3]*Main_3/(N-N_s) )
  }
  p <- Problem(obj, constraints = constraints)

  set.seed(seed)
  res <- solve(p, solver = solver, warm_start = FALSE)

  ## Checking Results
  try_again <- FALSE
  try_num <- 0
  if(res$status == "infeasible" | res$status == "infeasible_inaccurate"){
    if(is.null(stratify) == TRUE){
      try_again <- TRUE
    }else{
      message("\n Constraints are infeasible. Please check whether conditions in `stratify` are feasible. \n")
    }
  }
  if(try_again == TRUE){
    try_num   <- 1
    feastol_use <- 1e-08
    reltol_use  <- 1e-08
    abstol_use  <- 1e-08
    while(try_again == TRUE & try_num <= max_iter){
      set.seed(seed + try_num)
      feastol_use <- feastol_use*10
      reltol_use  <- reltol_use*10
      abstol_use <- abstol_use*10
      res <- solve(p,
                   solver  = solver,
                   feastol = feastol_use,
                   reltol  = reltol_use,
                   abstol  = abstol_use)

      if(res$status == "infeasible" | res$status == "infeasible_inaccurate"){
        if(is.null(stratify) == TRUE){
          try_again <- TRUE
        }
      }else{
        try_again <- FALSE
      }
      try_num <- try_num + 1
    }
  }

  ## Collecting Results
  ss0 <- res$getValue(S)
  Z_out <- NA
  W_out <- res$getValue(W)
  Q_out <- NA

  obj_1  <- sum(Z_out^2)/(L*(N-N_s))
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
                       "N" = N,
                       "W_internal" = W_out,
                       "try_num" = try_num)

  if(any(is.na(selected_sites))){
    message("\n The optimization fails. \n")
    if(res$status == "infeasible" | res$status == "infeasible_inaccurate"){
      if(is.null(stratify) == FALSE){
        message("\n Constraints are infeasible. Please check whether conditions in `stratify` are feasible. \n")
      }
    }else if(res$status == "unbounded" | res$status == "unbounded_inaccurate"){
      message("\n The objective function is unbounded. Some problems in X. \n")
    }else if(res$status == "solver_error"){
      message("\n Errors in the underlying package `CVXR`. \n")
    }
  }

  # compute W
  if(any(is.na(selected_sites)) == FALSE){
    out_W <- sps_weights(X = X, site = ss_out, site_name = rownames(X))
    W_out <- out_W$W
  }

  out <- list("selected_sites" = selected_sites,
              "W" = W_out,
              "obj" = c(obj_1, obj_2, obj_3),
              "internal" = internal_use)
  class(out) <- c(class(out), "sps")
  return(out)
}
