# basic Synthetic Site Weights
synth_base <- function(target_X, X_s){

  # weights
  w <- Variable(nrow(X_s))

  # constraints
  constraint_1 <- w >= 0
  constraint_2 <- sum(w) == 1
  constraints <- list(constraint_1, constraint_2)

  obj <- Minimize(cvxr_norm(target_X - t(X_s) %*% w))

  # Optimization
  p <- Problem(obj, constraints)
  res <- solve(p)
  w_sol <- res$getValue(w)

  RMSE <- sqrt(mean((target_X - t(X_s) %*% w_sol)^2))

  out <- list("w" = w_sol, "RMSE" = RMSE)

  return(out)
}

# between-site variance
between_var_sps <- function(estimate, se, X_selected, site_name = NULL){

  s <- nrow(X_selected)
  ind <- seq(1:s)

  e_res <- c()
  within_var <- c()
  W_sq <- c()
  for(i in 1:s){
    t_use <- i
    s_use <- setdiff(ind, i)
    X_s <- X_selected[s_use, , drop = FALSE]
    target_X_use <- as.numeric(X_selected[t_use,])
    sol <- synth_base(target_X = target_X_use, X_s = X_s)
    W0 <- sol$w

    estimate_target <- sum(W0 * estimate[s_use])
    e_res[i] <- estimate[t_use] - estimate_target
    within_var[i] <- se[t_use]^2 + sum((W0^2) * (se[s_use]^2))
    W_sq[i] <- sum((W0^2))
  }

  # Estimate between-site variance
  num  <- sum(e_res^2) - (sum(within_var))
  deno <- sum(1 + sum(W_sq))
  bet_var <- max(0, num/deno)

  return(bet_var)
}

## Stratify

#' @export
stratify_sps <- function(X, num_site = NULL, condition = NULL){

  col_use <- condition[1]

  if((condition[2] %in% c("larger than or equal to",
                          "smaller than or equal to",
                          "between")) == FALSE){
    stop("The second element of `condition` should be `larger than or equal to`, `smaller than or equal to`, or `between` ")
  }
  if((num_site[1] %in% c("at least",
                         "at most")) == FALSE){
    stop("The first element of `num_site` should be `larger than or equal to` or `smaller than or equal to` ")
  }

  type_use <- case_when(condition[2] == "larger than or equal to" ~ "geq",
                        condition[2] == "smaller than or equal to" ~ "leq",
                        condition[2] == "between" ~ "between")
  if(type_use == "between"){
    value <- c(condition[3], condition[4])
  }else{
    value <- c(condition[3])
  }

  # C_l <- list()
  # for(i in 1:length(columns)){
  #   col_use <- columns[i]
  #   C_l[[i]] <- stratify_base(X = X, columns = col_use, type = type_use, value = value)
  # }
  # C_use <- do.call("rbind", C_l)

  C_use <- stratify_base(X = X, columns = col_use, type = type_use, value = value)

  if(num_site[1] == "at least"){
    c0_value <- as.numeric(num_site[2])
  }else if(num_site[1] == "at most"){
    c0_value <- -1*as.numeric(num_site[2])
  }
  c0_use <- rep(c0_value, nrow(C_use))

  st <- list("C" = C_use, "c0" = c0_use)
  class(st) <- c(class(st), "st_direct")
  return(st)
}

stratify_base <- function(X, columns = NULL, type = "geq", value = 1){

  if(is.null(columns)){
    columns <- seq(1:ncol(X))
  }
  C_mat <- matrix(NA, nrow = length(columns), ncol = nrow(X))
  for(i in 1:length(columns)){
    if(type == "geq"){
      C_mat[i,1:ncol(C_mat)] <- as.numeric(X[, columns[i]] >= value)
    }else if(type == "leq"){
      C_mat[i,1:ncol(C_mat)] <- as.numeric(X[, columns[i]] <= value)
    }else if(type == "between"){
      C_mat[i, 1:ncol(C_mat)] <- as.numeric(X[, columns[i]] >= value[1])*as.numeric(X[, columns[i]] <= value[2])
    }
  }
  return(C_mat)
}

