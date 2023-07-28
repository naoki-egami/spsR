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
stratify <- function(X, columns = NULL, type = "geq", value = 1){

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

