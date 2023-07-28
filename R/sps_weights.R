sps_weights <- function(X, site, site_name = NULL){

  if(is.null(site_name)){
    site_name <- seq(1:length(site))
  }

  N <- nrow(X)
  s_site <- which(site == 1)
  X_s <- X[s_site, , drop = FALSE]
  t_site <- setdiff(seq(1:nrow(X)), s_site)
  X_t <- X[t_site, , drop = FALSE]
  s <- sum(site)

  # Estimate for each target site
  RMSE <- c()
  # Rows: Selected Sites
  # Cols: Target Sites
  W <- matrix(NA, nrow = s, ncol = length(t_site))
  for(j in 1:length(t_site)){
    target_X_use <- as.numeric(X_t[j,])
    sol <- synth_base(target_X = target_X_use, X_s = X_s)
    RMSE[j] <- sol$RMSE
    W[1:s,j] <- sol$w
  }
  names(RMSE) <- site_name[t_site]
  colnames(W) <- site_name[t_site]
  rownames(W) <- site_name[s_site]

  out <- list("W" = round(W, 2), "RMSE" = RMSE)
  return(out)
}
