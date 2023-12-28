#' Summary function
#' @param object Output from \code{sps_estimator()}.
#' @param ci (Default = \code{0.95}.) Coverage of the estimated confidence intervals.
#' @param ... Other arguments.
#' @export

summary.sps_estimator <- function(object, ci = 0.95, ...){

  out_main <- object$average_site_ATE
  alpha_h <- 1 - (1 - ci)/2
  ci_low  <- out_main[1] - qnorm(alpha_h)*out_main[2]
  ci_high <- out_main[1] + qnorm(alpha_h)*out_main[2]
  p_v <- 1 - pnorm(abs(out_main[1]/out_main[2]))
  out_tab <- c(out_main, ci_low, ci_high, p_v)

  # add significance
  sig <- rep("", length(out_tab[5]))
  sig[out_tab[5]  < 0.001] <- "***"
  sig[out_tab[5]  >= 0.001 & out_tab[5]  < 0.01] <- "**"
  sig[out_tab[5]  >= 0.01 & out_tab[5]  < 0.05] <- "*"
  sig[out_tab[5]  >= 0.05 & out_tab[5]  < 0.1] <- "."
  out_tab_print <- as.data.frame(matrix(out_tab, nrow = 1))
  out_tab_print$Sig <- sig
  names(out_tab_print) <- c("Estimate", "Std. Error", "CI Lower", "CI Upper", "p value", "")

  if(is.null(object$subgroup) == FALSE){
    sub_main <- object$subgroup_average_site_ATE
    alpha_h <- 1 - (1 - ci)/2
    ci_low  <- sub_main[, 1] - qnorm(alpha_h)*sub_main[, 2]
    ci_high <- sub_main[, 1] + qnorm(alpha_h)*sub_main[, 2]
    p_v <- 1 - pnorm(abs(sub_main[, 1]/sub_main[, 2]))
    out_tab_sub <- cbind(sub_main, ci_low, ci_high, p_v)

    # add significance
    sig <- rep("", length(out_tab_sub[, 5]))
    sig[out_tab_sub[,5]  < 0.001] <- "***"
    sig[out_tab_sub[,5]  >= 0.001 & out_tab_sub[,5]  < 0.01] <- "**"
    sig[out_tab_sub[,5]  >= 0.01 & out_tab_sub[,5]  < 0.05] <- "*"
    sig[out_tab_sub[,5]  >= 0.05 & out_tab_sub[,5]  < 0.1] <- "."
    out_tab_sub_print <- as.data.frame(as.matrix(out_tab_sub))
    out_tab_sub_print$Sig <- sig
    names(out_tab_sub_print) <- c("Estimate", "Std. Error", "CI Lower", "CI Upper", "p value", "")
  }

  cat("\n")
  cat("Average-Site ATE:\n")
  print(out_tab_print, row.names = FALSE)

  if(is.null(object$subgroup) == FALSE){
    cat("\nSubgroup Average-Site ATE:\n")
    print(out_tab_sub_print, row.names = TRUE)
  }
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

  invisible(out_tab_print)

}
