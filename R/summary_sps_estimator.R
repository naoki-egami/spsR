#' Summary function
#' @param object the output from sps_estimator()
#' @param ci digits used for a table of covariate profiles. Default = 2.
#' @param ... Other arguments.
#' @export

summary.sps_estimator <- function(object, ci = 0.95, ...){

  out_main <- object$overall
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

  cat("\n")
  # cat("Estimate:\n")
  print(out_tab_print, row.names = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

  invisible(out_tab_print)

}
