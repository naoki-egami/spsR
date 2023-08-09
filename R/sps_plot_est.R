#' Plot function for \code{sps_estimator}
#' @param estimates_selected data.frame with two columns: the first column represents estimates of the site-specific ATEs for the selected sites and the second column represents its corresponding standard error. The number of rows is equal to the number of the selected sites and \code{rownames(estimates_selected)} should be names of the selected sites.
#' @param sps_est Output from function \code{sps_estimator()}
#' @return \code{sps_plot_est} plots the estimates and 95% confidence interval bars for site-specific ATEs and an average-site ATE.
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export

sps_plot_est <- function(estimates_selected, sps_est){
  est <- c(estimates_selected[,1], sps_est$overall[1])
  se  <- c(estimates_selected[,2], sps_est$overall[2])
  group <- c(rep('Site-Specific ATEs', nrow(estimates_selected)), 'Average-Site ATE')
  site  <- c(row.names(estimates_selected), 'Average-Site')

  pdata <- data.frame(est, se, group, site)
  pdata$group <- factor(pdata$group, levels = c('Site-Specific ATEs', 'Average-Site ATE'))
  pdata$site  <- factor(pdata$site, levels = pdata[order(pdata$est, decreasing = T), 'site'])

  ggplot(data = pdata,
         aes(x = pdata$site,
             y = pdata$est,
             ymin = pdata$est - qnorm(1 - 0.05/2) * pdata$se,
             ymax = pdata$est + qnorm(1 - 0.05/2) * pdata$se,
             color = pdata$group)) +
    geom_point() +
    geom_errorbar(aes(width = 0.1)) +
    #geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
    scale_color_manual(values = c('red', 'black')) +
    facet_grid(. ~ group, scales = 'free_x', space = 'free') +
    ylab('Estimates') + xlab('Sites') +
    theme_bw() +
    theme(legend.position = 'none',
          panel.grid = element_blank()) -> g

  return(g)
}
