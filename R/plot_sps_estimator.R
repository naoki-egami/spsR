#' Plot function for \code{sps_estimator}
#' @param x Output from function \code{sps_estimator()}
#' @param ... Other arguments.
#' @return \code{sps_plot_est} plots estimates and 95% confidence intervals for the site-specific ATEs and average-site ATE.
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export

plot.sps_estimator <- function(x, ...){

  if(is.null(x$subgroup) == FALSE){
    subgroup <- TRUE
  }else{
    subgroup <- FALSE
  }

  if(subgroup == FALSE){
    estimates_selected <- x$estimates_selected

    est <- c(estimates_selected[,1], x$average_site_ATE[1])
    se  <- c(estimates_selected[,2], x$average_site_ATE[2])

    group <- c(rep('Site-Specific ATEs', nrow(estimates_selected)), 'Average-Site ATE')
    site  <- c(row.names(estimates_selected), 'Average-Site')

    pdata <- data.frame(est, se, group, site)
    pdata$group <- factor(pdata$group, levels = c('Site-Specific ATEs', 'Average-Site ATE'))
    pdata$site  <- factor(pdata$site, levels = pdata[order(pdata$est, decreasing = FALSE), 'site'])

    g <- ggplot(data = pdata,
                aes(x = pdata$site,
                    y = pdata$est,
                    ymin = pdata$est - qnorm(1 - 0.05/2) * pdata$se,
                    ymax = pdata$est + qnorm(1 - 0.05/2) * pdata$se,
                    color = pdata$group)) +
      geom_point(size = 4) +
      geom_errorbar(aes(width = 0.1), linewidth = 1.2) +
      #geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
      scale_color_manual(values = c('black', 'blue')) +
      facet_grid(. ~ group, scales = 'free_x', space = 'free') +
      ylab('Estimates') + xlab('Sites') +
      theme_bw() +
      theme(legend.position = 'none',
            panel.grid = element_blank())

    suppressWarnings(print(g))
  }else{
    estimates_selected <- x$estimates_selected

    est <- c(estimates_selected[,1], x$subgroup_average_site_ATE[, 1])
    se  <- c(estimates_selected[,2], x$subgroup_average_site_ATE[, 2])

    group <- c(rep('Site-Specific ATEs', nrow(estimates_selected)), rep('Subgroup Average-Site ATE', nrow(x$subgroup_average_site_ATE)))
    site  <- c(row.names(estimates_selected), rownames(x$subgroup_average_site_ATE))

    pdata <- data.frame(est, se, group, site)
    pdata$group <- factor(pdata$group, levels = c('Site-Specific ATEs', 'Subgroup Average-Site ATE'))
    pdata_specific <- pdata[1:nrow(estimates_selected), ]
    specific_order <- pdata_specific[order(pdata_specific$est, decreasing = FALSE), 'site']
    order_site <- c(specific_order, rownames(x$subgroup_average_site_ATE))
    pdata$site  <- factor(pdata$site, levels = order_site)

    g <- ggplot(data = pdata,
                aes(x = pdata$site,
                    y = pdata$est,
                    ymin = pdata$est - qnorm(1 - 0.05/2) * pdata$se,
                    ymax = pdata$est + qnorm(1 - 0.05/2) * pdata$se,
                    color = pdata$group)) +
      geom_point(size = 4) +
      geom_errorbar(aes(width = 0.1), linewidth = 1.2) +
      #geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
      scale_color_manual(values = c('black', 'blue')) +
      facet_grid(. ~ group, scales = 'free_x', space = 'free') +
      ylab('Estimates') + xlab('Sites') +
      theme_bw() +
      theme(legend.position = 'none',
            panel.grid = element_blank())
  }
  suppressWarnings(print(g))
}
