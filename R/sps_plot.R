#' Plot function for \code{sps}
#' @param out Output from function \code{sps()}. When \code{before_selection = TRUE}, set this argument to be \code{NULL}.
#' @param title Title of the plot
#' @param columns (Optional. Default = \code{NULL}) Names of columns users want to visualize. When \code{NULL}, the function plots every column.
#' @param before_selection Logical (\code{TRUE} or \code{FALSE}. Default = \code{FALSE}). When \code{FALSE}, the function compares selected and non-selected sites. When \code{TRUE}, the function simply visualizes the distribution of site-level variables in the target population of sites.
#' @param X (Optional. Use this only when \code{before_selection = TRUE}.) Site-level variables for the target population of sites.
#' @param histogram Logical (\code{TRUE} or \code{FALSE}. Default = \code{FALSE}). When \code{NULL}, the function plots the histograms of site-level variables separately for target population and selected sites. Only available when \code{before_selection = FALSE}.
#' @importFrom ggfittext geom_fit_text
#' @return \code{sps_plot} visualizes the distribution of site-level variables.
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
sps_plot <- function(out, title = NULL, columns = NULL, before_selection = FALSE, X = NULL, histogram = FALSE){

  if(before_selection == TRUE){

    if(is.null(X) == TRUE){
      stop(" When `before_selection = TRUE`, `X` should be specified. ")
    }

    if(before_selection == TRUE & histogram == TRUE){
      stop(" The histogram option is only available when `before_selection = FALSE`. ")
    }

    if(is.null(columns) == TRUE){
      columns <- colnames(X)
      X_use <- X
    }else{
      X_use <- X[, columns]
    }

    p <- sps_plot_base(X = X_use, title = title)

    suppressWarnings(print(p))

  }else{

    if(("sps" %in% class(out)) == FALSE){
      stop(" `out` should be an output from function `sps()` ")
    }

    N <- length(out$internal$ss)
    Xp <- as.data.frame(out$internal$X)

    if(is.null(columns) == TRUE){
      columns <- colnames(Xp)
    }else{
      Xp <- Xp[, columns]
    }

    if(histogram == TRUE){

      col_use <- rep("Not Selected", N)
      col_use[out$internal$ss == 1] <- "SPS"
      Xp$site <- rownames(Xp)
      Xp$summary_var <- col_use

      Xp_all <- Xp
      Xp_all$summary_var <- 'Target Population'
      Xp <- rbind(Xp_all, Xp[Xp$summary_var=='SPS',])

      Xp_long <- reshape(Xp,
                         varying = list(names(Xp)[1:(ncol(Xp)-2)]),
                         v.names = "value",
                         times = names(Xp)[1:(ncol(Xp)-2)],
                         timevar = "name",
                         direction = "long")
      Xp_long$summary_var <- factor(Xp_long$summary_var, levels = c('Target Population', 'SPS'))
      p <- ggplot(Xp_long, aes(x = value, fill = summary_var)) +
        geom_histogram(alpha=1, position="identity") +
        facet_grid(summary_var ~ name, scale = 'free') +
        scale_fill_manual(name = '', values = c('black', 'red')) +
        facetted_pos_scales(y = list(version == 'Original'  ~ scale_y_continuous(limits = c(0,2)),
                                     version == 'SPS' ~ scale_y_continuous(limits = c(0,2)),
                                     version == 'Target Pop' ~ scale_y_continuous(limits = c(0,22)))) +
        xlab('') + ylab(NULL) + theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        ggtitle(label = paste0(title))

      suppressWarnings(print(p))

    }else{

      col_use <- rep("Not", N)
      col_use[out$internal$ss == 1] <- "Selected"
      Xp$summary_var <- col_use
      N_s <- sum(out$internal$ss)
      N_r <- N - N_s
      colname_use <- colnames(Xp)
      colname_use[colname_use == "summary_var"] <- "Summary"

      p <- ggpairs(Xp, aes(color = Xp$summary_var, shape = Xp$summary_var),
                   columnLabels = NULL,   # NAOKI: Big Change
                   upper = "blank",
                   progress = FALSE,
                   lower = list(combo = wrap("facethist", binwidth = 0.2),
                                continuous = wrap(ggally_points, size = 2)),
                   axisLabels = c("show")) +  # Naoki
        scale_shape_manual(values = c(16, 17)) +
        # DIANA: theme_bw() prevents panel-specific axis change, so all elements of theme_bw() are manually inserted below.
        theme(panel.background = element_rect(fill = "white", colour = NA),
              panel.border = element_rect(fill = NA, colour = "grey20"),
              strip.background = element_rect(fill = "grey85", colour = "grey20"),
              legend.key = element_rect(fill = "white", colour = NA)) +
        ggtitle(label = paste0(title))

      # Change color manually.
      # Loop through each plot changing relevant scales
      for(i in 1:p$nrow) {
        for(j in 1:p$ncol){
          p[i,j] <- p[i,j] +
            scale_fill_manual(values = c(adjustcolor("black", 1),
                                         adjustcolor("red", 1))) +
            scale_color_manual(values = c("black", "red"))
        }
      }

      # add_text <- ggally_text(
      #   label = paste0("N = ",N,"\nSelected = ", N_s),
      #   color = "black",
      #   hjust = "left",
      #   xP = 0.1,
      #   fontface = "bold",
      #   size = size)
      #
      # # DIANA: Removing X-axis text and grid lines in the last panel (Summary X Summary)
      # p[ncol(Xp), ncol(Xp)] <- add_text + theme(axis.text.x  = element_blank(),
      #                                           axis.ticks.x = element_blank(),
      #                                           panel.grid.major   = element_blank(),
      #                                           panel.grid.minor   = element_blank())

      ## NEW PART
      d_t <- data.frame(x = 0.5, y = 0.5, colname = paste0("N = ",N,"\nSelect = ", N_s))
      p[ncol(Xp), ncol(Xp)] <- ggplot(data = d_t, aes(x, y, label = colname)) +
        geom_tile(fill = 'white') +
        ggfittext::geom_fit_text(reflow = F) +
        theme(axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid   = element_blank(),
              panel.background = element_rect(fill = 'white'))

      # DIANA: Removing Y-axis text in panels in the last row
      for (i in 1:ncol(Xp)) {
        p[ncol(Xp),i] <- p[ncol(Xp),i] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      }


      # DIANA: Adding grid lines to all panels except for the last one
      # Diana Update: removing gridlines
      # for (i in 1:ncol(Xp)) {
      #   for (j in 1:(ncol(Xp)-1)) {
      #     if(i != j){
      #       # p[i, j] <- p[i, j] + theme(panel.grid = element_line(colour = "grey92"),
      #       #                            panel.grid.minor = element_line(linewidth = rel(0.5)))
      #     }
      #   }
      # }

      ## Naoki Big Change: Adding Variable Name
      ncl <- ncol(Xp) - 1
      d_t_l <- list()
      for(i in 1:ncl){
        d_t_l[[i]] <- data.frame(x = 0.5, y = 0.5, colname = colnames(Xp)[i])
        p[i, i] <- ggplot(data = d_t_l[[i]], aes(x, y, label= colname)) +
          # geom_text(x = 0.5, y = 0.5, label= colname) +
          ggfittext::geom_fit_text(reflow = FALSE) +
          theme(axis.text.x  = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid   = element_blank())
        # geom_text(x = 0.5, y = 0.5, label= colnames(Xp)[i]) +
        # geom_tile(fill = 'white') +
        # ,panel.background = element_rect(fill = 'white')
      }

      # NAOKI: Removing Y-axis text in panels in the first row
      for (i in 1:ncol(Xp)) {
        p[1,i] <- p[1,i] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      }
      suppressWarnings(print(p))
    }
  }
}


sps_plot_base <- function(X, title = NULL){

  Xp <- as.data.frame(X)
  L <- ncol(Xp)
  Xp$selection <- "1"

  p <- ggpairs(Xp[, 1:L], aes(color = Xp$selection),
               columnLabels = NULL, # Diana Update
               upper = "blank",
               progress = FALSE,
               lower = list(combo = wrap("facethist", binwidth = 0.2))) +
    scale_shape_manual(values=c(16, 17)) +
    # DIANA: theme_bw() prevents panel-specific axis change, so all elements of theme_bw() are manually inserted below.
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          legend.key = element_rect(fill = "white", colour = NA)) +
    ggtitle(label = paste0(title))

  # theme_bw() +
  # ggtitle(label = paste0(title))

  # Change color manually.
  # Loop through each plot changing relevant scales
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] +
        scale_fill_manual(values = c(adjustcolor("black", 0.3))) +
        scale_color_manual(values = c("black"))
    }
  }

  for(el in 1:p$nrow){
    d_t <- data.frame(x = 0.5, y = 0.5)
    p[el, el] <- ggplot(data = d_t, aes(d_t$x, d_t$y, label = colnames(Xp)[el])) +
      geom_tile(fill = 'white') +
      ggfittext::geom_fit_text(reflow = F, grow = T) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid   = element_blank(),
            panel.background = element_rect(fill = 'white'))
  }

  ## Diana Update: Adding Variable Name
  ncl <- ncol(Xp) - 1
  d_t_l <- list()
  for(i in 1:ncl){
    d_t_l[[i]] <- data.frame(x = 0.5, y = 0.5, colname = colnames(Xp)[i])
    p[i, i] <- ggplot(data = d_t_l[[i]], aes(x, y, label= colname)) +
      ggfittext::geom_fit_text(reflow = FALSE) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major   = element_blank(),
            panel.grid.minor   = element_blank())
  }
  # suppressWarnings(print(p))
  return(p)
}
