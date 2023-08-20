#' Plot function for \code{sps}
#' @param out Output from function \code{sps()}
#' @param title Title of the plot
#' @param columns (Optional. Default = \code{NULL}) Names of columns users want to visualize. When \code{NULL}, the function plots every column.
#' @param size (Default = \code{2}) Size of texts on the summary panel.
#' @param before_selection Logical (\code{TRUE} or \code{FALSE}. Default = \code{FALSE}). When \code{FALSE}, the function will compare selected and non-selected sites. When \code{TRUE}, the function simply visualizes the distribution of site-level variables in the target population of sites.
#' @param X (Optional. Use this only when \code{before_selection = TRUE}.) Site-level variables for the target population of sites.
#' @importFrom ggfittext geom_fit_text
#' @return \code{sps_plot} visualizes the distribution of site-level variables.
#' @references Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
#' @export
sps_plot <- function(out, title = NULL, columns = NULL, size = 2, before_selection = FALSE, X = NULL){

  if(before_selection == TRUE){

    if(is.null(columns) == TRUE){
      columns <- colnames(X)
    }else{
      X_use <- X[, columns]
    }


    p <- sps_plot_base(X = X_use, title = title)

    suppressWarnings(print(p))

  }else{
    N <- length(out$internal$ss)
    Xp <- as.data.frame(out$internal$X)

    if(is.null(columns) == TRUE){
      columns <- colnames(Xp)
    }else{
      Xp <- Xp[, columns]
    }

    col_use <- rep("Not", N)

    col_use[out$internal$ss == 1] <- "Selected"
    Xp$summary_var <- col_use
    N_s <- sum(out$internal$ss)
    N_r <- N - N_s
    colname_use <- colnames(Xp)
    colname_use[colname_use == "summary_var"] <- "Summary"

    p <- ggpairs(Xp, aes(color = Xp$summary_var, shape = Xp$summary_var),
                 columnLabels = colname_use,
                 upper = "blank",
                 progress = FALSE,
                 lower = list(combo = wrap("facethist", binwidth = 0.2),
                              continuous = wrap(ggally_points, size = 2)),
                 axisLabels = c("show")) +
      scale_shape_manual(values=c(16, 17)) +
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
          scale_fill_manual(values = c(adjustcolor("black", 0.3),
                                       adjustcolor("red", 0.3))) +
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
    d_t <- data.frame(x = 0.5, y = 0.5)
    p[ncol(Xp), ncol(Xp)] <- ggplot(data = d_t, aes(d_t$x, d_t$y, label = paste0("N = ",N,"\nSelected = ", N_s))) +
      geom_tile(fill = 'white') +
      ggfittext::geom_fit_text(reflow = F, grow = T) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major   = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.background = element_rect(fill = 'white'))

    # DIANA: Removing Y-axis text in panels in the last row
    for (i in 1:ncol(Xp)) {
      p[ncol(Xp),i] <- p[ncol(Xp),i] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }

    # DIANA: Adding grid lines to all panels except for the last one
    for (i in 1:ncol(Xp)) {
      for (j in 1:(ncol(Xp)-1)) {
        p[i, j] <- p[i, j] + theme(panel.grid = element_line(colour = "grey92"),
                                   panel.grid.minor = element_line(linewidth = rel(0.5)))
      }
    }
    suppressWarnings(print(p))
  }
}

# sps_plot_custom <- function(X, selected, columns = NULL, title = NULL, size = 2){
#
#   N <- nrow(X)
#   Xp <- as.data.frame(X)
#
#   if(is.null(columns) == TRUE){
#     columns <- colnames(Xp)
#   }else{
#     Xp <- Xp[, columns]
#   }
#
#   col_use <- rep("Not", N)
#   col_use[selected == 1] <- "Selected"
#
#   Xp$Summary <- col_use
#   # levels(Xp$selection)
#   N_s <- sum(out$ss)
#   N_r <- N - N_s
#
#   p <- ggpairs(Xp, aes(color = Summary),
#                upper = "blank",
#                progress = FALSE,
#                lower = list(combo = wrap("facethist", binwidth = 0.2),
#                             continuous = wrap(ggally_points, size = 2))) +
#     # theme_bw() +
#     # DIANA: theme_bw() prevents panel-specific axis change, so all elements of theme_bw() are manually inserted below.
#     theme(panel.background = element_rect(fill = "white", colour = NA),
#           panel.border = element_rect(fill = NA, colour = "grey20"),
#           strip.background = element_rect(fill = "grey85", colour = "grey20"),
#           legend.key = element_rect(fill = "white", colour = NA)) +
#     ggtitle(label = paste0(title))
#   # Change color manually.
#   # Loop through each plot changing relevant scales
#   for(i in 1:p$nrow) {
#     for(j in 1:p$ncol){
#       p[i,j] <- p[i,j] +
#         scale_fill_manual(values = c(adjustcolor("black", 0.3),
#                                      adjustcolor("red", 0.3))) +
#         scale_color_manual(values = c("black", "red"))
#     }
#   }
#
#   add_text <- ggally_text(
#     label = paste0("N = ",N,"\nSelected = ", N_s),
#     color = "black",
#     hjust = "left",
#     xP = 0.1,
#     fontface = "bold",
#     size = size)
#
#   # DIANA: Removing X-axis text and grid lines in the last panel (Summary X Summary)
#   p[ncol(Xp), ncol(Xp)] <- add_text + theme(axis.text.x  = element_blank(),
#                                             axis.ticks.x = element_blank(),
#                                             panel.grid.major   = element_blank(),
#                                             panel.grid.minor   = element_blank())
#
#   # DIANA: Removing Y-axis text in panels in the last row
#   for (i in 1:ncol(Xp)) {
#     p[ncol(Xp),i] <- p[ncol(Xp),i] + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
#   }
#
#   # DIANA: Adding grid lines to all panels except for the last one
#   for (i in 1:ncol(Xp)) {
#     for (j in 1:(ncol(Xp)-1)) {
#       p[i, j] <- p[i, j] + theme(panel.grid = element_line(colour = "grey92"),
#                                  panel.grid.minor = element_line(linewidth = rel(0.5)))
#     }
#   }
#
#   p
# }

sps_plot_base <- function(X, title = NULL){

  Xp <- as.data.frame(X)
  L <- ncol(Xp)
  Xp$selection <- "1"

  p <- ggpairs(Xp[, 1:L], aes(color = Xp$selection),
               upper = "blank",
               progress = FALSE,
               lower = list(combo = wrap("facethist", binwidth = 0.2))) + theme_bw() +
    ggtitle(label = paste0(title))
  # Change color manually.
  # Loop through each plot changing relevant scales
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] +
        scale_fill_manual(values = c(adjustcolor("black", 0.3))) +
        scale_color_manual(values = c("black"))
    }
  }
  p
}
