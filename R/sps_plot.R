sps_plot <- function(out, columns = NULL, title = NULL, size = 2){

  N <- length(out$ss)
  Xp <- as.data.frame(out$X)

  if(is.null(columns) == TRUE){
    columns <- colnames(Xp)
  }else{
    Xp <- Xp[, columns]
  }

  col_use <- rep("Not", N)

  col_use[out$ss == 1] <- "Selected"
  Xp$summary_var <- col_use
  N_s <- sum(out$ss)
  N_r <- N - N_s
  colname_use <- colnames(Xp)
  colname_use[colname_use == "summary_var"] <- "Summary"

  p <- ggpairs(Xp, aes(color = Xp$summary_var),
               columnLabels = colname_use,
               upper = "blank",
               progress = FALSE,
               lower = list(combo = wrap("facethist", binwidth = 0.2),
                            continuous = wrap(ggally_points, size = 2)),
               axisLabels = c("show")) +
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

  add_text <- ggally_text(
    label = paste0("N = ",N,"\nSelected = ", N_s),
    color = "black",
    hjust = "left",
    xP = 0.1,
    fontface = "bold",
    size = size)

  # DIANA: Removing X-axis text and grid lines in the last panel (Summary X Summary)
  p[ncol(Xp), ncol(Xp)] <- add_text + theme(axis.text.x  = element_blank(),
                                            axis.ticks.x = element_blank(),
                                            panel.grid.major   = element_blank(),
                                            panel.grid.minor   = element_blank())

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
  p
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
