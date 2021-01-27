#' Graphical representation of the structure of disagreement among the panel of raters
#'
#' @description Draws the graphs representing the structure of disagreement among the set of ratings. Tunning parameters such as the color of the clusters or the axis of PCA to be plotted can be specified.
#' @usage plot_agreeclust(res, choice = "all", interact = FALSE, col_clust = NULL, axis = c(1, 2),
#' name_rater = "rater", ext_dev_Rstudio = FALSE, vignette = FALSE)

#' @param res An object of class agreeclust returned by the functions of the package.
#' @param choice A character element specifying the graphs to be plotted ('seg' for the representation of the segmentation process (dendrogram + partitioning if a consolidation process has been implemented), 'mul' for the multidimensional representation of the structure of disagreement (PCA), 'all' for both representations. By default, both representations are plotted.)
#' @param interact  A boolean specifying if the graphical outputs should be interactive (with plotly) or not.
#' @param col_clust A vector with as many color as clusters of raters. Colors can be specified by a character such as "blue" or by a hexadecimal code such as "#0D10D5". All hexadecimal codes can be found at \url{http://htmlcolorcodes.com}. By default, no colors are specified and default colors are used.
#' @param axis A length 2 numeric vector specifying the PCA components to be plotted. By default, the first 2 components are plotted.
#' @param name_rater A character element indicating how the raters should be named (e.g. "rater", "participant", "psychologists")
#' @param ext_dev_Rstudio A boolean specifying if the graphical outputs should be plotted in the Rstudio plot pane or not.
#' @param vignette A boolean specifying if the graphical outputs are plotted in a vignette or not.

#' @import ggplot2
#' @importFrom grid unit textGrob
#' @importFrom ggrepel geom_text_repel
#' @importFrom reshape2 melt
#' @importFrom doBy summaryBy
#' @importFrom manipulateWidget combineWidgets
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom shiny tags
#' @importFrom grDevices dev.cur dev.new dev.off

#' @return Returns the graphical representations
#' @export
#'
#' @examples
#' data(binary_data_for_example)
#' res_pedag <- get_agreeclust_bin(dta = binary_data_for_example,
#'                                 id_info_rater = 9 : nrow(binary_data_for_example),
#'                                 type_info_rater = c(rep("cat", 2), "cont"),
#'                                 id_info_stim = 21 : ncol(binary_data_for_example),
#'                                 type_info_stim = c(rep("cont", 4), "cat"),
#'                                 paral_null = FALSE
#'                                 )
#' plot_agreeclust(res_pedag, col_clust = c("blue", "green"))
#' plot_agreeclust(res_pedag, choice = "mul", interact = TRUE, col_clust = c("blue", "green"))
plot_agreeclust <- function(res, choice = "all", interact = FALSE, col_clust = NULL, axis = c(1, 2),
                            name_rater = "rater", ext_dev_Rstudio = FALSE, vignette = FALSE) {

  # check the format of the arguments
  if (!inherits(res, "agreeclust_binary") & !inherits(res, "agreeclust_continuous")) {
    stop("Non convenient data - res should be an agreeclust object")
  }
  choice <- match.arg(choice, c("all", "seg", "mul"))
  mat.partition <- cbind.data.frame(res$partition, names(res$partition))
  colnames(mat.partition) <- c("cluster", "rater")
  mat.partition$cluster <- as.factor(mat.partition$cluster)
  nb.clust <- nlevels(mat.partition$cluster)
  if (!is.null(col_clust)) {
    if (length(col_clust) < nb.clust) {
      stop("Non convenient specification of colors - col_clust should contain at least as many elements as clusters")
    }
  }
  if (length(axis) != 2 | length(unique(axis)) != 2 | class(axis) != "numeric") {
    stop("Non convenient specification of axis - axis should be a numeric vector of 2 different elements")
  }

  # default palette
  palette.col <- c("#90B08F", "#EA485C", "#FF8379", "#009193", "#FFCEA5", "#A9A9A9", "#B0983D", "#941751", "#333333", "#A8D9FF")

  # plot the segmentation
  if (choice == "all" | choice == "seg") {
    plot.dendro <- ggplot2::ggplot(NULL) +
      ggplot2::geom_segment(data = res$res_plot_segment$data_segments, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#444444") +
      ggplot2::geom_text(data = res$res_plot_segment$data_labels, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = cluster), size = 2.1) +
      ggplot2::ylim(-1, (max(res$res_plot_segment$data_segments$y) + (0.2 * max(res$res_plot_segment$dendrogram$height) / 10))) +
      ggplot2::theme(
        legend.key = ggplot2::element_rect(colour = "white", fill = "white"),
        legend.title = ggplot2::element_text(colour = "#444444"),
        legend.text = ggplot2::element_text(colour = "#444444"),
        panel.background = ggplot2::element_rect(fill = 'white', colour = "white"),
        panel.grid.major = ggplot2::element_line(colour = "white"),
        panel.grid.minor = ggplot2::element_line(colour = "white"),
        plot.title = ggplot2::element_text(hjust = 0.5, vjust = -1, size = 10, colour = "#444444"),
        plot.margin = grid::unit(c(0.5,0,0,0), "cm"),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "none")
    if (!is.null(col_clust)) {
      plot.dendro <- plot.dendro +
        ggplot2::scale_colour_manual(values = col_clust[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
    } else {
      plot.dendro <- plot.dendro +
        ggplot2::scale_colour_manual(values = palette.col[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
    }
    if (!is.null(res$res_plot_segment$data_labels_partitioning)) {
      plot.dendro <- plot.dendro +
        ggplot2::ggtitle("Before consolidation")
    }
    plot.dendro <- plot.dendro +
      ggplot2::geom_hline(data = res$res_plot_segment$res_test, ggplot2::aes(yintercept = y), colour = "#444444", linetype = 2) +
      ggplot2::geom_text(data = res$res_plot_segment$res_test, ggplot2::aes(x = 1, y = (y + (0.2 * max(res$res_plot_segment$dendrogram$height) / 10))), label = res$res_plot_segment$res_test[,"pval"], colour = "#444444", size = 2.5)
    if (!is.null(res$res_plot_segment$res_test_noLC)) {
      plot.dendro <- plot.dendro +
        ggplot2::ylim(-1, (max(res$res_plot_segment$data_segments$y) + (0.3 * max(res$res_plot_segment$dendrogram$height) / 10))) +
        ggplot2::geom_point(data = res$res_plot_segment$res_test_noLC, ggplot2::aes(x = x, y = y), colour = "#444444", size = 3, shape = 18) +
        ggplot2::geom_text(data = res$res_plot_segment$res_test_noLC, ggplot2::aes(x = x, y = (y + (0.3 * max(res$res_plot_segment$dendrogram$height) / 10))), label = res$res_plot_segment$res_test_noLC[, "pval"], colour = "#444444", size = 2.5) +
        ggplot2::guides(colour = guide_legend(override.aes = list(size=2.5)))
    }
    text.legend.test.height <- "p-value associated to the test of H0: this K-latent class structure is not significant (w.r.t. the (K-1)-latent class structure)"
    text.legend.test.noLC <- "p-value associated to the test of H0: the perfect agreement model well fits the data (displayed only if the number of clusters found equals 1)"
    plot.legend.dendro <- ggplot2::ggplot(NULL) +
      ggplot2::coord_fixed() +
      ggplot2::geom_point(data = res$res_plot_segment$coord_legend_dendro, ggplot2::aes(x = x, y = y), colour = "white", size = 2, shape = 18) +
      ggplot2::geom_segment(data = res$res_plot_segment$coord_test_height, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), colour = "#444444", linetype = 2) +
      ggplot2::geom_text(data = res$res_plot_segment$coord_test_height, ggplot2::aes(x = x, y = (y + 0.2)), label = "p-value", colour = "#444444", hjust = 0, size = 2) +
      ggplot2::geom_point(data = res$res_plot_segment$coord_test_noLC, ggplot2::aes(x = x, y = y), colour = "#444444", size = 2, shape = 18) +
      ggplot2::geom_text(data = res$res_plot_segment$coord_test_noLC, ggplot2::aes(x = x, y = (y - 0.2)), label = "p-value", colour = "#444444", size = 2) +
      ggplot2::geom_text(data = res$res_plot_segment$coord_legend_test_height, ggplot2::aes(x = x, y = y + 0.1), label = text.legend.test.height, hjust = 0, colour = "#444444", size = 2) +
      ggplot2::geom_text(data = res$res_plot_segment$coord_legend_test_noLC, ggplot2::aes(x = x, y = y - 0.1), label = text.legend.test.noLC, hjust = 0, colour = "#444444", size = 2) +
      ggplot2::theme(
        legend.key = ggplot2::element_rect(colour = "white", fill = "white"),
        legend.title = ggplot2::element_text(colour = "#444444"),
        legend.text = ggplot2::element_text(colour = "#444444"),
        panel.background = ggplot2::element_rect(fill = 'white', colour = "white"),
        panel.grid.major = ggplot2::element_line(colour = "white"),
        panel.grid.minor = ggplot2::element_line(colour = "white"),
        plot.margin = grid::unit(c(0.5, 0, 0, 0), "cm"),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank())
    if (!is.null(res$res_plot_segment$data_labels_partitioning)) {
      plot.partitioning <- ggplot2::ggplot(NULL) +
        ggplot2::geom_text(data = res$res_plot_segment$data_labels_partitioning, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = cluster), size = 2.1) +
        ggplot2::ylim(-0.3, 0) +
        ggplot2::ggtitle("After consolidation") +
        ggplot2::theme(
          legend.key = ggplot2::element_rect(colour = "white",fill = "white"),
          legend.title = ggplot2::element_text(colour = "#444444"),
          legend.text = ggplot2::element_text(colour = "#444444"),
          panel.background = ggplot2::element_rect(fill = 'white', colour = "white"),
          panel.grid.major = ggplot2::element_line(colour = "white"),
          panel.grid.minor = ggplot2::element_line(colour = "white"),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 10, colour = "#444444"),
          plot.margin = grid::unit(c(0.5, 0, 0, 0.45), "cm"),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          legend.position = "none")
      if (!is.null(col_clust)) {
        plot.partitioning <- plot.partitioning +
          ggplot2::scale_colour_manual(values = col_clust[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
      } else {
        plot.partitioning <- plot.partitioning +
          ggplot2::scale_colour_manual(values = palette.col[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
      }
    }
    plot.legend.clust <- ggplot2::ggplot(NULL) +
      ggplot2::geom_label(data = res$res_plot_segment$data_labels, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = cluster), colour = "transparent") +
      ggplot2::theme(
        plot.margin = grid::unit(c(0,0,0,0), "cm"),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size=8, colour = "#444444"),
        legend.text = ggplot2::element_text(size=8, colour = "#444444"),
        legend.margin = ggplot2::margin(t=0, unit='cm'),
        legend.key = ggplot2::element_rect(size=4),
        legend.key.size = grid::unit(0.4, "cm"))
    if (!is.null(col_clust)) {
      plot.legend.clust <- plot.legend.clust +
        ggplot2::scale_fill_manual(values = col_clust[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
    } else {
      plot.legend.clust <- plot.legend.clust +
        ggplot2::scale_fill_manual(values = palette.col[1 : nlevels(res$res_plot_segment$data_labels$cluster)])
    }
    empty.dev <- (dev.cur() == 1)
    legend.plot <- get.legend(plot.legend.clust)
    if (empty.dev == TRUE) {
      dev.off()
    }
    main.title <- grid::textGrob(paste(paste0(simpleCap(name_rater), "s"), "clustering"), gp = grid::gpar(fontsize = 12, font = 2, col = "#444444"))
    if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext_dev_Rstudio == TRUE) {
      dev.new(noRStudioGD = TRUE)
    }
    if (is.null(res$res_plot_segment$data_labels_partitioning)) {
      gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.dendro + ggplot2::theme(legend.position = "none"),
                               plot.legend.dendro + ggplot2::theme(legend.position = "none"),
                               ncol = 1, nrow = 2, heights = c(4, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    } else {
      gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + ggplot2::theme(legend.position = "none"),
                               plot.partitioning + ggplot2::theme(legend.position = "none"),
                               ncol = 1, nrow = 3, heights = c(4, 1, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  }

  # plot the multidimensional representation of disagreement
  if (choice == "all" | choice == "mul") {
    res.pca <- res$res_pca
    coord.raters <- res.pca$ind$coord[, axis]
    mat.coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
    colnames(mat.coord.raters) <- c("rater", "AxeA", "AxeB")
    mat.partition <- cbind.data.frame(names(res$partition), res$partition)
    colnames(mat.partition) <- c("rater", "cluster")
    coord.raters <- merge(mat.coord.raters, mat.partition, by = "rater")
    coord.raters$cluster <- as.factor(coord.raters$cluster)
    rownames(coord.raters) <- coord.raters[, "rater"]
    coord.raters <- coord.raters[, -1]
    amplitude.x <- max(coord.raters$AxeA) - min(coord.raters$AxeA)
    amplitude.y <- max(coord.raters$AxeB) - min(coord.raters$AxeB)
    if (amplitude.x == max(c(amplitude.x, amplitude.y))) {
      xlim <- c((min(coord.raters$AxeA)-0.1), (max(coord.raters$AxeA)+0.1))
      ylim <- c((-(abs(min(coord.raters$AxeB))/amplitude.y*amplitude.x)-0.1), (max(coord.raters$AxeB)/amplitude.y*amplitude.x+0.1))
    } else {
      xlim <- c((-(abs(min(coord.raters$AxeA))/amplitude.x*amplitude.y)-0.1), (max(coord.raters$AxeA)/amplitude.x*amplitude.y+0.1))
      ylim <- c((min(coord.raters$AxeB)-0.1), (max(coord.raters$AxeB)+0.1))
    }
    if (length(which(round(coord.raters$AxeB, 10) == rep(0, nrow(coord.raters)))) == nrow(coord.raters)) {
      ylim <- xlim
    }
    if (length(which(round(coord.raters$AxeA, 10) == rep(0, nrow(coord.raters)))) == nrow(coord.raters)) {
      xlim <- ylim
    }
    if (interact == FALSE) {
      plot.ind.pca <- ggplot2::ggplot(NULL) +
        ggplot2::labs(x = paste("Dim ", axis[1]," - ", round(res.pca$eig[axis[1], 2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %", sep = "")) +
        ggplot2::coord_fixed() +
        ggplot2::xlim(xlim[1], xlim[2]) +
        ggplot2::ylim(ylim[1], ylim[2]) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "black", size = 0.2) +
        ggplot2::geom_point(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, color = cluster)) +
        ggrepel::geom_text_repel(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, label = rownames(coord.raters), color = cluster), segment.color = "#444444", segment.size = 0.3, size = 2.3) +
        ggplot2::geom_point(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, color = cluster)) +
        ggplot2::ggtitle(paste("Representation of the", paste0(name_rater, "s"))) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 12, colour = "#444444"),
          plot.margin = grid::unit(c(1, 0.5, 0.5, 0.5), "cm"),
          panel.background = ggplot2::element_rect(fill = 'white', colour = "#444444"),
          panel.grid.major = ggplot2::element_line(colour = "white"),
          panel.grid.minor = ggplot2::element_line(colour = "white"),
          axis.text = ggplot2::element_text(colour = "#444444"),
          axis.ticks = ggplot2::element_line(colour = "#444444"),
          axis.title = ggplot2::element_text(colour = "#444444"),
          legend.position = "none")
      if (!is.null(col_clust)) {
        plot.ind.pca <- plot.ind.pca +
          ggplot2::scale_color_manual(values = col_clust[1 : nlevels(coord.raters$cluster)])
      } else {
        plot.ind.pca <- plot.ind.pca +
          ggplot2::scale_color_manual(values = palette.col[1 : nlevels(coord.raters$cluster)])
      }
    } else if (interact == TRUE) {
      coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
      colnames(coord.raters)[1] <- "rater"
      coord.raters$cluster <- paste("cluster", coord.raters$cluster)
      if (!is.null(col_clust)) {
        col <- col_clust[1 : nlevels(as.factor(coord.raters$cluster))]
      } else {
        col = palette.col[1 : nlevels(as.factor(coord.raters$cluster))]
      }
      text.tooltip <- paste(paste0(simpleCap(name_rater), ":"), coord.raters$rater, '<br>cluster:', gsub("cluster ", "", coord.raters$cluster))
      if (!is.null(res$call$id_info_rater)) {
        info.rater <- res$call$dta[res$call$id_info_rater, ]
        for (i in 1 : length(res$call$id_info_rater)) {
          info <- droplevels(info.rater[i, 1 : length(coord.raters$rater)][order(info.rater[i, 1 : length(coord.raters$rater)], coord.raters$rater)])
          text.tooltip <- paste(text.tooltip, paste0("<br>", rownames(info.rater)[i], ": ", as.character(unlist(info))))
        }
      }
      plot.ind.pca.interact <- plotly::plot_ly(coord.raters) %>%
        plotly::add_trace(coord.raters,
                  x = ~AxeA ,
                  y = ~AxeB,
                  color = ~cluster,
                  hoverinfo = 'text',
                  text = text.tooltip,
                  type = "scatter", mode = "markers",
                  marker = list(size = 6),
                  showlegend = TRUE,
                  colors = col) %>%
        plotly::layout(
          legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.25),
          margin = list(l = 50, r = 50, t = 50, b = 50),
          title = list(text = paste("Representation of the", paste0(name_rater, "s")), font = list(size = 14, color = "#444444")),
          xaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "y", showgrid = FALSE, title = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1),
          yaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "x", showgrid = FALSE, title = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1)
        )
    }
    coord.stimuli <- as.data.frame(res.pca$var$coord[, axis])
    colnames(coord.stimuli) <- c("AxeA", "AxeB")
    amplitude.x <- max(coord.stimuli$AxeA) - min(coord.stimuli$AxeA)
    amplitude.y <- max(coord.stimuli$AxeB) - min(coord.stimuli$AxeB)
    if (amplitude.x == max(c(amplitude.x, amplitude.y))) {
      xlim <- c((min(coord.stimuli$AxeA)-0.1), (max(coord.stimuli$AxeA)+0.1))
      ylim <- c((-(abs(min(coord.stimuli$AxeB))/amplitude.y*amplitude.x)-0.1), (max(coord.stimuli$AxeB)/amplitude.y*amplitude.x+0.1))
    } else {
      xlim <- c((-(abs(min(coord.stimuli$AxeA))/amplitude.x*amplitude.y)-0.1), (max(coord.stimuli$AxeA)/amplitude.x*amplitude.y+0.1))
      ylim <- c((min(coord.stimuli$AxeB)-0.1), (max(coord.stimuli$AxeB)+0.1))
    }
    if (length(which(round(coord.stimuli$AxeB, 10) == rep(0, nrow(coord.stimuli)))) == nrow(coord.stimuli)) {
      ylim <- xlim
    }
    if (length(which(round(coord.stimuli$AxeA, 10) == rep(0, nrow(coord.stimuli)))) == nrow(coord.stimuli)) {
      xlim <- ylim
    }
    if (interact == FALSE) {
      plot.var.pca <- ggplot2::ggplot(NULL) +
        ggplot2::labs(x = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %",sep = "")) +
        ggplot2::coord_fixed() +
        ggplot2::xlim(xlim[1], xlim[2]) +
        ggplot2::ylim(ylim[1], ylim[2]) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
        ggplot2::geom_vline(xintercept = 0,  linetype = 2, color = "black", size = 0.2) +
        ggplot2::geom_segment(data = coord.stimuli, ggplot2::aes(x = 0, y = 0, xend = AxeA, yend = AxeB), alpha = 1, color = "black", size = 0.3, arrow = arrow(length = unit(0.3, "cm"))) +
        ggrepel::geom_text_repel(data = coord.stimuli, ggplot2::aes(x = AxeA, y = AxeB, label = rownames(coord.stimuli)), segment.color = "transparent", segment.size = 0.3, size = 2.3) +
        ggplot2::ggtitle("Representation of the stimuli") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 12, colour = "#444444"),
          plot.margin = grid::unit(c(1, 0.5, 0.5, 0.5), "cm"),
          panel.background = ggplot2::element_rect(fill = 'white', colour = "#444444"),
          panel.grid.major = ggplot2::element_line(colour = "white"),
          panel.grid.minor = ggplot2::element_line(colour = "white"),
          axis.text = ggplot2::element_text(colour = "#444444"),
          axis.ticks = ggplot2::element_line(colour = "#444444"),
          axis.title = ggplot2::element_text(colour = "#444444"),
          legend.position = "none")
      plot.legend.clust <- ggplot2::ggplot(NULL) +
        geom_label(data = res$res_plot_segment$data_labels, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = cluster), colour = "transparent") +
        ggplot2::theme(
          plot.margin = grid::unit(c(0,0,0,0), "cm"),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(size=8),
          legend.text = ggplot2::element_text(size=8),
          legend.margin = ggplot2::margin(t=0, unit='cm'),
          legend.key = ggplot2::element_rect(size=4),
          legend.key.size = grid::unit(0.4, "cm"))
      if (!is.null(col_clust)) {
        plot.legend.clust <- plot.legend.clust +
          ggplot2::scale_fill_manual(values = col_clust[1 : nlevels(coord.raters$cluster)])
      } else {
        plot.legend.clust <- plot.legend.clust +
          ggplot2::scale_fill_manual(values = palette.col[1 : nlevels(coord.raters$cluster)])
      }
      empty.dev <- (dev.cur() == 1)
      legend.plot <- get.legend(plot.legend.clust)
      if (empty.dev == TRUE) {
        dev.off()
      }
      main.title <- grid::textGrob(paste("Multidimensional representation of the structure \n of disagreement among the panel of", paste0(name_rater, "s")), gp = gpar(fontsize = 12, font = 2, col = "#444444"))
      if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext_dev_Rstudio == TRUE) {
        dev.new(noRStudioGD = TRUE)
      }
      gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.ind.pca + ggplot2::theme(legend.position="none"),
                             plot.var.pca + ggplot2::theme(legend.position="none"),
                             ncol = 2, nrow = 1),
                 legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    } else if (interact == TRUE) {
      coord.stimuli <- cbind.data.frame(rownames(coord.stimuli), coord.stimuli, paste("cluster", c(1 : nlevels(as.factor(coord.raters$cluster)), sample(1 : nlevels(as.factor(coord.raters$cluster)), (nrow(coord.stimuli) - nlevels(as.factor(coord.raters$cluster))), replace = TRUE))))
      colnames(coord.stimuli)[c(1, 4)] <- c("stimulus", "FakeGroup")
      coord.stimuli$FakeGroup <- as.factor(coord.stimuli$FakeGroup)
      text.tooltip <- paste("stimulus:", coord.stimuli$stimulus)
      if (!is.null(res$call$id_info_rater)) {
        dta <- res$call$dta[-res$call$id_info_rater,]
        dta <- droplevels(dta)
      }
      if (!is.null(res$call$id_info_stim)) {
        dta <- dta[, -res$call$id_info_stim]
        dta <- droplevels(dta)
      }
      melted.data <- reshape2::melt(as.matrix(dta))
      colnames(melted.data) <- c("stimulus", "rater", "rating")
      mat.partition <- cbind.data.frame(names(res$partition), res$partition)
      colnames(mat.partition) <- c("rater", "cluster")
      melted.cluster <- merge(melted.data, mat.partition, by = "rater")
      melted.cluster$rating <- as.numeric(as.character(melted.cluster$rating))
      if (inherits(res, "agreeclust_binary")) {
        sum.pos.ratings <- doBy::summaryBy(rating ~ stimulus : cluster, data = melted.cluster, FUN = "sum")
        colnames(sum.pos.ratings) <- c("stimulus", "cluster", "NbPosrating")
        size.clust <- summary(as.factor(mat.partition$cluster))
        for (i in 1 : nlevels(as.factor(mat.partition$cluster))) {
          tab.clust <- sum.pos.ratings[which(sum.pos.ratings$cluster == levels(as.factor(mat.partition$cluster))[i]), ]
          tab.clust[, "NbPosrating"] <- tab.clust[, "NbPosrating"] / size.clust[which(names(size.clust) == levels(as.factor(mat.partition$cluster))[i])] * 100
          text.tooltip <- paste(text.tooltip,
                                paste(paste0("<br>Positive ratings in cluster ", levels(as.factor(mat.partition$cluster))[i], ":"),
                                      paste0(round(tab.clust[order(tab.clust$stimulus, coord.stimuli$stimulus), "NbPosrating"], 1), "%")))
        }
      }
      if (inherits(res, "agreeclust_continuous")) {
        sum.ratings <- doBy::summaryBy(rating ~ stimulus : cluster, data = melted.cluster, FUN = "mean")
        colnames(sum.ratings) <- c("stimulus", "cluster", "Meanrating")
        for (i in 1 : nlevels(as.factor(mat.partition$cluster))) {
          tab.clust <- sum.ratings[which(sum.ratings$cluster == levels(as.factor(mat.partition$cluster))[i]), ]
          text.tooltip <- paste(text.tooltip,
                                paste(paste0("<br>Average rating in cluster ", levels(as.factor(mat.partition$cluster))[i], ":"),
                                      paste0(round(tab.clust[order(tab.clust$stimulus, coord.stimuli$stimulus), "Meanrating"], 1), "")))
        }
      }
      if (!is.null(res$call$id_info_stim)) {
        info.stim <- res$call$dta[, res$call$id_info_stim]
        for (i in 1 : length(res$call$id_info_stim)) {
          info <- info.stim[1 : length(coord.stimuli$stimulus), i][order(info.stim[1 : length(coord.stimuli$stimulus), i], coord.stimuli$stimulus)]
          text.tooltip <- paste(text.tooltip, paste0("<br>", colnames(info.stim)[i], ": ", as.character(unlist(info))))
        }
      }
      plot.var.pca.interact <- plotly::plot_ly(coord.stimuli) %>%
        plotly::add_trace(coord.stimuli,
                  x = ~AxeA ,
                  y = ~AxeB,
                  color = ~FakeGroup,
                  hoverinfo = 'text',
                  hoverlabel = list(bgcolor = "black", font = list(color = "white")),
                  text = text.tooltip,
                  type = "scatter", mode = "markers",
                  marker = list(size = 6, color = "white"),
                  showlegend = TRUE) %>%
        plotly::add_annotations(axref = "x", ax = coord.stimuli$AxeA, aref = "x", x = 0,
                  ayref = "y", ay = coord.stimuli$AxeB, aref = "y", y = 0,
                  text = coord.stimuli$stimulus, arrowhead = 0, arrowcolor = "black") %>%
        plotly::layout(
          legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.25, font = list(color = "transparent")),
          margin = list(l = 50, r = 50, t = 50, b = 50),
          title = list(text = "Representation of the stimuli", font = list(size = 14, color = "#444444")),
          xaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "y", showgrid = FALSE, title = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1),
          yaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "x", showgrid = FALSE, title = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1)
        )
      combine.plot <- manipulateWidget::combineWidgets(plot.ind.pca.interact, plot.var.pca.interact, nrow = 1, ncol = 2,
                                                 title = paste("Multidimensional representation of the structure <br>of disagreement among the panel of", paste0(name_rater, "s")),
                                                 titleCSS = "font-size:16px",
                                                 header = shiny::tags$div(shiny::tags$span("header", style = "color:white;font-size:20px")),
                                                 leftCol = shiny::tags$div(shiny::tags$span("header", style = "color:white;font-size:10px")),
                                                 rightCol = shiny::tags$div(shiny::tags$span("header", style = "color:white;font-size:10px")))
      if (vignette == FALSE) {
        if (ext_dev_Rstudio == TRUE) {
          old.viewer <- options()$viewer
          options(viewer = NULL)
          print(combine.plot)
          options(viewer = old.viewer)
        } else {
          print(combine.plot)
        }
      }
    }
  }

  # end the function
  if (vignette == TRUE & interact == TRUE & (choice == "all" | choice == "mul")) {
    message("Representations plotted")
    return(combine.plot)
  } else {
    return(message("Representations plotted"))
  }

}
