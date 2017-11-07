plot.AgreeClust <- function(res, choice = "all", col.clust = NULL, axis = c(1, 2)) {

  options(warn = -1)

  # load packages
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  suppressPackageStartupMessages(require(ggrepel, quietly = TRUE))

  # check the format of the arguments
  if (!inherits(res, "AgreeClust")) {
    stop("Non convenient data - res should be an AgreeClust object")
  }
  choice <- match.arg(choice, c("all", "seg", "mul"))
  mat.partition <- cbind.data.frame(res[[5]], names(res[[5]]))
  colnames(mat.partition) <- c("Cluster", "Rater")
  mat.partition$Cluster <- as.factor(mat.partition$Cluster)
  nb.clust <- nlevels(mat.partition$Cluster)
  if (!is.null(col.clust)) {
    if (length(col.clust) != nb.clust) {
      stop("Non convenient specification of colors - col.clust should contain as many elements as clusters")
    }
  }
  if (length(axis) != 2 | length(unique(axis)) != 2 | class(axis) != "numeric") {
    stop("Non convenient specification of axis - axis should be a numeric vector of 2 different elements")
  }

  # function to extract legend
  get.legend <- function(plot){
    tmp <- ggplot_gtable(ggplot_build(plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  # plot the segmentation
  if (choice == "all" | choice == "seg") {
    if (length(res[[6]]) == 4) { # consolidation performed
      plot.dendro <- res[[6]][[1]]
      plot.legend.dendro <- res[[6]][[2]]
      plot.partitioning <- res[[6]][[3]]
      plot.legend.clust <- res[[6]][[4]]
    } else if (length(res[[6]]) == 3) { # consolidation not performed
      plot.dendro <- res[[6]][[1]]
      plot.legend.dendro <- res[[6]][[2]]
      plot.legend.clust <- res[[6]][[3]]
    }
    if (!is.null(col.clust)) {
      plot.dendro <- plot.dendro +
        scale_colour_manual(values = col.clust)
      plot.legend.clust <- plot.legend.clust +
        scale_fill_manual(values = col.clust)
      if (length(res[[6]]) == 4) {
        plot.partitioning <- plot.partitioning +
          scale_colour_manual(values = col.clust)
      }
    }
    legend.plot <- get.legend(plot.legend.clust)
    main.title <- textGrob("Raters clustering", gp = gpar(fontsize = 12, font = 2))
    if (length(res[[6]]) == 4) {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               plot.partitioning + theme(legend.position = "none"),
                               ncol = 1, nrow = 3, heights = c(4, 1, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    } else if (length(res[[6]]) == 3) {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               ncol = 1, nrow = 2, heights = c(4, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  }

  # plot the multidimensional representation of disagreement
  if (choice == "all" | choice == "mul") {
    res.pca <- res[[7]]
    coord.raters <- res.pca$ind$coord[, axis]
    mat.coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
    colnames(mat.coord.raters) <- c("Rater", "AxeA", "AxeB")
    coord.raters <- merge(mat.coord.raters, mat.partition, by = "Rater")
    coord.raters$Cluster <- as.factor(coord.raters$Cluster)
    rownames(coord.raters) <- coord.raters[, "Rater"]
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
    plot.ind.pca <- ggplot(NULL) +
      labs(x = paste("Dim ", axis[1]," - ", round(res.pca$eig[axis[1], 2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %", sep = "")) +
      coord_fixed() +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2]) +
      geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
      geom_vline(xintercept = 0, linetype = 2, color = "black", size = 0.2) +
      geom_point(data = coord.raters, aes(x = AxeA, y = AxeB, color = Cluster)) +
      geom_text_repel(data = coord.raters, aes(x = AxeA, y = AxeB, label = rownames(coord.raters), color = Cluster), segment.color = "grey", segment.size = 0.3, size = 2.3) +
      geom_point(data = coord.raters, aes(x = AxeA, y = AxeB, color = Cluster)) +
      ggtitle("Representation of the raters") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = 'white', colour = "black"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        legend.position = "none")
    if (!is.null(col.clust)) {
      plot.ind.pca <- plot.ind.pca +
        scale_color_manual(values = col.clust)
    } else {
      palette.col <- c("#90B08F", "#EA485C", "#FF8379", "#009193", "#FFCEA5", "#A9A9A9", "#B0983D", "#941751", "#333333", "#A8D9FF")
      plot.ind.pca <- plot.ind.pca +
        scale_color_manual(values = palette.col[1 : nlevels(coord.raters$Cluster)])
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
    plot.var.pca <- ggplot(NULL) +
      labs(x = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %",sep = "")) +
      coord_fixed() +
      xlim(xlim[1], xlim[2]) +
      ylim(ylim[1], ylim[2]) +
      geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
      geom_vline(xintercept = 0,  linetype = 2, color = "black", size = 0.2) +
      geom_segment(data = coord.stimuli, aes(x = 0, y = 0, xend = AxeA, yend = AxeB), alpha = 1, color = "black", size = 0.3, arrow = arrow(length = unit(0.3, "cm"))) +
      geom_text_repel(data = coord.stimuli, aes(x = AxeA, y = AxeB, label = rownames(coord.stimuli)), segment.color = "transparent", segment.size = 0.3, size = 2.3) +
      ggtitle("Representation of the stimuli") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        panel.background = element_rect(fill = 'white', colour = "black"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        legend.position = "none")
    plot.legend.clust <- res[[6]][[length(res[[6]])]]
    if (!is.null(col.clust)) {
      plot.legend.clust <- plot.legend.clust +
        scale_fill_manual(values = col.clust)
    }
    legend.plot <- get.legend(plot.legend.clust)
    main.title <- textGrob("Multidimensional representation of the structure \n of disagreement among the panel of raters", gp = gpar(fontsize = 12,font = 2))
    grid.arrange(arrangeGrob(plot.ind.pca + theme(legend.position="none"),
                           plot.var.pca + theme(legend.position="none"),
                           ncol = 2, nrow = 1),
               legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
  }

  # end the function
  options(warn = 0)
  return("Representations plotted")

}
