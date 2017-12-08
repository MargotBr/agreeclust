plot.AgreeClust <- function(res, choice = "all", interact = FALSE, col.clust = NULL, axis = c(1, 2), name.rater = "rater", ext.dev.Rstudio = FALSE, vignette = FALSE) {

  options(warn = -1)

  # load packages
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  suppressPackageStartupMessages(require(ggrepel, quietly = TRUE))
  suppressPackageStartupMessages(require(plotly, quietly = TRUE))
  suppressPackageStartupMessages(require(manipulateWidget, quietly = TRUE))
  suppressPackageStartupMessages(require(shiny, quietly = TRUE))
  suppressPackageStartupMessages(require(reshape2, quietly = TRUE))
  suppressPackageStartupMessages(require(doBy, quietly = TRUE))

  # check the format of the arguments
  if (!inherits(res, "AgreeClustBin") & !inherits(res, "AgreeClustCont")) {
    stop("Non convenient data - res should be an AgreeClust object")
  }
  choice <- match.arg(choice, c("all", "seg", "mul"))
  mat.partition <- cbind.data.frame(res$partition, names(res$partition))
  colnames(mat.partition) <- c("Cluster", "Rater")
  mat.partition$Cluster <- as.factor(mat.partition$Cluster)
  nb.clust <- nlevels(mat.partition$Cluster)
  if (!is.null(col.clust)) {
    if (length(col.clust) < nb.clust) {
      stop("Non convenient specification of colors - col.clust should contain at least as many elements as clusters")
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

  # function to put the first letter as capital letter
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
  }

  # default palette
  palette.col <- c("#90B08F", "#EA485C", "#FF8379", "#009193", "#FFCEA5", "#A9A9A9", "#B0983D", "#941751", "#333333", "#A8D9FF")

  # plot the segmentation
  if (choice == "all" | choice == "seg") {
    plot.dendro <- ggplot(NULL) +
      geom_segment(data = res$res.plot.segment$data.segments, aes(x = x, y = y, xend = xend, yend = yend), colour = "#444444") +
      geom_text(data = res$res.plot.segment$data.labels, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = Cluster), size = 2.1) +
      ylim(-1, (max(res$res.plot.segment$data.segments$y) + (0.2 * max(res$res.plot.segment$dendrogram$height) / 10))) +
      theme(
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.title = element_text(colour = "#444444"),
        legend.text = element_text(colour = "#444444"),
        panel.background = element_rect(fill = 'white', colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 10, colour = "#444444"),
        plot.margin = unit(c(0.5,0,0,0), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
    if (!is.null(col.clust)) {
      plot.dendro <- plot.dendro +
        scale_colour_manual(values = col.clust[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
    } else {
      plot.dendro <- plot.dendro +
        scale_colour_manual(values = palette.col[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
    }
    if (!is.null(res$res.plot.segment$data.labels.partitioning)) {
      plot.dendro <- plot.dendro +
        ggtitle("Before consolidation")
    }
    plot.dendro <- plot.dendro +
      geom_hline(data = res$res.plot.segment$res.test, aes(yintercept = y), colour = "#444444", linetype = 2) +
      geom_text(data = res$res.plot.segment$res.test, aes(x = 1, y = (y + (0.2 * max(res$res.plot.segment$dendrogram$height) / 10))), label = res$res.plot.segment$res.test[,"pval"], colour = "#444444", size = 2.5)
    if (!is.null(res$res.plot.segment$res.test.noLC)) {
      plot.dendro <- plot.dendro +
        ylim(-1, (max(res$res.plot.segment$data.segments$y) + (0.3 * max(res$res.plot.segment$dendrogram$height) / 10))) +
        geom_point(data = res$res.plot.segment$res.test.noLC, aes(x = x, y = y), colour = "#444444", size = 3, shape = 18) +
        geom_text(data = res$res.plot.segment$res.test.noLC, aes(x = x, y = (y + (0.3 * max(res$res.plot.segment$dendrogram$height) / 10))), label = res$res.plot.segment$res.test.noLC[, "pval"], colour = "#444444", size = 2.5) +
        guides(colour = guide_legend(override.aes = list(size=2.5)))
    }
    text.legend.test.height <- "p-value associated to the test of H0: this K-latent class structure is not significant (w.r.t. the (K-1)-latent class structure)"
    text.legend.test.noLC <- "p-value associated to the test of H0: the perfect agreement model well fits the data (displayed only if the number of clusters found equals 1)"
    plot.legend.dendro <- ggplot(NULL) +
      coord_fixed() +
      geom_point(data = res$res.plot.segment$coord.legend.dendro, aes(x = x, y = y), colour = "white", size = 2, shape = 18) +
      geom_segment(data = res$res.plot.segment$coord.test.height, aes(x = x, xend = xend, y = y, yend = yend), colour = "#444444", linetype = 2) +
      geom_text(data = res$res.plot.segment$coord.test.height, aes(x = x, y = (y + 0.2)), label = "p-value", colour = "#444444", hjust = 0, size = 2) +
      geom_point(data = res$res.plot.segment$coord.test.noLC, aes(x = x, y = y), colour = "#444444", size = 2, shape = 18) +
      geom_text(data = res$res.plot.segment$coord.test.noLC, aes(x = x, y = (y - 0.2)), label = "p-value", colour = "#444444", size = 2) +
      geom_text(data = res$res.plot.segment$coord.legend.test.height, aes(x = x, y = y + 0.1), label = text.legend.test.height, hjust = 0, colour = "#444444", size = 2) +
      geom_text(data = res$res.plot.segment$coord.legend.test.noLC, aes(x = x, y = y - 0.1), label = text.legend.test.noLC, hjust = 0, colour = "#444444", size = 2) +
      theme(
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.title = element_text(colour = "#444444"),
        legend.text = element_text(colour = "#444444"),
        panel.background = element_rect(fill = 'white', colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
    if (!is.null(res$res.plot.segment$data.labels.partitioning)) {
      plot.partitioning <- ggplot(NULL) +
        geom_text(data = res$res.plot.segment$data.labels.partitioning, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = Cluster), size = 2.1) +
        ylim(-0.3, 0) +
        ggtitle("After consolidation") +
        theme(
          legend.key = element_rect(colour = "white",fill = "white"),
          legend.title = element_text(colour = "#444444"),
          legend.text = element_text(colour = "#444444"),
          panel.background = element_rect(fill = 'white', colour = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          plot.title = element_text(hjust = 0.5, size = 10, colour = "#444444"),
          plot.margin = unit(c(0.5, 0, 0, 0.45), "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")
      if (!is.null(col.clust)) {
        plot.partitioning <- plot.partitioning +
          scale_colour_manual(values = col.clust[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
      } else {
        plot.partitioning <- plot.partitioning +
          scale_colour_manual(values = palette.col[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
      }
    }
    plot.legend.clust <- ggplot(NULL) +
      geom_label(data = res$res.plot.segment$data.labels, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = Cluster), colour = "transparent") +
      theme(
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "bottom",
        legend.title = element_text(size=8, colour = "#444444"),
        legend.text = element_text(size=8, colour = "#444444"),
        legend.margin = margin(t=0, unit='cm'),
        legend.key = element_rect(size=4),
        legend.key.size = unit(0.4, "cm"))
    if (!is.null(col.clust)) {
      plot.legend.clust <- plot.legend.clust +
        scale_fill_manual(values = col.clust[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
    } else {
      plot.legend.clust <- plot.legend.clust +
        scale_fill_manual(values = palette.col[1 : nlevels(res$res.plot.segment$data.labels$Cluster)])
    }
    empty.dev <- (dev.cur() == 1)
    legend.plot <- get.legend(plot.legend.clust)
    if (empty.dev == TRUE) {
      dev.off()
    }
    main.title <- textGrob(paste(paste0(simpleCap(name.rater), "s"), "clustering"), gp = gpar(fontsize = 12, font = 2, col = "#444444"))
    if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
      dev.new(noRStudioGD = TRUE)
    }
    if (is.null(res$res.plot.segment$data.labels.partitioning)) {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               ncol = 1, nrow = 2, heights = c(4, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    } else {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               plot.partitioning + theme(legend.position = "none"),
                               ncol = 1, nrow = 3, heights = c(4, 1, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  }

  # plot the multidimensional representation of disagreement
  if (choice == "all" | choice == "mul") {
    res.pca <- res$res.pca
    coord.raters <- res.pca$ind$coord[, axis]
    mat.coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
    colnames(mat.coord.raters) <- c("Rater", "AxeA", "AxeB")
    mat.partition <- cbind.data.frame(names(res$partition), res$partition)
    colnames(mat.partition) <- c("Rater", "Cluster")
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
    if (interact == FALSE) {
      plot.ind.pca <- ggplot(NULL) +
        labs(x = paste("Dim ", axis[1]," - ", round(res.pca$eig[axis[1], 2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %", sep = "")) +
        coord_fixed() +
        xlim(xlim[1], xlim[2]) +
        ylim(ylim[1], ylim[2]) +
        geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
        geom_vline(xintercept = 0, linetype = 2, color = "black", size = 0.2) +
        geom_point(data = coord.raters, aes(x = AxeA, y = AxeB, color = Cluster)) +
        geom_text_repel(data = coord.raters, aes(x = AxeA, y = AxeB, label = rownames(coord.raters), color = Cluster), segment.color = "#444444", segment.size = 0.3, size = 2.3) +
        geom_point(data = coord.raters, aes(x = AxeA, y = AxeB, color = Cluster)) +
        ggtitle(paste("Representation of the", paste0(name.rater, "s"))) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, colour = "#444444"),
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
          panel.background = element_rect(fill = 'white', colour = "#444444"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          axis.text = element_text(colour = "#444444"),
          axis.ticks = element_line(colour = "#444444"),
          axis.title = element_text(colour = "#444444"),
          legend.position = "none")
      if (!is.null(col.clust)) {
        plot.ind.pca <- plot.ind.pca +
          scale_color_manual(values = col.clust[1 : nlevels(coord.raters$Cluster)])
      } else {
        plot.ind.pca <- plot.ind.pca +
          scale_color_manual(values = palette.col[1 : nlevels(coord.raters$Cluster)])
      }
    } else if (interact == TRUE) {
      coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
      colnames(coord.raters)[1] <- "Rater"
      coord.raters$Cluster <- paste("Cluster", coord.raters$Cluster)
      if (!is.null(col.clust)) {
        col <- col.clust[1 : nlevels(as.factor(coord.raters$Cluster))]
      } else {
        col = palette.col[1 : nlevels(as.factor(coord.raters$Cluster))]
      }
      text.tooltip <- paste(paste0(simpleCap(name.rater), ":"), coord.raters$Rater, '<br>Cluster:', gsub("Cluster ", "", coord.raters$Cluster))
      if (!is.null(res$call$id.info.rater)) {
        info.rater <- res$call$dta[res$call$id.info.rater, ]
        for (i in 1 : length(res$call$id.info.rater)) {
          info <- droplevels(info.rater[i, 1 : length(coord.raters$Rater)][order(info.rater[i, 1 : length(coord.raters$Rater)], coord.raters$Rater)])
          text.tooltip <- paste(text.tooltip, paste0("<br>", rownames(info.rater)[i], ": ", as.character(unlist(info))))
        }
      }
      plot.ind.pca.interact <- plot_ly(coord.raters) %>%
        add_trace(coord.raters,
                  x = ~AxeA ,
                  y = ~AxeB,
                  color = ~Cluster,
                  hoverinfo = 'text',
                  text = text.tooltip,
                  type = "scatter", mode = "markers",
                  marker = list(size = 6),
                  showlegend = TRUE,
                  colors = col) %>%
        layout(
          legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.25),
          margin = list(l = 50, r = 50, t = 50, b = 50),
          titlefont = list(size = 14, color = "#444444"),
          title = paste("Representation of the", paste0(name.rater, "s")),
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
          plot.title = element_text(hjust = 0.5, size = 12, colour = "#444444"),
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
          panel.background = element_rect(fill = 'white', colour = "#444444"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"),
          axis.text = element_text(colour = "#444444"),
          axis.ticks = element_line(colour = "#444444"),
          axis.title = element_text(colour = "#444444"),
          legend.position = "none")
      plot.legend.clust <- ggplot(NULL) +
        geom_label(data = res$res.plot.segment$data.labels, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = Cluster), colour = "transparent") +
        theme(
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "bottom",
          legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.margin = margin(t=0, unit='cm'),
          legend.key = element_rect(size=4),
          legend.key.size = unit(0.4, "cm"))
      if (!is.null(col.clust)) {
        plot.legend.clust <- plot.legend.clust +
          scale_fill_manual(values = col.clust[1 : nlevels(coord.raters$Cluster)])
      } else {
        plot.legend.clust <- plot.legend.clust +
          scale_fill_manual(values = palette.col[1 : nlevels(coord.raters$Cluster)])
      }
      empty.dev <- (dev.cur() == 1)
      legend.plot <- get.legend(plot.legend.clust)
      if (empty.dev == TRUE) {
        dev.off()
      }
      main.title <- textGrob(paste("Multidimensional representation of the structure \n of disagreement among the panel of", paste0(name.rater, "s")), gp = gpar(fontsize = 12, font = 2, col = "#444444"))
      if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
        dev.new(noRStudioGD = TRUE)
      }
      grid.arrange(arrangeGrob(plot.ind.pca + theme(legend.position="none"),
                             plot.var.pca + theme(legend.position="none"),
                             ncol = 2, nrow = 1),
                 legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    } else if (interact == TRUE) {
      coord.stimuli <- cbind.data.frame(rownames(coord.stimuli), coord.stimuli, paste("Cluster", c(1 : nlevels(as.factor(coord.raters$Cluster)), sample(1 : nlevels(as.factor(coord.raters$Cluster)), (nrow(coord.stimuli) - nlevels(as.factor(coord.raters$Cluster))), replace = TRUE))))
      colnames(coord.stimuli)[c(1, 4)] <- c("Stimulus", "FakeGroup")
      coord.stimuli$FakeGroup <- as.factor(coord.stimuli$FakeGroup)
      text.tooltip <- paste("Stimulus:", coord.stimuli$Stimulus)
      if (!is.null(res$call$id.info.rater)) {
        dta <- res$call$dta[-res$call$id.info.rater,]
        dta <- droplevels(dta)
      }
      if (!is.null(res$call$id.info.stim)) {
        dta <- dta[, -res$call$id.info.stim]
        dta <- droplevels(dta)
      }
      melted.data <- melt(as.matrix(dta))
      colnames(melted.data) <- c("Stimulus", "Rater", "Rating")
      mat.partition <- cbind.data.frame(names(res$partition), res$partition)
      colnames(mat.partition) <- c("Rater", "Cluster")
      melted.cluster <- merge(melted.data, mat.partition, by = "Rater")
      melted.cluster$Rating <- as.numeric(as.character(melted.cluster$Rating))
      if (inherits(res, "AgreeClustBin")) {
        sum.pos.ratings <- summaryBy(Rating ~ Stimulus : Cluster, data = melted.cluster, FUN = "sum")
        colnames(sum.pos.ratings) <- c("Stimulus", "Cluster", "NbPosRating")
        size.clust <- summary(as.factor(mat.partition$Cluster))
        for (i in 1 : nlevels(as.factor(mat.partition$Cluster))) {
          tab.clust <- sum.pos.ratings[which(sum.pos.ratings$Cluster == levels(as.factor(mat.partition$Cluster))[i]), ]
          tab.clust[, "NbPosRating"] <- tab.clust[, "NbPosRating"] / size.clust[which(names(size.clust) == levels(as.factor(mat.partition$Cluster))[i])] * 100
          text.tooltip <- paste(text.tooltip,
                                paste(paste0("<br>Positive ratings in Cluster ", levels(as.factor(mat.partition$Cluster))[i], ":"),
                                      paste0(tab.clust[order(tab.clust$Stimulus, coord.stimuli$Stimulus), "NbPosRating"], "%")))
        }
      }
      if (inherits(res, "AgreeClustCont")) {
        sum.ratings <- summaryBy(Rating ~ Stimulus : Cluster, data = melted.cluster, FUN = "mean")
        colnames(sum.ratings) <- c("Stimulus", "Cluster", "MeanRating")
        for (i in 1 : nlevels(as.factor(mat.partition$Cluster))) {
          tab.clust <- sum.ratings[which(sum.ratings$Cluster == levels(as.factor(mat.partition$Cluster))[i]), ]
          text.tooltip <- paste(text.tooltip,
                                paste(paste0("<br>Average rating in Cluster ", levels(as.factor(mat.partition$Cluster))[i], ":"),
                                      paste0(tab.clust[order(tab.clust$Stimulus, coord.stimuli$Stimulus), "MeanRating"], "")))
        }
      }
      if (!is.null(res$call$id.info.stim)) {
        info.stim <- res$call$dta[, res$call$id.info.stim]
        for (i in 1 : length(res$call$id.info.stim)) {
          info <- info.stim[1 : length(coord.stimuli$Stimulus), i][order(info.stim[1 : length(coord.stimuli$Stimulus), i], coord.stimuli$Stimulus)]
          text.tooltip <- paste(text.tooltip, paste0("<br>", colnames(info.stim)[i], ": ", as.character(unlist(info))))
        }
      }
      plot.var.pca.interact <- plot_ly(coord.stimuli) %>%
        add_trace(coord.stimuli,
                  x = ~AxeA ,
                  y = ~AxeB,
                  color = ~FakeGroup,
                  hoverinfo = 'text',
                  hoverlabel = list(bgcolor = "black", font = list(color = "white")),
                  text = text.tooltip,
                  type = "scatter", mode = "markers",
                  marker = list(size = 6, color = "white"),
                  showlegend = TRUE) %>%
        add_annotations(axref = "x", ax = coord.stimuli$AxeA, aref = "x", x = 0,
                  ayref = "y", ay = coord.stimuli$AxeB, aref = "y", y = 0,
                  text = coord.stimuli$Stimulus, arrowhead = 0, arrowcolor = "black") %>%
        layout(
          legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.25, font = list(color = "transparent")),
          margin = list(l = 50, r = 50, t = 50, b = 50),
          title = "Representation of the stimuli",
          titlefont = list(size = 14, color = "#444444"),
          xaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "y", showgrid = FALSE, title = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1),
          yaxis = list(zerolinecolor = "#D6D5D5", scaleanchor = "x", showgrid = FALSE, title = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2],2],2), "%", sep=""), titlefont = list(color = "#444444", size = 13), tickfont = list(size = 10, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1)
        )
      combine.plot <- combineWidgets(plot.ind.pca.interact, plot.var.pca.interact, nrow = 1, ncol = 2,
                                     title = paste("Multidimensional representation of the structure <br>of disagreement among the panel of", paste0(name.rater, "s")),
                                     titleCSS = "font-size:16px",
                                     header = tags$div(tags$span("header", style = "color:white;font-size:20px")),
                                     leftCol = tags$div(tags$span("header", style = "color:white;font-size:10px")),
                                     rightCol = tags$div(tags$span("header", style = "color:white;font-size:10px")))
      if (vignette == FALSE) {
        if (ext.dev.Rstudio == TRUE) {
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
  options(warn = 0)
  if (vignette == TRUE & interact == TRUE & (choice == "all" | choice == "mul")) {
    message("Representations plotted")
    return(combine.plot)
  } else {
    return(message("Representations plotted"))
  }

}
