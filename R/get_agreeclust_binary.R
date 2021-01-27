#' Agreement-based clustering of binary ratings
#'
#' @param dta A binary matrix or data frame of dimensions S x R (S=number of stimuli, R=number of raters).
#' @param model A formula for the model. Either 'rating ~ rater + stimulus' (default), 'rating ~ rater', or 'rating ~ 1'.
#' @param max_clust An integer specifying the maximum number of clusters of raters. By default, this number is fixed to 10.
#' @param approx_null A boolean indicating if the null LRT distribution should be approximated using Satterthwaite's approximation. By default, the null LRT distribution is approximated.
#' @param paral_null A boolean indicating if the computation of the null LRT distribution should be parallelized. By default, the computation of the null LRT distribution is parallelized on nb.cores-1 cores. During the process, a text file 'TestDendrogram_processing.txt' is created.
#' @param consol A boolean indicating if a k-means consolidation of the partition of raters should be performed. By default, the partition is consolidated.
#' @param id_info_rater A vector of integer elements composed of the identification of the lines containing the supplementary information (i.e. covariates) about the raters. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the raters.
#' @param type_info_rater A vector of character elements composed of the type of the covariates about the raters. This vector must be of the same length that id.info.rater. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the raters.
#' @param id_info_stim A vector of integer elements composed of the identification of the columns containing the supplementary information (i.e. covariates) about the stimuli. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.
#' @param type_info_stim A vector of character elements composed of the type of the covariates about the stimuli. This vector must be of the same length that id.info.stim. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.
#' @param graph A boolean specifying if the graphical outputs should be plotted or not. By default, they are plotted.
#' @param ext_dev_Rstudio A boolean specifying if the graphical outputs should be plotted in the Rstudio plot pane or not.

#' @importFrom stats hclust
#' @importFrom reshape2 melt
#' @importFrom parallel detectCores makeCluster clusterExport clusterApply stopCluster
#' @importFrom ggdendro dendro_data
#' @import ggplot2
#' @importFrom grid unit textGrob gpar arrow
#' @importFrom FactoMineR PCA catdes
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom stats as.dendrogram as.dist binomial coef cutree deviance dist glm kmeans model.matrix na.omit pchisq pnorm predict quantile rchisq residuals runif summary.glm var vcov
#'
#' @return list
#' \itemize{
#' \item{profiles_residuals}{A matrix of dimensions S x R (S=number of stimuli, R=number of raters) containing the residuals profiles of the raters obtained through the modelling of the set of binary ratings.}
#' \item{mat_disag}{A matrix of dimensions R x R (R=number of raters) corresponding to the dissimilarity matrix between the raters.}
#' \item{pval_dendro}{A vector containing the probabilities associated to the statistical test realized at each level of the dendrogram.}
#' \item{nb_clust_found}{An integer corresponding to the number of clusters found among the panel.}
#' \item{partition}{A vector representing the partition of the raters (consolidated partition if consol = TRUE).}
#' \item{res_plot_segment}{All the graphical results of the segmentation.}
#' \item{res_pca}{All the results of the PCA.}
#' \item{charact_clust}{The results of the description of the clusters by information describing the raters and/or the stimuli.}
#' }
#'
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
#' res_pedag
#'
#' \dontrun{
#' data(goodgesture)
#' res_goodgesture <- get_agreeclust_bin(dta = goodgesture,
#'                                       model = "rating ~ rater + stimulus",
#'                                       id_info_rater = 40 : nrow(goodgesture),
#'                                       type_info_rater = rep("cat", 4),
#'                                       id_info_stim = 73 : ncol(goodgesture),
#'                                       type_info_stim = c(rep("cat", 3), rep("cont", 11)),
#'                                       paral_null = FALSE
#'                                       )
#' }
get_agreeclust_bin <- function(dta, model = "rating ~ rater + stimulus", max_clust = 10, approx_null = TRUE, paral_null = TRUE, consol = TRUE, id_info_rater = NULL, type_info_rater = NULL, id_info_stim = NULL, type_info_stim = NULL, graph = TRUE, ext_dev_Rstudio = FALSE) {

  # save the data set
  dta.sauv <- dta

  # remove external information about raters and stimuli
  if (!is.null(id_info_rater)) {
    dta <- dta[-id_info_rater,]
    dta <- droplevels(dta)
  }
  if (!is.null(id_info_stim)) {
    dta <- dta[, -id_info_stim]
    dta <- droplevels(dta)
  }

  # calculate the numbers of raters and stimuli
  nbrater <- ncol(dta)
  nbstim <- nrow(dta)

  # create a res object to save the results
  res <- list()

  # return the important arguments
  res[[1]] <- list(dta.sauv, id_info_rater, type_info_rater, id_info_stim, type_info_stim)
  names(res[[1]]) <- c("dta", "id_info_rater", "type_info_rater", "id_info_stim", "type_info_stim")

  # create the melted data set
  melted.data <- reshape2::melt(as.matrix(dta))
  melted.data <- melted.data[, c(2,1,3)]
  colnames(melted.data) <- c("rater", "stimulus", "rating")
  melted.data$rater <- as.factor(melted.data$rater)
  melted.data$stimulus <- as.factor(melted.data$stimulus)
  melted.data$rating <- as.factor(melted.data$rating)

  # adjust the no-latent class model
  mod.noLC <- glm(model, data = melted.data, family = binomial)

  # compute the disagreement matrix
  mat.resids <- matrix(residuals(mod.noLC, type = "deviance"), nrow(dta), ncol(dta))
  dimnames(mat.resids) <- list(rownames(dta), colnames(dta))
  mat.resids <- t(as.data.frame(mat.resids))
  disagmat <- dist(mat.resids, method = "euclidean")
  disagmat <- as.dist(disagmat)
  res[[2]] <- mat.resids
  res[[3]] <- as.matrix(disagmat)

  # construct the dendrogram
  dendrogram <- stats::hclust(disagmat, method = "ward.D2")

  # compute p-values for each level of the dendrogram
  message("Computation of the dendrogram testing in progress")
  message("If the computation is parallelized: Processing status can be checked in the 'TestDendrogram_processing.txt' file")
  message("If the computation is not parallelized: Processing status can be checked in the console")
  model2 <- paste(model, "+ cluster + stimulus:cluster")
  cut.model.rater <- strsplit(model2, "rater")[[1]]
  if (length(cut.model.rater) == 2) {
    model2 <- paste0(cut.model.rater[1], "rater%in%cluster", cut.model.rater[2])
  }

  control.break <<- TRUE
  list.pval <- lapply(1 : max_clust, compute_pval, dendrogram = dendrogram, melted.data = melted.data, model = model, model2 = model2, approx_null = approx_null,
                      dta = dta, paral_null = paral_null)
  pval <- unlist(list.pval)
  if (all(pval <= 0.05)) {
    nb.found <- max_clust
  } else {
    nb.found <- which(pval > 0.05)
  }
  partition.noconsol <- cutree(dendrogram, k = nb.found)
  mat.partition.noconsol <- cbind.data.frame(partition.noconsol, names(partition.noconsol))
  colnames(mat.partition.noconsol) <- c("cluster", "rater")
  res[[4]] <- pval
  res[[5]] <- nb.found

  # test the goodness of fit of the no-latent class model if nb.found = 1
  if (nb.found == 1) {
    res.test.noLC <- as.data.frame(matrix(NA, 1, 3))
    colnames(res.test.noLC) <- c("x", "y", "pval")
    res.test.noLC[1, "pval"] <- round(1 - pchisq(mod.noLC$deviance, mod.noLC$df.residual), 2)
  }

  # implement a partitioning algorithm to consolidate the partition
  if (consol == TRUE) {
    centers <- by(mat.resids, partition.noconsol, colMeans)
    centers <- matrix(unlist(centers), ncol = ncol(mat.resids), byrow = TRUE)
    res.consol <- kmeans(mat.resids, centers = centers, iter.max = 10)
    partition.consol <- res.consol$cluster
    mat.partition.consol <- cbind.data.frame(partition.consol, names(partition.consol))
    colnames(mat.partition.consol) <- c("cluster", "rater")
    res[[6]] <- partition.consol
  } else {
    res[[6]] <- partition.noconsol
  }

  # plot the basic dendrogram
  palette.col <- c("#90B08F", "#EA485C", "#FF8379", "#009193", "#FFCEA5", "#A9A9A9", "#B0983D", "#941751", "#333333", "#A8D9FF")
  dendrogram.info <- as.dendrogram(dendrogram)
  dendrogram.data <- ggdendro::dendro_data(dendrogram.info)
  data.labels <- dendrogram.data$labels
  colnames(data.labels)[3] <- "rater"
  data.labels <- merge(data.labels, mat.partition.noconsol, by = "rater")
  data.labels$cluster <- as.factor(data.labels$cluster)
  data.segments <- dendrogram.data$segments
  plot.dendro <- ggplot2::ggplot(NULL) +
    ggplot2::geom_segment(data = data.segments, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#444444") +
    ggplot2::geom_text(data = data.labels, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = cluster), size = 2.1) +
    ggplot2::scale_colour_manual(values = palette.col[1 : nlevels(data.labels$cluster)]) +
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
  if (consol == TRUE) {
    plot.dendro <- plot.dendro +
      ggplot2::ggtitle("Before consolidation")
  }

  # add p-values on the dendrogram
  res.test <- cbind.data.frame(dendrogram$height[length(dendrogram$height) : 1], rep(NA, length(dendrogram$height)), rep(NA, length(dendrogram$height)))
  colnames(res.test) <- c("height", "pval", "y")
  res.test[1 : length(pval), "pval"] <- pval
  res.test[(length(pval) + 1) : length(dendrogram$height), "pval"] <- NA
  for (i in 1 : nrow(res.test)) {
    if (i != nrow(res.test)) {
      res.test[i, "y"] <- res.test[i, "height"] - ((res.test[i, "height"] - res.test[(i + 1), "height"]) / 2)
    } else {
      res.test[i, "y"] <- res.test[i, "height"] - ((res.test[i, "height"] - 0) / 2)
    }
  }
  res.test <- res.test[-which(is.na(res.test[,"pval"])|round(res.test[,"height"], 5)==0),]
  plot.dendro <- plot.dendro +
    ggplot2::geom_hline(data = res.test, ggplot2::aes(yintercept = y), colour = "#444444", linetype = 2) +
    ggplot2::geom_text(data = res.test, ggplot2::aes(x = 1, y = (y + (0.2 * max(dendrogram$height) / 10))), label = res.test[,"pval"], colour = "#444444", size = 2.5)

  # add the p-value corresponding to the no-latent class model on the dendrogram
  if (nb.found == 1) {
    coord.all.nodes <- get_nodes_xy(dendrogram.info)
    res.test.noLC[, c(1,2)] <- coord.all.nodes[which(coord.all.nodes[, 2] == max(coord.all.nodes[, 2])), ]
    plot.dendro <- plot.dendro +
      ggplot2::ylim(-1, (max(data.segments$y) + (0.4 * max(dendrogram$height) / 10))) +
      ggplot2::geom_point(data = res.test.noLC, aes(x = x, y = y), colour = "#444444", size = 3, shape = 18) +
      ggplot2::geom_text(data = res.test.noLC, aes(x = x, y = (y + (0.4 * max(dendrogram$height) / 10))), label = res.test.noLC[, "pval"], colour = "#444444", size = 2.5) +
      ggplot2::guides(colour = guide_legend(override.aes = list(size=2.5)))
  } else {
    plot.dendro <- plot.dendro +
      ggplot2::ylim(-1, (max(data.segments$y) + (0.2 * max(dendrogram$height) / 10)))
  }

  # add legend to the dendrogram
  coord.legend.dendro <- as.data.frame(matrix(NA, 4, 2))
  coord.legend.dendro[, 1] <- c(0, 0, 10, 10)
  coord.legend.dendro[, 2] <- c(0, 1, 0, 1)
  colnames(coord.legend.dendro) <- c("x", "y")
  coord.test.height <- data.frame(x = 0, y = 0.8, xend = 0.6, yend = 0.8)
  coord.test.noLC <- data.frame(x = 0.3, y = 0.4)
  coord.legend.test.height <- data.frame(x = 1, y = 0.8)
  coord.legend.test.noLC <- data.frame(x = 1, y = 0.4)
  text.legend.test.height <- "p-value associated to the test of H0: this K-latent class structure is not significant (w.r.t. the (K-1)-latent class structure)"
  text.legend.test.noLC <- "p-value associated to the test of H0: the perfect agreement model well fits the data (displayed only if the number of clusters found equals 1)"
  plot.legend.dendro <- ggplot2::ggplot(NULL) +
    ggplot2::coord_fixed() +
    ggplot2::geom_point(data = coord.legend.dendro, ggplot2::aes(x = x, y = y), colour = "white", size = 2, shape = 18) +
    ggplot2::geom_segment(data = coord.test.height, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), colour = "#444444", linetype = 2) +
    ggplot2::geom_text(data = coord.test.height, ggplot2::aes(x = x, y = (y + 0.2)), label = "p-value", colour = "#444444", hjust = 0, size = 2) +
    ggplot2::geom_point(data = coord.test.noLC, ggplot2::aes(x = x, y = y), colour = "#444444", size = 2, shape = 18) +
    ggplot2::geom_text(data = coord.test.noLC, ggplot2::aes(x = x, y = (y - 0.2)), label = "p-value", colour = "#444444", size = 2) +
    ggplot2::geom_text(data = coord.legend.test.height, ggplot2::aes(x = x, y = y + 0.1), label = text.legend.test.height, hjust = 0, colour = "#444444", size = 2) +
    ggplot2::geom_text(data = coord.legend.test.noLC, ggplot2::aes(x = x, y = y - 0.1), label = text.legend.test.noLC, hjust = 0, colour = "#444444", size = 2) +
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

  # plot the partition of raters after consolidation if consol = TRUE
  if (consol == TRUE) {
    data.labels.partitioning <- merge(data.labels[, -4], mat.partition.consol, by = "rater")
    data.labels.partitioning$cluster <- as.factor(data.labels.partitioning$cluster)
    plot.partitioning <- ggplot2::ggplot(NULL) +
      ggplot2::geom_text(data = data.labels.partitioning, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = cluster), size = 2.1) +
      ggplot2::scale_colour_manual(values = palette.col[1 : nlevels(data.labels$cluster)]) +
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
  }

  # save the global legend of the clusters
  plot.legend.clust <- ggplot2::ggplot(NULL) +
    ggplot2::geom_label(data = data.labels, ggplot2::aes(label = rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = cluster), colour = "transparent") +
    ggplot2::scale_fill_manual(values = palette.col[1 : nlevels(data.labels$cluster)]) +
    ggplot2::theme(
      plot.margin = grid::unit(c(0,0,0,0), "cm"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size=8, colour = "#444444"),
      legend.text = ggplot2::element_text(size=8, colour = "#444444"),
      legend.margin = ggplot2::margin(t=0, unit='cm'),
      legend.key = ggplot2::element_rect(size=4),
      legend.key.size = grid::unit(0.4, "cm"))

  if ((Sys.getenv("RSTUDIO") == "1") == FALSE) {
    empty.dev <- (dev.cur() == 1)
  }
  legend.plot <- get.legend(plot.legend.clust)
  if ((Sys.getenv("RSTUDIO") == "1") == FALSE) {
    if (empty.dev == TRUE) {
      dev.off()
    }
  }

  # combine all plots
  main.title <- grid::textGrob("Raters clustering", gp = grid::gpar(fontsize = 12, font = 2, col = "#444444"))
  if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext_dev_Rstudio == TRUE) {
    dev.new(noRStudioGD = TRUE)
  }
  res[[7]] <- list(data.segments, data.labels, dendrogram, res.test, coord.legend.dendro, coord.test.height, coord.test.noLC, coord.legend.test.height, coord.legend.test.noLC)
  names(res[[7]]) <- c("data_segments", "data_labels", "dendrogram", "res_test", "coord_legend_dendro", "coord_test_height", "coord_test_noLC", "coord_legend_test_height", "coord_legend_test_noLC")
  if (nb.found == 1) {
    res[[7]][[length(res[[7]]) + 1]] <- res.test.noLC
    names(res[[7]])[length(res[[7]])] <- "res_test_noLC"
  } else {
    res[[7]][length(res[[7]]) + 1] <- list(NULL)
    names(res[[7]])[length(res[[7]])] <- "res_test_noLC"
  }
  if (consol == FALSE) {
    res[[7]][length(res[[7]]) + 1] <- list(NULL)
    names(res[[7]])[length(res[[7]])] <- "data_labels_partitioning"
    if (graph == TRUE) {
      gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.dendro + ggplot2::theme(legend.position = "none"),
                               plot.legend.dendro + ggplot2::theme(legend.position = "none"),
                               ncol = 1, nrow = 2, heights = c(4, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  } else if (consol == TRUE) {
    res[[7]][[length(res[[7]]) + 1]] <- data.labels.partitioning
    names(res[[7]])[length(res[[7]])] <- "data_labels_partitioning"
    if (graph == TRUE) {
      gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.dendro + ggplot2::theme(legend.position = "none"),
                               plot.legend.dendro + ggplot2::theme(legend.position = "none"),
                               plot.partitioning + ggplot2::theme(legend.position = "none"),
                               ncol = 1, nrow = 3, heights = c(4, 1, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  }

  # plot the multidimensional representation of the disagreement
  if (consol == TRUE) {
    mat.partition <- mat.partition.consol
  } else if (consol == FALSE) {
    mat.partition <- mat.partition.noconsol
  }
  res.pca <- FactoMineR::PCA(mat.resids, scale.unit = FALSE, ncp = Inf, graph = FALSE)
  res[[8]] <- res.pca
  axis = c(1, 2)
  coord.raters <- res.pca$ind$coord[, axis]
  mat.coord.raters <- cbind.data.frame(rownames(coord.raters), coord.raters)
  colnames(mat.coord.raters) <- c("rater", "AxeA", "AxeB")
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
  plot.ind.pca <- ggplot2::ggplot(NULL) +
    ggplot2::labs(x = paste("Dim ", axis[1]," - ", round(res.pca$eig[axis[1], 2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %", sep = "")) +
    ggplot2::coord_fixed() +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "black", size = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "black", size = 0.2) +
    ggplot2::geom_point(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, color = cluster)) +
    ggplot2::scale_color_manual(values = palette.col[1 : nlevels(data.labels$cluster)]) +
    ggrepel::geom_text_repel(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, label = rownames(coord.raters), color = cluster), segment.color = "#444444", segment.size = 0.3, size = 2.3) +
    ggplot2::geom_point(data = coord.raters, ggplot2::aes(x = AxeA, y = AxeB, color = cluster)) +
    ggplot2::ggtitle("Representation of the raters") +
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
  plot.var.pca <- ggplot2::ggplot(NULL) +
    ggplot2::labs(x = paste("Dim ", axis[1], " - ", round(res.pca$eig[axis[1],2], 2) , " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.pca$eig[axis[2], 2], 2), " %",sep = "")) +
    ggplot2::coord_fixed() +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#444444", size = 0.2) +
    ggplot2::geom_vline(xintercept = 0,  linetype = 2, color = "#444444", size = 0.2) +
    ggplot2::geom_segment(data = coord.stimuli, ggplot2::aes(x = 0, y = 0, xend = AxeA, yend = AxeB), alpha = 1, color = "black", size = 0.3, arrow = grid::arrow(length = grid::unit(0.3, "cm"))) +
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
  main.title <- grid::textGrob("Multidimensional representation of the structure \n of disagreement among the panel of raters", gp = grid::gpar(fontsize = 12, font = 2, col = "#444444"))
  if (graph == TRUE) {
    if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext_dev_Rstudio == TRUE) {
      dev.new(noRStudioGD = TRUE)
    }
    gridExtra::grid.arrange(gridExtra::arrangeGrob(plot.ind.pca + ggplot2::theme(legend.position="none"),
                             plot.var.pca + ggplot2::theme(legend.position="none"),
                             ncol = 2, nrow = 1),
                 legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
  }

  # interpret the clusters
  if (nb.found > 1) {
    res[[9]] <- list()
    mat.partition$cluster <- as.factor(mat.partition$cluster)
    # supplement the interpretation of the clusters with information about the raters
    list.charact_cluster_rater <- lapply(1 : nlevels(mat.partition$cluster), charact_cluster_rater, mat.partition = mat.partition, res.pca = res.pca, id_info_rater = id_info_rater, dta.sauv = dta.sauv,
                                         mat.resids = mat.resids, type_info_rater = type_info_rater, id_info_stim = id_info_stim)

    # supplement the interpretation of the clusters with information about the stimuli
    if (!is.null(id_info_stim)) {
      melted.data.clusters <- merge(melted.data, mat.partition, by = "rater")
      melted.data.clusters$rater <- as.factor(melted.data.clusters$rater)
      melted.data.clusters$stimulus <- as.factor(melted.data.clusters$stimulus)
      melted.data.clusters$rating <- as.factor(melted.data.clusters$rating)
      melted.data.clusters$cluster <- as.factor(melted.data.clusters$cluster)
      list.charact_cluster_stim <- lapply(1 : nlevels(mat.partition$cluster), charact_cluster_stim, mat.partition = mat.partition, id_info_rater = id_info_rater, id_info_stim = id_info_stim,
                                          dta.sauv = dta.sauv, dta = dta, type_info_stim = type_info_stim, melted.data.clusters = melted.data.clusters)
    }

    # combine the results
    for (i in 1 : nlevels(mat.partition$cluster)) {
      res[[9]][[i]] <- list.charact_cluster_rater[[i]]
      if (!is.null(id_info_stim)) {
        res[[9]][[i]][[length(res[[9]][[i]]) + 1]] <- list.charact_cluster_stim[[i]]
        if (length(list.charact_cluster_stim[[i]]) == 1 & is.null(res[[9]][[i]]$info.stim[[1]])) {
          res[[9]][[i]][length(res[[9]][[i]])] <- list(NULL)
        }
        names(res[[9]][[i]])[length(res[[9]][[i]])] <- "info.stim"
      }
    }
    names(res[[9]]) <- 1 : nlevels(mat.partition$cluster)
  } else {
    res[9] <- list(NULL)
  }

  # return the results
  names(res) <- c("call", "profiles_residuals", "mat_disag", "pval_dendro", "nb_clust_found", "partition", "res_plot_segment", "res_pca", "charact_clust")
  message("Clustering performed")

  class(res) <- c("agreeclust_binary", "list")
  return(res)

}
