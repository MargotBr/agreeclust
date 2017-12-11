AgreeClustCont <- function(dta, model = "Rating ~ Rater + Stimulus", max.clust = 10, paral.null = TRUE, consol = TRUE, id.info.rater = NULL, type.info.rater = NULL, id.info.stim = NULL, type.info.stim = NULL, graph = TRUE, ext.dev.Rstudio = FALSE) {

  options(warn = -1)

  # load packages
  suppressPackageStartupMessages(require(reshape2, quietly = TRUE))
  suppressPackageStartupMessages(require(parallel, quietly = TRUE))
  suppressPackageStartupMessages(require(ggdendro, quietly = TRUE))
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  suppressPackageStartupMessages(require(dendextend, quietly = TRUE))
  suppressPackageStartupMessages(require(grid, quietly = TRUE))
  suppressPackageStartupMessages(require(gridExtra, quietly = TRUE))
  suppressPackageStartupMessages(require(FactoMineR, quietly = TRUE))
  suppressPackageStartupMessages(require(ggrepel, quietly = TRUE))

  # save the data set
  dta.sauv <- dta

  # remove external information about raters and stimuli
  if (!is.null(id.info.rater)) {
    dta <- dta[-id.info.rater,]
    dta <- droplevels(dta)
  }
  if (!is.null(id.info.stim)) {
    dta <- dta[, -id.info.stim]
    dta <- droplevels(dta)
  }

  # calculate the numbers of raters and stimuli
  nbrater <- ncol(dta)
  nbstim <- nrow(dta)

  # create a res object to save the results
  res <- list()

  # return the important arguments
  res[[1]] <- list(dta.sauv, id.info.rater, type.info.rater, id.info.stim, type.info.stim)
  names(res[[1]]) <- c("dta", "id.info.rater", "type.info.rater", "id.info.stim", "type.info.stim")

  # create the melted data set
  melted.data <- melt(as.matrix(dta))
  melted.data <- melted.data[, c(2,1,3)]
  colnames(melted.data) <- c("Rater", "Stimulus", "Rating")
  melted.data$Rater <- as.factor(melted.data$Rater)
  melted.data$Stimulus <- as.factor(melted.data$Stimulus)
  melted.data$Rating <- as.numeric(as.character(melted.data$Rating))

  # adjust the no-latent class model
  mod.noLC <- lm(model, data = melted.data)

  # compute the disagreement matrix
  mat.resids <- matrix(residuals(mod.noLC), nrow(dta), ncol(dta))
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
  model2 <- paste(model, "+ Cluster + Stimulus:Cluster")
  cut.model.Rater <- strsplit(model2, "Rater")[[1]]
  if (length(cut.model.Rater) == 2) {
    model2 <- paste0(cut.model.Rater[1], "Rater%in%Cluster", cut.model.Rater[2])
  }
  #remove_outliers <- function(x, na.rm = TRUE, ...) {
  #  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  #  H <- 10 * IQR(x, na.rm = na.rm)
  #  y <- x
  #  y[x < (qnt[1] - H)] <- NA
  #  y[x > (qnt[2] + H)] <- NA
  #  y
  #}
  compute.dist.null <- function(j, list.data.null, K) {

    # Null dendrogram
    message(paste0("Comparison test: ", (K + 1), "-latent class model VS " , (K), "-latent class model | Progress: ", round(j / length(list.data.null) * 100, 1), "%"))
    data.null <- as.matrix(list.data.null[[j]])
    melted.data.null <- melt(data.null)
    melted.data.null <- melted.data.null[, c(2,1,3)]
    colnames(melted.data.null) <- c("Rater", "Stimulus", "Rating")
    melted.data.null$Rater <- as.factor(melted.data.null$Rater)
    melted.data.null$Stimulus <- as.factor(melted.data.null$Stimulus)
    melted.data.null$Rating <- as.numeric(as.character(melted.data.null$Rating))
    mod.noLC.null <- lm(model, data = melted.data.null)
    mat.resids.null <- matrix(residuals(mod.noLC.null), nrow(data.null), ncol(data.null))
    dimnames(mat.resids.null) <- list(rownames(data.null), colnames(data.null))
    mat.resids.null <- t(as.data.frame(mat.resids.null))
    disagmat.null <- dist(mat.resids.null, method = "euclidean")
    disagmat.null <- as.dist(disagmat.null)
    dendrogram.null <- stats::hclust(disagmat.null, method = "ward.D2")

    # Null K-latent class model
    partition.KLC.null <- cutree(dendrogram.null, k = K)
    partition.KLC.null <- cbind.data.frame(names(partition.KLC.null), partition.KLC.null)
    colnames(partition.KLC.null) <- c("Rater", "Cluster")
    clustered.data.KLC.null <- merge(melted.data.null, partition.KLC.null, by = "Rater")
    clustered.data.KLC.null$Cluster <- as.factor(clustered.data.KLC.null$Cluster)
    if (K == 1) {
      mod.KLC.null <- lm(model, data = clustered.data.KLC.null)
    } else {
      mod.KLC.null <- lm(model2, data = clustered.data.KLC.null)
      X.KLC.null <- model.matrix(mod.KLC.null)
      X.KLC.null <- X.KLC.null[, !is.na(coef(mod.KLC.null))]
      mod.KLC.null <- lm(Rating ~ ., data = data.frame(Rating = clustered.data.KLC.null$Rating, X.KLC.null[, -1]))
    }

    # Null (K+1)-latent class model
    partition.Kplus1LC.null <- cutree(dendrogram.null, k = K + 1)
    partition.Kplus1LC.null <- cbind.data.frame(names(partition.Kplus1LC.null), partition.Kplus1LC.null)
    colnames(partition.Kplus1LC.null) <- c("Rater", "Cluster")
    clustered.data.Kplus1LC.null <- merge(melted.data.null, partition.Kplus1LC.null, by = "Rater")
    clustered.data.Kplus1LC.null$Cluster <- as.factor(clustered.data.Kplus1LC.null$Cluster)
    mod.Kplus1LC.null <- lm(model2, data = clustered.data.Kplus1LC.null)
    X.Kplus1LC.null <- model.matrix(mod.Kplus1LC.null)
    X.Kplus1LC.null <- X.Kplus1LC.null[, !is.na(coef(mod.Kplus1LC.null))]
    mod.Kplus1LC.null <- lm(Rating ~ ., data = data.frame(Rating = clustered.data.Kplus1LC.null$Rating, X.Kplus1LC.null[, -1]))

    # Compute the null F statistic
    f.H0 <- anova(mod.KLC.null, mod.Kplus1LC.null, test = "F")$F[2]
    return(f.H0)

  }
  compute.pval <- function(K) {

    if(control.break) {

      # K-latent class model
      partition.KLC <- cutree(dendrogram, k = K)
      partition.KLC <- cbind.data.frame(names(partition.KLC), partition.KLC)
      colnames(partition.KLC) <- c("Rater", "Cluster")
      clustered.data.KLC <- merge(melted.data, partition.KLC, by = "Rater")
      clustered.data.KLC$Cluster <- as.factor(clustered.data.KLC$Cluster)
      if (K == 1) {
        mod.KLC <- lm(model, data = clustered.data.KLC)
      } else {
        mod.KLC <- lm(model2, data = clustered.data.KLC)
        X.KLC <- model.matrix(mod.KLC)
        X.KLC <- X.KLC[, !is.na(coef(mod.KLC))]
        mod.KLC <- lm(Rating ~ ., data = data.frame(Rating = clustered.data.KLC$Rating, X.KLC[, -1]))
      }

      # (K+1)-latent class model
      partition.Kplus1LC <- cutree(dendrogram, k = K + 1)
      partition.Kplus1LC <- cbind.data.frame(names(partition.Kplus1LC), partition.Kplus1LC)
      colnames(partition.Kplus1LC) <- c("Rater", "Cluster")
      clustered.data.Kplus1LC <- merge(melted.data, partition.Kplus1LC, by = "Rater")
      clustered.data.Kplus1LC$Cluster <- as.factor(clustered.data.Kplus1LC$Cluster)
      mod.Kplus1LC <- lm(model2, data = clustered.data.Kplus1LC)
      X.Kplus1LC <- model.matrix(mod.Kplus1LC)
      X.Kplus1LC <- X.Kplus1LC[,!is.na(coef(mod.Kplus1LC))]
      mod.Kplus1LC <- lm(Rating~., data=data.frame(Rating = clustered.data.Kplus1LC$Rating, X.Kplus1LC[, -1]))

      # Compute the F statistic
      f.obs <- anova(mod.KLC, mod.Kplus1LC, test = "F")$F[2]

      # Generate the null data sets
      value.null <- predict(mod.KLC, type = "response")
      #if (approx.null == FALSE) {
        nb.simul.null <- 1000
      #} else if (approx.null == TRUE) {
      #  nb.simul.null <- 250
      #}
      list.data.null <- lapply(1:nb.simul.null, function(l, value.null) {
        as.numeric(value.null + runif(length(value.null)))}, value.null = value.null)
      list.data.null <- lapply(list.data.null, function(vec, k, dnames) matrix(vec, nrow = k, dimnames = dnames), k = nlevels(melted.data$Stimulus), dnames = dimnames(dta))

      # Compute p-value
      if (paral.null == FALSE) {
        list.f.H0 <- lapply(1 : nb.simul.null, compute.dist.null, list.data.null = list.data.null, K = K)
      } else if (paral.null == TRUE) {
        environment(compute.dist.null) <- .GlobalEnv
        nb.cores <- detectCores()
        cl <- makeCluster(nb.cores - 1, outfile = "TestDendrogram_processing.txt")
        clusterExport(cl, varlist = c("melt", "model", "model2", "list.data.null", "K"), environment())
        list.f.H0 <- clusterApply(cl, 1 : nb.simul.null, compute.dist.null, list.data.null = list.data.null, K = K)
        stopCluster(cl)
        if ("TestDendrogram_processing.txt" %in% list.files()) {
          file.remove("TestDendrogram_processing.txt")
        }
      }
      f.H0 <- unlist(list.f.H0)
      #if (approx.null == TRUE) {
      #  if (mean(lrt.H0) == 0 & var(lrt.H0) == 0) {
      #    lrt.H0 <- rep(0, 1000)
      #  } else {
      #    lrt.H0 <- remove_outliers(lrt.H0)
      #    if (length(which(is.na(lrt.H0))) != 0) {
      #      lrt.H0 <- lrt.H0[-which(is.na(lrt.H0))]
      #    }
      #    nu <- (mean(lrt.H0)^2) / ((1/2) * var(lrt.H0))
      #    c <- mean(lrt.H0) / nu
      #    lrt.H0 <- rep(0, 1000)
      #    for (i in 1 : length(lrt.H0)) {
      #      lrt.H0[i] = c * rchisq(1, nu, ncp = 0)
      #   }
      #  }
      #}
      pval <- length(which(f.obs <= f.H0)) / 1000
      if (pval >= 0.05) {
        control.break <<- FALSE
      }
      return(pval)

    } else {
      return(NULL)
    }
  }
  control.break <- TRUE
  list.pval <- lapply(1 : max.clust, compute.pval)
  pval <- unlist(list.pval)
  if (all(pval <= 0.05)) {
    nb.found <- max.clust
  } else {
    nb.found <- which(pval > 0.05)
  }
  partition.noconsol <- cutree(dendrogram, k = nb.found)
  mat.partition.noconsol <- cbind.data.frame(partition.noconsol, names(partition.noconsol))
  colnames(mat.partition.noconsol) <- c("Cluster", "Rater")
  res[[4]] <- pval
  res[[5]] <- nb.found

  # test the goodness of fit of the no-latent class model if nb.found = 1
  if (nb.found == 1) {
    res.test.noLC <- as.data.frame(matrix(NA, 1, 3))
    colnames(res.test.noLC) <- c("x", "y", "pval")
    res.test.noLC[1, "pval"] <- round(1 - pf(summary(mod.noLC)$fstatistic[1], summary(mod.noLC)$fstatistic[2], summary(mod.noLC)$fstatistic[3]), 2)
  }

  # implement a partitioning algorithm to consolidate the partition
  if (consol == TRUE) {
    centers <- by(mat.resids, partition.noconsol, colMeans)
    centers <- matrix(unlist(centers), ncol = ncol(mat.resids), byrow = TRUE)
    res.consol <- kmeans(mat.resids, centers = centers, iter.max = 10)
    partition.consol <- res.consol$cluster
    mat.partition.consol <- cbind.data.frame(partition.consol, names(partition.consol))
    colnames(mat.partition.consol) <- c("Cluster", "Rater")
    res[[6]] <- partition.consol
  } else {
    res[[6]] <- partition.noconsol
  }

  # plot the basic dendrogram
  palette.col <- c("#90B08F", "#EA485C", "#FF8379", "#009193", "#FFCEA5", "#A9A9A9", "#B0983D", "#941751", "#333333", "#A8D9FF")
  dendrogram.info <- as.dendrogram(dendrogram)
  dendrogram.data <- dendro_data(dendrogram.info)
  data.labels <- label(dendrogram.data)
  colnames(data.labels)[3] <- "Rater"
  data.labels <- merge(data.labels, mat.partition.noconsol, by = "Rater")
  data.labels$Cluster <- as.factor(data.labels$Cluster)
  data.segments <- segment(dendrogram.data)
  plot.dendro <- ggplot(NULL) +
    geom_segment(data = data.segments, aes(x = x, y = y, xend = xend, yend = yend), colour = "#444444") +
    geom_text(data = data.labels, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = Cluster), size = 2.1) +
    scale_colour_manual(values = palette.col[1 : nlevels(data.labels$Cluster)]) +
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
  if (consol == TRUE) {
    plot.dendro <- plot.dendro +
      ggtitle("Before consolidation")
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
    geom_hline(data = res.test, aes(yintercept = y), colour = "#444444", linetype = 2) +
    geom_text(data = res.test, aes(x = 1, y = (y + (0.2 * max(dendrogram$height) / 10))), label = res.test[,"pval"], colour = "#444444", size = 2.5)

  # add the p-value corresponding to the no-latent class model on the dendrogram
  if (nb.found == 1) {
    coord.all.nodes <- get_nodes_xy(dendrogram.info)
    res.test.noLC[, c(1,2)] <- coord.all.nodes[which(coord.all.nodes[, 2] == max(coord.all.nodes[, 2])), ]
    plot.dendro <- plot.dendro +
      ylim(-1, (max(data.segments$y) + (0.4 * max(dendrogram$height) / 10))) +
      geom_point(data = res.test.noLC, aes(x = x, y = y), colour = "#444444", size = 3, shape = 18) +
      geom_text(data = res.test.noLC, aes(x = x, y = (y + (0.4 * max(dendrogram$height) / 10))), label = res.test.noLC[, "pval"], colour = "#444444", size = 2.5) +
      guides(colour = guide_legend(override.aes = list(size=2.5)))
  } else {
    plot.dendro <- plot.dendro +
      ylim(-1, (max(data.segments$y) + (0.2 * max(dendrogram$height) / 10)))
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
  plot.legend.dendro <- ggplot(NULL) +
    coord_fixed() +
    geom_point(data = coord.legend.dendro, aes(x = x, y = y), colour = "white", size = 2, shape = 18) +
    geom_segment(data = coord.test.height, aes(x = x, xend = xend, y = y, yend = yend), colour = "#444444", linetype = 2) +
    geom_text(data = coord.test.height, aes(x = x, y = (y + 0.2)), label = "p-value", colour = "#444444", hjust = 0, size = 2) +
    geom_point(data = coord.test.noLC, aes(x = x, y = y), colour = "#444444", size = 2, shape = 18) +
    geom_text(data = coord.test.noLC, aes(x = x, y = (y - 0.2)), label = "p-value", colour = "#444444", size = 2) +
    geom_text(data = coord.legend.test.height, aes(x = x, y = y + 0.1), label = text.legend.test.height, hjust = 0, colour = "#444444", size = 2) +
    geom_text(data = coord.legend.test.noLC, aes(x = x, y = y - 0.1), label = text.legend.test.noLC, hjust = 0, colour = "#444444", size = 2) +
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

  # plot the partition of raters after consolidation if consol = TRUE
  if (consol == TRUE) {
    data.labels.partitioning <- merge(data.labels[, -4], mat.partition.consol, by = "Rater")
    data.labels.partitioning$Cluster <- as.factor(data.labels.partitioning$Cluster)
    plot.partitioning <- ggplot(NULL) +
      geom_text(data = data.labels.partitioning, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, colour = Cluster), size = 2.1) +
      scale_colour_manual(values = palette.col[1 : nlevels(data.labels$Cluster)]) +
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
  }

  # save the global legend of the clusters
  plot.legend.clust <- ggplot(NULL) +
    geom_label(data = data.labels, aes(label = Rater, x = x, y = -0.1, angle = 90, hjust = 1, fill = Cluster), colour = "transparent") +
    scale_fill_manual(values = palette.col[1 : nlevels(data.labels$Cluster)]) +
    theme(
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = "bottom",
      legend.title = element_text(size=8, colour = "#444444"),
      legend.text = element_text(size=8, colour = "#444444"),
      legend.margin = margin(t=0, unit='cm'),
      legend.key = element_rect(size=4),
      legend.key.size = unit(0.4, "cm"))
  get.legend <- function(plot){
    tmp <- ggplot_gtable(ggplot_build(plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
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
  main.title <- textGrob("Raters clustering", gp = gpar(fontsize = 12, font = 2, col = "#444444"))
  if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
    dev.new(noRStudioGD = TRUE)
  }
  res[[7]] <- list(data.segments, data.labels, dendrogram, res.test, coord.legend.dendro, coord.test.height, coord.test.noLC, coord.legend.test.height, coord.legend.test.noLC)
  names(res[[7]]) <- c("data.segments", "data.labels", "dendrogram", "res.test", "coord.legend.dendro", "coord.test.height", "coord.test.noLC", "coord.legend.test.height", "coord.legend.test.noLC")
  if (nb.found == 1) {
    res[[7]][[length(res[[7]]) + 1]] <- res.test.noLC
    names(res[[7]])[length(res[[7]])] <- "res.test.noLC"
  } else {
    res[[7]][length(res[[7]]) + 1] <- list(NULL)
    names(res[[7]])[length(res[[7]])] <- "res.test.noLC"
  }
  if (consol == FALSE) {
    res[[7]][length(res[[7]]) + 1] <- list(NULL)
    names(res[[7]])[length(res[[7]])] <- "data.labels.partitioning"
    if (graph == TRUE) {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               ncol = 1, nrow = 2, heights = c(4, 1)),
                   legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
    }
  } else if (consol == TRUE) {
    res[[7]][[length(res[[7]]) + 1]] <- data.labels.partitioning
    names(res[[7]])[length(res[[7]])] <- "data.labels.partitioning"
    if (graph == TRUE) {
      grid.arrange(arrangeGrob(plot.dendro + theme(legend.position = "none"),
                               plot.legend.dendro + theme(legend.position = "none"),
                               plot.partitioning + theme(legend.position = "none"),
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
  res.pca <- PCA(mat.resids, scale.unit = FALSE, ncp = Inf, graph = FALSE)
  res[[8]] <- res.pca
  axis = c(1, 2)
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
    scale_color_manual(values = palette.col[1 : nlevels(data.labels$Cluster)]) +
    geom_text_repel(data = coord.raters, aes(x = AxeA, y = AxeB, label = rownames(coord.raters), color = Cluster), segment.color = "#444444", segment.size = 0.3, size = 2.3) +
    geom_point(data = coord.raters, aes(x = AxeA, y = AxeB, color = Cluster)) +
    ggtitle("Representation of the raters") +
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
      plot.title = element_text(hjust = 0.5, size = 12, colour = "#444444"),
      plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
      panel.background = element_rect(fill = 'white', colour = "#444444"),
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_line(colour = "white"),
      axis.text = element_text(colour = "#444444"),
      axis.ticks = element_line(colour = "#444444"),
      axis.title = element_text(colour = "#444444"),
      legend.position = "none")
  main.title <- textGrob("Multidimensional representation of the structure \n of disagreement among the panel of raters", gp = gpar(fontsize = 12, font = 2, col = "#444444"))
  if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
    dev.new(noRStudioGD = TRUE)
  }
  if (graph == TRUE) {
    grid.arrange(arrangeGrob(plot.ind.pca + theme(legend.position="none"),
                           plot.var.pca + theme(legend.position="none"),
                           ncol = 2, nrow = 1),
               legend.plot, nrow = 2, top = main.title, heights = c(8, 1))
  }

  # interpret the clusters
  if (nb.found > 1) {
    res[[9]] <- list()
    mat.partition$Cluster <- as.factor(mat.partition$Cluster)
    # supplement the interpretation of the clusters with information about the raters
    charact.cluster.rater <- function (i) {

      res.clust.rater <- list()
      clust <- levels(mat.partition$Cluster)[i]

      # calculate the number of raters
      nb.rater.clust <- length(which(mat.partition$Cluster == clust))
      res.clust.rater[[1]] <- nb.rater.clust

      # calculate the percentage of raters
      perc.rater.clust <- round(nb.rater.clust / nrow(mat.partition) * 100, 0)
      res.clust.rater[[2]] <- perc.rater.clust

      # find the parangon
      coord.bary.clust <- apply(res.pca$ind$coord[which(rownames(res.pca$ind$coord)%in%mat.partition[which(mat.partition$Cluster == clust), "Rater"]),], 2, mean)
      coord.rater.bary <- rbind.data.frame(res.pca$ind$coord, coord.bary.clust)
      rownames(coord.rater.bary)[nrow(coord.rater.bary)] <- "barycentre"
      mat.dist.bary.clust <- as.data.frame(as.matrix(dist(coord.rater.bary, "euclidean")))
      dist.bary.clust <- mat.dist.bary.clust[nrow(mat.dist.bary.clust), -ncol(mat.dist.bary.clust)]
      parangon.clust <- names(which.min(dist.bary.clust))
      res.clust.rater[[3]] <- parangon.clust

      # interpret the cluster with external information about the raters
      if (!is.null(id.info.rater)) {
        info.rater <- dta.sauv[id.info.rater, -id.info.stim]
        dta.info.rater <- cbind.data.frame(rownames(mat.resids), t(info.rater))
        for (j in 2:ncol(dta.info.rater)) {
          if (type.info.rater[j-1] == "cat") {
            dta.info.rater[,j] <- as.factor(dta.info.rater[, j])
          }
          if (type.info.rater[j-1] == "cont") {
            dta.info.rater[,j] <- as.numeric(dta.info.rater[, j])
          }
        }
        colnames(dta.info.rater)[1] <- "Rater"
        dta.info.rater <- merge(dta.info.rater, mat.partition, by = "Rater")
        dta.info.rater <- dta.info.rater[, -1]
        dta.info.rater$Cluster <- as.factor(dta.info.rater$Cluster)
        res.info.rater <- catdes(dta.info.rater, ncol(dta.info.rater), 0.05)
        res.clust.rater[[4]] <- list()
        if (is.null(res.info.rater$quanti[[i]])) {
          nbquanti <- 0
        } else {
          nbquanti <- nrow(res.info.rater$quanti[[i]])
        }
        if (is.null(res.info.rater$category[[i]])) {
          nbquali <- 0
        } else {
          nbquali <- nrow(res.info.rater$category[[i]])
        }
        if (nbquanti == 0 & nbquali == 0) {
          res.clust.rater[4] <- list(NULL)
        } else {
          info.rater.sup.clust <- as.data.frame(matrix(NA, nbquanti + nbquali, 3))
          info.rater.sup.clust[, 1] <- c(rownames(res.info.rater$quanti[[i]]), rownames(res.info.rater$category[[i]]))
          colnames(info.rater.sup.clust) <- c("information", "sign statistic test", "pvalue")
          for (j in 1 : nrow(info.rater.sup.clust)) {
            if (length(which(rownames(res.info.rater$quanti[[i]]) == info.rater.sup.clust[j, 1])) != 0) {
              info.rater.sup.clust[j, 3] <- res.info.rater$quanti[[i]][which(rownames(res.info.rater$quanti[[i]]) == info.rater.sup.clust[j, 1]), 6]
              if (res.info.rater$quanti[[i]][which(rownames(res.info.rater$quanti[[i]]) == info.rater.sup.clust[j, 1]), 1] > 0) {
                info.rater.sup.clust[j,2] <- "+"
              } else {
                info.rater.sup.clust[j,2] <- "-"
              }
            }
            if (length(which(rownames(res.info.rater$category[[i]]) == info.rater.sup.clust[j,1])) != 0) {
              info.rater.sup.clust[j, 3] <- res.info.rater$category[[i]][which(rownames(res.info.rater$category[[i]]) == info.rater.sup.clust[j, 1]), 4]
              if (res.info.rater$category[[i]][which(rownames(res.info.rater$category[[i]]) == info.rater.sup.clust[j, 1]), 5] > 0) {
                info.rater.sup.clust[j, 2] <- "+"
              } else {
                info.rater.sup.clust[j, 2] <- "-"
              }
            }
          }
          info.rater.sup.clust.plus <- info.rater.sup.clust[which(info.rater.sup.clust[, 2] == "+"), ]
          info.rater.sup.clust.plus <- info.rater.sup.clust.plus[order(info.rater.sup.clust.plus[, 3], decreasing = FALSE), ]
          info.rater.sup.clust.moins <- info.rater.sup.clust[which(info.rater.sup.clust[, 2] == "-"), ]
          info.rater.sup.clust.moins <- info.rater.sup.clust.moins[order(info.rater.sup.clust.moins[, 3], decreasing = TRUE), ]
          info.rater.sup.clust <- rbind.data.frame(info.rater.sup.clust.plus, info.rater.sup.clust.moins)
          res.clust.rater[[4]] <- info.rater.sup.clust
        }
      }

      names(res.clust.rater)[1 : 3] <- c("nb.raters", "percent.of.panel", "parangon")
      if (!is.null(id.info.rater)) {
        names(res.clust.rater)[4] <- "info.raters"
      }

      return(res.clust.rater)
    }
    list.charact.cluster.rater <- lapply(1 : nlevels(mat.partition$Cluster), charact.cluster.rater)
    # supplement the interpretation of the clusters with information about the stimuli
    if (!is.null(id.info.stim)) {
      melted.data.clusters <- merge(melted.data, mat.partition, by = "Rater")
      melted.data.clusters$Rater <- as.factor(melted.data.clusters$Rater)
      melted.data.clusters$Stimulus <- as.factor(melted.data.clusters$Stimulus)
      melted.data.clusters$Rating <- as.numeric(as.character(melted.data.clusters$Rating))
      melted.data.clusters$Cluster <- as.factor(melted.data.clusters$Cluster)
      charact.cluster.stim <- function (i) {

        # interpret the cluster with external information about the stimuli
        if (!is.null(id.info.rater)) {
          info.stim <- as.data.frame(dta.sauv[-id.info.rater, id.info.stim])
        } else {
          info.stim <- as.data.frame(dta.sauv[, id.info.stim])
        }
        colnames(info.stim) <- colnames(dta.sauv)[id.info.stim]
        dta.info.stim <- cbind.data.frame(rownames(dta), info.stim)
        for (j in 2 : ncol(dta.info.stim)) {
          if (type.info.stim[j - 1] == "cat") {
            dta.info.stim[,j] <- as.factor(dta.info.stim[,j])
          }
          if (type.info.stim[j - 1] == "cont") {
            dta.info.stim[,j] <- as.numeric(dta.info.stim[,j])
          }
        }
        colnames(dta.info.stim)[1] <- "Stimulus"
        melted.data.clusters.info.sup.stim <- merge(dta.info.stim, melted.data.clusters, by = "Stimulus")
        melted.data.clusters.info.sup.stim <- melted.data.clusters.info.sup.stim[, -which(colnames(melted.data.clusters.info.sup.stim)%in%c("Stimulus","Rater"))]

        info.stim.sup <- as.data.frame(matrix(NA, 1, 3))
        colnames(info.stim.sup) <- c("information", "sign statistic test", "pvalue")
        for (j in 1 : length(id.info.stim)) {
          name.Supp <- colnames(dta.sauv)[id.info.stim[j]]
          melted.data.clusters.info.sup.stim.j <- melted.data.clusters.info.sup.stim[, which(colnames(melted.data.clusters.info.sup.stim) %in% c("Rating", name.Supp, "Cluster"))]
          colnames(melted.data.clusters.info.sup.stim.j)[1] <- "Supp"
          res.AovSum <- AovSum(Rating ~ Cluster * Supp, data = melted.data.clusters.info.sup.stim.j)
          res.AovSum <- list(res.AovSum[[1]], res.AovSum[[2]])
          names(res.AovSum) <- c("GlobTest", "LocTest")
          if (res.AovSum$GlobTest["Cluster:Supp", 5] < 0.05) {
            int.clust <- as.data.frame(res.AovSum$LocTest[grep(paste("Cluster - ", i, " :", sep = ""), rownames(res.AovSum$LocTest)), ])
            if ((nrow(int.clust) == 4) & (ncol(int.clust) == 1)) {
              int.clust <- t(int.clust)
            }
            if (length(which(int.clust[, 4] < 0.05)) != 0) {
              for (k in 1 : length(which(int.clust[, 4] < 0.05))) {
                pos <- which(int.clust[, 4] < 0.05)[k]
                if (type.info.stim[j] == "cat") {
                  name.info <- paste(name.Supp, gsub(paste("Cluster - ", i, " : Supp - ", sep = ""),"", rownames(int.clust)[pos]), sep = " - ")
                }
                if (type.info.stim[j] == "cont") {
                  name.info <- name.Supp
                }
                charac.signif <- c(name.info, int.clust[pos, 1], int.clust[pos, 4])
                info.stim.sup <- rbind.data.frame(info.stim.sup, charac.signif)
              }
            }
          }
        }
        info.stim.sup <- info.stim.sup[-1, ]
        if (nrow(info.stim.sup) == 0) {
          res.clust.stim <- list(NULL)
        } else {
          info.stim.sup[which(info.stim.sup[, 2] < 0), 2] <- "-"
          info.stim.sup[which(info.stim.sup[, 2] > 0), 2] <- "+"
          info.stim.sup[, 3] <- as.numeric(info.stim.sup[, 3])
          info.stim.sup.plus <- info.stim.sup[which(info.stim.sup[, 2] == "+"), ]
          info.stim.sup.plus <- info.stim.sup.plus[order(info.stim.sup.plus[, 3], decreasing = FALSE), ]
          info.stim.sup.moins <- info.stim.sup[which(info.stim.sup[, 2] == "-"), ]
          info.stim.sup.moins <- info.stim.sup.moins[order(info.stim.sup.moins[, 3], decreasing = TRUE), ]
          info.stim.sup <- rbind.data.frame(info.stim.sup.plus, info.stim.sup.moins)
          res.clust.stim <- info.stim.sup
        }

        return(res.clust.stim)

      }
      list.charact.cluster.stim <- lapply(1 : nlevels(mat.partition$Cluster), charact.cluster.stim)
    }
    # combine the results
    for (i in 1 : nlevels(mat.partition$Cluster)) {
      res[[9]][[i]] <- list.charact.cluster.rater[[i]]
      if (!is.null(id.info.stim)) {
        res[[9]][[i]][[length(res[[9]][[i]]) + 1]] <- list.charact.cluster.stim[[i]]
        if (length(list.charact.cluster.stim[[i]]) == 1 & is.null(res[[9]][[i]]$info.stim[[1]])) {
          res[[9]][[i]][length(res[[9]][[i]])] <- list(NULL)
        }
        names(res[[9]][[i]])[length(res[[9]][[i]])] <- "info.stim"
      }
    }
    names(res[[9]]) <- 1 : nlevels(mat.partition$Cluster)
  } else {
    res[9] <- list(NULL)
  }

  # return the results
  names(res) <- c("call", "profiles.residuals", "mat.disag", "pval.dendro", "nb.clust.found", "partition", "res.plot.segment", "res.pca", "charact.clust")
  message("Clustering performed")
  options(warn = 0)
  class(res) <- c("AgreeClustCont", "list ")
  return(res)

}
