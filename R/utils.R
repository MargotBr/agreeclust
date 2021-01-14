#' Function remove_outliers
GlmSum <- function(formula, data, na.action = na.omit) {

  old.contr <- options()$contrasts
  on.exit(options(contrasts = old.contr))
  options(contrasts = c("contr.sum", "contr.sum"))
  don <- data

  modele <- glm(formula, data = don, family = binomial)
  test.global <- car::Anova(modele)
  test.local <- summary.glm(modele)$coef

  cov.mat <- vcov(modele)
  facteurs <- rownames(attr(modele$terms, "factors"))[-1]
  interact <- NULL
  if (length(colnames(attr(modele$terms, "factors"))) > length(facteurs)) {
    interact <- colnames(attr(modele$terms, "factors"))[-(1:length(facteurs))]
  }
  niveau <- list()
  for (i in 1:length(facteurs)) {
    if (is.factor(don[, facteurs[i]])) {
      niveau[[i]] <- paste(facteurs[i], levels(don[, facteurs[i]]), sep = " - ")
    } else {
      niveau[[i]] <- facteurs[i]
    }
  }
  res <- test.local[c(1, 1), ]

  iinit <- 2
  for (i in 1:length(facteurs)) {
    old.rownames <- rownames(res)
    if (is.factor(don[, facteurs[i]])) {
      indices <- iinit:(iinit + nlevels(don[, facteurs[i]]) - 2)
      coeff.dern.mod <- -sum(test.local[indices, 1])
      std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
      zval.dern.mod <- coeff.dern.mod/std.dern.mod
      pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
      dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
      res <- rbind(res, test.local[indices, ], dern.mod)
      rownames(res) <- c(old.rownames, niveau[[i]])
      iinit <- iinit + nlevels(don[, facteurs[i]]) - 1
    } else {
      indices <- iinit
      res <- rbind(res, test.local[indices, ])
      rownames(res) <- c(old.rownames, niveau[[i]])
      iinit <- iinit + 1
    }
  }
  res <- res[-1, ]

  if (!is.null(interact)) {
    for (k in 1:length(interact)) {
      fact.int <- rownames(attr(modele$terms, "factors"))[which(attr(modele$terms, "factors")[, interact[k]] == 1)]
      old.rownames <- rownames(res)
      fact1 <- fact.int[1]
      fact2 <- fact.int[2]
      iinit0 <- iinit
      if ((is.factor(don[, fact1])) & (is.factor(don[,fact2]))) {
        for (l in 1:(nlevels(don[, fact2]) - 1)) {
          indices <- iinit:(iinit + (nlevels(don[, fact1]) -2))
          coeff.dern.mod <- -sum(test.local[indices, 1])
          std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
          zval.dern.mod <- coeff.dern.mod/std.dern.mod
          pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
          dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
          res <- rbind(res, test.local[indices, ], dern.mod)
          iinit <- iinit + (nlevels(don[, fact1]) - 1)
        }
        iinit = iinit0
        for (l in 1:(nlevels(don[, fact1]) - 1)) {
          indices <- iinit + (nlevels(don[, fact1]) - 1) * (0:(nlevels(don[, fact2]) - 2))
          coeff.dern.mod <- -sum(test.local[indices, 1])
          std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
          zval.dern.mod <- coeff.dern.mod/std.dern.mod
          pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
          dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
          res <- rbind(res, dern.mod)
          iinit <- iinit + 1
        }
        indices <- iinit0:(iinit0 + (nlevels(don[, fact1]) - 1) * (nlevels(don[, fact2]) - 1) - 1)
        coeff.dern.mod <- sum(test.local[indices, 1])
        std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
        zval.dern.mod <- coeff.dern.mod/std.dern.mod
        pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
        dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
        res <- rbind(res, dern.mod)
        iinit <- iinit0 + (nlevels(don[, fact1]) - 1) * (nlevels(don[, fact2]) - 1)
        nom <- old.rownames
        aa <- paste(fact2, levels(don[, fact2]), sep = " - ")
        for (i in 1:length(aa)) nom <- c(nom, paste(paste(fact1, levels(don[, fact1]), sep = " - "), aa[i], sep = " : "))
      }
      if ((is.factor(don[, fact1])) & (!is.factor(don[,fact2]))) {
        indices <- iinit:(iinit + (nlevels(don[, fact1]) -  2))
        coeff.dern.mod <- -sum(test.local[indices, 1])
        std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
        zval.dern.mod <- coeff.dern.mod/std.dern.mod
        pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
        dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
        res <- rbind(res, test.local[indices, ], dern.mod)
        iinit <- iinit + (nlevels(don[, fact1]) - 1)
        nom <- c(old.rownames, paste(paste(fact1, levels(don[, fact1]), sep = " - "), fact2, sep = " : "))
      }
      if ((!is.factor(don[, fact1])) & (is.factor(don[,fact2]))) {
        indices <- iinit:(iinit + (nlevels(don[, fact2]) - 2))
        coeff.dern.mod <- -sum(test.local[indices, 1])
        std.dern.mod <- sqrt(sum(cov.mat[indices,indices]))
        zval.dern.mod <- coeff.dern.mod/std.dern.mod
        pvalue.dern.mod <- 2 * pnorm(-abs(zval.dern.mod))
        dern.mod <- c(coeff.dern.mod, std.dern.mod, zval.dern.mod, pvalue.dern.mod)
        res <- rbind(res, test.local[indices, ], dern.mod)
        iinit <- iinit + (nlevels(don[, fact2]) - 1)
        nom <- c(old.rownames, paste(paste(fact2, levels(don[, fact2]), sep = " - "), fact1, sep = " : "))
      }
      if ((!is.factor(don[, fact1])) & (!is.factor(don[,fact2]))) {
        indices <- iinit
        res <- rbind(res, test.local[indices, ])
        iinit <- iinit + 1
        nom <- c(old.rownames, paste(fact1, fact2, sep = " : "))
      }
      rownames(res) <- nom
    }
  }

  result <- list(GlobTest = test.global, LocTest = res)
  class(result) <- "GlmSum"
  options(contrasts = old.contr)
  return(result)
}

#' Function remove_outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 10 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#' Function compute_dist_null
compute_dist_null <- function(j, list.data.null, K, model, model2) {

  # Null dendrogram
  message(paste0("Comparison test: ", (K + 1), "-latent class model VS " , (K), "-latent class model | Progress: ", round(j / length(list.data.null) * 100, 1), "%"))
  data.null <- as.matrix(list.data.null[[j]])
  melted.data.null <- reshape2::melt(data.null)
  melted.data.null <- melted.data.null[, c(2,1,3)]
  colnames(melted.data.null) <- c("rater", "stimulus", "rating")
  melted.data.null$rater <- as.factor(melted.data.null$rater)
  melted.data.null$stimulus <- as.factor(melted.data.null$stimulus)
  melted.data.null$rating <- as.factor(melted.data.null$rating)
  mod.noLC.null <- glm(model, data = melted.data.null, family = binomial)
  mat.resids.null <- matrix(residuals(mod.noLC.null, type = "deviance"), nrow(data.null), ncol(data.null))
  dimnames(mat.resids.null) <- list(rownames(data.null), colnames(data.null))
  mat.resids.null <- t(as.data.frame(mat.resids.null))
  disagmat.null <- dist(mat.resids.null, method = "euclidean")
  disagmat.null <- as.dist(disagmat.null)
  dendrogram.null <- stats::hclust(disagmat.null, method = "ward.D2")

  # Null K-latent class model
  partition.KLC.null <- cutree(dendrogram.null, k = K)
  partition.KLC.null <- cbind.data.frame(names(partition.KLC.null), partition.KLC.null)
  colnames(partition.KLC.null) <- c("rater", "cluster")
  clustered.data.KLC.null <- merge(melted.data.null, partition.KLC.null, by = "rater")
  clustered.data.KLC.null$cluster <- as.factor(clustered.data.KLC.null$cluster)
  if (K == 1) {
    mod.KLC.null <- glm(model, data = clustered.data.KLC.null, family = binomial)
  } else {
    mod.KLC.null <- glm(model2, data = clustered.data.KLC.null, family = binomial)
    X.KLC.null <- model.matrix(mod.KLC.null)
    X.KLC.null <- X.KLC.null[, !is.na(coef(mod.KLC.null))]
    mod.KLC.null <- glm(rating ~ ., data = data.frame(rating = clustered.data.KLC.null$rating, X.KLC.null[, -1]), family = binomial)
  }

  # Null (K+1)-latent class model
  partition.Kplus1LC.null <- cutree(dendrogram.null, k = K + 1)
  partition.Kplus1LC.null <- cbind.data.frame(names(partition.Kplus1LC.null), partition.Kplus1LC.null)
  colnames(partition.Kplus1LC.null) <- c("rater", "cluster")
  clustered.data.Kplus1LC.null <- merge(melted.data.null, partition.Kplus1LC.null, by = "rater")
  clustered.data.Kplus1LC.null$cluster <- as.factor(clustered.data.Kplus1LC.null$cluster)
  mod.Kplus1LC.null <- glm(model2, data = clustered.data.Kplus1LC.null, family = binomial)
  X.Kplus1LC.null <- model.matrix(mod.Kplus1LC.null)
  X.Kplus1LC.null <- X.Kplus1LC.null[, !is.na(coef(mod.Kplus1LC.null))]
  mod.Kplus1LC.null <- glm(rating ~ ., data = data.frame(rating = clustered.data.Kplus1LC.null$rating, X.Kplus1LC.null[, -1]), family = binomial)

  # Compute the null LRT statistic
  lrt.H0 <- deviance(mod.KLC.null) - deviance(mod.Kplus1LC.null)
  return(lrt.H0)

}

#' Function compute_pval
compute_pval <- function(K, dendrogram, melted.data, model, model2, approx_null, dta, paral_null) {

  if(control.break) {

    # K-latent class model
    partition.KLC <- cutree(dendrogram, k = K)
    partition.KLC <- cbind.data.frame(names(partition.KLC), partition.KLC)
    colnames(partition.KLC) <- c("rater", "cluster")
    clustered.data.KLC <- merge(melted.data, partition.KLC, by = "rater")
    clustered.data.KLC$cluster <- as.factor(clustered.data.KLC$cluster)
    if (K == 1) {
      mod.KLC <- glm(model, data = clustered.data.KLC, family = binomial)
    } else {
      mod.KLC <- glm(model2, data = clustered.data.KLC, family = binomial)
      X.KLC <- model.matrix(mod.KLC)
      X.KLC <- X.KLC[, !is.na(coef(mod.KLC))]
      mod.KLC <- glm(rating ~ ., data = data.frame(rating = clustered.data.KLC$rating, X.KLC[, -1]), family = binomial)
    }

    # (K+1)-latent class model
    partition.Kplus1LC <- cutree(dendrogram, k = K + 1)
    partition.Kplus1LC <- cbind.data.frame(names(partition.Kplus1LC), partition.Kplus1LC)
    colnames(partition.Kplus1LC) <- c("rater", "cluster")
    clustered.data.Kplus1LC <- merge(melted.data, partition.Kplus1LC, by = "rater")
    clustered.data.Kplus1LC$cluster <- as.factor(clustered.data.Kplus1LC$cluster)
    mod.Kplus1LC <- glm(model2, data = clustered.data.Kplus1LC, family = binomial)
    X.Kplus1LC <- model.matrix(mod.Kplus1LC)
    X.Kplus1LC <- X.Kplus1LC[,!is.na(coef(mod.Kplus1LC))]
    mod.Kplus1LC <- glm(rating~., data=data.frame(rating = clustered.data.Kplus1LC$rating, X.Kplus1LC[, -1]), family = binomial)

    # Compute the LRT statistic
    lrt.obs <- deviance(mod.KLC) - deviance(mod.Kplus1LC)

    # Generate the null data sets
    proba.null <- predict(mod.KLC, type = "response")
    if (approx_null == FALSE) {
      nb.simul.null <- 1000
    } else if (approx_null == TRUE) {
      nb.simul.null <- 250
    }
    list.data.null <- lapply(1:nb.simul.null, function(l, proba.null) {
      as.numeric(runif(length(proba.null)) <= proba.null) }, proba.null = proba.null)
    list.data.null <- lapply(list.data.null, function(vec, k, dnames) matrix(vec, nrow = k, dimnames = dnames), k = nlevels(melted.data$stimulus), dnames = dimnames(dta))

    # Compute p-value
    if (paral_null == FALSE) {
      list.lrt.H0 <- lapply(1 : nb.simul.null, compute_dist_null, list.data.null = list.data.null, K = K, model = model, model2 = model2)
    } else if (paral_null == TRUE) {
      environment(compute_dist_null) <- .GlobalEnv
      nb.cores <- parallel::detectCores()
      cl <- parallel::makeCluster(nb.cores - 1, outfile = "TestDendrogram_processing.txt", setup_strategy = "sequential")
      melt_copy <- reshape2::melt
      parallel::clusterExport(cl, varlist = c("melt_copy", "model", "model2", "list.data.null", "K"), environment())
      list.lrt.H0 <- parallel::clusterApply(cl, 1 : nb.simul.null, compute_dist_null, list.data.null = list.data.null, K = K)
      parallel::stopCluster(cl)
      if ("TestDendrogram_processing.txt" %in% list.files()) {
        file.remove("TestDendrogram_processing.txt")
      }
    }
    lrt.H0 <- unlist(list.lrt.H0)
    if (approx_null == TRUE) {
      if (mean(lrt.H0) == 0 & var(lrt.H0) == 0) {
        lrt.H0 <- rep(0, 1000)
      } else {
        lrt.H0 <- remove_outliers(lrt.H0)
        if (length(which(is.na(lrt.H0))) != 0) {
          lrt.H0 <- lrt.H0[-which(is.na(lrt.H0))]
        }
        nu <- (mean(lrt.H0)^2) / ((1/2) * var(lrt.H0))
        c <- mean(lrt.H0) / nu
        lrt.H0 <- rep(0, 1000)
        for (i in 1 : length(lrt.H0)) {
          lrt.H0[i] = c * rchisq(1, nu, ncp = 0)
        }
      }
    }
    pval <- length(which(lrt.obs <= lrt.H0)) / 1000
    if (pval >= 0.05) {
      control.break <<- FALSE
    }
    return(pval)

  } else {
    return(NULL)
  }
}

#' Function charact_cluster_rater
charact_cluster_rater <- function(i, mat.partition, res.pca, id_info_rater, dta.sauv, mat.resids, type_info_rater, id_info_stim) {

  res.clust.rater <- list()
  clust <- levels(mat.partition$cluster)[i]

  # calculate the number of raters
  nb.rater.clust <- length(which(mat.partition$cluster == clust))
  res.clust.rater[[1]] <- nb.rater.clust

  # calculate the percentage of raters
  perc.rater.clust <- round(nb.rater.clust / nrow(mat.partition) * 100, 0)
  res.clust.rater[[2]] <- perc.rater.clust

  # find the parangon
  coord.bary.clust <- apply(res.pca$ind$coord[which(rownames(res.pca$ind$coord)%in%mat.partition[which(mat.partition$cluster == clust), "rater"]),], 2, mean)
  coord.rater.bary <- rbind.data.frame(res.pca$ind$coord, coord.bary.clust)
  rownames(coord.rater.bary)[nrow(coord.rater.bary)] <- "barycentre"
  mat.dist.bary.clust <- as.data.frame(as.matrix(dist(coord.rater.bary, "euclidean")))
  dist.bary.clust <- mat.dist.bary.clust[nrow(mat.dist.bary.clust), -ncol(mat.dist.bary.clust)]
  parangon.clust <- names(which.min(dist.bary.clust))
  res.clust.rater[[3]] <- parangon.clust

  # interpret the cluster with external information about the raters
  if (!is.null(id_info_rater)) {
    info.rater <- dta.sauv[id_info_rater, -id_info_stim]
    dta.info.rater <- cbind.data.frame(rownames(mat.resids), t(info.rater))
    for (j in 2:ncol(dta.info.rater)) {
      if (type_info_rater[j-1] == "cat") {
        dta.info.rater[,j] <- as.factor(dta.info.rater[, j])
      }
      if (type_info_rater[j-1] == "cont") {
        dta.info.rater[,j] <- as.numeric(dta.info.rater[, j])
      }
    }
    colnames(dta.info.rater)[1] <- "rater"
    dta.info.rater <- merge(dta.info.rater, mat.partition, by = "rater")
    dta.info.rater <- dta.info.rater[, -1]
    dta.info.rater$cluster <- as.factor(dta.info.rater$cluster)
    res.info.rater <- FactoMineR::catdes(dta.info.rater, ncol(dta.info.rater), 0.05)
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
  if (!is.null(id_info_rater)) {
    names(res.clust.rater)[4] <- "info.raters"
  }

  return(res.clust.rater)
}

#' Function charact_cluster_stim
charact_cluster_stim <- function(i, mat.partition, id_info_rater, id_info_stim, dta.sauv, dta, type_info_stim, melted.data.clusters) {

  # interpret the cluster with external information about the stimuli
  if (!is.null(id_info_rater)) {
    info.stim <- as.data.frame(dta.sauv[-id_info_rater, id_info_stim])
  } else {
    info.stim <- as.data.frame(dta.sauv[, id_info_stim])
  }
  colnames(info.stim) <- colnames(dta.sauv)[id_info_stim]
  dta.info.stim <- cbind.data.frame(rownames(dta), info.stim)
  for (j in 2 : ncol(dta.info.stim)) {
    if (type_info_stim[j - 1] == "cat") {
      dta.info.stim[,j] <- as.factor(dta.info.stim[,j])
    }
    if (type_info_stim[j - 1] == "cont") {
      dta.info.stim[,j] <- as.numeric(dta.info.stim[,j])
    }
  }
  colnames(dta.info.stim)[1] <- "stimulus"
  melted.data.clusters.info.sup.stim <- merge(dta.info.stim, melted.data.clusters, by = "stimulus")
  melted.data.clusters.info.sup.stim <- melted.data.clusters.info.sup.stim[, -which(colnames(melted.data.clusters.info.sup.stim)%in%c("stimulus","rater"))]

  info.stim.sup <- as.data.frame(matrix(NA, 1, 3))
  colnames(info.stim.sup) <- c("information", "sign statistic test", "pvalue")
  for (j in 1 : length(id_info_stim)) {
    name.Supp <- colnames(dta.sauv)[id_info_stim[j]]
    melted.data.clusters.info.sup.stim.j <- melted.data.clusters.info.sup.stim[, which(colnames(melted.data.clusters.info.sup.stim) %in% c("rating", name.Supp, "cluster"))]
    colnames(melted.data.clusters.info.sup.stim.j)[1] <- "Supp"
    res.GlmSum <- GlmSum(rating ~ cluster * Supp, data = melted.data.clusters.info.sup.stim.j)
    if (res.GlmSum$GlobTest["cluster:Supp", 3] < 0.05) {
      int.clust <- as.data.frame(res.GlmSum$LocTest[grep(paste("cluster - ", i, " :", sep = ""), rownames(res.GlmSum$LocTest)), ])
      if ((nrow(int.clust) == 4) & (ncol(int.clust) == 1)) {
        int.clust <- t(int.clust)
      }
      if (length(which(int.clust[, 4] < 0.05)) != 0) {
        for (k in 1 : length(which(int.clust[, 4] < 0.05))) {
          pos <- which(int.clust[, 4] < 0.05)[k]
          if (type_info_stim[j] == "cat") {
            name.info <- paste(name.Supp, gsub(paste("cluster - ", i, " : Supp - ", sep = ""),"", rownames(int.clust)[pos]), sep = " - ")
          }
          if (type_info_stim[j] == "cont") {
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
    info.stim.sup.plus <- info.stim.sup[which(info.stim.sup[, 2] == "+"), ]
    info.stim.sup.plus <- info.stim.sup.plus[order(info.stim.sup.plus[, 3], decreasing = FALSE), ]
    info.stim.sup.moins <- info.stim.sup[which(info.stim.sup[, 2] == "-"), ]
    info.stim.sup.moins <- info.stim.sup.moins[order(info.stim.sup.moins[, 3], decreasing = TRUE), ]
    info.stim.sup <- rbind.data.frame(info.stim.sup.plus, info.stim.sup.moins)
    res.clust.stim <- info.stim.sup
  }

  return(res.clust.stim)

}

#' Function get.legend
get.legend <- function(plot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Function simpleCap
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}
