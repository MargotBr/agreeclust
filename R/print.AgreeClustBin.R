print.AgreeClustBin <- function (res){

  if (!inherits(res, "AgreeClustBin")){
    stop("Non convenient data - res should be an AgreeClustBin object")
  }

  cat("** Results for the agreement-based clustering (AgreeClustBin) **\n")
  cat("\n")
  cat("The analysis was performed on", (ncol(res$call$dta) - length(res$call$id.info.stim)),
      "raters who assessed", (nrow(res$call$dta) - length(res$call$id.info.rater)), "stimuli\n")
  cat("The results are available in the following objects:\n\n")
  res.desc <- array("", c(9, 2), list(1 : 9, c("name", "description")))
  res.desc[1, ] <- c("$call", "arguments used in the AgreeClust function")
  res.desc[2, ] <- c("$profiles.residuals", "matrix of profiles of deviance residuals")
  res.desc[3, ] <- c("$mat.disag", "disagreement matrix")
  res.desc[4, ] <- c("$pval.dendro", "p-values in the dendrogram")
  res.desc[5, ] <- c("$nb.clust.found", "number of clusters of raters found")
  res.desc[6, ] <- c("$partition", "partition of raters found (consolidated or not)")
  res.desc[7, ] <- c("$res.plot.segment", "graphical results of the clustering (not needed)")
  res.desc[8, ] <- c("$res.pca", "PCA results of the multidimensional analysis of the structure of disagreement")
  res.desc[9, ] <- c("$charact.clust", "description of the clusters of raters")
  print(res.desc)

}
