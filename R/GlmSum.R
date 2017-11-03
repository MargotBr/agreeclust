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
