\name{AgreeClust-package}
\alias{AgreeClust-package}
\alias{AgreeClust}
\docType{package}
\title{
\packageTitle{AgreeClust}
}
\description{
\packageDescription{AgreeClust}
}

\author{
\packageAuthor{AgreeClust}

Maintainer: \packageMaintainer{AgreeClust}
}

\keyword{ package }
\seealso{
\code{\link{AgreeClustBin}}
\code{\link{plot.AgreeClust}}
}
\examples{
# Small pedagogic data set
data(pedagdataBin)
res.pedag <- AgreeClustBin(dta = pedagdataBin, model = "Rating ~ Rater + Stimulus", id.info.rater = 9 : nrow(pedagdataBin), type.info.rater = rep("cat", 3), id.info.stim = 21 : ncol(pedagdataBin), type.info.stim = c(rep("cont", 4), "cat"))
plot.AgreeClust(res.pedag, col.clust = c("blue", "green"))
}
