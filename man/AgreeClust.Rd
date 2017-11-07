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
Margot Brard, Sebastien Le, David Causeur

Maintainer: Margot Brard \email{margot.brard@agrocampus-ouest.fr}
}

\keyword{ package }
\seealso{
\code{\link{AgreeClustBin}}
\code{\link{AgreeClustCont}}
\code{\link{plot.AgreeClust}}
}
\examples{
# Small binary pedagogic data set
data(pedagdataBin)
res.pedagBin <- AgreeClustBin(dta = pedagdataBin, model = "Rating ~ Rater + Stimulus", id.info.rater = 9 : nrow(pedagdataBin), type.info.rater = rep("cat", 3), id.info.stim = 21 : ncol(pedagdataBin), type.info.stim = c(rep("cont", 4), "cat"))
plot.AgreeClust(res.pedagBin, col.clust = c("blue", "green"))

# Small continuous pedagogic data set
data(pedagdataCont)
res.pedagCont <- AgreeClustCont(dta = pedagdataBin, model = "Rating ~ Rater + Stimulus", id.info.rater = 9 : nrow(pedagdataCont), type.info.rater = rep("cat", 3), id.info.stim = 21 : ncol(pedagdataCont), type.info.stim = c(rep("cont", 4), "cat"))
plot.AgreeClust(res.pedagCont, col.clust = c("blue", "green"))
}
