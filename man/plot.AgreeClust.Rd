\name{plot.AgreeClust}
\alias{plot.AgreeClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of the structure of disagreement among the panel of raters}
\description{Draws the graphs representing the structure of disagreement among the set of ratings. Tunning parameters such as the color of the clusters or the axis of PCA to be plotted can be specified.}
\usage{
plot.AgreeClust(res, choice = "all", interact = FALSE, col.clust = NULL, axis = c(1, 2), name.rater = "rater", ext.dev.Rstudio = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
An object of class AgreeClust returned by the functions of the package.}
  \item{choice}{
A character element specifying the graphs to be plotted ('seg' for the representation of the segmentation process (dendrogram + partitioning if a consolidation process has been implemented), 'mul' for the multidimensional representation of the structure of disagreement (PCA), 'all' for both representations. By default, both representations are plotted.)}
  \item{interact}{
A boolean specifying if the graphical outputs should be interactive (with plotly) or not.}
  \item{col.clust}{
A vector with as many color as clusters of raters. Colors can be specified by a character such as "blue" or by a hexadecimal code such as "#0D10D5". All hexadecimal codes can be found at \url{http://htmlcolorcodes.com}. By default, no colors are specified and default colors are used.}
  \item{axis}{
A length 2 numeric vector specifying the PCA components to be plotted. By default, the first 2 components are plotted.}
 \item{name.rater}{
A character element indicating how the raters should be named (e.g. "rater", "participant", "psychologists")}
\item{ext.dev.Rstudio}{
A boolean specifying if the graphical outputs should be plotted in the Rstudio plot pane or not.}
}
\value{
Returns the graphical representations.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Margot Brard \email{margot.brard@agrocampus-ouest.fr}}
\seealso{
\code{\link{AgreeClustBin} \link{AgreeClustCont}}
}
\examples{
# Small pedagogic data set
data(pedagdataBin)
res.pedag <- AgreeClustBin(dta = pedagdataBin, model = "Rating ~ Rater + Stimulus", id.info.rater = 9 : nrow(pedagdataBin), type.info.rater = rep("cat", 3), id.info.stim = 21 : ncol(pedagdataBin), type.info.stim = c(rep("cont", 4), "cat"))
plot.AgreeClust(res.pedag, col.clust = c("blue", "green"))
plot.AgreeClust(res.pedag, choice = "mul", interact = TRUE, col.clust = c("blue", "green"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
