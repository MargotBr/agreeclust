\name{AgreeClustBin}
\alias{AgreeClustBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Agreement-based unsupervised clustering of binary ratings}
\description{Performs an agreement-based clustering of binary ratings to obtain a segmentation of the panel of raters. The method is able to manage covariates providing informations about the raters and the stimuli (i.e. objects assessed) to facilitate the interpretation of the clusters of raters.}
\usage{
AgreeClustBin(dta, model = "Rating ~ Rater + Stimulus", max.clust = 10, approx.null = TRUE, paral = TRUE, consol = TRUE, id.info.rater = NULL, type.info.rater = NULL, id.info.stim = NULL, type.info.stim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{
A binary matrix or data frame of dimensions S x R (S=number of stimuli, R=number of raters).}
  \item{model}{
A formula for the model. Either 'Rating~Rater+Stimulus' (default), 'Rating~Rater', or 'Rating~1'.}
  \item{max.clust}{
An integer specifying the maximum number of clusters of raters. By default, this number is fixed to 10.}
  \item{approx.null}{
A boolean indicating if the null LRT distribution should be approximated using Satterthwaite's approximation. By default, the null LRT distribution is approximated.}
  \item{paral}{
A boolean indicating if the computation of the null LRT distribution should be parallelized. By default, the computation of the null LRT distribution is parallelized on nb.cores-1 cores. During the process, a text file 'TestDendrogram_processing.txt' is created.}
  \item{consol}{
A boolean indicating if a k-means consolidation of the partition of raters should be performed. By default, the partition is consolidated.}
  \item{id.info.rater}{
A vector of integer elements composed of the identification of the lines containing the supplementary information (i.e. covariates) about the raters. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the raters.}
  \item{type.info.rater}{
A vector of character elements composed of the type of the covariates about the raters. This vector must be of the same length that id.info.rater. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the raters.}
  \item{id.info.stim}{
A vector of integer elements composed of the identification of the columns containing the supplementary information (i.e. covariates) about the stimuli. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.}
  \item{type.info.stim}{
A vector of character elements composed of the type of the covariates about the stimuli. This vector must be of the same length that id.info.stim. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{profiles.residuals}{A matrix of dimensions S x R (S=number of stimuli, R=number of raters) containing the residuals profiles of the raters obtained through the modelling of the set of binary ratings.}
\item{mat.disag}{A matrix of dimensions R x R (R=number of raters) corresponding to the dissimilarity matrix between the raters.}
\item{pval.dendro}{A vector containing the probabilities associated to the statistical test realized at each level of the dendrogram.}
\item{nb.clust.found}{An integer corresponding to the number of clusters found among the panel.}
\item{partition}{A vector representing the partition of the raters (consolidated partition if consol = TRUE).}
\item{res.plot.segment}{All the graphical results of the segmentation.}
\item{res.pca}{All the results of the PCA.}
\item{charact.clust }{The results of the description of the clusters by information describing the raters and/or the stimuli.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{Margot Brard \email{margot.brard@agrocampus-ouest.fr}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.AgreeClust}}
}
\examples{
# Small pedagogic data set
data(pedagdata)
res.pedag <- AgreeClustBin(dta = pedagdata, model = "Rating ~ Rater + Stimulus", id.info.rater = 9 : nrow(pedagdata), type.info.rater = rep("cat", 3), id.info.stim = 21 : ncol(pedagdata), type.info.stim = c(rep("cont", 4), "cat"))
res.pedag

# 'Good gesture' data set
data(ggdata)
res.gg <- AgreeClustBin(dta = ggdata, model = "Rating ~ Rater + Stimulus", id.info.rater = 40 : nrow(ggdata), type.info.rater = rep("cat", 4), id.info.stim = 73 : ncol(ggdata), type.info.stim = c(rep("cat", 3), rep("cont", 11)))
res.gg
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
