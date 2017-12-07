\name{goodgesture}
\alias{goodgesture}
\docType{data}
\title{
 Good gesture data set
}
\description{
In the 'good gesture' study, a panel of 72 participants was recruited and asked to sort each of 39 videos of culinary gestures into one of two predefined categories: 'videos representative of the concept of good gesture' and 'videos not representative of the concept of good gesture'. The experiment led to 72 profiles of 39-dimensional binary ratings. Additionally, external information about the raters were collected, such as their gender (F or M), their age bracket (18-20, 21-30, 31-40, 41-50, 51-60, or 71-80), their frequency of eating oysters (never, occasionally, or regularly) and their level of expertise for opening oysters (beginner, intermediate, or expert). Data were also supplemented by fourteen covariates providing information about the stimuli: the gender of the person who realized the gesture (F or M), the frequency of oysters opening of the person who realized the gesture (never, once a year, once a month, and once a week), the level of expertise for opening oysters of the person who realized the gesture (beginner, intermediate, or expert) and eleven numeric evaluations of several defects observed in the video, like the damages on the shell, the dangerousness of the use of the knife, the imprecision of the gesture, and so forth.
}
\usage{data("goodgesture")}
\format{
An initial data frame of dimensions 39 x 72 containing the binary ratings on the 39 stimuli provided by the panel of 72 raters. The initial data frame is supplemented by the 4 covariates providing information about the raters (all categorical); and by the 14 covariates providing information about the stimuli (the first 4 categorical and the last 11 continuous).
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\author{Margot Brard \email{margot.brard@agrocampus-ouest.fr}}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(goodgesture)
}
\keyword{datasets}
