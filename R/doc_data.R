#' Pedagogic data set with binary ratings
#' @name binary_data_for_example
#' @docType data
#' @author Margot Brard \email{margot.brard@@yahoo.fr}
#' @description This data set is a small artificial pedagogic data set.
#' @format An initial data frame of dimensions 8 x 20 containing the binary ratings on 8 stimuli provided by a panel of 20 raters. The initial data frame is supplemented by 3 covariates providing information about the raters (Gender, Age, Freq.use; the first 2 categorical and the last one continuous); and by 5 covariates providing information about the stimuli (Citrus.fruits, Vanilla, Wood, Lotus, Packaging; the first 4 continuous and the last one categorical).
#' @examples
#' data("continuous_data_for_example")
#'
"binary_data_for_example"

#' Good gesture data set
#' @name goodgesture
#' @docType data
#' @author Margot Brard \email{margot.brard@@yahoo.fr}
#' @description In the 'good gesture' study, a panel of 72 participants was recruited and asked to sort each of 39 videos of culinary gestures into one of two predefined categories: 'videos representative of the concept of good gesture' and 'videos not representative of the concept of good gesture'. The experiment led to 72 profiles of 39-dimensional binary ratings. Additionally, external information about the raters were collected, such as their gender (F or M), their age bracket (18-20, 21-30, 31-40, 41-50, 51-60, or 71-80), their frequency of eating oysters (never, occasionally, or regularly) and their level of expertise for opening oysters (beginner, intermediate, or expert). Data were also supplemented by fourteen covariates providing information about the stimuli: the gender of the person who realized the gesture (F or M), the frequency of oysters opening of the person who realized the gesture (never, once a year, once a month, and once a week), the level of expertise for opening oysters of the person who realized the gesture (beginner, intermediate, or expert) and eleven numeric evaluations of several defects observed in the video, like the damages on the shell, the dangerousness of the use of the knife, the imprecision of the gesture, and so forth.
#' @format An initial data frame of dimensions 39 x 72 containing the binary ratings on the 39 stimuli provided by the panel of 72 raters. The initial data frame is supplemented by the 4 covariates providing information about the raters (all categorical); and by the 14 covariates providing information about the stimuli (the first 4 categorical and the last 11 continuous).
#' @examples
#' data("goodgesture")
"goodgesture"
