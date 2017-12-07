The `AgreeClust` package considers a latent class regression modeling framework for highlighting the structure of disagreement among panels of raters involved in an inquiry. On the contrary to popular approaches, the present method considers the ratings data provided by all raters when studying the structure of disagreement among the panel. More precisely, the structure of disagreement is captured through the profiles of residuals of a no-latent class regression model adjusted on the entire set of binary ratings, and can be visualized by using exploratory data analysis tools. The disagreement between two raters is then quantify in a concise way through the Euclidean distance between their respective profiles of residuals, this disagreement index being used as a basis to construct a dendrogram representing the structure of disagreement among the panel. The proper number of disagreed clusters among the panel of raters is then chosen by implementing a sequential strategy to test the significance of each K-clusters structure of disagreement.

# <span style="color: #EA485C">Installing AgreeClust</span>

To get the current development version from GitHub:

  ```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("MargotBr/AgreeClust", build_vignettes = TRUE)
library(AgreeClust)
```

# <span style="color: #EA485C">Using AgreeClust</span>

To get an overview of the functionalities of the package, read the corresponding vignette:

  ```{r eval=FALSE}
vignette("AgreeClust")
```
