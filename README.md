
<!-- README.md is generated from README.Rmd. Please edit that file -->

# agreeclust <img src='man/figures/hex-agreeclust-white.png' align="right" height="139" />

The `{agreeclust}` package considers a latent class regression modeling
framework for highlighting the structure of disagreement among panels of
raters involved in an inquiry. On the contrary to popular approaches,
the present method considers the ratings data provided by all raters
when studying the structure of disagreement among the panel. More
precisely, the structure of disagreement is captured through the
profiles of residuals of a no-latent class regression model adjusted on
the entire set of binary ratings, and can be visualized by using
exploratory data analysis tools. The disagreement between two raters is
then quantify in a concise way through the Euclidean distance between
their respective profiles of residuals, this disagreement index being
used as a basis to construct a dendrogram representing the structure of
disagreement among the panel. The proper number of disagreed clusters
among the panel of raters is then chosen by implementing a sequential
strategy to test the significance of each \(K\)-clusters structure of
disagreement.

To get the current version from GitHub:

``` r
if(!requireNamespace("devtools")){install.packages("devtools")}
devtools::install_github("MargotBr/agreeclust", build_vignettes = TRUE)
library(agreeclust)
```

# Repex

``` r
library(agreeclust)
data(binary_data_for_example)
res_pedag <- get_agreeclust_bin(dta = binary_data_for_example,
                                id_info_rater = 9 : nrow(binary_data_for_example),
                                type_info_rater = c(rep("cat", 2), "cont"),
                                id_info_stim = 21 : ncol(binary_data_for_example),
                                type_info_stim = c(rep("cont", 4), "cat"),
                                paral_null = FALSE,
                                graph = FALSE)
names(res_pedag)
#> [1] "call"               "profiles_residuals" "mat_disag"         
#> [4] "pval_dendro"        "nb_clust_found"     "partition"         
#> [7] "res_plot_segment"   "res_pca"            "charact_clust"
plot_agreeclust(res_pedag)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r
plot_agreeclust(res_pedag, interact = TRUE)
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />
