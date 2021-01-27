## dev history
usethis::use_build_ignore("dev")
usethis::use_pkgdown()
usethis::use_mit_license('ThinkR')
usethis::use_namespace(roxygen = TRUE)
usethis::use_pipe()

## Vignettes
usethis::use_vignette("a-overview-methodology")
usethis::use_vignette("b-data-needed")
usethis::use_vignette("c-binary-ratings")
usethis::use_vignette("d-continuous-ratings")

# Use pkgdown
usethis::use_pkgdown()

## Description
attachment::att_amend_desc()

## Pkgdown
pkgdown::build_site()

## Check
devtools::check()
