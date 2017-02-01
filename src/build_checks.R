### From a great vignette
### https://cran.r-project.org/web/packages/Rd2roxygen/vignettes/Rd2roxygen.html

library(Rd2roxygen)
library(roxygen2)

## This converts the original package to roxygen2 format
## Rd2roxygen('enaR')

## Process the package to generate NAMESPACE and help files
roxygenize('../../enaR')

## This builds the package from Roxygen
## roxygen_and_build('enaR')
## test <- rab('../../enaR')

## This builds and checks enaR

