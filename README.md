# metaFlexB
This is a README file of the R package _metaFlexB_. In our paper, we develop new Bayesian procedures for estimating and testing the overall treatment effect and inter-study heterogeneity in random-effects meta-analysis of rare binary events.

## Installation of the package

To install our package, you may execute the following codes:

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("chriszhangm/metaFlexB")
library(metaFlexB)
```
For Mac Users who cannot compile the code, please refer [this answer](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).

## A example of using _metaFlexB_

We show a toy example to apply the function `main_draw`, which produces random posterior draws of global parameters of interest.

## Notes
Alternative way to access the function:

Step1: Download src folder entirely

Step2: Install Rcpp,RcppArmadillo and RcppDist packages.

Step3: Let your R locate to src folder.

Step4: Use the following R code

Library(Rcpp)
SourceCpp('META.cpp')

Then, you should see the function named main_draw() in the section of "Environment".
