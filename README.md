# metaFlexB
An initial version of the R package metaFlexB. 

For Mac User:

-Install Xcode

-Check gcc version, and install gfortran binary (gfortran 8.2) if your device compiles the package fails.

-More details: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/


Alternative way to access the function:

Step1: Download src folder entirely

Step2: Install Rcpp,RcppArmadillo and RcppDist packages.

Step3: Let your R locate to src folder.

Step4: Use the following R code

Library(Rcpp)
SourceCpp('META.cpp')

Then, you should see the function named main_draw() in the section of "Environment".
