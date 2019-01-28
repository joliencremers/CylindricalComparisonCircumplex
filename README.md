# CylindricalComparisonCircumplex

This repository contain the data archive for the paper:

*Regression models for Cylindrical data in Psychology*

The folder `Manuscript` contains information to reproduce the paper manuscript and analyses. 

The file `ManuscriptMarkdownCV.Rmd` contains code to reproduce all elements of the manuscript in `ManuscriptMarkdownCV.pdf`. Supporting functions to perform the analyses are found in the folder R-code. This folder contains code for the posterior samplers and optimization routines of the four cylindrical models in the R-scripts:

* `Posterior Sampling CL-PN.R`
* `Posterior Sampling CL-GPN.R`
* `Abe-Ley optimization.R`
* `Posterior Sampling Joint GPN-SSN.R`

These scripts are annotated to describe their use. Note that these functions are written for the specific example in the current paper (1 linear outcome, 1 circular outcome, 1 covariate). Examples of their use for the analyses in the paper are found in the file `ManuscriptMarkdownCV.Rmd`.

The file `Sample Abe-Ley.R` contains code for sampling from the Abe-Ley distribution.
