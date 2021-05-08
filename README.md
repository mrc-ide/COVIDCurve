
# COVIDCurve
<!-- badges: start -->
[![Travis build status](https://travis-ci.org/mrc-ide/COVIDCurve.svg?branch=master)](https://travis-ci.org/mrc-ide/COVIDCurve)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/COVIDCurve/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/COVIDCurve?branch=master)  
<!-- badges: end -->
  
<description>

## IMPORTANT NOTES

:warning: This package depends on the R-package [DrJacoby](https://github.com/mrc-ide/drjacoby) version 1.2.0.  


## Installation & Use
The purpose of this package is to provide a framework to fit age-specific infection fatality ratios (IFRs) in an epidemic that is using serologic data to estimate prevalence.  The model was originally designed in response to the COVID-19 pandemic (hence the name) but is broadly applicable.

Accurate inference of the IFR based on serological data is challenging due to a number of factors that can bias estimates away from the true value, including: (1) the delay between infection and death, (2) the dynamical process of seroconversion and *seroreversion*, (3) potential differences in age-specific attack rates, and (4) serological test characteristics. These delay-distributions and serologic test characteristics are expected to cause both shifts in when deaths are observed as well as when and how many true infections are observed, which will then affect the IFR. This framework is summarized in the conceptual diagram below.

The model can be briefly summarized as having to pieces: (1) an _infection curve shape_ that is informed by the shape the of the observed death curve but "thrown"" backwards in time; (2) a _pin_ or a point on the y-axis that the _shape_ must go through (_i.e._ a pin to inform the area under the curve, or cumulative infections). This pin is informed by the observed serologic data, accounting for test characteristics (_e.g._ sensitivity and specificity) as well as the time from infection to seroconversion and *seroreversion*. Accounting for seroreversion, or the waning of antibodies over time that then produces a false negative test, becomes increasingly important as epidemics progress and/or as daily infections start to contract.  

<p align="center">
<img src="https://raw.githubusercontent.com/mrc-ide/COVIDCurve/master/R_ignore/images/Fig_concept_diagram.png" width="400" height="500">
</p>

## Installation, Dependencies, and Use
Please see the vignettes for installation (and dependencies) as well as use instructions. Thank you for your interest in our package! 

R package dependencies are listed in the `DESCRIPTION` file and tested operating systems are indicated in the `travis.yml` file. 
