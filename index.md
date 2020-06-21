# COVIDCurve 

 <!-- badges: start -->
[![Travis build status](https://travis-ci.org/mrc-ide/COVIDCurve.svg?branch=master)](https://travis-ci.org/mrc-ide/COVIDCurve)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/COVIDCurve/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/COVIDCurve?branch=master)  
<!-- badges: end -->
<br/>


## Introduction
The purpose of this "package" is to provide a framework to fit age-specific infection fatality rates (IFRs) or the number of individuals infected in an epidemic under several key assumptions (see below). 

### Overview of the `COVIDCurve` Model 
We assume that the only data avaialble is the cumulative number of deaths on a given day during an epidemic. In addition, we assume that the epidemic incidence curve can be caputred with a cubic spline function. We then "cast" our incidence curve on to delay of disease-onset to death curve that we assume is gamma distributed. 

<br/>
<img src="https://raw.githubusercontent.com/mrc-ide/COVIDCurve/master/R_ignore/images/curve_example.png" height="200px" width="300px" align="center" />
<br/>

As a result, the expected number of deaths on a given day is a function of the total number of individuals infected on that day accounting for the offset of time infection to (potential) death.   
  
Given that there are two (joint) features to this function, one can either put a strong prior on the: (1) IFR and infer the incidence curve; or (2) Incidence curve and infer the IFR. 
