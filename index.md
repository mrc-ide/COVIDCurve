# Curve Aware

 <!-- badges: start -->
[![Travis build status](https://travis-ci.org/nickbrazeau/CurveAware.svg?branch=master)](https://travis-ci.org/nickbrazeau/CurveAware)
[![Codecov test coverage](https://codecov.io/gh/nickbrazeau/CurveAware/branch/master/graph/badge.svg)](https://codecov.io/gh/nickbrazeau/CurveAware?branch=master)  
<!-- badges: end -->
<br/>


## Introduction
The purpose of this "package" is to provide a framework to fit age-specific infection fatality rates (IFRs) or the initial number of individuals infected (*I0*) in an epidemic under several key assumptions (see below). We fit our model using Metropolis-Coupled MCMC with the [`drjacoby`](https://mrc-ide.github.io/drjacoby/) package (so much so that this "package" is more appropriately called a wrapper of `drjacoby`). 

### Overview of the `CurveAware` Model 
We assume that the only data avaialble is the cumulative number of deaths on a given day during an epidemic. In addition, we assume that the epidemic is growing exponentially (_i.e._ saturation of the fraction susceptible has not been reached). We then "cast" our incidence curve on to delay of disease-onset to death curve that we assume is gamma distributed. 

<br/>
<img src="https://raw.githubusercontent.com/nickbrazeau/CurveAware/master/R_ignore/images/curve_example.png" height="200px" width="300px" align="center" />
<br/>

As a result, the expected number of deaths on a given day is a function of the total number of individuals infected on that day accounting for the offset of time infection to (potential) death.   
  
Given that there are two (joint) features to this function, one can either put a strong prior on the: (1) IFR and infer the _I0_ (and the "exponential curve" of the epidemic); or (2) _I0_ and infer the IFR. 
