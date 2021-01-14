## Installation

You can install the released version of sumSome with:

``` r
devtools::install_github("annavesely/sumSome")
```


## Simulation

Generate matrix of permutation p-values
``` r
library(sumSome)

seed <- 42
alpha <- 0.4 # significance level

G <- simData(prop = 0.6, m = 5, B = 10, alpha = alpha, seed)
```
Determine the subset of interest

``` r
S <- c(1,2)
```
Create object of class sumSome

``` r
res <- sumPvals(G, S, alpha = alpha, type = "vovk.wang", r = 0)

res
summary(res)
```
Algorithm converged after 2 iterations. Get 0.6 confidence bounds for the number of true discoveries, true discovery proportion (TDP) and false discovery proportion (FDP).

``` r
discoveries(res) # 1
tdp(res) # 0.5
fdp(res) # 0.5
```
With 0.6 confidence, the true discoveries are at least 1, the TDP is at least 0.5, and the FDP is at most 0.5.


## fMRI Data

If needed, install the package fMRIdata from Github

``` r
devtools::install_github("angeella/fMRIdata")
```
Data
``` r
library(fMRIdata)
data("Auditory_copes")
data("Auditory_mask")
data("Auditory_clusterTH3_2")
```
Create object of class sumBrain

``` r
res <- brainPvals(copes = Auditory_copes, mask = Auditory_mask, B = 200, type = "cauchy", seed = 42)

res
summary(res)
```
Cluster analysis. Produce confidence bound for the number of true discoveries and the true discovery proportion within clusters (may require some minutes)

``` r
out <- clusterAnalysis(res, clusters = Auditory_clusterTH3_2, nMax=50, silent=FALSE)
```
Write the TDP map as Nifti file
``` r
library(RNifti)
RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz")
```
