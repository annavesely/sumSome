# sumSome
[![DOI](https://zenodo.org/badge/324800427.svg)](https://zenodo.org/badge/latestdoi/324800427)

sumSome is the package developed to quickly perform permutation-based closed testing by sum tests. The procedure applies to any global test which is sum-based, i.e., such that a group statistic may be written as a sum of contributions per feature (sum of t-scores, p-value combinations etc.). It adapts to the unknown joint distribution of the data through random permutations.

The method allows to construct lower confidence bounds for the proportion of true discoveries (TDP), simultaneously over all subsets of hypotheses. Simultaneity ensures control of the TDP even when the subset of interest is selected post hoc, after seeing the data.

As a main feature, the package produces simultaneous lower confidence bounds for the proportion of active voxels in different clusters for fMRI cluster analysis. Moreover, it allows to analyze gene expression data and generic permutation statistics.


## Installation

The latest version of the package can be installed with:

``` r
devtools::install_github("annavesely/sumSome")
```


## fMRI Data
To study active voxels in different clusters, we start from the list of copes (constrast maps for each subject) and the mask. Different fMRI datasets may be found in the package [fMRIdata](https://github.com/angeella/fMRIdata). As an example, here we use data from the [Auditory dataset](https://openneuro.org/datasets/ds000116/versions/00003).

``` r
require(fMRIdata) # install from https://github.com/angeella/fMRIdata
data("Auditory_copes") # list of copes
data("Auditory_mask") # mask
```

First, we compute permutation test statistics for each voxel inside the brain, through one-sample t tests. For the usual significance levels (e.g., 5%), 200 permutations are generally suitable, since they provide sufficient power without requiring much computation time. Time may be further reduced by truncating the statistics, so that each value less extreme than a parameter ```truncFrom``` is set to a given value ```truncTo```. We store information on the statistics and other analysis parameters in a ```sumBrain``` object. There are two options, as following.

**1.** The function ```brainScores``` computes t-statistics:

``` r
res <- brainScores(copes = Auditory_copes, mask = Auditory_mask, alternative = "two.sided",
                   alpha = 0.05, B = 200, seed = 42, truncFrom = 3.2, truncTo = 0)
res
summary(res)
```

**2.** The function ```brainPvals``` computes p-value combinations (Fisher, Pearson, Liptak, Edgington, Cauchy, Vovk and Wang with parameter ```r```):

``` r
res <- brainPvals(copes = Auditory_copes, mask = Auditory_mask, alternative = "two.sided",
                  alpha = 0.05, B = 200, seed = 42, truncFrom = 0.05, truncTo = 0.5,
                  type = "vovk.wang", r = 0)
res
summary(res)
```

Subsequently, we construct lower confidence bounds for the proportion of active voxels (TDP) inside different clusters. While the algorithm may require many iterations to converge to full closed testing results, in most cases 50 iterations are sufficient to obtain a TDP very close to that given by closed testing.

``` r
data("Auditory_clusterTH3_2") # cluster map
out <- brainAnalysis(sumBrain = res, clusters = Auditory_clusterTH3_2, nMax = 50, silent = FALSE)
```


Finally, we may write the TDP map as a Nifti file. In this case, we need the Nifti file for the mask.

``` r
require(RNifti)
maskNifti <- "mask.nii.gz" # name of mask Nifti file
RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz", template = maskNifti)
```


## Gene Expression Data
To study differences in gene expression between two populations, we use the expression values of different samples. Here we take and pre-process the ```montpick``` dataset from [ReCount](http://bowtie-bio.sourceforge.net/recount/index.shtml).

``` r
require(BiocManager)
require(dynamicTreeCut)
require(Biobase)

load(file=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData"))
pheno <- Biobase::phenoData(montpick.eset)
labels <- as.factor(pheno$population) # labels for two populations
expr <- Biobase::exprs(montpick.eset) # expression data

expr <- log(expr + 1) # log transform
expr <- ex[rowMeans(expr) > 0, ] # genes with non-null expression
```

Analogously to fMRI data analysis, we compute permutation test statistics for each gene, using two-sample t tests, and we store information on the analysis in a ```sumGene``` object. There are two options.

**1.** The function ```geneScores``` computes t-statistics:

``` r
res <- geneScores(expr = expr, labels = labels, alternative = "two.sided",
                   alpha = 0.05, B = 200, seed = 42)
res
summary(res)
```

**2.** The function ```genePvals``` computes p-value combinations:

``` r
res <- genePvals(expr = expr, labels = labels, alternative = "two.sided",
                  alpha = 0.05, B = 200, seed = 42, truncFrom = 0.05, truncTo = 0.5,
                  type = "vovk.wang", r = 0)
res
summary(res)
```

Subsequently, we compute lower confidence bounds for the proportion of differentially expressed genes (TDP) inside clusters.

``` r
eDist <- dist(expr) # distance between genes
tree <- hclust(eDist, method = "ward.D2")
clusters <- unname(dynamicTreeCut::cutreeDynamic(tree, distM = as.matrix(eDist)))

out <- geneAnalysis(sumGene=res, clusters=clusters, nMax = 50, silent = FALSE)
```



## General Setting
The analysis employs a matrix of statistics, where columns correspond to hypotheses, and rows to data transformations (the first is the identity). Such a matrix may be simulated with the function ```simData```. Here, we are generating p-values corresponding to 5 hypotheses and 10 permutations, where 60% of the null hypotheses are false.

``` r 
G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, p = TRUE, seed = 42)
```

Then we may analyze any subset of hypotheses, storing the results into a ```sumObj``` object. There are two options, as follows.

**1.** The function ```sumStats``` analyzes generic statistics:

``` r
S <- c(1,2) # subset of interest
res <- sumStats(G = G, S = S, alternative = "lower", alpha = 0.4, truncFrom = 0.4, truncTo = 0.5)
res
summary(res)
```

**2.** The function ```sumPvals``` analyzes p-value combinations (Fisher, Pearson, Liptak, Edgington, Cauchy, Vovk and Wang with parameter ```r```):

``` r
S <- c(1,2) # subset of interest
res <- sumPvals(G = G, S = S, alpha = 0.4, truncFrom = 0.4, truncTo = 0.5, type = "vovk.wang", r = 0)
res
summary(res)
```

The resulting ```sumObj``` object contains lower confidence bounds for the number of true discoveries and the TDP, as well as upper confidence bounds for the false discovery proportion (FDP). 

``` r
discoveries(res) # lower confidence bound for the number of true discoveries
tdp(res) # lower confidence bound for the TDP
fdp(res) # upper confidence bound for the FDP
```

# References
Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.

Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.

Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.

# Did you find some bugs?

Please write to anna.vesely@phd.unipd.it.

