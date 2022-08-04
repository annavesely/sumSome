<img src="sticker.png" align="right" alt="" width="200" />

# sumSome
[![DOI](https://zenodo.org/badge/324800427.svg)](https://zenodo.org/badge/latestdoi/324800427)

sumSome is the package developed to quickly perform closed testing by sum tests. The procedure applies to any global test which is sum-based, i.e., such that a group statistic may be written as a sum of contributions per feature (sum of t scores, p-value combinations etc.). The method allows to construct lower confidence bounds for the proportion of true discoveries (TDP), simultaneously over all subsets of hypotheses. Simultaneity ensures control of the TDP even when the subset of interest is selected post hoc, after seeing the data.

As main features, the package uses permutations to produce simultaneous lower confidence bounds for the proportion of active voxels in clusters for fMRI data, and for the proportion of differentially expressed genes in pathways for gene expression data. Moreover, it allows to analyze generic statistics using both permutations and a parametric approach.


## Installation

The latest version of the package can be installed with:

``` r
devtools::install_github("annavesely/sumSome")
```


## fMRI data
To study active voxels in different clusters, we start from the list of copes (constrast maps for each subject) and the mask. Different fMRI datasets may be found in the package [fMRIdata](https://github.com/angeella/fMRIdata). As an example, here we use data from the [Auditory dataset](https://openneuro.org/datasets/ds000116/versions/00003).

``` r
require(fMRIdata) # install from https://github.com/angeella/fMRIdata
data("Auditory_copes") # list of copes
data("Auditory_mask") # mask
```

First, we compute permutation test statistics for each voxel inside the brain, through one-sample t tests. For the usual significance levels (e.g., 5%), 200 permutations are generally suitable, since they provide sufficient power without requiring much computation time. Time may be further reduced by truncating the statistics, so that each value less extreme than a parameter ```truncFrom``` is set to a given value ```truncTo```. We store information on the statistics and other analysis parameters in a ```sumBrain``` object. There are two options, as follows.

**1.** ```brainScores``` computes t-statistics:

``` r
res <- brainScores(copes = Auditory_copes, mask = Auditory_mask, seed = 42)
```

**2.** ```brainPvals``` computes p-value combinations (Fisher, Pearson, Liptak, Edgington, Cauchy, harmonic mean, Vovk and Wang with parameter ```r```):

``` r
res <- brainPvals(copes = Auditory_copes, mask = Auditory_mask, seed = 42, type = "vovk.wang", r = 0)
```

Subsequently, we construct lower confidence bounds for the proportion of active voxels (TDP) inside different clusters. While the algorithm may require many iterations to converge to full closed testing results, in most cases 50 iterations are sufficient to obtain adequate results.

``` r
summary(res)
data("Auditory_clusterTH3_2") # cluster map
out <- brainAnalysis(sumBrain = res, clusters = Auditory_clusterTH3_2)
```


Finally, we may write the TDP map as a Nifti file. In this case, we need the Nifti file for the mask (e.g., ```mask.nii.gz``` saved in the working directory).

``` r
require(RNifti)
maskNifti <- "mask.nii.gz" # name of mask Nifti file
RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz", template = maskNifti)
```


## Gene expression data
To study differences in gene expression between two populations, we use the expression values of different samples. Here we take and pre-process the ```BRCA``` (Breast Invasive Carcinoma) dataset from [TGCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) to compare two types of primary solid tumor: infiltrating lobular carcinoma and infiltrating ductal carcinoma.

``` r
require(curatedTCGAData)
require(TCGAutils)
require(EnrichmentBrowser)

d <- curatedTCGAData::curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2Gene*", dry.run = FALSE)
d <- TCGAutils::splitAssays(d, sampleCodes = "01", exclusive = TRUE) # primary solid tumor
d <- d[ ,d$histological_type %in% c("infiltrating lobular carcinoma", "infiltrating ductal carcinoma")]
labels <- d$histological_type

d <- d[rowMeans(assay(d)) >= 10, ] # select genes with mean expression >= 10
d <- EnrichmentBrowser::idMap(d[[1]], org = "hsa", from = "SYMBOL", to = "ENTREZID") # map gene ID types
expr <- assay(d)
```

Analogously to fMRI data analysis, we compute permutation test statistics for each gene, using two-sample t tests, and we store information on the analysis in a ```sumGene``` object.

**1.** ```geneScores``` computes t-statistics:

``` r
res <- geneScores(expr = expr, labels = labels, seed = 42)
```

**2.** ```genePvals``` computes p-value combinations:

``` r
res <- genePvals(expr = expr, labels = labels, seed = 42, type = "vovk.wang", r = -1)
```

Subsequently, we compute lower confidence bounds for the proportion of differentially expressed genes (TDP) inside pathways.

``` r
summary(res)
pathways <- EnrichmentBrowser::getGenesets(org = "hsa", db = "kegg")
out <- geneAnalysis(sumGene = res, pathways = pathways)
```


## General setting

### Permutation approach
In the general setting, we start with a matrix of statistics, where columns correspond to hypotheses, and rows to data transformations (the first is the identity). Such a matrix may be simulated with the function ```simData```. Here, we are generating p-values corresponding to 5 hypotheses and 10 permutations, where 60% of the hypotheses are false.

``` r 
G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, p = TRUE, seed = 42)
S <- c(1,2) # subset of interest
```

Then we may analyze any subset of hypotheses, storing results into a ```sumObj``` object.

**1.** ```sumStats``` analyzes generic statistics:

``` r
res <- sumStats(G = G, S = S, alternative = "lower", alpha = 0.4, truncFrom = 0.4, truncTo = 0.5)
```

**2.** ```sumPvals``` analyzes p-value combinations:

``` r
res <- sumPvals(G = G, S = S, alpha = 0.4, truncFrom = 0.4, truncTo = 0.5, type = "vovk.wang", r = 0)
```

The resulting ```sumObj``` object contains confidence bounds for the number of true discoveries, the TDP and the false discovery proportion (FDP). 

``` r
summary(res)
discoveries(res) # lower confidence bound for the number of true discoveries
tdp(res) # lower confidence bound for the TDP
fdp(res) # upper confidence bound for the FDP
```



### Parametric approach
The analysis in the parametric framework is similar to the previous one. We start with a vector of statistics, each corresponding to a hypothesis. Here we generate p-values corresponding to 5 hypotheses.

``` r 
g <- as.vector(simData(prop = 0.6, m = 5, B = 1, alpha = 0.4, p = TRUE, seed = 42))
S <- c(1,2) # subset of interest
```

Then we may analyze any subset of hypotheses, obtaining a ```sumObj``` object.

**1.** ```sumStatsPar``` analyzes generic statistics, relying on a user-defined vector of critical values (e.g., here we consider Fisher combination of p-values):

``` r
g <- -2 * log(g) # statistics
cvs <- qchisq(p = 0.4, df = 2 * seq(5), lower.tail = FALSE) # critical values
res <- sumStatsPar(g = g, S = S, alpha = 0.4, cvs = cvs)
```

**2.** ```sumPvalsPar``` analyzes p-value combinations (Fisher, Pearson, Liptak, Cauchy, harmonic mean, Vovk and Wang with parameter ```r```), either under the assumption of independence or for general dependence structures:

``` r
res <- sumPvalsPar(g = g, S = S, alpha = 0.4, type = "harmonic", independence = FALSE)
```

The ```sumObj``` object can be accessed as in the previous section.

``` r
summary(res)
discoveries(res) # lower confidence bound for the number of true discoveries
tdp(res) # lower confidence bound for the TDP
fdp(res) # upper confidence bound for the FDP
```


# References
Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.

Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.

Tian, J., Chen, X., Katsevich, E., Goeman, J. J. and Ramdas, A. (2021). Large-scale simultaneous inference under dependence. Scandinavian Journal of Statistics, to appear. (Pre-print arXiv:2102.11253)

Vesely, A., Finos, L., and Goeman, J. J. (2021). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.


# Did you find some bugs?

Please write to anna.vesely@unipd.it.

