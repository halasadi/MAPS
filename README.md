recent Migration and Population Surface estimation (MAPs)
=====================================

Long streches of sequence similarity are indicative of recent ancestry. When long enough, these long streches of similarity (aka "IBD" segments) corresponds to segments with no intervening recombination. Using coalescent theory, we model the number of these segments between individuals across space and infer a population size surface and migration surface.

** we define long to be > L for some sufficiently long user defined L, e.g. 6cM

**Note: the software is located in the "src" folder.**

## data input

The MAPs software is adapted from the eems software and as a result, the usage is similar. **Here, we only highlight the differences**. Please see https://github.com/dipetkov/eems for all the information not listed here. 

Instead of an disssimilarity matrix as in eems, an similarity matrix based on IBD sharing must be used as input for MAPs. To create a IBD similarity matrix, first apply IBD calling softwares such as refinedIBD to the data. We use refinedIBD with default settings. It is important to discard IBD segments in sparse regions (e.g. regions with less than a SNP every centimorgan).

Afterwards, a matrix must be constructed such that the i,jth entry of the matrix is number of IBD segments greater than L between individuals i and j. Individuals must be haploid so the user must first phase the data (fortunately refinedIBD automatically phases the data.) The similarity matrix must end with with the prefix .sims, e.g. `eems_4_Inf.sims`. 

The `L` (on **cM** scale) parameter is specified in the params file withe the parameter `lowerBound`

For example,
```
datapath = ../data/europe_4_Inf.sims
lowerBound = 4
```

We also allow the capability to visualize IBD segments between a threshold. For example betweem 4cM and 8cM, which can be specified like this:
```
datapath = ../data/europe_4_8.sims
lowerBound = 4
upperBound = 8
```
## parameter configuration

Almost all the parameters are the same as in eems. The differences are

* `genomeSize` (optional defaults to 3000cM)

* `nIndiv` (required, number of haplods)

* `lowerBound` (required, lower bound in cM)

* `upperBound` (optional defaults to infinity, upper bound in cM)

* `qrateMuProposalS2` (optional, acts in the same was as in mrateMuProposalS2 except for the coalescent rates)

```
./runeems2 --params params.ini --seed 123
```

## plotting

Finally, the MAPs results can be visualized with the function `eems.plots` defined in the R package `rEEMSplot2`. The instructions are the same as in eems. 