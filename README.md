Migration and Population Surface estimation (MAPS)
=====================================

The MAPs software is adapted from the eems software and the software usage is very similar. **Here, we only highlight the differences between the usage of MAPs and EEMS**. Please see https://github.com/dipetkov/eems for the usage in EEMS.

## data input

Maps requires the

* IBD sharing matrix (ends in .sims) such
* coordinate file (ends in .coord). See EEMS documentation for more information
* outer file (ends in .outer). See EEMS documentation for more information.

An IBD sharing matrix is required for MAPs (instead of an disssimilarity matrix as in eems). The IBD sharing matrix ${X}$ is defined such that $X_{i,j}$ is the number of IBD segments shared in a length bin R between haploid $i$ and haploid $j$, the length bin or range R is described below. In the MAPs paper we use the software refinedIBD to call and phase diploid data. The sharing matrix must end with with the prefix .sims, e.g. `popressard_2_Inf.sims`. 

The length bin R is defined by a lowerbound and an upperbound on the **cM** scale, and can be specified in the params file withe the parameter `lowerBound` and `upperBound`. If `upperBound` entry is blank, it is assumed to be infinity. 

For example,
```
datapath = popressard_2_Inf.sims
lowerBound = 2
```

We also allow the capability to visualize IBD segments in a length region, for example betweem 2cM and 8cM.
```
datapath = popressard_2_Inf.sims
lowerBound = 4
upperBound = 8
```
## parameter configuration

As mentioned above, the parameters in MAPs are nearly identical. The differences are

* `genomeSize` (optional defaults to 3000cM)

* `nIndiv` (required, number of haploid individuals)

* `lowerBound` (required, lower bound in cM)

* `upperBound` (optional defaults to infinity, upper bound in cM)

* `qrateMuProposalS2` (optional, acts in the same was as in mrateMuProposalS2 except for the coalescent rates)

* `usebootstrap` (optional, with options = 0 or 1, defaults to 1. MAPS uses the bootstrap to estimate the effective sample size in the data. For very long segments >10cM, the bootstrap is not accurate and we suggest setting the option to 0)

```
./runeems2 --params params-test.ini --seed 123
```

## plotting

Finally, the MAPs results can be visualized with the function `eems.plots` defined in the R package `rEEMSplot2`. The instructions are the same as in eems. However in MAPs, one can specify constants for the inferred migration rates and population sizes to be approximately independent of the grid-density. (see plot-MAPS.R for instructions and an example)

## comparing MAPS runs for different IBD length bins

Somtimes, it can be difficult to interpret differences of MAPS runs between two IBD lengths bins. For example, 2-6cM versus >6cM. There is an option in MAPS to allow users to more easily interpret differences. For example, if I want to visualize the differences in migration surfaces between 2-6cM and >6cM.

The procedure is as follow.

1. Run MAPS on 2-6cM and plot using `eems.plots` function in R
2. Run MAPS on >6cM and plot using `eems.plots` function in R
3. Use the `eems.plot.difference` function in R to visualize difference between surfaces


## examples

Please see the `data` folder for an example.

## confused?

The MAPS software is still being tested. I will appreciate any bugs/comments with the documentation and improve the software.