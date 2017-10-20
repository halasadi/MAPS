Migration and Population Surface estimation (MAPS)
=====================================

The MAPs software is built from the eems software and the software usage is very similar. **Here, we only highlight the differences between the usage of MAPs and EEMS**. Please see https://github.com/dipetkov/eems for the usage in EEMS.

## installing MAPS

I recommend using *conda* to install MAPS

*  install conda, see: https://conda.io/docs/
* ```conda create -n MAPS``` (this creates a conda enviornment for MAPS)
* ```source activate MAPS```
* ```conda install boost=1.5.7``` (MAPS only works with this version of boost)
* ```conda install eigen```
* clone the repository and in the ```src``` directory, type ```make``` 


## preparing data for MAPS

MAPS requires for input

* IBD sharing matrix (ends in .sims) (see below)
* coordinate file (ends in .coord). see eems documentation
* outer file (ends in .outer). see eems documentation

An IBD sharing matrix is required for MAPs (instead of an disssimilarity matrix as in eems). The IBD sharing matrix ${X}$ is defined such that $X_{i,j}$ is the number of IBD segments shared in a length bin R between haploid $i$ and haploid $j$, the length bin or range R is described below. In the MAPs paper we use the software refinedIBD to call and phase diploid data. The sharing matrix must end with with the prefix .sims, e.g. `popressard_2_Inf.sims`. 

The length bin R is defined by a lowerbound and an upperbound on the **cM** scale, and can be specified in the `params` file withe the parameter `lowerBound` and `upperBound`. If `upperBound` entry is blank, it is assumed to be infinity. 

For example,
```
datapath = popressard_2_Inf.sims
lowerBound = 2
```

We also allow the capability to visualize IBD segments in a length region, for example betweem 2cM and 8cM.
```
datapath = popressard_2_8.sims
lowerBound = 2
upperBound = 8
```
## parameter configuration

As mentioned above, the parameters in MAPs are nearly identical. However, there are a few additional arguments in the `params` file, 

* `genomeSize` (optional defaults to 3000cM)

* `nIndiv` (required, number of haploid individuals)

* `lowerBound` (required, lower bound in cM)

* `upperBound` (optional defaults to infinity, upper bound in cM)

* `qrateMuProposalS2` (optional, acts in the same was as in mrateMuProposalS2 except for the coalescent rates)

* `usebootstrap` (optional, with options = 0 or 1, defaults to 1. MAPS uses the bootstrap to estimate the effective sample size in the data. For very long segments >10cM, the bootstrap is not accurate and we suggest setting the option to 0)

## running MAPS

You can run MAPS with the command such as this

```
./runeems2 --params params.ini --seed 123
```

## plotting

Finally, the MAPs results can be visualized with the function `eems.plots` defined in the R package `rEEMSplot2`. The MAPS R plotting scripts are built upon the eems plotting scripts and therfore the usage is very similar. You must install the R plotting scripts from source, 

```
## Part 1: Install rEEMSplots2
## Check that the current directory contains the rEEMSplots source directory
if (file.exists("./rEEMSplots2")) {
  install.packages("rEEMSplots2", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```

See `examples` on how to plot.

### comparing MAPS runs for different IBD length bins

Sometimes, it can be difficult to interpret differences of MAPS runs between two IBD lengths bins. For example, 2-6cM versus >6cM. There is an option in MAPS to plot the differences between MAPS runs. For example, if I want to visualize the differences in migration surfaces between 2-6cM and >6cM.

The procedure is as follow.

1. Run MAPS on 2-6cM and plot using `eems.plots` function in R
2. Run MAPS on >6cM and plot using `eems.plots` function in R
3. Use the `eems.plot.difference` function in R to visualize difference between surfaces


## examples

Please see the `examples` folder for examples data-input, how to run MAPS, and how to plot the results.

## confused?

The MAPS software is still being tested. I will appreciate any bugs/comments with the code/documentation. Please post an issue and I will get to it.

## credits

The MAPS software was developed by Hussein Al-Asadi with Desislava Petkova at the University of Chicago.  