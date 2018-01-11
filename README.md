Migration and Population Surface estimation (MAPS)
=====================================

The MAPs software is built from the eems software and the software
usage is very similar. **Here, we only highlight the differences
between the usage of MAPs and EEMS**. Please see
https://github.com/dipetkov/eems for the usage in EEMS.

## installing MAPS

I recommend using *conda* to install MAPS

*  install conda, see: https://conda.io/docs/
* ```conda create -n MAPS``` (this creates a conda enviornment for MAPS)
* ```source activate MAPS```
* ```conda install boost=1.57.0``` (MAPS only works with this version of boost)
* ```conda install eigen```
* clone the repository 
* change the paths accordingly the Makefile
   For example, the paths to Eigen and Boost in my Makefile are
   ```
   EIGEN_INC = /Users/halasadi/anaconda/envs/MAPS/include/eigen3
   BOOST_LIB = /Users/halasadi/anaconda/envs/MAPS/lib
   BOOST_INC = /Users/halasadi/anaconda/envs/MAPS/include/boost
   ```
* in the ```src``` directory, type ```make``` 


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

* `olderpath` (optional, path to a run with a older time period, MAPS will only estimate the difference between rates
               from the older time period)

## running MAPS

You can run MAPS with the command such as this

```
./runeems2 --params params.ini --seed 123
```

## plotting

Please use the plotmaps package, https://github.com/halasadi/plotmaps

## examples

Please see the `examples` folder for examples data-input, how to run MAPS, and how to plot the results.

## confused?

The MAPS software is still being tested. I will appreciate any bugs/comments with the code/documentation. Please post an issue and I will get to it.
