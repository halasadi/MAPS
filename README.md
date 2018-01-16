# Migration and Population Surface estimation (MAPS)

*Add brief (1--3 sentence) description of MAPS here. Even better,
include a cool plot generated using your software; see
[here](https://github.com/genetics-statistics/GEMMA) for an example.*

The MAPs software is built from the eems software and the software
usage is very similar. **Here, we only highlight the differences
between the usage of MAPs and EEMS**. Please see
https://github.com/dipetkov/eems for the usage in EEMS.

*The MAPS program has been tested with ...*

Please post bugs, questions and feature requests or suspected bugs to
[Github issues](https://github.com/halasadi/MAPS/issues).

## Citing this work

If you find the MAPS program, or any source code contained in this
repository, useful for your work, please cite:

> Add citation here.

## License

Copyright (c) 2017-2018, Hussein Al-Asadi.

The software and example programs in this repository are made
available under the terms of the
[MIT license](https://opensource.org/licenses/mit-license.html).
See [LICENSE](LICENSE) for the full text of the license.

## Quick Start

Follow these steps to quickly get started using MAPS.

For installing the software dependencies, we provide detailed
instructions using the
[conda package manager](https://conda.io/docs). This is only a
recommended aproach---conda is not required to use MAPS, and other
approaches can be taken to install the dependencies (e.g., by directly
downloading the source code, or by using [Homebrew](http://brew.sh) for
Mac).

1. Install [GNU Make](https://www.gnu.org/software/make).

2. Install a [Conda](https://conda.io/docs) distribution such as
   [Anaconda](https://www.anaconda.com/download) or
   [Miniconda](https://conda.io/miniconda.html) (*optional*).

3. Create a new conda environment for MAPS (*optional*).

   ```bash
   conda create -n MAPS
   source activate MAPS
   ```

3. Install the [Boost](http://www.boost.org) C++ source libraries
   version 1.57.0. If using conda, run `conda install boost=1.57.0`.

4. Install the [Eigen](http://eigen.tuxfamily.org) C++ matrix algebra
   library. If using conda, run `conda install eigen`.

*More instructions will go here.*

## installing MAPS

I recommend using *conda* to install MAPS

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

## Preparing data for MAPS

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

## Credits

This project was developed by
[Hussein Al-Hasadi](https://github.com/halasadi) at the University of
Chicago.

Thanks to [Matthew Stephens](stephenslab.uchicago.edu) for his support
and mentorship.
