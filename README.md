recent Migration and Population Surface estimation (MAPs)
=====================================

Long streches of sequence similarity are indicative of recent ancestry. When long enough, these long streches of similarity (aka "IBD" segments) corresponds to segments with no intervening recombination. Using coalescent theory, we model the number of these segments between individuals across space and infer a population size surface and migration surface.

The MAPs software is adapted from the eems software and the software usage is very similar. **Here, we only highlight the differences between the usage of MAPs and EEMS**. Please see https://github.com/dipetkov/eems for the usage in EEMS.

## data input

* IBD sharing matrix (ends in .sims)
* coordinate file (ends in .coord). See EEMS documentation for more information
* outer file (ends in .outer). See EEMS documentation for more information.

Instead of an disssimilarity matrix as input (as in in eems), a IBD sharing matrix is required for MAPs. The IBD sharing matrix is such that $N_{i,j}$ is the number of segments greater than $L$ cM shared between haploid $i$ and haploid $j$. Fortunately, a popular IBD calling software refinedIBD automatically phases diplid data.The sharing matrix must end with with the prefix .sims, e.g. `eems_4_Inf.sims`. 

The `L` (on **cM** scale) parameter is specified in the params file withe the parameter `lowerBound`

For example,
```
datapath = ../data/europe_4_Inf.sims
lowerBound = 4
```

We also allow the capability to visualize IBD segments in a length region, for example betweem 4cM and 8cM.
```
datapath = ../data/europe_4_8.sims
lowerBound = 4
upperBound = 8
```
## parameter configuration

As mentioned above, the parameters in MAPs are nearly identical. The differences are

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