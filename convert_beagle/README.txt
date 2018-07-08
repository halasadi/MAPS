Here is an example of how to convert Beagle output to a .sims and .coord file format required by MAPS

Files included:

* POPRES_CHR1.merged.nosparse_cm.finalqc.ibd: the output of the Snakemake pipeline: https://github.com/halasadi/ibd_data_pipeline
* novembre2008.txt: the geographic coordinates of a subset of individuals from the POPRES individuals (see https://github.com/NovembreLab/Novembre_etal_2008_misc)
* example_beagle_to_sims.R: an Rscript for converting the .ibd file (generated from the Snakemake pipeline) to a .sims and .coord file


