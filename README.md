# Growth rates of 2023 H3 influenza viruses versus sequencing-based neutralization assay titers
Analysis by Caroline Kikawa, [Jesse Bloom](https://jbloomlab.org/), John Huddleston, [Trevor Bedford](https://bedford.io/).

## Overview

The overall goal of this analysis is to determine if the growth rates of different human H3N2 influenza strains correlate with human neutralizing antibody titers against these strains as measured in a high-throughput [sequencing-based neutralization assay](https://www.biorxiv.org/content/10.1101/2024.03.08.584176v1).
These titers were measured by Caroline Kikawa using a library of influenza viruses with HAs primarily from strains circulating in late 2023 against a variety of children and adult sera, with the neutralization data at [https://github.com/jbloomlab/flu_seqneut_H3N2_2023-2024](https://github.com/jbloomlab/flu_seqneut_H3N2_2023-2024).

This analysis uses multinomial logistic regression as implemented in [evofr](https://github.com/blab/evofr) to estimate the growth advantages of all of the strains in the library with sufficient sequences.
It then compares these growth advantages to the neutralization titers.

Doing this analysis requires making various choices about how to count natural sequences as corresponding to sequencing-based neutralization assay library strains, how many sequences to require to make a growth advantage estimate for a strain, what timeframe to analyze, and what sera to use.
This pipeline is structured to allow exploration of many different values for these choices.
For a good summary of what we deem to be the most "reasonable" choice, see the following plots:
  - [interactive chart comparing growth advantages to neutralization titers](https://jbloomlab.github.io/flu_H3_2023_seqneut_vs_growth/growth_vs_titers_gisaid-ha1-within1_2023-mincounts80_child-and-adultprevax-sera.html)
  - [interactive chart showing multinomial logistic regression fits for estimating growth advantages](https://jbloomlab.github.io/flu_H3_2023_seqneut_vs_growth/mlr_gisaid-ha1-within1_2023-mincounts80.html)

The key result shown by above plots is that the growth advantages correlate extremely well with the neutralization titers, with strains with lower titers having higher growth advantages.
To explore similar plots for other values for the various choices, see [https://jbloomlab.github.io/flu_H3_2023_seqneut_vs_growth](https://jbloomlab.github.io/flu_H3_2023_seqneut_vs_growth); the key finding of a strong correlation turns out to be quite robust to the choices.

## Running the pipeline
Build the `conda` environment in [environment.yml](environment.yml), then activate it and run the `snakemake` pipeline in [Snakefile](Snakefile), eg:

    conda activate flu_H3_2023_seqneut_vs_growth
    snakemake -j 8 --software-deployment-method conda

If you are using the Hutch cluster, you can also just run via the [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script, which can be submitted via `sbatch`.

## Input data and configuration
The configuration is specified in [config.yaml](config.yaml), and should be largely self-explanatory from the comments in that file.

Input data are in [./data/](data):

 - [data/H3_library_constructs_protein.fasta](data/H3_library_constructs_protein.fasta): file from Caroline with the strains in the library. Note that this file has the full chimeric HAs used in the sequencing-based neutralization assay, which have non-native signal peptide and endodomain.

 - [data/ncbi_flu_h3_prots.fa](data/ncbi_flu_h3_prots.fa): set of sequences downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) on Sept-4-2024 selecting the following options:
   - *Protein* sequences
   - *Virus/Taxonomy* of *Influenza A virus, taxid:11320*
   - *Genotype* of *H3N2*
   - *Host* of *Homo sapiens (human), taxid:9606*
   - *Collection Date* from *Jan 1, 2022* to *Sept 1, 2024*
   - *Segment* set to *4* (hemagglutinin)
   - *Lab Passaged* set to *Exclude*

 - [data/gisaid_flu_h3_prots.fa](data/gisaid_flu_h3_prots.fa): (**due to GISAID data sharing rules, this file is not tracked in the GitHub repository**) set of sequences downloaded from [GISAID](https://gisaid.org/) EpiFlu on Sept-5-2024 selecting the following options:
   - Influenza A, H3N2
   - Host of human
   - *Collection Date* from *Jan 1, 2022* to *Sept 1, 2024*
   - Only keeping *original* sequences (excluding lab passaged)
   - Downloading just HA protein sequences.

 - [data/H3_seqneut_titers.csv](data/H3_seqneut_titers.csv): titers from Caroline's 2023 H3N2 sequencing-based neutrlization experiments as taken from [https://github.com/jbloomlab/flu_seqneut_H3N2_2023-2024](https://github.com/jbloomlab/flu_seqneut_H3N2_2023-2024) on Sept-11-2024.

## Workflow and results
The results are placed in [./results/](results), some of which may not be tracked in this GitHub repository.
Essentially, the workflow proceeds as follows:

#### Get counts of sequences matching each strain
This step is executed by the rule `strain_counts` in [Snakefile](Snakefile).

For each strain in the library (specified under `strain_prots` in [config.yaml](config.yaml)), we find how many naturally occurring sequences for each date match to that sequence.
We look for matches in the sets of naturally occurring sequence specified under `protsets` in [config.yaml](config.yaml).
To determine if a sequence is a match to a protein in the library, we first trim the library proteins according to the boundaries specified in `trim` under `protsets`.
At least some trimming is important because the library strains have chimeric signal peptides and endodomains (transmembrane domains and cytoplasmic tails).
We then go through the natural sequences and see which ones match (contain a substring) that matches the trimmed library protein, allowing the number of differences (single amino-acid substitutions or indels) specified under `maxdiff` under `protesets` in [config.yaml](config.yaml).
If there are multiple matches, we take the best one (fewest differences).
If there are multiple matches with the same number of differences, we assign fractional weights to the matching (eg, a sequence that matrches two library strains is assigned a weight of 0.5 to each).
Note that we also require these sequences to have a complete date in the `YYYY-MM-DD` format.

The results of this rule are in [./results/strain_counts/](results/strain_counts) with files named as:
  - `<protset>_counts_overall.csv`: overall number of sequences that match to each library protein (variant) for each `protset`.
  - `<protset>_counts_by_date.csv`: number of sequences that match to each library protein for each date. This file can serve as input to [evofr](https://github.com/blab/evofr).
  - `<protset>_strain_matches.csv`: the matches for each library protein for each sequence

Note that sequences that do not match any library proteins are assigned to the variant *other*.

#### Fit MLR fitness estimates to the strain counts
This step is executed by the rule `mlr` in [Snakefile](Snakefile).

There are different settings for the MLR defined under `mlrfits` in [config.yaml](config.yaml).
First each set of sequences matched to strains (eg, the `protset`s defined above), we first figure out how many strains have enough sequences matching them to make a reasonable fitness estimate, restricting to the date ranges specified for that `mlrfit` and requiring the specified minimum counts over the date range.
We then fit the `mlr` models, either dropping or keeping as pseudo strains the sequences not in the library or assigned to library strains with insufficient counts, as specified for that `mlrfit`.
The fitting is done using [evofr](https://github.com/blab/evofr), with the settings specified by `mlr_settings` in [config.yaml](config.yaml).

The results of each `mlrfit` for each `protset` are saved in [./results/mlr/](results/mlr) and are as follows:
  - `growth_advantages_<protset>_<mlrfit>.csv`: the growth advantages for each strain fit by MLR.
  - `mlr_<proset>_<mlrfit>.html`: interactive [altair](https://altair-viz.github.io) chart displaying key information about the protein set and fit.
  - `counts_to_fit_<protset>_<mlrfit>.csv`: the counts for strains being fit in MLR.
  - `mlr_<protset>_<mlrfit>.ipynb`: Jupyter notebook that does fitting (not tracked in repo).

#### Compare MLR fits
This step is executed by the rule `compare_mlr_fits` in [Snakefile](Snakefile).

This rule simply compares the growth advantages estimated for the different `protset` and `mlrfit` settings.
The results are saved in [./results/compare_mlr_fits](results/compare_mlr_fits):
  - [results/compare_mlr_fits/mlrfits_corr.html](results/compare_mlr_fits/mlrfits_corr.html): plot of correlation coefficients of growth advantages for different fits.
  - [results/compare_mlr_fits/mlrfits_scatter.html](results/compare_mlr_fits/mlrfits_scatter.html): scatter plots comparing the per-strain growth advantages for different fits.
  - `results/compare_mlr_fits/compare_mlr_fits.ipynb`: Jupyter notebook that makes the plots (not tracked in repo).

#### Compare titers to growth rates
This step is executed by the rule `growth_vs_titers` in [Snakefile](Snakefile).

This rule compares the growth advantages of strains to their titers as measured in the sequencing-based neutralization assays, using the sera titers specified under `sera` in [config.yaml](config.yaml).

For the comparison, we first get the set of strains that have both a growth advantage estimate and titer data.
We then compare (correlate) the growth advantages to the titers summarized across sera in three different ways:
  - the **geometric** mean titer
  - the median titer
  - the fraction of sera with a titer less than some cutoff, with the cutoff chosen to maximize the correlation

For the mean and median titers, we correlate growth advantage with the log of these quantities; for the fraction below the titer we correlate with that quantity directly.

We also estimate P-values for each correlation by randomizing the growth data among strains and re-computing the correlation for these randomized data (the P-value is the fraction of randomizations with as good of a correlation as the actual data).
The P-values are one sided, and for the mean and median titer represent the fraction of randomizations with Pearson R less than the actual value, while for fraction below a cutoff they represent the fraction of randomizations with a Pearson R greater than the actual value.
The number of randomizations is determined by the `nrandom` parameter under `growth_vs_titer_params` in [config.yaml](config.yaml).
For the fraction below a cutoff titer, both for the actual data and each randomization we test a range of titers to pick the cutoff that maximizes the correlation.
We do this by testing `corr_titer_cutoff_points` cutoffs spaced logarithmically in `corr_titer_cutoff_range`, where these are parameters specified under `growth_vs_titer_params` in [config.yaml](config.yaml).

The results of this analysis are placed in [results/growth_vs_titers](results/growth_vs_titers/) as follows:
 - `growth_vs_titers_<protset>_<mlrfit>_<sera>.html`: chart showing correlations, and randomizations for the cutoff method.
 - `growth_vs_titers_<protset>_<mlrfit>_<sera>.ipynb`: Jupyter notebook doing the analysis.

#### Make `./docs/` folder with HTML plots for rendering with GitHub Pages
The rule `charts_to_docs` makes a [./docs/](docs) folder that contains the interactive HTML charts alongside a rudimentary index that can be rendered via GitHub Pages to allow easy inspection of the chart.
