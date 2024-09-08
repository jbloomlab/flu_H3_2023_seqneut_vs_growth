# Growth rates of influenza viruses with HAs in H3 2023-2024 sequencing-based neutralization libraries

The overall goal of this analysis is to estimate a growth rate in the 2023-2024 timeframe (probably centered around mid to late 2023) for influenza viruses with the HA proteins in a [sequencing-based neutralization assay](https://www.biorxiv.org/content/10.1101/2024.03.08.584176v1) library to examine how neutralization titers against the strains are related to their growth rates.

In order to do that, it is necessary to first figure out sequence counts over time for sequences that "match" (by some reasonable definition) the strains in the library, and then estimate their growth rates such as using the approaches implemented in [evofr](https://github.com/blab/evofr).

## Running the pipeline
Build the `conda` environment in [environment.yml](environment.yml), then activate it and run the `snakemake` pipeline in [Snakefile](Snakefile), eg:

    conda activate flu_H3_2023-2024_strains_MLR
    snakemake -j 8 --software-deployment-method conda

If you are using the Hutch cluster, you can also just run via the [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script, which can be submitted via `sbatch`.

## Input data and configuration
The configuration is specified in [config.yaml](config.yaml)

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

## Workflow and results
The results are placed in [./results/](results), some of which may not be tracked in this GitHub repository.
Essentially, the workflow proceeds as follows:

#### Get counts of sequences matching each strain
This step is executed by the rule `strain_counts` in [Snakefile](Snakefile).

For each strain in the library (specified under `strain_prots` in [config.yaml](config.yaml)), we find how many naturally occurring sequences for each date match to that sequence.
We look for matches in the sets of naturally occurring sequence specified under `protsets` in [config.yaml](config.yaml).
To determine if a sequence is a match to a protein in the library, we first trim the library proteins according to the boundaries specified in `trim` under `protsets`.
At least some trimming is important because the library strains have chimeric signal peptides and endodomains (transmembrane domains and cytoplasmic tails).
We then go through the natural sequences and see which ones match (contain a substring) that matches the trimmed library protein, allowing the number of differences (single amino-acid substitutions or indels) specified under `protset_maxdiffs` in [config.yaml](config.yaml).
If there are multiple matches, we take the best one (fewest differences).
If there are multiple matches with the same number of differences, we assign fractional weights to the matching (eg, a sequence that matrches two library strains is assigned a weight of 0.5 to each).
Note that we also require these sequences to have a complete date in the `YYYY-MM-DD` format.

The results of this rule are in [./results/strain_counts/](results/strain_counts) with files named as:
  - `<protset>_<maxdiff>_counts_overall.csv`: overall number of sequences that match to each library protein (variant) for each `protset` and `maxdiff`.
  - `<protset>_<maxdiff>_counts_by_date.csv`: number of sequences that match to each library protein for each date. This file can serve as input to [evofr](https://github.com/blab/evofr).
  - `<protset>_<maxdiff>_strain_matches.csv`: the matches for each library protein for each sequence

Note that sequences that do not match any library proteins are assigned to the variant *other*.

#### Fit MLR fitness estimates to the strain counts
This step is executed by the rule `mlr` in [Snakefile](Snakefile).

For each set of sequences matched to strains (based on the different `protset` and `maxdiff` values), we first figure out how many strains have enough sequences matching them to make a reasonable fitness estimate.
We look at sequences in each date range specified under `dateranges` in [config.yaml](config.yaml), and consider a strain to have enough sequences to fit if it has at least `min_counts` sequences as specified in [config.yaml](config.yaml) over the total date range.
First, we make some plots determining how many librarystrains have enough sequences, and we group the others together as a new set of library strains with insufficient counts.
Interactive [altair](https://altair-viz.github.io) charts are created that show the counts and frequencies of each library strain.
These charts are in [./results/mlr/](results/mlr) as files named as `results/mlr/strain_counts_<protset>_<maxdiff>_<daterange>.html`.

In addition, the Jupyter notebooks doing the entire analysis for this step are saved as files named as `results/mlr/mlr_<protset>_<maxdiff>_<daterange>.ipynb`.

**The next step is the MLR, which is not yet implemented**

**Note that the charts and notebooks created by this rule are currently not tracked in this GitHub repo.**

