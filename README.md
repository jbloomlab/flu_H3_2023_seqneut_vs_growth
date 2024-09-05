# Analysis of growth rates of influenza viruses with HAs in H3 2023-2024 sequencing-based neutralization libraries

## Input data
Input data are in [./data/](data):

 - [data/H3_library_constructs_protein.fasta](data/H3_library_constructs_protein.fasta): file from Caroline with the strains in the library. Note that this file has the full chimeric HAs used in the sequencing-based neutralization assay, which have non-native signal peptide and endodomain.

 - [data/ncbi_flu_h3_prots.fa](data/ncbi_flu_h3_prots.fa): set of sequences downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) on Sept-4-2024 selecting the following options:
   - *Protein* sequences
   - *Virus/Taxonomy* of *Influenza A virus, taxid:11320*
   - *Genotype* of *H3N2*
   - *Host* of *Homo sapiens (human), taxid:9606*
   - *Collection Date* from *Jan 1, 2023* to *Jun 30, 2024*
   - *Segment* set to *4* (hemagglutinin)
   - *Lab Passaged* set to *Exclude*

 - [data/gisaid_flu_h3_prots.fa](data/gisaid_flu_h3_prots.fa): (**due to GISAID data sharing rules, this file is not tracked in the GitHub repository**) set of sequences downloaded from [GISAID](https://gisaid.org/) EpiFlu on Sept-5-2024 selecting the following options:
   - Influenza A, H3N2
   - Host of human
   - *Collection Date* from *Jan 1, 2023* to *Jun 30, 2024*
   - Only keeping *original* sequences (excluding lab passaged)
   - Downloading just HA protein sequences.
