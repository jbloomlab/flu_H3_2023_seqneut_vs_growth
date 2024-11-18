# Getting number of mutations from MRCA for the 2022 / 2023 library strains

The file [h3n2_ha.json](h3n2_ha.json) is a Nextstrain JSON of all of the strains at the time of library design provided by Caroline.

The notebook [get_mutations.ipynb](get_mutations.ipynb) processes that tree to get the mutations relative to the MRCA into the file [HA_mutations_from_MRCA.csv](HA_mutations_from_MRCA.csv) for all of the strains with library titers.

Note that the notebook does some adjustments to account for the fact that in the tree the strain names lack underscores and that in the tree `A/SOUTHAFRICA/R07876/2023` is misnamed as `A/SOUTHAFRICA/R07876/202023`.
