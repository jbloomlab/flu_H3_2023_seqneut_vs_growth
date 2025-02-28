# Multiple regression of neutralization titers and HA1 mutations with growth rates
In the main pipeline, we performed an analysis that first fit MLR models to strain frequencies over time, and then used those model fits to estiamte strain-specific growth rates. We found that the fraction of individuals with low neutralization titers and number of HA1 mutations correlated strongly with these MLR-estimated strain growth rates. However, the fraction of individuals with low neutralization titers and number of HA1 mutations were also collinear. Here, using multiple regression, we try to determine which of the predictors (neutralization titers or HA1 mutations) more fully explain the dependent outcome variable (growth rate). 

## Build `conda` environment
First build and activate the conda environment in [environment.yml](environment.yml) with

    conda env create -f environment.yml
    conda activate multiple_regression_titers_and_muts

*Note that if not on the Hutch cluster, the environment.yml file should be edited according to the comments in that file.*

## Run the analysis
All input data are pulled directly from the `flu_H3_2023_seqneut_vs_growth` output, and are specified in the notebook. The analyses themselves are are performed interactively in Jupyter notebook [multiple_regression_titers_and_muts.ipynb](multiple_regression_titers_and_muts.ipynb) which includes additional notes and commentary on the analysis.