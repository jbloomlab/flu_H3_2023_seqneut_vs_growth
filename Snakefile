"""``snakemake`` file that runs analysis."""


import itertools
import re


configfile: "config.yaml"


# output charts
charts = {
    "Strain titers versus growth advantage": {
        f"{protset=}": {
            f"{mlrfit=}": {
                f"{sera=}": f"results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}.html"
                for sera in config["sera"]
            }
            for mlrfit in config["mlrfits"]
        }
        for protset in config["protsets"] 
    },
    "MLR fits of growth advantage": {
        f"{protset=}": {
            f"{mlrfit=}": f"results/mlr/mlr_{protset}_{mlrfit}.html"
            for mlrfit in config["mlrfits"]
        }
        for protset in config["protsets"]
    },
    "Comparison of MLR fits": {
        "correlation matrix": "results/compare_mlr_fits/mlrfits_corr.html",
        "scatter plot": "results/compare_mlr_fits/mlrfits_scatter.html",
    },
}


def extract_final_values(d):
    """Get list of final values in nested dict."""
    values = []
    for key, value in d.items():
        if isinstance(value, dict):
            values.extend(extract_final_values(value))
        else:
            values.append(value)
    return values


rule all:
    """Target rule."""
    input:
        expand(
            [
                "results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_scatter.csv",
                "results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_corr.csv",
            ],
            protset=config["protsets"],
            mlrfit=config["mlrfits"],
            sera=config["sera"],
        ),
        "docs",


rule strain_counts:
    """Get counts of each strain."""
    input:
        strain_prots=config["strain_prots"],
        protset=lambda wc: config["protsets"][wc.protset]["protset"],
    output:
        counts_overall="results/strain_counts/{protset}_counts_overall.csv",
        counts_by_date="results/strain_counts/{protset}_counts_by_date.csv",
        strain_matches="results/strain_counts/{protset}_strain_matches.csv",
    params:
        trim=lambda wc: config["protsets"][wc.protset]["trim"],
        maxdiff=lambda wc: config["protsets"][wc.protset]["maxdiff"],
    log:
        "results/logs/strain_counts_{protset}.txt",
    conda:
        "environment.yml"
    script:
        "scripts/strain_counts.py"


rule mlr:
    """Fit MLR estimates of fitness (growth) advantages."""
    input:
        counts_by_date="results/strain_counts/{protset}_counts_by_date.csv",
    output:
        chart="results/mlr/mlr_{protset}_{mlrfit}.html",
        counts_to_fit="results/mlr/counts_to_fit_{protset}_{mlrfit}.csv",
        growth_advantages="results/mlr/growth_advantages_{protset}_{mlrfit}.csv",
    params:
        date_start=lambda wc: config["mlrfits"][wc.mlrfit]["date_start"],
        date_end=lambda wc: config["mlrfits"][wc.mlrfit]["date_end"],
        min_counts=lambda wc: config["mlrfits"][wc.mlrfit]["min_counts"],
        keep_not_in_library=lambda wc: config["mlrfits"][wc.mlrfit]["keep_not_in_library"],
        keep_insufficient_counts=lambda wc: config["mlrfits"][wc.mlrfit]["keep_insufficient_counts"],
        **config["mlr_settings"],
    log:
        notebook="results/mlr/mlr_{protset}_{mlrfit}.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/mlr.py.ipynb"


rule growth_vs_titers:
    """Compare MLR estimated growth advantages to measured titers."""
    input:
        growth="results/mlr/growth_advantages_{protset}_{mlrfit}.csv",
        titers=lambda wc: config["sera"][wc.sera]["csv"],
        muts_from_mrca=config["muts_from_mrca"],
    output:
        chart="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}.html",
        simple_corr_chart="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_corr.html",
        simple_cutoff_chart="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_cutoff.html",
        muts_vs_titers_chart="results/growth_vs_titers/muts_vs_titers_{protset}_{mlrfit}_{sera}.html",
        scatter_csv="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_scatter.csv",
        corr_csv="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_corr.csv",
    params:
        **config["growth_vs_titer_params"],
        sera_regex=lambda wc: config["sera"][wc.sera]["sera_regex"],
        pool=lambda wc: config["sera"][wc.sera]["pool"],
        simple_plot_scale=config["simple_plot_scale"],
    log:
        notebook="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/growth_vs_titers.py.ipynb"


rule compare_mlr_fits:
    """Compare the different MLR fit results."""
    input:
        growth_advantages=expand(
            "results/mlr/growth_advantages_{protset}_{mlrfit}.csv",
            protset=config["protsets"],
            mlrfit=config["mlrfits"],
        ),
    output:
        corr_chart="results/compare_mlr_fits/mlrfits_corr.html",
        scatter_chart="results/compare_mlr_fits/mlrfits_scatter.html",
    params:
        protsets_mlrfits=list(itertools.product(config["protsets"], config["mlrfits"])),
    log:
        notebook="results/compare_mlr_fits/compare_mlr_fits.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/compare_mlr_fits.py.ipynb"


rule charts_to_docs:
    """Copy and write all the charts to a `./docs/` subdirectory for GitHub Pages."""
    input:
        extract_final_values(charts),
        expand(
            [
                "results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_corr.html",
                "results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}_cutoff.html",
                "results/growth_vs_titers/muts_vs_titers_{protset}_{mlrfit}_{sera}.html",
            ],
            protset=config["protsets"],
            mlrfit=config["mlrfits"],
            sera=config["sera"],
        ),
    output:
        docsdir=directory("docs"),
    params:
        charts=charts,
        title=config["docs_title"],
        heading=config["docs_heading"],
    log:
        "results/charts_to_docs.txt",
    conda:
        "environment.yml"
    script:
        "scripts/charts_to_docs.py"
