"""``snakemake`` file that runs analysis."""


import itertools


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        "results/compare_mlr_fits/mlrfits_corr.html",
        "results/compare_mlr_fits/mlrfits_scatter.html",
        expand(
            "results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}.html",
            protset=config["protsets"],
            mlrfit=config["mlrfits"],
            sera=config["sera"],
        ),


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
    output:
        chart="results/growth_vs_titers/growth_vs_titers_{protset}_{mlrfit}_{sera}.html",
    params:
        **config["growth_vs_titer_params"],
        sera_regex=lambda wc: config["sera"][wc.sera]["sera_regex"],
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
