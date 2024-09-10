"""``snakemake`` file that runs analysis."""


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        expand(
            "results/mlr/counts_{protset}_{mlrfit}.html",
            protset=config["protsets"],
            mlrfit=config["mlrfits"],
        )


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
        counts_chart="results/mlr/counts_{protset}_{mlrfit}.html",
        counts_to_fit="results/mlr/counts_to_fit_{protset}_{mlrfit}.csv",
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
