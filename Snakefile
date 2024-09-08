"""``snakemake`` file that runs analysis."""


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        expand(
            "results/mlr/strain_counts_{protset}_{maxdiff}.html",
            protset=config["protsets"],
            maxdiff=config["protset_maxdiffs"],
        )


rule strain_counts:
    """Get counts of each strain."""
    input:
        strain_prots=config["strain_prots"],
        protset=lambda wc: config["protsets"][wc.protset]["protset"],
    output:
        counts_overall="results/strain_counts/{protset}_{maxdiff}_counts_overall.csv",
        counts_by_date="results/strain_counts/{protset}_{maxdiff}_counts_by_date.csv",
        strain_matches="results/strain_counts/{protset}_{maxdiff}_strain_matches.csv",
    params:
        trim_strain_prots=lambda wc: config["protsets"][wc.protset]["trim"],
    log:
        "results/logs/strain_counts_{protset}_{maxdiff}.txt",
    conda:
        "environment.yml"
    script:
        "scripts/strain_counts.py"


rule mlr:
    """Fit MLR estimates of fitness (growth) advantages."""
    input:
        counts_by_date="results/strain_counts/{protset}_{maxdiff}_counts_by_date.csv",
    output:
        counts_chart="results/mlr/strain_counts_{protset}_{maxdiff}.html",
    params:
        min_counts=config["min_counts"],
        plot_window_frame_days=config["plot_window_frame_days"],
    log:
        notebook="results/mlr/mlr_{protset}_{maxdiff}.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/mlr.py.ipynb"
