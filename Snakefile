"""``snakemake`` file that runs analysis."""


configfile: "config.yaml"


rule all:
    """Target rule."""
    input:
        expand(
            "results/strain_counts/{prot_set}_counts_by_date.csv",
            prot_set=config["prot_sets"],
        )


rule strain_counts:
    """Get counts of each strain."""
    input:
        strain_prots=config["strain_prots"],
        prot_set=lambda wc: config["prot_sets"][wc.prot_set]["prot_set"],
    output:
        counts_overall="results/strain_counts/{prot_set}_counts_overall.csv",
        counts_by_date="results/strain_counts/{prot_set}_counts_by_date.csv",
        strain_matches="results/strain_counts/{prot_set}_strain_matches.csv",
    params:
        trim_strain_prots=lambda wc: config["prot_sets"][wc.prot_set]["trim"],
    log:
        "results/logs/strain_counts_{prot_set}.txt",
    conda:
        "environment.yml"
    script:
        "scripts/strain_counts.py"
