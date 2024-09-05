"""Get counts of sequences for each strain."""


import re
import sys

import Bio.SeqIO

import pandas as pd


sys.stdout = sys.stderr = log = open(snakemake.log[0], "w")

aas = "ACDEFGHIKLMNPQRSTVWY"

trim_start, trim_end = snakemake.params.trim_strain_prots
assert 1 <= trim_start < trim_end

print(f"Reading strain_prots from {snakemake.input.strain_prots}")
print(f"Trimming to {trim_start} to {trim_end}")
strain_prots = {}
for seq in Bio.SeqIO.parse(snakemake.input.strain_prots, "fasta"):
    name = seq.id.split("|")[0]
    if name in strain_prots:
        raise ValueError(f"Duplicate {name=} in {snakemake.input.strain_prots}")
    s = str(seq.seq).upper()
    while s.endswith("*"):  # remove trailing stop codons
        s = s[: -1]
    assert re.fullmatch(f"[{aas}]+", s), f"{name=} has invalid residues:\n{s}\n{seq}"
    assert len(s) >= trim_end, f"{trim_end=}, {len(s)=}"
    strain_prots[name] = s[trim_start - 1: trim_end]
print(f"Read {len(strain_prots)=} strain proteins\n")
assert len(strain_prots) == len(set(strain_prots.values()))

print(f"Reading protein set from {snakemake.input.prot_set}")
strain_match_records = []
assert "other" not in strain_prots
accessions = set()
for prot in Bio.SeqIO.parse(snakemake.input.prot_set, "fasta"):

    p = str(prot.seq).upper()

    try:
        accession, _, _, date = (t.strip() for t in prot.description.split("|"))
    except ValueError:
        raise ValueError(f"Problem parsing header {prot.description=}")

    assert pd.notnull(accession), prot.description
    if accession in accessions:
        continue  # duplicate accession
    accessions.add(accession)

    for strain_name, strain_seq in strain_prots.items():
        if strain_seq in p:
            strain_match_records.append((strain_name, accession, date))
            break
    else:
        strain_match_records.append(("other", accession, date))

strain_matches = pd.DataFrame(
    strain_match_records, columns=["variant", "accession", "date"]
)

n_acc = strain_matches["accession"].nunique()
assert len(strain_matches) == n_acc, f"{len(strain_matches)=}, {n_acc=}" 

print(f"Read {len(strain_matches)} sequences.")

overall_counts = (
    strain_matches
    .groupby("variant", as_index=False)
    .aggregate(n_sequences=pd.NamedAgg("accession", "nunique"))
    .merge(
        pd.Series(strain_prots).rename_axis("variant").rename("sequence").reset_index(),
        on="variant",
        how="outer",
        validate="one_to_one",
    )
    .assign(n_sequences=lambda x: x["n_sequences"].fillna(0).astype(int))
    .sort_values("n_sequences", ascending=False)
)

print("Here are overall number matching each strain:")
print(overall_counts)
print(f"Writing to {snakemake.output.counts_overall}")
overall_counts.to_csv(snakemake.output.counts_overall, index=False)

print(f"\nWriting all strain matches to {snakemake.output.strain_matches}")
strain_matches.to_csv(snakemake.output.strain_matches, index=False)

print(f"\nWriting strain counts by date to {snakemake.output.counts_by_date}")
(
    strain_matches
    .groupby(["variant", "date"], as_index=False)
    .aggregate(sequences=pd.NamedAgg("accession", "nunique"))
    .to_csv(snakemake.output.counts_by_date, index=False)
)
