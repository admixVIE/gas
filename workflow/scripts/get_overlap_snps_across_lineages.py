# Copyright 2026 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import os
import pandas as pd


def load_positions(paths):
    """Load candidate SNP positions, grouped by species/population."""
    positions = {}

    for p in paths:
        species = os.path.basename(p).split(".")[0]

        try:
            df = pd.read_csv(
                p,
                sep="\t",
                dtype={"Chr": int, "Start": int},
            )

            species_positions = set(zip(df["Chr"], df["Start"]))
            positions[species] = positions.get(species, set()) | species_positions

        except pd.errors.EmptyDataError:
            positions.setdefault(species, set())

    return positions


# Load overlapping genes
gene_overlap = pd.read_csv(snakemake.input.gene_overlap, sep="\t")
overlap_genes = set(gene_overlap.iloc[:, 0].tolist())

# Load primary and secondary candidate positions, grouped by species/population
primary_label = snakemake.params.primary_label
secondary_label = snakemake.params.secondary_label

primary_positions = load_positions(snakemake.input.primary_candidates)
secondary_positions = (
    load_positions(snakemake.input.secondary_candidates) if secondary_label else {}
)

species_list = list(primary_positions.keys())
rows = []

# For each species/population, collect SNPs belonging to overlap genes only
for focal_species in species_list:
    for p in snakemake.input.primary_candidates:
        species = os.path.basename(p).split(".")[0]

        if species != focal_species:
            continue

        try:
            df = pd.read_csv(
                p,
                sep="\t",
                dtype={"Chr": int, "Start": int},
            )

            df = df[df["Gene.refGene"].isin(overlap_genes)]

            for _, row in df.iterrows():
                pos = (row["Chr"], row["Start"])

                result = {
                    "Chr": pos[0],
                    "Start": pos[1],
                    "Gene": row["Gene.refGene"],
                }

                for species in species_list:
                    result[f"{primary_label}_{species}"] = pos in primary_positions.get(
                        species, set()
                    )

                    if secondary_label:
                        result[f"{secondary_label}_{species}"] = (
                            pos in secondary_positions.get(species, set())
                        )

                rows.append(result)

        except pd.errors.EmptyDataError:
            continue


df = pd.DataFrame(rows).drop_duplicates()

df[f"{primary_label}_count"] = df[
    [f"{primary_label}_{species}" for species in species_list]
].sum(axis=1)

if secondary_label:
    df[f"{secondary_label}_count"] = df[
        [f"{secondary_label}_{species}" for species in species_list]
    ].sum(axis=1)

bool_cols = df.select_dtypes(include=["bool", "boolean"]).columns
df[bool_cols] = df[bool_cols].replace({True: "Y", False: "N"})

df.sort_values(["Chr", "Start", "Gene"]).to_csv(
    snakemake.output.snps,
    sep="\t",
    index=False,
)
