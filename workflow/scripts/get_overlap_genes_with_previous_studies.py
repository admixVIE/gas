# Copyright 2025 Xin Huang
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


import pandas as pd


def get_sublineages(lineage):
    mapping = {
        "GB": ["GB", "GBB", "GBG"],
        "PT": ["PT", "PTE", "PTS", "PTT", "PTV"],
    }

    return mapping.get(lineage, [lineage])


try:
    this_study = pd.read_csv(snakemake.input.genes, sep="\t")
    this_genes = this_study.iloc[:, 0].tolist()
except pd.errors.EmptyDataError:
    this_genes = []

this_genes_upper = {g: g.upper() for g in this_genes}

sublineages = get_sublineages(snakemake.wildcards.lineage)
prev_genes = {}
for study_name, base_dir in snakemake.params.items():
    union = set()
    for l in sublineages:
        path = f"{base_dir}/{l}.candidate.genes"
        try:
            prev_study = pd.read_csv(path, sep="\t", header=None)
            if prev_study.shape[0] > 0:
                genes = prev_study.iloc[:, 0].tolist()
                union.update(g.upper() for g in genes if g)
        except pd.errors.EmptyDataError:
            pass
    prev_genes[study_name] = union

study_cols = list(snakemake.params.keys())
rows = []
for g in this_genes:
    gU = this_genes_upper[g]
    row = {"gene": g}
    for s in study_cols:
        if not prev_genes.get(s):
            row[s] = "NA"
        else:
            row[s] = "Y" if gU in prev_genes[s] else "N"
    rows.append(row)

out_df = pd.DataFrame(rows, columns=["gene"] + study_cols)
out_df["count"] = (out_df[study_cols] == "Y").sum(axis=1)
out_df = out_df.sort_values(["count", "gene"], ascending=[False, True])
out_df.to_csv(snakemake.output.overlap, sep="\t", index=False)
