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


import os
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet


gene_path = snakemake.input.genes

genes = {}
all_genes = set()

for p in gene_path:
    filename = os.path.basename(p)
    lineage = filename.split(".")[0]
    try:
        df = pd.read_csv(p, sep="\t")
        genes[lineage] = df.iloc[:, 0].tolist()
    except pd.errors.EmptyDataError:
        genes[lineage] = []

    all_genes.update(g for g in genes[lineage])

lineages = list(genes.keys())
rows = []

for g in all_genes:
    row = {"gene": g}
    for l in lineages:
        row[l] = True if g in genes[l] else False
    rows.append(row)

df = pd.DataFrame(rows).set_index(lineages)

plt.figure(figsize=(8,6))
UpSet(df, subset_size="count", sort_categories_by="input").plot()
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
