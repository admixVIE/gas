# Copyright 2025 Xin Huang and Simon Chen
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
import matplotlib.pyplot as plt
import numpy as np

enrichment_data = pd.read_csv(
    snakemake.input["enrichment"],
    sep="\t",
    header=None,
    names=[
        "GO_ID",
        "avg_genes_sim",
        "genes_found",
        "p_value",
        "p_adjusted",
        "genes_uniq",
        "genes_max",
        "genes_total",
        "description",
        "gene_list",
    ],
)

plot_data = enrichment_data.sort_values("p_adjusted").head(20).copy()

plot_data["short_description"] = plot_data["description"].apply(
    lambda x: x[:47] + "..." if len(x) > 50 else x
)


def get_color(p_val):
    if p_val < 0.001:
        return "#0000FF"
    elif p_val < 0.01:
        return "#800080"
    elif p_val < 0.05:
        return "#FF0000"
    else:
        return "#808080"


plot_data["color"] = plot_data["p_adjusted"].apply(get_color)

plot_data = plot_data.sort_values("genes_found", ascending=True)

fig, ax = plt.subplots(figsize=(12, max(8, len(plot_data) * 0.4)))

bars = ax.barh(
    range(len(plot_data)), plot_data["genes_found"], color=plot_data["color"]
)

ax.set_yticks(range(len(plot_data)))
ax.set_yticklabels(plot_data["short_description"], fontsize=9)
ax.set_xlabel("Count", fontsize=11)

plt.suptitle(
    f"GO Enrichment Analysis (Gowinda) - Top {len(plot_data)} terms",
    fontsize=14,
    fontweight="bold",
    x=0.45,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", alpha=0.3)

from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor="#0000FF", label="< 0.001"),
    Patch(facecolor="#800080", label="< 0.01"),
    Patch(facecolor="#FF0000", label="< 0.05"),
    Patch(facecolor="#808080", label=">= 0.05"),
]
ax.legend(handles=legend_elements, title="p.adjust", loc="lower right")

plt.tight_layout()
plt.savefig(snakemake.output["plot"], dpi=300, bbox_inches="tight")
plt.close()
