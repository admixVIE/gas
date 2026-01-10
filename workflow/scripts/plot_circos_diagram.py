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


import re
import pandas as pd
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Function to load and filter BED data
def load_bed_without_sex_chr(bed_file):
    df = pd.read_csv(bed_file, sep="\t", header=None, dtype={0: str})
    df.columns = ["chr", "start", "end", "gene"]
    df["chr"] = df["chr"].astype(str).str.strip()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[~df["chr"].isin(["chrX", "chrY"])]
    return df


title = snakemake.params.title
title = title.replace("xpehh", "XP-EHH").replace("xpnsl", "XP-nSL")
title = title.replace("top 0.0005", "top 0.05%").replace("top 0.00005", "top 0.005%")


# Read gene BED files from each species
pan_df = load_bed_without_sex_chr(snakemake.input.pan_bed)
pongo_df = load_bed_without_sex_chr(snakemake.input.pongo_bed)
gorilla_df = load_bed_without_sex_chr(snakemake.input.gorilla_bed)

# Initialize circos sectors from chromosome BED
chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")

# Read and filter to only autosomes (chr1-22)
chr_bed_df = pd.read_csv(chr_bed_file, sep="\t", header=None)
autosomes = [f"chr{i}" for i in range(1, 23)]
chr_bed_df = chr_bed_df[chr_bed_df[0].isin(autosomes)]

# Create temporary filtered file
filtered_chr_bed = snakemake.output.plot.replace(".svg", "_chr_filtered.bed").replace(".png", "_chr_filtered.bed")
chr_bed_df.to_csv(filtered_chr_bed, sep="\t", header=False, index=False)

# Initialize circos with filtered chromosomes
circos = Circos.initialize_from_bed(filtered_chr_bed, space=2)

# Add cytoband tracks only for autosomes
circos.add_cytoband_tracks((81, 85), cytoband_file)

# Track positions (inner to outer)
track_configs = [
    ("Pan", pan_df, (65, 80), "#1f77b4"),    # Blue
    ("Pongo", pongo_df, (50, 65), "#ff7f0e"),  # Orange
    ("Gorilla", gorilla_df, (35, 50), "#2ca02c"),  # Green
]

for sector in circos.sectors:
    # Outer chromosome axis track
    axis_track = sector.add_track((81, 85), r_pad_ratio=0.0)
    axis_track.axis(lw=0.6)
    sector.text(sector.name, r=88, size=8)

    # Add tracks for each species
    for species_name, df, (r_inner, r_outer), color in track_configs:
        gene_track = sector.add_track((r_inner, r_outer), r_pad_ratio=0.0)
        gene_track.axis(lw=0.3)

        # Filter genes on this chromosome
        sub = df[df["chr"] == sector.name]

        # Draw gene blocks
        for s, e in zip(sub["start"].to_numpy(), sub["end"].to_numpy()):
            gene_track.rect(float(s), float(e), fc=color, ec="none", lw=1)

circos.text(title, size=10, r=0)

# Create figure and add legend
fig = circos.plotfig()

# Create legend with unique gene counts
legend_elements = [
    mpatches.Patch(color="#1f77b4", label=f"Pan (n={pan_df['gene'].nunique()})"),
    mpatches.Patch(color="#ff7f0e", label=f"Pongo (n={pongo_df['gene'].nunique()})"),
    mpatches.Patch(color="#2ca02c", label=f"Gorilla (n={gorilla_df['gene'].nunique()})"),
]

# Position legend outside the plot area
legend = fig.legend(
    handles=legend_elements,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=10,
    frameon=True,
    fancybox=True,
    shadow=True,
    title="Species",
    title_fontsize=10
)

# Save figure with extra space for legend
fig.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")
plt.close()
