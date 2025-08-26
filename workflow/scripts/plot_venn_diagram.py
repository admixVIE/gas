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


from matplotlib_venn import venn3
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rcParams["font.size"] = 8


def load_gene_set(filepath):
    with open(filepath, "r") as f:
        genes = {line.strip() for line in f if line.strip()}
    return genes


set_pan = load_gene_set(snakemake.input.pan_genes)
set_pongo = load_gene_set(snakemake.input.pongo_genes)
set_gorilla = load_gene_set(snakemake.input.gorilla_genes)

plt.figure(figsize=(6, 6))
venn3([set_pan, set_pongo, set_gorilla], set_labels=("Pan", "Pongo", "Gorilla"))

plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight", dpi=300)
plt.close()
