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


LINEAGES = {
    "Pan": ["PPA", "PTE", "PTS", "PTT", "PTV", "PT"],
    "Pongo": ["PA", "PP"],
    "Gorilla": ["GBB", "GBG", "GGG", "GB"],
}

POPULATIONS = {
    "Pan": ["PPA", "PTE", "PTS", "PTT", "PTV"],
    "Pongo": ["PA", "PP"],
    "Gorilla": ["GBB", "GBG", "GGG"],
}


rule all:
    input:
        expand(
            "results/selscan/all/all.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/betascan/all/all.m_{core_freq}.b1.top.{cutoff}.candidate.genes.overlap",
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/plots/selscan/all/all.{method}_{maf}.top.{cutoff}.candidate.genes.svg",
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/plots/betascan/all/all.m_{core_freq}.b1.top.{cutoff}.candidate.genes.svg",
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
        ),


rule selscan_overlap_across_species:
    input:
        genes=sum(
            [
                expand(
                    "results/selscan/{sp}/lineages/{method}_{maf}/candidates/{ln}.{method}_{maf}.top.{cutoff}.candidate.genes",
                    sp=[sp],
                    ln=LINEAGES[sp],
                    allow_missing=True,
                )
                for sp in LINEAGES
            ],
            [],
        ),
    output:
        genes="results/selscan/all/all.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
    script:
        "scripts/get_overlap_genes_across_lineages.py"


rule betascan_overlap_across_species:
    input:
        genes=sum(
            [
                expand(
                    "results/betascan/{sp}/{pp}/m_{core_freq}/candidates/{pp}.b1.top.{cutoff}.candidate.genes",
                    sp=[sp],
                    pp=POPULATIONS[sp],
                    allow_missing=True,
                )
                for sp in LINEAGES
            ],
            [],
        ),
    output:
        genes="results/betascan/all/all.m_{core_freq}.b1.top.{cutoff}.candidate.genes.overlap",
    script:
        "scripts/get_overlap_genes_across_lineages.py"


rule plot_selscan_overlap_across_species:
    input:
        pan_genes="results/selscan/Pan/lineages/{method}_{maf}/all/Pan.{method}_{maf}.top.{cutoff}.candidate.genes",
        pongo_genes="results/selscan/Pongo/lineages/{method}_{maf}/all/Pongo.{method}_{maf}.top.{cutoff}.candidate.genes",
        gorilla_genes="results/selscan/Gorilla/lineages/{method}_{maf}/all/Gorilla.{method}_{maf}.top.{cutoff}.candidate.genes",
    output:
        plot="results/plots/selscan/all/all.{method}_{maf}.top.{cutoff}.candidate.genes.svg",
    script:
        "scripts/plot_venn_diagram.py"


rule plot_betascan_overlap_across_species:
    input:
        pan_genes="results/betascan/Pan/all/Pan.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
        pongo_genes="results/betascan/Pongo/all/Pongo.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
        gorilla_genes="results/betascan/Gorilla/all/Gorilla.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
    output:
        plot="results/plots/betascan/all/all.m_{core_freq}.b1.top.{cutoff}.candidate.genes.svg",
    script:
        "scripts/plot_venn_diagram.py"
