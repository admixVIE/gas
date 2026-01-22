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



# Merge population-level annotated.candidates to species-level
def get_selscan_annotated_candidates(wildcards):
    lineages = LINEAGES[wildcards.species]
    return expand(
        "results/selscan/{{species}}/lineages/{{method}}_{{maf}}/candidates/{lineage}.{{method}}_{{maf}}.top.{{cutoff}}.annotated.candidates",
        lineage=lineages
    )


def get_betascan_annotated_candidates(wildcards):
    populations = POPULATIONS[wildcards.species]
    return expand(
        "results/betascan/{{species}}/{ppl}/m_{{core_freq}}/candidates/{ppl}.b1.top.{{cutoff}}.annotated.candidates",
        ppl=populations
    )


rule merge_selscan_annotated_candidates:
    input:
        candidates=get_selscan_annotated_candidates,
    output:
        merged="results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.annotated.candidates",
    shell:
        """
        head -1 {input.candidates[0]} > {output.merged}
        for file in {input.candidates}; do
            tail -n +2 $file >> {output.merged}
        done
        """


rule merge_betascan_annotated_candidates:
    input:
        candidates=get_betascan_annotated_candidates,
    output:
        merged="results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.annotated.candidates",
    shell:
        """
        head -1 {input.candidates[0]} > {output.merged}
        for file in {input.candidates}; do
            tail -n +2 $file >> {output.merged}
        done
        """


# Extract and validate gene regions for circos
rule extract_selscan_gene_regions_for_circos:
    input:
        refgene="resources/tools/annovar/humandb/hg38_refGene.txt",
        genes="results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.candidate.genes",
        candidates=rules.merge_selscan_annotated_candidates.output.merged,
    output:
        bed="results/plots/selscan/{species}/circos/{species}.{method}_{maf}.top.{cutoff}.genes.bed",
    script:
        "scripts/extract_gene_regions.py"


rule extract_betascan_gene_regions_for_circos:
    input:
        refgene="resources/tools/annovar/humandb/hg38_refGene.txt",
        genes="results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
        candidates=rules.merge_betascan_annotated_candidates.output.merged,
    output:
        bed="results/plots/betascan/{species}/circos/{species}.m_{core_freq}.b1.top.{cutoff}.genes.bed",
    script:
        "scripts/extract_gene_regions.py"


# Plot circos
rule plot_selscan_overlap_across_species:
    input:
        pan_bed="results/plots/selscan/Pan/circos/Pan.{method}_{maf}.top.{cutoff}.genes.bed",
        pongo_bed="results/plots/selscan/Pongo/circos/Pongo.{method}_{maf}.top.{cutoff}.genes.bed",
        gorilla_bed="results/plots/selscan/Gorilla/circos/Gorilla.{method}_{maf}.top.{cutoff}.genes.bed",
    output:
        plot="results/plots/selscan/all/all.{method}_{maf}.top.{cutoff}.candidate.genes.svg",
    params:
        title="{method}\n(top {cutoff})",
    script:
        "scripts/plot_circos_diagram.py"


rule plot_betascan_overlap_across_species:
    input:
        pan_bed="results/plots/betascan/Pan/circos/Pan.m_{core_freq}.b1.top.{cutoff}.genes.bed",
        pongo_bed="results/plots/betascan/Pongo/circos/Pongo.m_{core_freq}.b1.top.{cutoff}.genes.bed",
        gorilla_bed="results/plots/betascan/Gorilla/circos/Gorilla.m_{core_freq}.b1.top.{cutoff}.genes.bed",
    output:
        plot="results/plots/betascan/all/all.m_{core_freq}.b1.top.{cutoff}.candidate.genes.svg",
    params:
        title="B1\n(top {cutoff})",
    script:
        "scripts/plot_circos_diagram.py"
