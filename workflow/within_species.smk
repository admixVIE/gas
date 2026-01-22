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


import numpy as np


species = config["species"]
ppl = config["populations"]
pair = config["population_pairs"]


def get_species_lineages(species):
    lineages = {
        "Pan": [
            "PPA",
            "PT",
            "PTE",
            "PTS",
            "PTT",
            "PTV",
        ],
        "Gorilla": [
            "GB",
            "GGG",
            "GBB",
            "GBG",
        ],
        "Pongo": ["PA", "PP"],
    }
    return lineages.get(species, [])


lineages_for_species = get_species_lineages(species)


rule all:
    input:
        expand(
            "results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.{cutoff}.candidate.genes",
            species=species,
            ppl=ppl,
            core_freq=[0.15],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.biv_{dfe}.InferDFE.bestfits",
            species=species,
            pair=pair,
            demog=["split_mig"],
            dfe=["lognormal"],
        ),
        expand(
            "results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.png",
            species=species,
            ppl=ppl,
            demog=["two_epoch"],
        ),
        expand(
            "results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.{dfe}.png",
            species=species,
            ppl=ppl,
            demog=["two_epoch"],
            dfe=["lognormal"],
        ),
        expand(
            "results/dadi/{species}/dfe/{ppl}/StatDFE/{ppl}.{demog}.{dfe}.godambe.ci",
            species=species,
            ppl=ppl,
            demog=["two_epoch"],
            dfe=["lognormal"],
        ),
        expand(
            "results/dadi/{species}/jdfe/{pair}/StatDFE/{pair}.{demog}.biv_{dfe}.godambe.ci",
            species=species,
            pair=pair,
            demog=["split_mig"],
            dfe=["lognormal"],
        ),
        expand(
            "results/plots/dadi/{species}/dfe/{ppl}/{ppl}.{demog}.{dfe}.mut.prop.png",
            species=species,
            ppl=ppl,
            demog=["two_epoch"],
            dfe=["lognormal"],
        ),
        expand(
            "results/dadi/{species}/jdfe/{pair}/InferDFE/{pair}.{demog}.constants",
            species=species,
            pair=pair,
            demog=["split_mig"],
            dfe=["lognormal"],
        ),
        expand(
            "results/plots/enrichment/selscan/{species}/lineages/{method}_{maf}/{lineage}.{method}_{maf}.top.{cutoff}.gowinda.enrichment.png",
            species=species,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
            lineage=lineages_for_species,
        ),
        expand(
            "results/enrichment/betascan/{species}/{lineage}/m_{core_freq}/{lineage}.b1.top.{cutoff}.gowinda.enrichment.tsv",
            species=species,
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
            lineage=ppl,
        ),
        expand(
            "results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
            species=species,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
            lineage=lineages_for_species,
        ),
        expand(
            "results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.candidate.genes.overlap",
            species=species,
            lineage=ppl,
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
            species=species,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
            lineage=lineages_for_species,
        ),
        expand(
            "results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.candidate.genes.overlap",
            species=species,
            lineage=ppl,
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/plots/selscan/lineages/{species}.{method}_{maf}.top.{cutoff}.candidate.genes.png",
            species=species,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.candidate.genes",
            species=species,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
            species=species,
            core_freq=["0.15"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/plots/selscan/{species}/2pop/{pair}/{method}_{maf}/{pair}.normalized.{method}.scores.png",
            species=species,
            pair=pair,
            method=["xpehh", "xpnsl"],
            maf=["0.05"],
            cutoff=["0.0005", "0.00005"],
        ),
        expand(
            "results/plots/betascan/{species}/{ppl}/m_{core_freq}/{ppl}.b1.scores.png",
            species=species,
            ppl=ppl,
            core_freq=["0.15"],
        ),
        expand(
            "results/plots/dadi/{species}/jdfe/{pair}/{pair}.{demog}.png",
            species=species,
            pair=pair,
            demog=["split_mig"],
        ),
        expand(
            "results/plots/dadi/{species}/jdfe/{pair}/{pair}.{demog}.biv_{dfe}.png",
            species=species,
            pair=pair,
            demog=["split_mig"],
            dfe=["lognormal"],
        ),


include: "rules/download.smk"
include: "rules/preprocess.smk"
include: "rules/run_selscan_xp.smk"
include: "rules/run_betascan.smk"
include: "rules/run_dadi_dfe.smk"
include: "rules/run_dadi_jdfe.smk"
include: "rules/annotate_candidates.smk"
include: "rules/enrichment_analysis.smk"
include: "rules/make_plots.smk"
