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


ruleorder: selscan_overlap_across_lineages > selscan_overlap_with_previous_studies


rule annotate_selscan_xp_candidates:
    input:
        outliers="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.candidates.scores",
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.hg38_multianno.txt",
            i=list(range(1, 23)),
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    script:
        "../scripts/get_annotated_candidates.py"


rule split_xp_candidates_by_sign:
    input:
        annotated=rules.annotate_selscan_xp_candidates.output.annotated_candidates,
    output:
        positive_candidates="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.positive.annotated.candidates",
        negative_candidates="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.negative.annotated.candidates",
        positive_genes="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.positive.candidate.genes",
        negative_genes="results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.negative.candidate.genes",
    shell:
        """
        head -1 {input.annotated} > {output.positive_candidates}
        awk 'NR > 1 && $NF > 0' {input.annotated} >> {output.positive_candidates} || touch {output.positive_candidates}        
        head -1 {input.annotated} > {output.negative_candidates}
        awk 'NR > 1 && $NF < 0' {input.annotated} >> {output.negative_candidates} || touch {output.negative_candidates}
        
        awk 'NR > 1 && $NF > 0 {{print $7}}' {input.annotated} | grep -v ";" | sort | uniq > {output.positive_genes} || touch {output.positive_genes}
        awk 'NR > 1 && $NF < 0 {{print $7}}' {input.annotated} | grep -v ";" | sort | uniq > {output.negative_genes} || touch {output.negative_genes}
        """


rule consolidate_lineages:
    input:
        expand(
            "results/selscan/{species}/2pop/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.{sign}.annotated.candidates",
            pair=pair,
            sign=["positive", "negative"],
            allow_missing=True,
        ),
    params:
        base_dir="results/selscan/{species}/2pop",
    output:
        candidates="results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.annotated.candidates",
        genes="results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes",
    script:
        "../scripts/consolidate_lineages.py"


rule selscan_overlap_with_previous_studies:
    input:
        genes=rules.consolidate_lineages.output.genes,
    output:
        overlap="results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
    params:
        Cagan2016="resources/candidates/Cagan2016/positive_selection",
        Mattle_Greminger2018="resources/candidates/Mattle-Greminger2018/positive_selection",
        Nye2020="resources/candidates/Nye2020/positive_selection",
        Zhao2023="resources/candidates/Zhao2023/positive_selection",
        vanderValk2024="resources/candidates/vanderValk2024/positive_selection",
        Yoo2025="resources/candidates/Yoo2025/positive_selection",
    script:
        "../scripts/get_overlap_genes_with_previous_studies.py"


rule selscan_overlap_across_lineages:
    input:
        genes=expand(
            "results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes",
            species=species,
            lineage=lineages_for_species,
            allow_missing=True,
        ),
    output:
        genes="results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.candidate.genes.overlap",
    script:
        "../scripts/get_overlap_genes_across_lineages.py"


rule merge_all_selscan_candidate_genes:
    input:
        genes=expand(
            "results/selscan/{species}/lineages/{method}_{maf}/candidates/{lineage}.{method}_{maf}.top.{cutoff}.candidate.genes",
            species=species,
            lineage=lineages_for_species,
            allow_missing=True,
        ),
    output:
        genes="results/selscan/{species}/lineages/{method}_{maf}/all/{species}.{method}_{maf}.top.{cutoff}.candidate.genes",
    shell:
        """
        cat {input.genes} | grep -v gene | awk '{{print $1}}' | sort | uniq > {output.genes}
        """


rule annotate_betascan_candidates:
    input:
        outliers="results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.{cutoff}.candidates.scores",
        annotation=expand(
            "results/annotated_data/{species}/all/chr{i}.biallelic.snps.hg38_multianno.txt",
            i=list(range(1, 23)),
            allow_missing=True,
        ),
    output:
        annotated_candidates="results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.{cutoff}.annotated.candidates",
    resources:
        mem_gb=32,
    script:
        "../scripts/get_annotated_candidates.py"


rule get_betascan_candidate_genes:
    input:
        betascan_candidates=rules.annotate_betascan_candidates.output.annotated_candidates,
    output:
        betascan_genes="results/betascan/{species}/{ppl}/m_{core_freq}/candidates/{ppl}.b1.top.{cutoff}.candidate.genes",
    shell:
        """
        sed '1d' {input.betascan_candidates} | grep -v ";" | awk '{{print $7}}' | sort | uniq -c | awk '{{print $2"\\t"$1}}' | sed '1igene\\tsnp_count' > {output.betascan_genes} || true
        """


rule betascan_overlap_with_previous_studies:
    input:
        genes="results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.candidate.genes",
    output:
        overlap="results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.candidate.genes.overlap",
    params:
        Cagan2016="resources/candidates/Cagan2016/balancing_selection",
        Zhao2023="resources/candidates/Zhao2023/balancing_selection",
    script:
        "../scripts/get_overlap_genes_with_previous_studies.py"


rule betascan_overlap_across_lineages:
    input:
        genes=expand(
            "results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.candidate.genes",
            species=species,
            lineage=ppl,
            allow_missing=True,
        ),
    output:
        genes="results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.candidate.genes.overlap",
    script:
        "../scripts/get_overlap_genes_across_lineages.py"


rule merge_all_betascan_candidate_genes:
    input:
        genes=expand(
            "results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.candidate.genes",
            species=species,
            lineage=ppl,
            allow_missing=True,
        ),
    output:
        genes="results/betascan/{species}/all/{species}.m_{core_freq}.b1.top.{cutoff}.candidate.genes",
    shell:
        """
        cat {input.genes} | grep -v gene | awk '{{print $1}}' | sort | uniq > {output.genes}
        """
