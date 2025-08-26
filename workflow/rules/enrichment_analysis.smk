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


rule convert_ncbi_gtf:
    input:
        gtf=rules.download_ncbi_gene_data.output.gtf,
    output:
        gtf="results/annotated_data/hg38.gowinda.gtf",
    script:
        "../scripts/ncbi_gtf2gowinda.py"


rule convert_ncbi_go:
    input:
        gtf=rules.download_ncbi_gene_data.output.gtf,
        gene2go=rules.download_ncbi_gene_data.output.gene2go,
    output:
        go2gene="results/annotated_data/hg38.gowinda.go2gene",
    params:
        species="9606",
    script:
        "../scripts/ncbi_go2gowinda.py"


rule selscan_functional_enrichment:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates=rules.consolidate_lineages.output.candidates,
        total=expand(
            "results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz",
            i=np.arange(1, 23),
            pair=pair,
            allow_missing=True,
        ),
    output:
        candidate_snps="results/enrichment/selscan/{species}/lineages/{method}_{maf}/{lineage}.{method}_{maf}.top.{cutoff}.candidate.snps.tsv",
        total_snps="results/enrichment/selscan/{species}/lineages/{method}_{maf}/{lineage}.{method}_{maf}.top.{cutoff}.total.snps.tsv",
        enrichment="results/enrichment/selscan/{species}/lineages/{method}_{maf}/{lineage}.{method}_{maf}.top.{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=64,
        cpus=32,
    shell:
        """
        sed '1d' {input.candidates} | awk '{{print "chr"$1"\\t"$2}}' > {output.candidate_snps}
        for i in {input.total}; do
            bcftools query -f "%CHROM\\t%POS\\n" $i
        done > {output.total_snps}

        java -Xmx{resources.mem_gb}g -jar {input.gowinda} \
            --snp-file {output.total_snps} \
            --candidate-snp-file {output.candidate_snps} \
            --gene-set-file {input.go2gene} \
            --annotation-file {input.gtf} \
            --simulations 1000000 \
            --min-significance 1 \
            --gene-definition gene \
            --threads {resources.cpus} \
            --output-file {output.enrichment} \
            --mode gene \
            --min-genes 1 || true
        """


rule betascan_functional_enrichment:
    input:
        gowinda=rules.download_gowinda.output.gowinda,
        go2gene=rules.convert_ncbi_go.output.go2gene,
        gtf=rules.convert_ncbi_gtf.output.gtf,
        candidates="results/betascan/{species}/{lineage}/m_{core_freq}/candidates/{lineage}.b1.top.{cutoff}.annotated.candidates",
        total=expand(
            "results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz",
            i=np.arange(1, 23),
            ppl=ppl,
            allow_missing=True,
        ),
    output:
        candidate_snps="results/enrichment/betascan/{species}/{lineage}/m_{core_freq}/{lineage}.b1.top.{cutoff}.candidate.snps.tsv",
        total_snps="results/enrichment/betascan/{species}/{lineage}/m_{core_freq}/{lineage}.b1.top.{cutoff}.total.snps.tsv",
        enrichment="results/enrichment/betascan/{species}/{lineage}/m_{core_freq}/{lineage}.b1.top.{cutoff}.gowinda.enrichment.tsv",
    resources:
        mem_gb=64,
        cpus=32,
    shell:
        """
        sed '1d' {input.candidates} | awk '{{print "chr"$1"\\t"$2}}' > {output.candidate_snps}
        for i in {input.total}; do
            bcftools query -f "%CHROM\\t%POS\\n" $i
        done > {output.total_snps}

        java -Xmx{resources.mem_gb}g -jar {input.gowinda} \
            --snp-file {output.total_snps} \
            --candidate-snp-file {output.candidate_snps} \
            --gene-set-file {input.go2gene} \
            --annotation-file {input.gtf} \
            --simulations 1000000 \
            --min-significance 1 \
            --gene-definition gene \
            --threads {resources.cpus} \
            --output-file {output.enrichment} \
            --mode gene \
            --min-genes 1 || true
        """
