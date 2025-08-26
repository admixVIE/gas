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


import numpy as np


rule extract_pair_snps:
    input:
        vcf=rules.polarize_2pop.output.vcf,
        pair_info=rules.create_pair_info.output.pair_info,
    output:
        vcf1="results/polarized_data/{species}/2pop/{pair}/pop1.chr{i}.biallelic.snps.vcf.gz",
        vcf2="results/polarized_data/{species}/2pop/{pair}/pop2.chr{i}.biallelic.snps.vcf.gz",
        idx1="results/polarized_data/{species}/2pop/{pair}/pop1.chr{i}.biallelic.snps.vcf.gz.tbi",
        idx2="results/polarized_data/{species}/2pop/{pair}/pop2.chr{i}.biallelic.snps.vcf.gz.tbi",
        map="results/polarized_data/{species}/2pop/{pair}/chr{i}.biallelic.snps.map",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)

        bcftools view {input.vcf} -S <(grep -w $pop1 {input.pair_info} | awk '{{print $1}}') --force-samples | bgzip -c > {output.vcf1}
        bcftools view {input.vcf} -S <(grep -w $pop2 {input.pair_info} | awk '{{print $1}}') --force-samples | bgzip -c > {output.vcf2}
        bcftools query -f "%CHROM\\t%CHROM:%POS:%REF:%ALT\\t%POS\\t%POS\\n" {input.vcf} > {output.map}
        tabix -p vcf {output.vcf1}
        tabix -p vcf {output.vcf2}
        """


rule estimate_selscan_xp_scores:
    input:
        selscan=rules.download_selscan.output.selscan,
        vcf1=rules.extract_pair_snps.output.vcf1,
        vcf2=rules.extract_pair_snps.output.vcf2,
        map=rules.extract_pair_snps.output.map,
    output:
        out="results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.chr{i}.{method}.out",
        formatted_out="results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.chr{i}.{method}.formatted.out",
    params:
        output_prefix="results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.chr{i}",
    resources:
        cpus=8,
    shell:
        """
        {input.selscan} --vcf {input.vcf1} --vcf-ref {input.vcf2} --map {input.map} --{wildcards.method} --out {params.output_prefix} --threads {resources.cpus} --maf {wildcards.maf} --unphased
        head -1 {output.out} > {output.formatted_out}
        sed '1d' {output.out} | awk -v chr={wildcards.i} 'BEGIN{{OFS="\\t"}}{{print chr,$2,$3,$4,$5,$6,$7,$8}}' >> {output.formatted_out}
        """


rule normalize_selscan_xp_scores:
    input:
        norm=rules.download_selscan.output.norm,
        scores=expand(
            "results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.chr{i}.{method}.formatted.out",
            i=np.arange(1, 23),
            allow_missing=True,
        ),
    output:
        scores=expand(
            "results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.chr{i}.{method}.formatted.out.norm",
            i=np.arange(1, 23),
            allow_missing=True,
        ),
        log="results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.normalized.{method}.log",
    shell:
        """
        {input.norm} --files {input.scores} --log {output.log} --{wildcards.method}
        """


rule merge_selscan_xp_scores:
    input:
        scores=rules.normalize_selscan_xp_scores.output.scores,
    output:
        merged_scores="results/selscan/{species}/2pop/{pair}/{method}_{maf}/scores/{pair}.normalized.{method}.scores",
    shell:
        """
        cat {input.scores} | grep -v id | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$9}}' | sed '1iSNP\\tCHR\\tBP\\tnormalized_{wildcards.method}' > {output.merged_scores}
        """
