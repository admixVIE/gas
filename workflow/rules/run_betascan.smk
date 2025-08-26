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


rule test_hwe:
    input:
        vcf=rules.polarize_1pop.output.vcf,
    output:
        hwe_outliers="results/polarized_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.hwe.outliers",
    params:
        output_prefix="results/polarized_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps",
    shell:
        """
        plink --vcf {input.vcf} --hardy --out {params.output_prefix} --set-missing-var-ids @:#
        awk '$7>$8' {params.output_prefix}.hwe | \
        sed '1d' | \
        awk '$9<0.001{{print $2}}' | \
        awk -F ":" '{{print $1"\\t"$2}}' > {output.hwe_outliers}
        """


rule convert_repeat_files:
    input:
        rmsk=rules.download_repeats.output.rmsk,
        segdup=rules.download_repeats.output.segdup,
        simrep=rules.download_repeats.output.simrep,
    output:
        rmsk="results/processed_data/repeats/hg38.rmsk.autosomes.bed",
        segdup="results/processed_data/repeats/hg38.seg.dups.autosomes.bed",
        simrep="results/processed_data/repeats/hg38.simple.repeats.autosomes.bed",
    shell:
        """
        zcat {input.rmsk} | awk 'BEGIN{{OFS="\\t"}}$6!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $6,$7,$8,$11,$2,$10}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.rmsk}
        zcat {input.segdup} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$6,$7}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.segdup}
        zcat {input.simrep} | awk 'BEGIN{{OFS="\\t"}}$2!~/chr(X|Y|Un|M|[0-9]_|[0-9][0-9]_)/{{print $2,$3,$4,$5,$11}}' | sed 's/^chr//' | sort -k1,1n -k2,2n -k3,3n > {output.simrep}
        """


rule get_allele_counts:
    input:
        vcf=rules.polarize_1pop.output.vcf,
        hwe_outliers=rules.test_hwe.output.hwe_outliers,
        rmsk=rules.convert_repeat_files.output.rmsk,
        segdup=rules.convert_repeat_files.output.segdup,
        simrep=rules.convert_repeat_files.output.simrep,
    output:
        ac="results/betascan/{species}/{ppl}/ac/chr{i}.ac",
    shell:
        """
        if [ -s {input.hwe_outliers} ]; then \
            bcftools view -T ^{input.rmsk} {input.vcf} | \
            bcftools view -T ^{input.simrep} | \
            bcftools view -T ^{input.segdup} | \
            bcftools view -T ^{input.hwe_outliers} | \
            bcftools view -i "(INFO/AC>2*N_SAMPLES*0.05) && (INFO/AC<2*N_SAMPLES*0.95)" | \
            bcftools query -f "%POS\\t%INFO/AC\\t%INFO/AN\\n" > {output.ac}; \
        else \
            bcftools view -T ^{input.rmsk} {input.vcf} | \
            bcftools view -T ^{input.simrep} | \
            bcftools view -T ^{input.segdup} | \
            bcftools view -i "(INFO/AC>2*N_SAMPLES*0.05) && (INFO/AC<2*N_SAMPLES*0.95)" | \
            bcftools query -f "%POS\\t%INFO/AC\\t%INFO/AN\\n" > {output.ac}; \
        fi
        """


rule estimate_b1_scores:
    input:
        ac=rules.get_allele_counts.output.ac,
        betascan=rules.download_betascan.output.betascan,
    output:
        scores="results/betascan/{species}/{ppl}/m_{core_freq}/scores/chr{i}.b1.scores",
    shell:
        """
        python {input.betascan} -i {input.ac} -m {wildcards.core_freq} | grep -v Position | awk -v chr={wildcards.i} '{{print chr"\\t"$0}}' > {output.scores}
        """


rule merge_b1_scores:
    input:
        scores=expand(
            "results/betascan/{species}/{ppl}/m_{core_freq}/scores/chr{i}.b1.scores",
            i=np.arange(1, 23),
            allow_missing=True,
        ),
    output:
        merged_scores="results/betascan/{species}/{ppl}/m_{core_freq}/scores/{ppl}.b1.scores",
    shell:
        """
        cat {input.scores} | awk '{{print $1":"$2"\\t"$1"\\t"$2"\\t"$3}}' | sed '1iSNP\\tCHR\\tBP\\tB1' > {output.merged_scores}
        """
