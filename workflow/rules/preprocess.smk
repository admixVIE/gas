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


rule extract_biallelic_snps:
    input:
        vcf=rules.download_greatapes_genomes.output.vcf,
    output:
        vcf="results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz",
        index="results/processed_data/{species}/all/chr{i}.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        if [ "{wildcards.species}" != "Pan" ]; then
            bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bcftools view -e 'COUNT(GT="het") == N_SAMPLES' | bgzip -c > {output.vcf}
        else
            bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        fi
        tabix -p vcf {output.vcf}
        """


rule create_pair_info:
    input:
        metadata=rules.download_greatapes_metadata.output.txt,
    output:
        pair_info="results/dadi/{species}/jdfe/{pair}/pop.pair.list",
    shell:
        """
        pair="{wildcards.pair}"
        pop1=$(echo $pair | cut -d'_' -f1)
        pop2=$(echo $pair | cut -d'_' -f2)
        grep -v captive {input.metadata} | awk -v pop=$pop1 '$2 ~ ("^" pop "$"){{print $4"\\t"$2}}' > {output.pair_info}
        grep -v captive {input.metadata} | awk -v pop=$pop2 '$2 ~ ("^" pop "$"){{print $4"\\t"$2}}' >> {output.pair_info}
        """


rule annotate_biallelic_snps:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        refgene=rules.download_annovar_db.output.refgene,
        avsnp150=rules.download_annovar_db.output.avsnp150,
        dbnsfp42c=rules.download_annovar_db.output.dbnsfp42c,
    output:
        avinput="results/annotated_data/{species}/all/chr{i}.biallelic.snps.avinput",
        txt="results/annotated_data/{species}/all/chr{i}.biallelic.snps.hg38_multianno.txt",
    resources:
        cpus=8,
        mem_gb=32,
    params:
        output_prefix="results/annotated_data/{species}/all/chr{i}.biallelic.snps",
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\n" {input.vcf} > {output.avinput}
        resources/tools/annovar/table_annovar.pl {output.avinput} resources/tools/annovar/humandb/ -buildver hg38 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . --thread {resources.cpus}
        """


rule extract_pop_data:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        metadata=rules.download_greatapes_metadata.output.txt,
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.vcf.gz",
        index="results/processed_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        if [ "{wildcards.species}" != "Pan" ]; then
            bcftools view {input.vcf} -S <(sed '1d' {input.metadata} |\
            grep -w {wildcards.ppl} |\
            awk '{{print $4}}') --force-samples |\
            bcftools view -e "INFO/AC==0 || INFO/AC==2*N_SAMPLES" |\
            bcftools view -e 'COUNT(GT="het") == N_SAMPLES' |\
            bcftools annotate -x ^FORMAT/GT | bgzip -c > {output.vcf}
        else
            bcftools view {input.vcf} -S <(sed '1d' {input.metadata} |\
            grep -w {wildcards.ppl} |\
            awk '{{print $4}}') --force-samples |\
            bcftools view -e "INFO/AC==0 || INFO/AC==2*N_SAMPLES" |\
            bcftools annotate -x ^FORMAT/GT | bgzip -c > {output.vcf}
        fi
        tabix -p vcf {output.vcf}
        """


rule extract_pair_data:
    input:
        vcf=rules.extract_biallelic_snps.output.vcf,
        samples=rules.create_pair_info.output.pair_info,
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/chr{i}.biallelic.snps.vcf.gz",
        index="results/processed_data/{species}/2pop/{pair}/chr{i}.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        if [ "{wildcards.species}" != "Pan" ]; then
            bcftools view {input.vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples |\
            bcftools view -e "INFO/AC==0 || INFO/AC==2*N_SAMPLES" |\
            bcftools view -e 'COUNT(GT="het") == N_SAMPLES' |\
            bcftools annotate -x ^FORMAT/GT | bgzip -c > {output.vcf}
        else
            bcftools view {input.vcf} -S <(awk '{{print $1}}' {input.samples}) --force-samples |\
            bcftools view -e "INFO/AC==0 || INFO/AC==2*N_SAMPLES" |\
            bcftools annotate -x ^FORMAT/GT | bgzip -c > {output.vcf}
        fi
        tabix -p vcf {output.vcf}
        """


rule extract_anc_info:
    input:
        hg38=rules.download_reference_genomes.output.hg38,
        rheMac10=rules.download_reference_genomes.output.rheMac10,
        chain=rules.download_reference_genomes.output.chain,
    output:
        anc_alleles="results/polarized_data/Human/anc_info/hg38.chr{i}.anc.alleles.bed.gz",
        index="results/polarized_data/Human/anc_info/hg38.chr{i}.anc.alleles.bed.gz.tbi",
    params:
        out="results/polarized_data/Human/anc_info/hg38.chr{i}.anc.alleles.bed",
    resources:
        mem_gb=128,
    shell:
        """
        python workflow/scripts/get_ancestral_info.py --src-fasta {input.hg38} --tgt-fasta {input.rheMac10} --liftover-chain {input.chain} --chr-name chr{wildcards.i} --output {params.out}
        bgzip -c {params.out} > {output.anc_alleles}
        tabix -p bed {output.anc_alleles}
        """


rule polarize_1pop:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        anc_alleles=rules.extract_anc_info.output.anc_alleles,
    output:
        vcf="results/polarized_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.vcf.gz",
        index="results/polarized_data/{species}/1pop/{ppl}/chr{i}.biallelic.snps.vcf.gz.tbi",
    resources:
        mem_gb=32,
    script:
        "../scripts/polarize_snps.py"


rule polarize_2pop:
    input:
        vcf=rules.extract_pair_data.output.vcf,
        anc_alleles=rules.extract_anc_info.output.anc_alleles,
    output:
        vcf="results/polarized_data/{species}/2pop/{pair}/chr{i}.biallelic.snps.vcf.gz",
        index="results/polarized_data/{species}/2pop/{pair}/chr{i}.biallelic.snps.vcf.gz.tbi",
    resources:
        mem_gb=32,
    script:
        "../scripts/polarize_snps.py"


rule extract_1pop_exonic_data:
    input:
        vcf=rules.extract_pop_data.output.vcf,
        anno=rules.annotate_biallelic_snps.output.txt,
        anc_alleles=rules.extract_anc_info.output.anc_alleles,
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/chr{i}.biallelic.{mut_type}.snps.vcf.gz",
        index="results/processed_data/{species}/1pop/{ppl}/chr{i}.biallelic.{mut_type}.snps.vcf.gz.tbi",
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_gb=32,
    shell:
        """
        bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bcftools annotate -a {input.anc_alleles} -h <(echo "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">") -c CHROM,FROM,TO,AA |\
            bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule concat_1pop_exonic_data:
    input:
        vcfs=expand(
            "results/processed_data/{species}/1pop/{ppl}/chr{i}.biallelic.{mut_type}.snps.vcf.gz",
            i=np.arange(1, 23),
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.vcf.gz",
        index="results/processed_data/{species}/1pop/{ppl}/{ppl}.biallelic.{mut_type}.snps.vcf.gz.tbi",
    shell:
        """
        bcftools concat {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_2pop_exonic_data:
    input:
        vcf=rules.extract_pair_data.output.vcf,
        anno=rules.annotate_biallelic_snps.output.txt,
        anc_alleles=rules.extract_anc_info.output.anc_alleles,
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/chr{i}.biallelic.{mut_type}.snps.vcf.gz",
        index="results/processed_data/{species}/2pop/{pair}/chr{i}.biallelic.{mut_type}.snps.vcf.gz.tbi",
    params:
        condition=lambda wildcards: (
            "$9~/^synonymous/"
            if wildcards.mut_type == "syn"
            else "$9~/^nonsynonymous/"
        ),
    resources:
        mem_gb=32,
    shell:
        """
        bcftools view {input.vcf} -R <(awk '{params.condition}{{print $1"\\t"$2}}' {input.anno}) |\
            bcftools annotate -a {input.anc_alleles} -h <(echo "##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">") -c CHROM,FROM,TO,AA |\
            bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule concat_2pop_exonic_data:
    input:
        vcfs=expand(
            "results/processed_data/{species}/2pop/{pair}/chr{i}.biallelic.{mut_type}.snps.vcf.gz",
            i=np.arange(1, 23),
            allow_missing=True,
        ),
    output:
        vcf="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.vcf.gz",
        index="results/processed_data/{species}/2pop/{pair}/{pair}.biallelic.{mut_type}.snps.vcf.gz.tbi",
    shell:
        """
        bcftools concat {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
