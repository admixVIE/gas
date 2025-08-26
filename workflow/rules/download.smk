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


rule download_greatapes_genomes:
    output:
        vcf="resources/data/{species}/chr{i}.filteranno.vcf.gz",
        index="resources/data/{species}/chr{i}.filteranno.vcf.gz.csi",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/{wildcards.species}/{wildcards.species}_wild_filtered/chr{wildcards.i}.filteranno.vcf.gz -O {output.vcf}
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/{wildcards.species}/{wildcards.species}_wild_filtered/chr{wildcards.i}.filteranno.vcf.gz.csi -O {output.index}
        """


rule download_greatapes_metadata:
    output:
        txt="resources/data/metadata_full.txt",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/metadata_full.txt -O {output.txt}
        """


rule download_reference_genomes:
    output:
        hg38="resources/data/refgenomes/hg38.fa",
        rheMac10="resources/data/refgenomes/rheMac10.fa",
        chain="resources/data/refgenomes/hg38ToRheMac10.over.chain",
    params:
        dir="resources/data/refgenomes",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O {params.dir}/hg38.fa.gz
        gzip -d {params.dir}/hg38.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz -O {params.dir}/rheMac10.fa.gz
        gzip -d {params.dir}/rheMac10.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToRheMac10.over.chain.gz -O {params.dir}/hg38ToRheMac10.over.chain.gz
        gzip -d {params.dir}/hg38ToRheMac10.over.chain.gz
        """


rule download_annovar_db:
    output:
        refgene="resources/tools/annovar/humandb/hg38_refGene.txt",
        avsnp150="resources/tools/annovar/humandb/hg38_avsnp150.txt",
        dbnsfp42c="resources/tools/annovar/humandb/hg38_dbnsfp42c.txt",
    shell:
        """
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar avsnp150 resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp42c resources/tools/annovar/humandb/
        """


rule download_ncbi_gene_data:
    output:
        gtf="resources/data/refgenomes/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
        gene2go="resources/data/refgenomes/gene2go.gz",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -O {output.gtf}
        wget -c https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz -O {output.gene2go}
        """


rule download_repeats:
    output:
        rmsk="resources/data/repeats/rmsk.txt.gz",
        segdup="resources/data/repeats/genomicSuperDups.txt.gz",
        simrep="resources/data/repeats/simpleRepeat.txt.gz",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -O {output.rmsk}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz -O {output.segdup}
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz -O {output.simrep}
        """


rule download_selscan:
    output:
        selscan="resources/tools/selscan/selscan",
        norm="resources/tools/selscan/norm",
    shell:
        """
        wget -c https://github.com/szpiech/selscan/archive/refs/tags/v2.0.3.tar.gz
        tar -xvzf v2.0.3.tar.gz
        mv selscan-2.0.3/bin/linux/selscan-2.0.3 {output.selscan}
        mv selscan-2.0.3/bin/linux/norm {output.norm}
        rm v2.0.3.tar.gz
        rm -rf selscan-2.0.3/
        """


rule download_betascan:
    output:
        betascan="resources/tools/betascan/BetaScan.py",
    shell:
        """
        git clone https://github.com/ksiewert/BetaScan
        mv BetaScan/BetaScan.py {output.betascan}
        rm -rf BetaScan
        """


rule download_gowinda:
    output:
        gowinda="resources/tools/gowinda/Gowinda-1.12.jar",
    shell:
        """
        wget -c https://sourceforge.net/projects/gowinda/files/Gowinda-1.12.jar/download -O {output.gowinda}
        """
