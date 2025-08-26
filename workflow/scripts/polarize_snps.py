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


import cyvcf2
import gzip
import subprocess
from cyvcf2 import VCF, Writer


anc_alleles = dict()
with gzip.open(snakemake.input.anc_alleles, "rt") as f:
    for l in f:
        e = l.rstrip().split("\t")
        anc_alleles[f"{e[0]}:{e[2]}"] = e[3]

vcf = VCF(snakemake.input.vcf)
w = Writer(snakemake.output.vcf, vcf, mode="wz")
for v in vcf:
    key = f"{v.CHROM}:{v.POS}"
    if key in anc_alleles.keys():
        if (v.REF != anc_alleles[key]) and (v.ALT[0] != anc_alleles[key]):
            continue
        if v.REF != anc_alleles[key]:
            tmp = v.ALT[0]
            v.ALT = [v.REF]
            v.REF = tmp
            v.genotypes = [[int(not g[0]), int(not g[1]), True] for g in v.genotypes]
        w.write_record(v)

w.close()
vcf.close()

subprocess.run(["tabix", "-p", "vcf", snakemake.output.vcf], check=True)
