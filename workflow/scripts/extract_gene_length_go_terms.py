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


import gzip
import re
from collections import defaultdict


gtf_file = snakemake.input.genomic_gtf
gene2go_file = snakemake.input.gene2go
tax_id_target = snakemake.params.species

geneid2go = defaultdict(set)

with gzip.open(gene2go_file, "rt") as f:
    header = next(f)
    for line in f:
        fields = line.strip().split("\t")
        tax_id, geneid, go_id = fields[0], fields[1], fields[2]
        if tax_id == tax_id_target:
            geneid2go[geneid].add(go_id)

gene_info = {}

with gzip.open(gtf_file, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue
        start = int(fields[3])
        end = int(fields[4])
        gene_length = end - start + 1

        gene_id_match = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_id = gene_id_match.group(1) if gene_id_match else "NA"

        geneid_match = re.search(r'db_xref "GeneID:(\d+)"', fields[8])
        geneid = geneid_match.group(1) if geneid_match else "NA"

        if geneid in geneid2go:
            go_terms = ";".join(sorted(geneid2go[geneid]))
        else:
            go_terms = "NA"

        gene_info[gene_id] = (geneid, gene_length, go_terms)

with open(snakemake.output.txt, "w") as out:
    out.write("gene_id\tGeneID\tgene_length\tGO_terms\n")
    for gene_id, (geneid, gene_length, go_terms) in gene_info.items():
        out.write(f"{gene_id}\t{geneid}\t{gene_length}\t{go_terms}\n")
