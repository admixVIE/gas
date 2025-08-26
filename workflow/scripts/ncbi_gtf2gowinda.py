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


import re
import gzip


def nc_to_chr(nc_accession):
    mapping = {
        "NC_000001": "chr1",
        "NC_000002": "chr2",
        "NC_000003": "chr3",
        "NC_000004": "chr4",
        "NC_000005": "chr5",
        "NC_000006": "chr6",
        "NC_000007": "chr7",
        "NC_000008": "chr8",
        "NC_000009": "chr9",
        "NC_000010": "chr10",
        "NC_000011": "chr11",
        "NC_000012": "chr12",
        "NC_000013": "chr13",
        "NC_000014": "chr14",
        "NC_000015": "chr15",
        "NC_000016": "chr16",
        "NC_000017": "chr17",
        "NC_000018": "chr18",
        "NC_000019": "chr19",
        "NC_000020": "chr20",
        "NC_000021": "chr21",
        "NC_000022": "chr22",
    }
    key = nc_accession.split(".")[0]
    return mapping.get(key, None)


input_file = snakemake.input.gtf
output_file = snakemake.output.gtf

with gzip.open(input_file, "rt") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom_accession = fields[0]
        feature_type = fields[2]
        attr_field = fields[8]

        chr_name = nc_to_chr(chrom_accession)
        if chr_name is None:
            continue
        if feature_type != "gene":
            continue

        match = re.search(r'gene_id "[^"]+";', attr_field)
        if match:
            gene_id_str = match.group(0)
            fields[0] = chr_name
            fields[2] = "exon"
            out_line = "\t".join(fields[:8]) + "\t" + gene_id_str + "\n"
            outfile.write(out_line)
