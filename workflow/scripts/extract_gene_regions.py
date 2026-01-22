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


import pandas as pd

# Inputs files from Snakemake
refgene_txt = snakemake.input.refgene
gene_list = snakemake.input.genes
candidates = snakemake.input.candidates
out_bed = snakemake.output.bed

# Read gene list
genes = pd.read_csv(gene_list, sep="\t", header=None)[0].astype(str)
gene_set = set(g.strip() for g in genes.tolist() if g.strip())

# Read outlier SNPs from annotated candidates
snps = pd.read_csv(candidates, sep="\t")

# Filter: keep only gene-associated SNPs (exclude intergenic)
snps = snps[snps["Func.refGene"] != "intergenic"].copy()

# Prepare SNP data
snps["Chr"] = snps["Chr"].astype(str).str.strip()
snps["Start"] = snps["Start"].astype(int)
snps["End"] = snps["End"].astype(int)

# Read ANNOVAR refGene (genePred format)
rg = pd.read_csv(
    refgene_txt,
    sep="\t",
    header=None,
    usecols=[2, 4, 5, 12],  # chr, txStart, txEnd, gene_name
    names=["chr", "start", "end", "gene"],
    dtype={"chr": str, "start": int, "end": int, "gene": str},
)

rg["gene"] = rg["gene"].astype(str).str.strip()
rg["chr"] = rg["chr"].astype(str).str.strip()

# Keep only chr1..chr22 (autosomes only for circos)
autosomes = {f"chr{i}" for i in range(1, 23)}
rg = rg[rg["chr"].isin(autosomes)].copy()

# Filter to candidate genes only
rg = rg[rg["gene"].isin(gene_set)].copy()

# Validate each gene region: keep only if it contains outlier SNPs
validated_regions = []

for _, region in rg.iterrows():
    chr_str = region["chr"]
    chr_num = chr_str.replace("chr", "")
    start = region["start"]
    end = region["end"]
    gene = region["gene"]
    
    # Check if ANY outlier SNPs fall within this gene region
#    snps_in_region = snps[
#        (snps["Chr"] == chr_num) &
#        (snps["Start"] >= start) &
#        (snps["End"] <= end)
#    ]
    # Check SNPs: inside gene OR upstream/downstream of its gene
    snps_in_region = snps[
        (snps["Chr"] == chr_num) &
        (
            # Case 1: Inside gene boundaries (intronic, exonic, UTR, etc.)
            ((snps["Start"] >= start) & (snps["End"] <= end)) |
            
            # Case 2: Upstream/downstream
            ((snps["Func.refGene"].isin(["upstream", "downstream"])) &
             (snps["Gene.refGene"] == gene))
        )
    ]
    
    # Only keep regions with at least 1 outlier SNP
    if len(snps_in_region) > 0:
        validated_regions.append([chr_str, start, end, gene])

# Create BED dataframe
bed = pd.DataFrame(validated_regions, columns=["chr", "start", "end", "gene"])

# Merge overlapping regions for the same gene on same chromosome
bed = (
    bed.groupby(["chr", "gene"], as_index=False)
    .agg(start=("start", "min"), end=("end", "max"))
    .sort_values(["chr", "start"])
)

bed = bed[["chr", "start", "end", "gene"]]

# Save BED file (no header)
bed.to_csv(out_bed, sep="\t", header=False, index=False)
