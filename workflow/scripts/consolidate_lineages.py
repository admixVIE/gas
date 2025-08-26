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


import os
import pandas as pd


def get_lineage_mapping(species):
    """Get species-specific lineage mappings"""
    mappings = {
        "Pan": {
            "PPA": [
                ("PPA_PTE", "positive"),
                ("PPA_PTS", "positive"),
                ("PPA_PTV", "positive"),
                ("PPA_PTT", "positive"),
            ],
            "PT": [
                ("PPA_PTE", "negative"),
                ("PPA_PTS", "negative"),
                ("PPA_PTV", "negative"),
                ("PPA_PTT", "negative"),
            ],
            "PTE": [
                ("PTE_PTS", "positive"),
                ("PTE_PTV", "positive"),
                ("PTE_PTT", "positive"),
            ],
            "PTS": [
                ("PTS_PTV", "positive"),
                ("PTS_PTT", "positive"),
                ("PTE_PTS", "negative"),
            ],
            "PTV": [
                ("PTV_PTT", "positive"),
                ("PTE_PTV", "negative"),
                ("PTS_PTV", "negative"),
            ],
            "PTT": [
                ("PTE_PTT", "negative"),
                ("PTS_PTT", "negative"),
                ("PTV_PTT", "negative"),
            ],
        },
        "Gorilla": {
            "GB": [("GBB_GGG", "positive"), ("GBG_GGG", "positive")],
            "GGG": [
                ("GBB_GGG", "negative"),
                ("GBG_GGG", "negative"),
            ],
            "GBB": [("GBB_GBG", "positive")],
            "GBG": [("GBB_GBG", "negative")],
        },
        "Pongo": {
            "PA": [("PA_PP", "positive")],
            "PP": [("PA_PP", "negative")],
        },
    }
    return mappings.get(species, {})


def main():
    species = snakemake.wildcards.species
    method = snakemake.wildcards.method
    maf = snakemake.wildcards.maf
    cutoff = snakemake.wildcards.cutoff
    lineage = snakemake.wildcards.lineage
    base_dir = snakemake.params.base_dir

    # Outputs
    output_candidates = snakemake.output.candidates  # single TSV of SNPs (Cf)
    output_genes = snakemake.output.genes  # aggregated per-gene TSV

    # Oriented sources for this lineage: list of (pair, orient_sign)
    lineage_mapping = get_lineage_mapping(species)
    sources = lineage_mapping[lineage]  # e.g., [("PPA_PTE","positive"), ...]
    pair_names = [pair for pair, _ in sources]

    # Containers
    P_sets, N_sets, P_rows = {}, {}, {}
    all_cols = None

    def read_candidates(pair, sign):
        """Read candidates for a given pair/sign; return DataFrame or None."""
        path = f"{base_dir}/{pair}/{method}_{maf}/candidates/{pair}.normalized.{method}.top.{cutoff}.{sign}.annotated.candidates"
        if os.path.exists(path) and os.path.getsize(path) > 0:
            try:
                df = pd.read_csv(path, sep="\t")
                if not df.empty:
                    return df
            except Exception as e:
                print(f"Error reading {path}: {e}")
        return None

    def key_set_from_df(df):
        """Return set of (Chr, Start, End)."""
        return set(
            map(tuple, df[["Chr", "Start", "End"]].itertuples(index=False, name=None))
        )

    # Read per-pair positives/negatives with orientation
    for pair, orient_sign in sources:
        # Orientation: which file corresponds to "A stronger" (P_i) vs the opposite (N_i)
        sign_P = orient_sign
        sign_N = "negative" if orient_sign == "positive" else "positive"

        df_P = read_candidates(pair, sign_P)
        df_N = read_candidates(pair, sign_N)

        # Initialize columns template
        for df in (df_P, df_N):
            if df is not None and all_cols is None:
                all_cols = df.columns.tolist()

        # Normalize P_i
        if df_P is not None:
            if not {"Chr", "Start", "End"}.issubset(df_P.columns):
                raise ValueError(
                    f"{pair} {sign_P} is missing required columns Chr/Start/End"
                )
            df_P = df_P.drop_duplicates(subset=["Chr", "Start", "End"])
            P_rows[pair] = df_P
            P_sets[pair] = key_set_from_df(df_P)
        else:
            P_rows[pair] = pd.DataFrame(columns=all_cols or ["Chr", "Start", "End"])
            P_sets[pair] = set()

        # Normalize N_i
        if df_N is not None:
            if not {"Chr", "Start", "End"}.issubset(df_N.columns):
                raise ValueError(
                    f"{pair} {sign_N} is missing required columns Chr/Start/End"
                )
            df_N = df_N.drop_duplicates(subset=["Chr", "Start", "End"])
            N_sets[pair] = key_set_from_df(df_N)
        else:
            N_sets[pair] = set()

    union_P = set().union(*P_sets.values()) if P_sets else set()
    union_N = set().union(*N_sets.values()) if N_sets else set()
    Cf_keys = union_P.difference(union_N)

    # Helper: merge all positive rows for annotations; de-duplicate by locus
    def merge_positive_rows():
        nonempty = [df for df in P_rows.values() if not df.empty]
        if not nonempty:
            return pd.DataFrame(columns=all_cols or ["Chr", "Start", "End"])
        merged = pd.concat(nonempty, ignore_index=True)
        merged = merged.drop_duplicates(subset=["Chr", "Start", "End"])
        return merged

    merged_P = merge_positive_rows()
    if merged_P.empty and Cf_keys:
        merged_P = pd.DataFrame(columns=all_cols or ["Chr", "Start", "End"])

    support_map_pairs = {}
    support_map_count = {}

    for key in Cf_keys:
        supports = [pair for pair in pair_names if key in P_sets[pair]]
        support_map_pairs[key] = ",".join(supports) if supports else ""
        support_map_count[key] = len(supports)

    if not merged_P.empty:
        merged_P["_key"] = list(
            map(
                tuple,
                merged_P[["Chr", "Start", "End"]].itertuples(index=False, name=None),
            )
        )
        cand = merged_P[merged_P["_key"].isin(Cf_keys)].copy()
        cand["support_count"] = (
            cand["_key"].map(support_map_count).fillna(0).astype(int)
        )
        cand["support_pairs"] = cand["_key"].map(support_map_pairs).fillna("")
        cand["total_comparisons"] = len(pair_names)
        cand.drop(columns=["_key"], inplace=True)
        if {"Chr", "Start", "End"}.issubset(cand.columns):
            cand = cand.sort_values(["Chr", "Start", "End"])
    else:
        cand = pd.DataFrame(
            columns=(all_cols or ["Chr", "Start", "End"])
            + ["support_count", "support_pairs", "total_comparisons"]
        )

    cand.to_csv(output_candidates, sep="\t", index=False)

    if len(cand.columns) > 6:
        gene_col = cand.columns[6]
        genes_series = cand[gene_col].dropna()
        genes_series = genes_series[genes_series != "."]
        genes_series = genes_series[~genes_series.str.contains(";", na=False)]

        if not genes_series.empty:
            # Map locus key -> support info (from cand)
            cand["_key"] = list(
                map(
                    tuple,
                    cand[["Chr", "Start", "End"]].itertuples(index=False, name=None),
                )
            )
            key_to_pairs = dict(
                zip(cand["_key"], cand["support_pairs"].fillna("").tolist())
            )

            # Build per-gene aggregation:
            # - snp_count: number of Cf SNPs mapped to the gene
            # - support_pairs: union of pairs supporting any SNP under the gene
            # - support_count: size of that union
            records = []

            for _, row in cand.iterrows():
                gene = row[gene_col]
                if pd.isna(gene) or gene == "." or (";" in str(gene)):
                    continue
                key = (row["Chr"], row["Start"], row["End"])
                pairs_str = key_to_pairs.get(key, "")
                pairs = set(filter(None, pairs_str.split(",")))
                records.append((gene, key, pairs))

            # Aggregate
            from collections import defaultdict

            gene_snp_keys = defaultdict(list)
            gene_pair_sets = defaultdict(set)

            for gene, key, pairs in records:
                gene_snp_keys[gene].append(key)
                gene_pair_sets[gene].update(pairs)

            out_gene_rows = []
            total_cmp = len(pair_names)
            for gene, keys in gene_snp_keys.items():
                pair_set = gene_pair_sets[gene]
                out_gene_rows.append(
                    {
                        "gene": gene,
                        "snp_count": len(set(keys)),
                        "support_count": len(pair_set),
                        "support_pairs": ",".join(sorted(pair_set)),
                        "total_comparisons": total_cmp,
                    }
                )

            df_genes = pd.DataFrame(out_gene_rows).sort_values(
                ["support_count", "gene"], ascending=[False, True]
            )
            df_genes.to_csv(output_genes, sep="\t", index=False)
        else:
            # No valid genes after filtering; write empty file with header
            pd.DataFrame(
                columns=[
                    "gene",
                    "snp_count",
                    "support_count",
                    "support_pairs",
                    "total_comparisons",
                ]
            ).to_csv(output_genes, sep="\t", index=False)
    else:
        # No gene column; emit empty file
        pd.DataFrame(
            columns=[
                "gene",
                "snp_count",
                "support_count",
                "support_pairs",
                "total_comparisons",
            ]
        ).to_csv(output_genes, sep="\t", index=False)

    # Log summary
    print(
        f"{species} - {lineage}: candidates={len(cand)} SNPs; genes file written ({output_genes}); comparisons={len(pair_names)}"
    )


if __name__ == "__main__":
    main()
