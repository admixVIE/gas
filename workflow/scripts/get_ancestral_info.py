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


import os
import argparse
import pysam
from pyliftover import LiftOver


def get_fasta_positions(fasta_path: str, chrom: str) -> list:
    """
    Extracts all positions from a reference genome FASTA file without filtering any bases.
    """
    positions = []
    try:
        fasta = pysam.FastaFile(fasta_path)
        seq_length = fasta.get_reference_length(chrom)
        positions = list(range(seq_length))
        fasta.close()
    except Exception as e:
        print(f"Error reading FASTA for {chrom}: {e}")
    return positions


def liftover_positions(chrom: str, positions: list, lo: LiftOver) -> dict:
    """
    Performs LiftOver coordinate conversion while preserving source chromosome names.
    Returns a dictionary mapping (source chrom, source pos) -> (target chrom, target pos, strand).
    """
    lifted = {}
    for pos in positions:
        result = lo.convert_coordinate(chrom, pos)  # 0-based
        if not result:
            continue  # Skip unmapped positions

        new_chrom, new_pos, strand = result[0][
            :3
        ]  # Extract (target chrom, position, strand)
        if new_pos < 1:
            continue  # Skip invalid positions

        lifted[(chrom, pos)] = (new_chrom, new_pos, strand)
    return lifted


def complement_base(base: str) -> str:
    """Returns the complementary DNA base for reverse strand correction."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complement.get(base, "N")  # Default to 'N' if an unexpected base appears


def get_fasta_sequences(lifted_positions: dict, tgt_fasta: str) -> dict:
    """
    Extracts target genome bases at the LiftOver-mapped positions.
    Filters out positions where the mapped base is not A, C, G, or T.
    Converts bases if mapped to the reverse strand.
    """
    sequences = {}
    try:
        fasta = pysam.FastaFile(tgt_fasta)
        if not lifted_positions:
            return sequences

        for (src_chrom, src_pos), (
            tgt_chrom,
            tgt_pos,
            strand,
        ) in lifted_positions.items():
            if tgt_chrom not in fasta.references:
                print(
                    f"Warning: Target chromosome {tgt_chrom} not found in {tgt_fasta}, skipping"
                )
                continue

            seq_length = fasta.get_reference_length(tgt_chrom)
            if tgt_pos >= seq_length:
                print(
                    f"Warning: Target position {tgt_pos} is out of range for {tgt_chrom}, skipping"
                )
                continue

            # Fetch single base (0-based start, 1-based end)
            base = fasta.fetch(tgt_chrom, tgt_pos, tgt_pos + 1).upper()

            if base not in "ACGT":
                continue  # Filter out positions with N or other characters

            if strand == "-":
                base = complement_base(base)

            sequences[(src_chrom, src_pos)] = base

        fasta.close()
    except Exception as e:
        print(f"Error reading FASTA: {e}")

    return sequences


def process_chromosome(
    chrom: str, src_fasta: str, tgt_fasta: str, output_file: str, liftover_chain: str
):
    """
    Processes a single chromosome: extracts positions from the source genome,
    maps them to the target genome, retrieves target genome sequences,
    and saves results in the specified output file using source coordinates.
    """
    lo = LiftOver(liftover_chain)

    positions = get_fasta_positions(src_fasta, chrom)
    if not positions:
        print(f"Skipping {chrom} (no valid positions).")
        return

    lifted_positions = liftover_positions(chrom, positions, lo)
    if not lifted_positions:
        print(f"Skipping {chrom} (no successful LiftOver).")
        return

    sequences = get_fasta_sequences(lifted_positions, tgt_fasta)

    with open(output_file, "w") as f:
        for (src_chrom, src_pos), base in sequences.items():
            f.write(f"{src_chrom}\t{src_pos}\t{src_pos+1}\t{base}\n")  # BED format

    print(f"Done processing {chrom}, output saved to {output_file}")


def main():
    """
    Main function to convert genomic coordinates using LiftOver.

    Command-line Usage:
    -------------------
    python get_ancestral_info.py --src-fasta source.fa --tgt-fasta target.fa \
        --liftover-chain sourceToTarget.over.chain --chr-name 10 \
        --output /data/output/chr10.anc.alleles.bed
    """
    parser = argparse.ArgumentParser(
        description="Convert genomic coordinates using LiftOver while preserving source coordinates."
    )
    parser.add_argument(
        "--src-fasta",
        type=str,
        required=True,
        help="Path to the source reference genome (e.g., source.fa).",
    )
    parser.add_argument(
        "--tgt-fasta",
        type=str,
        required=True,
        help="Path to the target genome (e.g., target.fa).",
    )
    parser.add_argument(
        "--liftover-chain",
        type=str,
        required=True,
        help="Path to the LiftOver chain file.",
    )
    parser.add_argument(
        "--chr-name",
        type=str,
        required=True,
        help="Chromosome to process (e.g., 10 or X).",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to the output BED file."
    )

    args = parser.parse_args()

    process_chromosome(
        args.chr_name, args.src_fasta, args.tgt_fasta, args.output, args.liftover_chain
    )


if __name__ == "__main__":
    main()
