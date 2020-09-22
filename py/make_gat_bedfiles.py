import argparse
import os
import sys

import pandas as pd


def create_files(foreground_df,
                 background_df,
                 colnames=["chrom_bed", "start", "end"]):

    foreground_file = "foreground.bed"
    background_file = "background.bed"

    foreground_df[colnames].to_csv(
        foreground_file, sep='\t', index=False, header=False)
    background_df[colnames].to_csv(
        background_file, sep='\t', index=False, header=False)


def prepare_for_gat(asbs):
    dataframe = asbs.copy()
    dataframe["chrom_bed"] = dataframe.apply(
        lambda row: row["chrom"][3:], axis=1, result_type='reduce')
    # BED files are zero-indexed, whereas all things BaalChIP are
    # 1-indexed, since based on VCF/SAM
    dataframe["start"] = dataframe["pos"] - 1
    dataframe["end"] = dataframe["pos"]
    dataframe = dataframe.drop_duplicates()
    dataframe = dataframe.sort_values(by=["chrom_bed", "start"])
    return dataframe


def fetch_asb_set(input_file):
    input_snps = pd.read_csv(input_file)
    input_snps.rename(columns={"CHROM": "chrom", "POS": "pos"}, inplace=True)
    asbs = input_snps[(input_snps["isASB"]) & (input_snps["peak"])]
    return prepare_for_gat(asbs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=('Run GAT enrichment analysis for '
                     'each set of ASB SNPs individually'))
    parser.add_argument(
        'baal_file', help=('output file from overlap_beds.py,'
                           ' containing snps detected as ASB'))
    parser.add_argument(
        'snp_file', help='file containing all het SNPs in the cell line')
    args = parser.parse_args()

    all_snps = prepare_for_gat(pd.read_csv(args.snp_file,
                                           skiprows=1,
                                           sep="\t",
                                           names=["rsid", "chrom",
                                                  "pos", "ref", "alt", "raf"],
                                           dtype={"pos": int}))
    sig_snps = fetch_asb_set(args.baal_file)

    # No significant SNPs were found. This is not necessarily an error
    # but should produce no GAT analysis file.
    if len(sig_snps) == 0:
        sys.exit(os.EX_DATAERR)

    create_files(sig_snps, all_snps)
