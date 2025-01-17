from argparse import ArgumentParser
from pathlib import Path
import dask.dataframe as dd
import pandas as pd
import numpy as np
import sys
import os
import motifutils as mu

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--tf", type = str, help = "Gene symbol for TF with ASB events"
    )
    parser.add_argument(
        "--mode", type = str, nargs = '?', default = "allMotifs", help = "Mode to compile under, depending on which motifs were found for the TF-of-interest"
    )
    parser.add_argument(
        "--out_dir", type = str, nargs = '?', default = ".", help = "Output directory"
    ) 

    return parser.parse_args()

if __name__ == '__main__':
    tf = get_input().tf
    mode = get_input().mode
    out_dir = get_input().out_dir

    snp_peak_data, snp_data = mu.get_snp_data(".", tf)
    nopeak_data = mu.get_nopeak(".", tf, mode)
    jaspar_data = mu.get_jaspar(".", tf, mode)

    # Filter for concordant & discordant bQTLs
    concordant_jaspar, concordant_nopeak = mu.subset_asbs(jaspar_data, nopeak_data, tf, concordant = True, isasb = True)
    discordant_jaspar, discordant_nopeak = mu.subset_asbs(jaspar_data, nopeak_data, tf, concordant = False, isasb = True)

    # Identify SNPs that map to high-quality motifs but are not ASB sites
    nonasbs_jaspar, nonasbs_nopeak = mu.subset_asbs(jaspar_data, nopeak_data, tf, isasb = False)

    # Get all concordant ASBs by row
    is_concordant = snp_data.apply(
            lambda row: mu.map_snp_by_row(
                row,
                concordant_jaspar,
                concordant_nopeak,
                discordant_jaspar,
                discordant_nopeak,
                nonasbs_jaspar,
                nonasbs_nopeak,
                snp_data
            ),
            axis=1
        )

    # Join with original full asb dataframe
    final = pd.concat(is_concordant.tolist(), ignore_index = True, sort = False)
    final['ASB_quality'] = final.apply(mu.classify_asb_quality, axis=1)
    final = final.convert_dtypes() # Convert NAs to be consistent

    # Order column names the same
    colnames = snp_data.columns.tolist() + ['High_quality_motif','Concordant','ASB_quality','Motifs','Motif_group']
    final = final[colnames]

    os.makedirs(out_dir, exist_ok = True)
    final.to_csv(f"{out_dir}/{tf}.withPeaks.withMotifs.csv", index=False)