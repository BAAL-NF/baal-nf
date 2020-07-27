from argparse import ArgumentParser

import pyranges as pr
import pandas as pd

if __name__ == '__main__':
    parser = ArgumentParser(description='identify SNPs in BaalChIP output that overlap with peaks from a BED file')
    parser.add_argument('asb_file', help='output file from BaalChIP')
    parser.add_argument('bed_file', help='bed file with peak data')
    parser.add_argument('output', help=('File to store output in. The output format will be as asb_file, '
                                        'but with an additional boolean column indicating whether a SNP '
                                        'was found to overlap with one or more peaks in the input BED file'))
    args = parser.parse_args()
    
    asb_df = pd.read_csv(args.asb_file, index_col=0)
    
    input_cols = list(asb_df.columns)
    
    asb_df["Start"] = asb_df["POS"]
    asb_df["End"] = asb_df["POS"]+1
    asb_df["Chromosome"] = asb_df["CHROM"]

    asb = pr.PyRanges(asb_df)

    peaks = pr.read_bed(args.bed_file)
    peaks = peaks.insert(pd.Series(name="peak", data=[1]*len(peaks)))

    intersection = asb.join(peaks, how='left')
    result_df = intersection.as_df()
   
    result_df["peak"] = result_df.apply(lambda row: row["peak"] == 1, axis=1)
    result_df[input_cols + ["peak"]].to_csv(args.output)



