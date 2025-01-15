from pathlib import Path
from argparse import ArgumentParser
import os
import pandas as pd

# Custom libraries
import motifutils as mu
from tfomics import ReferenceGenome

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--tf_motifs", type = str, help="Gene symbol for TF motifs-of-interest"
    )
    parser.add_argument(
        "--tf_asb", type = str, help = "Gene symbol for TF with ASB events"
    )
    parser.add_argument(
        "--genomepy_dir", type = Path, help="Directory where you have installed genomepy index for your genome assembly"
    )
    parser.add_argument(
        "--assembly", type = str, help="Genome assembly, ex: hg19"
    )
    parser.add_argument(
        "--reference_offset", nargs='?', default=25, type = str, help="bp window around heterozygous SNP for determining motif instance"
    )
    parser.add_argument(
        "--out_dir", type = str,  nargs='?', default=".", help="Number of CPUs used to cluster motifs"
    )
    return parser.parse_args()

if __name__ == '__main__':
    assembly = get_input().assembly
    genomepy_dir = get_input().genomepy_dir
    tf_motifs = get_input().tf_motifs
    tf_asb = get_input().tf_asb
    out_dir = get_input().out_dir
    fasta = f"{genomepy_dir}/{assembly}.fa"

    ReferenceGenome.offset = 25
    genome = ReferenceGenome(fasta) 

    # JASPAR motifs
    print("Reading in JASPAR motifs...")
    jaspar_motifs = mu.read_jaspar_motifs(".", tf_motifs)

    # Pull IC for each motif 
    ic = mu.pull_ic(jaspar_motifs)

    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    peak_data = mu.get_snps(".", tf_asb)

    print("Scoring each SNP against each motif-of-interest...")
    scores_peaks = mu.score_snps(peak_data, jaspar_motifs, tf_motifs, genome)

    print("Computing spearman's rho correlation coefficient for each motif...")
    correlations = mu.compute_correlations(scores_peaks)

    print("Filter for motifs with a significant and positive correlation coefficient rho...")
    high_quality_motifs, sig_mappings = mu.find_high_quality_motifs(scores_peaks, correlations, jaspar_motifs, tf_motifs)

    print("Reading in SNP data...")
    snp_data = mu.get_snps(".", tf_asb, filter_peak=False)

    print("Scoring all SNPs against high-quality motifs")
    scores = mu.score_snps(snp_data, high_quality_motifs, tf_motifs, genome)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = mu.compute_snp_counts(scores, snp_data, tf_motifs, tf_asb)

    # Get dataframe of all SNPs mapping to all motifs & add column information
    out_snps = mu.compile_scores_data(scores, tf_motifs)

    # Add concordant/discordant information
    concordant = mu.filter_asbs(out_snps, concordant = True)
    discordant = mu.filter_asbs(out_snps, concordant = False)
    concordant_bqtls = len(concordant.ID.unique())
    total_snps = len(out_snps.ID.unique())
    print(
        f"ASBs factor: {tf_asb}, Motifs factor: {tf_motifs}, Significant SNPs: {total_snps}, High-quality bQTLs: {concordant_bqtls}"
    ) 

    # Add Concordance information into output table
    asbs = pd.concat([concordant, discordant])
    final = pd.merge(out_snps, asbs, how = "left")
    final = final.convert_dtypes()

    # Compile motif information
    motif_inf = mu.compile_motif_data(correlations, ic, snp_counts, tf_motifs)

    # Write all CSV files
    motif_inf.to_csv(f"{out_dir}/JASPAR_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: JASPAR_Motif_metadata_{tf_motifs}.csv")

    # Export all CSV files
    print("Saving files to csv...")
    outname = "_".join([tf_motifs, "ASBs", "JASPAR", "Motifs", "AR", "score_diff.csv"])
    final.to_csv(f"{out_dir}/{outname}", index = False)
    print("Written file:", outname)
