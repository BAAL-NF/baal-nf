from collections import defaultdict
from pathlib import Path
from argparse import ArgumentParser
import pandas as pd
import numpy as np

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
        "--n_cpus", type = int,  nargs='?', default=1, help="Number of CPUs used to cluster motifs"
    )
    parser.add_argument(
        "--out_dir", type = str,  nargs='?', default=".", help="Number of CPUs used to cluster motifs"
    )
    parser.add_argument(
        "--sub_dir", type = bool,  nargs='?', default=False, help="Organize output into subdirectories"
    )
    parser.add_argument(
        "--save_logo_plots", type = bool,  nargs='?', default=False, help="Whether to save logo plots or not"
    )

    return parser.parse_args()

if __name__ == '__main__':
    assembly = get_input().assembly
    genomepy_dir = get_input().genomepy_dir
    tf_motifs = get_input().tf_motifs
    tf_asb = get_input().tf_asb
    n_cpus = get_input().n_cpus
    out_dir = get_input().out_dir
    sub_dir = get_input().sub_dir
    save_logo_plots = get_input().save_logo_plots

    fasta = f"{genomepy_dir}/{assembly}.fa"
    ReferenceGenome.offset = 25
    genome = ReferenceGenome(fasta) 

    # Read in NoPeak motifs
    print("Reading in NoPeak motifs...")
    nopeak_motifs = mu.read_nopeak_motifs(input_path = ".", tf = tf_motifs)

    # Pull IC and kmer count for each motif and filter out low kmer motifs
    motif_md = mu.pull_motif_metadata(motifs = nopeak_motifs, kmer_count = True)
    motif_md, nopeak_motifs = mu.filter_low_kmer_motifs(md = motif_md, motifs = nopeak_motifs, tf = tf_motifs, 
                                                        out_dir = out_dir, sub_dir = sub_dir)

    # Cluster nopeak_motifs if there are more than one
    clustered_motifs = mu.cluster_high_kmer_motifs(motifs = nopeak_motifs, tf = tf_motifs, ncpus = n_cpus,
        out_dir = out_dir, sub_dir = sub_dir)

    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    peak_data = mu.get_snps(".",tf_asb)

    print("Scoring each SNP against each motif-of-interest...")
    scores_peaks = mu.score_snps(peak_data, clustered_motifs, tf_motifs, genome)

    print("Computing spearman's rho correlation coefficient for each motif...")
    correlations = mu.compute_correlations(scores_peaks)

    print("Filter for motifs with a significant and positive correlation coefficient rho...")
    high_quality_motifs, sig_mappings = mu.find_high_quality_motifs(scores_peaks, correlations, clustered_motifs, tf_motifs)

    print("Computing information content across clustered motifs...")
    ic = mu.pull_motif_metadata(motifs = clustered_motifs)

    print("Reading in SNP data...")
    snp_data = mu.get_snps(".", tf_asb, filter_peak=False)

    print("Scoring all SNPs against high-quality motifs")
    scores = mu.score_snps(snp_data, high_quality_motifs, tf_motifs, genome)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = mu.compute_snp_counts(scores, snp_data, tf_motifs, tf_asb)

    # Get dataframe of all SNPs mapping to all motifs & add column information
    out_snps = mu.compile_scores_data(scores, tf_motifs)

    # Compare against JASPAR motifs
    print("Classifying NoPeak motifs by mapping to known JASPAR motifs...")
    redundant_motifs, accessory_motifs, denovo_motifs = mu.find_motif_matches(motifs = high_quality_motifs, 
                                                                                        tf = tf_motifs,
                                                                                        out_dir = out_dir)

    mu.save_logo_plot(motifs = redundant_motifs, motif_list = high_quality_motifs, tf = tf_motifs, 
                    motif_group = "redundant", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)
    mu.save_logo_plot(motifs = accessory_motifs, motif_list = high_quality_motifs, tf = tf_motifs, 
                    motif_group = "accessory", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)
    mu.save_logo_plot(motifs = denovo_motifs, motif_list = high_quality_motifs, tf = tf_motifs, 
                    motif_group = "denovo", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)

    # Add NoPeak motif information into SNP dataframe
    out_snps = mu.add_nopeak_information(out_snps, correlations, redundant_motifs, accessory_motifs, denovo_motifs)

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

    # Add NoPeak motif information into Motif metadata dataframe
    motif_inf = mu.compile_motif_data(correlations,ic,snp_counts,tf_motifs)
    motif_inf["motif"] = motif_inf.index
    motif_inf = mu.add_nopeak_information(motif_inf, correlations, redundant_motifs, accessory_motifs, denovo_motifs)
    motif_inf = motif_inf.drop(['motif'],axis=1)

    # Write all files to CSV
    print("Saving files to csv...")
    motif_inf.to_csv(f"{out_dir}/NoPeak_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: NoPeak_Motif_metadata_{tf_motifs}.csv")
    outname = "_".join([tf_motifs, "ASBs", "NoPeak", "Motifs", "AR_score_diff.csv"])
    final.to_csv(f"{out_dir}/{outname}", index = False)
    print("Written file:", outname)
