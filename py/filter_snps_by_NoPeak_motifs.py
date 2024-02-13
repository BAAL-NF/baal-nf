from collections import defaultdict
from pathlib import Path
from argparse import ArgumentParser

import pandas as pd
import numpy as np
import motifutils as mu

from scores import ScoreSet
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
        "--out_dir", type = str,  nargs='?', default=".", help="Output directory"
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
    if (len(nopeak_motifs[tf_motifs]) > 1):
        clustered_motifs = mu.cluster_high_kmer_motifs(motifs = nopeak_motifs, tf = tf_motifs, ncpus = n_cpus,
            out_dir = out_dir, sub_dir = sub_dir)
    elif (len(nopeak_motifs[tf_motifs]) == 1):
        clustered_motifs = nopeak_motifs
    elif (len(nopeak_motifs[tf_motifs]) == 0):
        print(f"No NoPeak motifs found with KMER count > 10 for {tf_motifs}. Exiting...")
        exit(0)
    
    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    asb_data = mu.get_asb(".",tf_asb)

    # If there are no ASBs, just exit without creating a new file
    if (asb_data.shape[0]==0):
        print(f"No ASB sites found for {tf_asb}. Exiting...")
        exit(0)
    
    # Create score set object through gimmemotifs for all facors and motifs within ASB sites for genome hg19
    print("Creating score set object through gimmemotifs for all factors and motifs in ASB sites...")
    score_sets = { factor : ScoreSet(asb_data, motif_list, genome) for factor, motif_list in clustered_motifs.items()}

    # Get best motif score with FDR < 0.05
    print("Pulling best motif score with FDR < 0.05 for all motifs in this set...")
    all_scores = { factor : set.get_best_scores(fpr=0.05) for factor, set in score_sets.items() }

    if all_scores[tf_motifs].empty:
        print(f"No sequences scored against motifs for {tf_motifs} found with FDR < 0.05")
        exit(0)
    
    # Remove score_diff = 0
    scores = { factor : set[set.score_diff != 0] for factor, set in all_scores.items()}

    print("Computing spearman's rho correlation coefficient for each motif...")
    correlations = mu.compute_correlations(scores)

    # Filter for motifs with PCC pval <= 0.05
    print("Filter for motifs with a significant correlation coefficient rho...")
    out = scores[tf_motifs].copy()
    filtered_motifs, sig = mu.find_high_quality_motifs(out, correlations, clustered_motifs, tf_motifs)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = mu.compute_snp_counts(out, asb_data, tf_motifs, tf_asb)
    snp_counts = snp_counts[snp_counts.motif.isin(correlations.motif)]

    # Write csv with score_diff and AR information
    out["tf_motif"] = tf_motifs
    out["High_quality_motif"] = out["motif"].isin(sig["motif"].unique())

    # If no significant motifs, do not match to JASPAR database
    if sig.shape[0]==0:
        print("No High-quality motifs")
        out["NoPeak_mapping"] = np.nan
    
    # Otherwise - match against JASPAR
    else:
        # Compare against JASPAR motifs
        print("Classifying NoPeak motifs by mapping to known JASPAR motifs...")
        redundant_motifs, accessory_motifs, denovo_motifs = mu.find_motif_matches(motifs = filtered_motifs, 
                                                                                            tf = tf_motifs,
                                                                                            out_dir = out_dir)

        mu.save_logo_plot(motifs = redundant_motifs, motif_list = filtered_motifs, tf = tf_motifs, 
                        motif_group = "redundant", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)
        mu.save_logo_plot(motifs = accessory_motifs, motif_list = filtered_motifs, tf = tf_motifs, 
                        motif_group = "accessory", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)
        mu.save_logo_plot(motifs = denovo_motifs, motif_list = filtered_motifs, tf = tf_motifs, 
                        motif_group = "denovo", out_dir = out_dir, sub_dir = sub_dir, save = save_logo_plots)

        # Add NoPeak motif information into dataframe
        out = mu.add_nopeak_information(out, correlations, redundant_motifs, accessory_motifs, denovo_motifs)
    
    # Filter bQTLs and add information into output file
    total_snps = len(out.ID.unique())

    # Identify concordant/discordant bQTLs on high quality motifs
    out_high_quality = mu.filter_for_high_quality_motifs(out) 
    concordant = mu.filter_bqtls(out_high_quality, concordant = True)
    discordant = mu.filter_bqtls(out_high_quality, concordant = False)
    concordant_bqtls = len(concordant.ID.unique())
    print(
        f"ASBs factor: {tf_asb}, Motifs factor: {tf_motifs}, Significant SNPs: {total_snps}, High-quality bQTLs: {concordant_bqtls}"
    ) 

    # Add Concordance information into output table
    concordant["Concordant"] = True
    discordant["Concordant"] = False
    bqtls = pd.concat([concordant, discordant])
    final = pd.merge(out, bqtls, how = "left")
    final = final.convert_dtypes()

    # Regenerate motif metadata and then concatenate all motif dataframe to export as one file
    ic = mu.pull_motif_metadata(motifs = clustered_motifs) 
    motif_inf = snp_counts.set_index('motif').join(correlations.set_index('motif')).join(ic.set_index('motif'))
    motif_inf = motif_inf.drop(['factor'],axis=1)
    motif_inf["high_quality"] = (motif_inf.pval <= 0.05) & (motif_inf.spearmans_rho > 0)

    # Export all CSV files
    print("Saving files to csv...")
    outname = "_".join([tf_motifs, "ASBs", "NoPeak", "Motifs", "AR_score_diff.csv"])
    final.to_csv(f"{out_dir}/{outname}", index = False)
    print("Written file:", outname)
    pd.DataFrame(motif_inf).to_csv(f"{out_dir}/NoPeak_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: NoPeak_Motif_metadata_{tf_motifs}.csv")
