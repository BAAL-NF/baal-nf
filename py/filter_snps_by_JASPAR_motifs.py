from pathlib import Path
from argparse import ArgumentParser
import os

import pandas as pd
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
        "--reference_offset", nargs='?', default=25, type = str, help="bp window around heterozygous SNP for determining motif instance"
    )
    parser.add_argument(
        "--out_dir", type = str,  nargs='?', default=".", help="Output directory"
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
    asb_data = mu.get_asb(".",tf_asb)

    # If there are no ASBs, just exit without creating a new file
    if (asb_data.shape[0]==0):
        print(f"No ASB sites found for {tf_asb}. Exiting...")
        exit(0)

    # Create score set object through gimmemotifs for all factors and motifs within ASB sites for genome hg19
    print("Creating score set object through gimmemotifs for all factors and motifs in ASB sites...")
    score_sets = { factor : ScoreSet(asb_data, motif_list, genome) for factor, motif_list in jaspar_motifs.items()}

    # Get best motif score with FDR < 0.05
    print("Pulling best motif score with FDR < 0.05 for all motifs in this set...")
    all_scores = { factor : set.get_best_scores(fpr=0.05) for factor, set in score_sets.items() }

    if all_scores[tf_motifs].empty:
        print(f"No sequences scored against motifs for {tf_motifs} found with FDR < 0.05")
        exit(0)
    
    # Remove score_diff = 0
    scores = { factor : set[set.score_diff != 0] for factor, set in all_scores.items()}
    correlations = mu.compute_correlations(scores)

    # Filter for motifs with PCC <= 0.05
    out = scores[tf_motifs].copy()
    filtered_motifs, sig = mu.find_high_quality_motifs(out, correlations, jaspar_motifs, tf_motifs)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = mu.compute_snp_counts(out, asb_data, tf_motifs, tf_asb)
    snp_counts = snp_counts[snp_counts.motif.isin(correlations.motif)]

    # Write csv with score_diff and AR information
    out["tf_motif"] = tf_motifs
    out["High_quality_motif"] = out["motif"].isin(sig["motif"].unique())

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

    # Concatenate all motif dataframe to export as one file
    motif_inf = snp_counts.set_index('motif').join(correlations.set_index('motif')).join(ic.set_index('motif'))
    motif_inf = motif_inf.drop(['factor'],axis=1)
    motif_inf["high_quality"] = (motif_inf.pval <= 0.05) & (motif_inf.spearmans_rho > 0)

    # Save all CSV files
    os.makedirs(out_dir, exist_ok=True)
    outname = "_".join([tf_motifs, "ASBs", "JASPAR", "Motifs", "AR", "score_diff.csv"])
    final.to_csv(f"{out_dir}/{outname}", index = False)
    print("Written file:", outname)
    print("Saving files to csv...")
    pd.DataFrame(motif_inf).to_csv(f"{out_dir}/JASPAR_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: JASPAR_Motif_metadata_{tf_motifs}.csv")
