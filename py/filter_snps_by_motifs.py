import itertools
import glob
from collections import defaultdict
from pathlib import Path
import os
from argparse import ArgumentParser

import dask.dataframe as dd
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from gimmemotifs.motif import Motif, read_motifs
from gimmemotifs.fasta import Fasta
from gimmemotifs.scanner import Scanner
from nopeak import NoPeakMotif
from scores import ScoreSet
from tfomics import ReferenceGenome

import pickle
from scipy.stats import spearmanr

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
        "--reference_offset", nargs='?', const=25, type = str, help="bp window around heterozygous SNP for determining motif instance"
    )
    return parser.parse_args()

def filter_bqtls(df):
    return df[
        df.isASB
        & (
            ((df["Corrected.AR"] > 0.5) & (df.score_diff > 0))
            | ((df["Corrected.AR"] < 0.5) & (df.score_diff < 0))
        )
    ].copy()

if __name__ == '__main__':
    assembly = get_input().assembly
    genomepy_idx = get_input().genomepy_dir
    tf_motifs = get_input().tf_motifs
    tf_asb = get_input().tf_asb
    fasta = f"{genomepy_idx}/{assembly}.fa"

    ReferenceGenome.offset = 25
    genome = ReferenceGenome(fasta) # works to here
    motifs = defaultdict(list)

    # JASPAR motifs
    print("Reading in JASPAR motifs...")
    jaspar_motifs = defaultdict(list)

    for file in Path(f"./").glob("*.jaspar"):
        jaspar_motifs[str(tf_motifs)] += read_motifs(str(file), fmt="jaspar")

    print(jaspar_motifs)

    # Pull IC for each motif and write to dataframe
    # Pull out information content for all motifs
    ic = pd.DataFrame(columns = ['motif','information_content'])
    for factor, motif_list in jaspar_motifs.items():
        ic = pd.DataFrame({'motif':[str(motif) for motif in motif_list], 
                           'information_content': [str(motif.information_content) for motif in motif_list]})

    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    asb_data = dd.read_csv(f"./*{tf_asb}*.csv", include_path_column=True)
    asb_data = asb_data[asb_data.peak]
    asb_data["cell_line"] = asb_data.apply(lambda row: Path(row.path).name.split("_")[0],  axis=1)
    asb_data["tf"] = asb_data.apply(lambda row: Path(row.path).name.split(".")[0].split("_")[1],  axis=1)
    asb_data = asb_data.compute().reset_index(drop=True)

    # Create score set object through gimmemotifs for all facors and motifs within ASB sites for genome hg19
    # This is configured for having 1 factor associated with ASB sites provided and 1 factor for motifs
    print("Creating score set object through gimmemotifs for all factors and motifs in ASB sites...")
    score_sets = { factor : ScoreSet(asb_data, motif_list, genome) for factor, motif_list in jaspar_motifs.items()}

    # Get best motif score with FDR < 0.05
    print("Pulling best motif score with FDR < 0.05 for all motifs in this set...")
    scores = { factor : set.get_best_scores(fpr=0.05) for factor, set in score_sets.items() }

    rsid_counts = {factor : len(set.ID.unique()) for factor, set in asb_data.groupby("tf")}

    print("Computing spearman's rho correlation coefficient for each motif...")
    correlations = defaultdict(dict)

    for factor, score_set in scores.items():
        for motif, group in score_set.groupby("motif"):
            correlation = spearmanr(group["Corrected.AR"],group.score_diff)
            correlations[factor][motif] = [correlation[0],correlation[1]]
            print(f"{factor}, {motif}: {correlation}") 

    # Export dataframe with motif/factor/correlation coefficient mapping
    corrs = pd.DataFrame(columns = ['factor','motif','spearmans_rho','pval'])

    for factor, vals in correlations.items():
        roes = [t[0] for t in vals.values()]
        pvals = [p[1] for p in vals.values()]
        tmp = pd.DataFrame({'factor': factor, 'motif' : vals.keys(), 'spearmans_roe' : roes, 'pval' : pvals})
        corrs = corrs.append(tmp)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = []
    for factor, group in scores.items():
        # SNP counts
        id_counts = pd.DataFrame(group.pivot(columns=("motif"), values=("ID")).apply(lambda series: len(series.dropna().unique()), axis=0), columns=["counts"]).sort_values(by="counts", ascending=False).reset_index()
        id_counts["fraction_of_snps"] = id_counts["counts"]/rsid_counts[tf_asb]
        id_counts["tf_motif"] = factor
        id_counts["tf_asb"] = tf_asb
        snp_counts.append(id_counts)
        snp_counts = pd.concat(snp_counts)

        # Write csv with score_diff and AR information
        out_df = pd.DataFrame(scores[factor])
        out_df["tf_motif"] = factor
        outname = "_".join([factor, "Motifs", tf_asb, "ASBs", "AR_score_diff.csv"])
        out_df.to_csv(outname, index = False)
        print("Written file:", outname)

        # Filter bQTLs and write to file
        total_snps = len(group.ID.unique())
        filtered = filter_bqtls(group)
        filtered_bqtls = len(filtered.ID.unique())
        with open(f"rsids_{tf_asb}_AR_concordant_with_{factor}_motif.txt", 'w') as outfile:
            outfile.writelines(map(lambda line: f"{line}\n", filtered.ID.unique()))
        print(
            f"ASBs factor: {tf_asb}, Motifs factor: {factor}, Significant SNPs: {total_snps}, Filtered bQTLs: {filtered_bqtls}"
        ) 

    # Export all CSV files
    print("Saving files to csv...")
    pd.DataFrame(snp_counts).to_csv(f"nSNPs_{tf_asb}_in_{tf_motifs}_motifs.csv", index = False)
    print(f"Written file: nSNPs_{tf_asb}_in_{tf_motifs}_motifs.csv")
    pd.DataFrame(corrs).to_csv(f"Spearmans_rho_pval_{tf_motifs}_motif_{tf_asb}_asb.csv", index = False)
    print(f"Written file: Spearmans_rho_pval_{tf_motifs}_motif_{tf_asb}_asb.csv")
    ic.to_csv(f"Information_content_{tf_motifs}_motifs.csv", index = False)
    print(f"Written file: Information_content_{tf_motifs}_motifs.csv")
