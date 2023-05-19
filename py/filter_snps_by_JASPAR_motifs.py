import glob
from collections import defaultdict
from pathlib import Path
from argparse import ArgumentParser
import sys

import dask.dataframe as dd
import pandas as pd

from gimmemotifs.motif import read_motifs
from scores import ScoreSet
from tfomics import ReferenceGenome

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

def read_jaspar_motifs(input_path, tf):
    motifs = defaultdict(list)

    for file in Path(f"{input_path}/").glob("*.jaspar"):
        motifs[str(tf)] += read_motifs(str(file), fmt="jaspar")

    return motifs

def pull_ic(motifs):
    df = pd.DataFrame(columns = ['motif','information_content'])
    for factor, motif_list in motifs.items():
        df = pd.DataFrame({'motif':[str(motif) for motif in motif_list], 
                           'information_content': [str(motif.information_content) for motif in motif_list]})

    return df

def get_asb(input_path, tf):
    data = dd.read_csv(f"{input_path}/*{tf}*.csv", include_path_column=True)
    data = data[data.peak]
    data["cell_line"] = data.apply(lambda row: Path(row.path).name.split("_")[0],  axis=1)
    data["tf"] = data.apply(lambda row: Path(row.path).name.split(".")[0].split("_")[1],  axis=1)
    data = data.compute().reset_index(drop=True)

    return data

def compute_correlations(scores_dict):
    corr_dict = defaultdict(dict)

    for factor, score_set in scores_dict.items():
        for motif, group in score_set.groupby("motif"):
            group = group[group.score_diff != 0]
            correlation = spearmanr(group["Corrected.AR"],group.score_diff)
            corr_dict[factor][motif] = [correlation[0],correlation[1]]
            print(f"{factor}, {motif}: {correlation}") 

    # Export dataframe with motif/factor/correlation coefficient mapping
    corrs = pd.DataFrame(columns = ['factor','motif','spearmans_rho','pval'])

    for factor, vals in corr_dict.items():
        roes = [t[0] for t in vals.values()]
        pvals = [p[1] for p in vals.values()]
        tmp = pd.DataFrame({'factor': factor, 'motif' : vals.keys(), 'spearmans_rho' : roes, 'pval' : pvals})
        corrs = corrs.append(tmp)

    return corrs

def find_significant_motifs(df, correlations, motifs, tf):
    filt = []
    significant_motifs = correlations[correlations.pval <= 0.05].motif.values
    table_filt = df[df['motif'].isin(significant_motifs)] 
    for motif in motifs[tf]:
        if str(motif) in significant_motifs:
            filt.append(motif)

    if table_filt.shape[0]==0:
        sys.stderr.write(f"No significant motifs found for {tf_motifs}.\n")

    return filt, table_filt

def compute_snp_counts(df, asb, motif_tf, asb_tf):
    rsid_counts = {factor : len(set.ID.unique()) for factor, set in asb.groupby("tf")}
    snp_cts = []
    id_counts = pd.DataFrame(df.pivot(columns=("motif"), values=("ID")).apply(lambda series: len(series.dropna().unique()), axis=0), columns=["counts"]).sort_values(by="counts", ascending=False).reset_index()
    id_counts["fraction_of_snps"] = id_counts["counts"]/rsid_counts[asb_tf]
    id_counts["tf_motif"] = motif_tf
    id_counts["tf_asb"] = asb_tf
    snp_cts.append(id_counts)
    snp_cts = pd.concat(snp_cts)

    return snp_cts

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
    genomepy_dir = get_input().genomepy_dir
    tf_motifs = get_input().tf_motifs
    tf_asb = get_input().tf_asb
    fasta = f"{genomepy_dir}/{assembly}.fa"

    ReferenceGenome.offset = 25
    genome = ReferenceGenome(fasta) 

    # JASPAR motifs
    print("Reading in JASPAR motifs...")
    jaspar_motifs = read_jaspar_motifs(".", tf_motifs)

    # Pull IC for each motif and write to dataframe
    ic = pull_ic(jaspar_motifs)

    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    asb_data = get_asb(".",tf_asb)

    # Create score set object through gimmemotifs for all facors and motifs within ASB sites for genome hg19
    # This is configured for having 1 factor associated with ASB sites provided and 1 factor for motifs
    print("Creating score set object through gimmemotifs for all factors and motifs in ASB sites...")
    score_sets = { factor : ScoreSet(asb_data, motif_list, genome) for factor, motif_list in jaspar_motifs.items()}

    # Get best motif score with FDR < 0.05
    print("Pulling best motif score with FDR < 0.05 for all motifs in this set...")
    scores = { factor : set.get_best_scores(fpr=0.05) for factor, set in score_sets.items() }

    correlations = compute_correlations(scores)

    # Filter for motifs with PCC <= 0.05
    out = scores[tf_motifs].copy()
    filtered_motifs, sig = find_significant_motifs(out, correlations, jaspar_motifs, tf_motifs)

    print("Computing the number of SNPs that fall within each motif...")
    snp_counts = compute_snp_counts(out, asb_data, tf_motifs, tf_asb)

    # Write csv with score_diff and AR information
    out["tf_motif"] = tf_motifs
    out["Significant_rho"] = out["motif"].isin(sig["motif"].unique())

    # Filter bQTLs and add information into output file
    total_snps = len(out.ID.unique())
    filtered = filter_bqtls(out)
    filtered_bqtls = len(filtered.ID.unique())
    print(
        f"ASBs factor: {tf_asb}, Motifs factor: {tf_motifs}, Significant SNPs: {total_snps}, Filtered bQTLs: {filtered_bqtls}"
    ) 

    # Add Concordance information into output table
    out["Concordant"] = out["ID"].isin(filtered["ID"].unique())

    outname = "_".join([tf_motifs, "ASBs", "JASPAR", "Motifs", "AR_score_diff.csv"])
    out.to_csv(f"{outname}", index = False)
    print("Written file:", outname)

    # Concatenate all motif dataframe to export as one file
    motif_inf = snp_counts.set_index('motif').join(correlations.set_index('motif')).join(ic.set_index('motif'))
    motif_inf = motif_inf.drop(['factor'],axis=1)

    # Export all CSV files
    print("Saving files to csv...")
    pd.DataFrame(motif_inf).to_csv(f"JASPAR_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: JASPAR_Motif_metadata_{tf_motifs}.csv")
