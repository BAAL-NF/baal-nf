from collections import defaultdict
from pathlib import Path
import sys
from argparse import ArgumentParser

import dask.dataframe as dd
import pandas as pd
import numpy as np

from gimmemotifs.motif import read_motifs
from gimmemotifs.comparison import MotifComparer
from nopeak import NoPeakMotif
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
        "--k", type = str, help="length of kmer for which NoPeak was run"
    )

    return parser.parse_args()

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

def read_motifs_nopeak(input_path):
    motifs = defaultdict(list)
    for file in Path(f"{input_path}/").glob("*_dedup.bam.motifs.txt"):
        motifs[str(tf_motifs)] += NoPeakMotif.from_file(
                str(file)
            )

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

def find_motif_matches(motifs, tf):
    mc = MotifComparer()
    compare_motifs = mc.get_all_scores(motifs=motifs, dbmotifs=read_motifs("JASPAR2022_vertebrates"), match="partial", metric="seqcor", combine="mean")

    a=[]
    for i in compare_motifs.keys():
        b = pd.DataFrame.from_dict(compare_motifs[i]).T
        b["query"] = i
        a.append(b)

    summary = pd.concat(a,axis=0)
    summary.columns = ['score','position','strand','NoPeak_motif']
    summary = summary.assign(TF = [tf for idx1, idx2, idx3, tf in summary.index.str.split('.')])
    summary_human = summary[summary['TF'].str.isupper()]

    # Now identify when you get a match to the expected TF
    matches_expected = np.unique(
        summary_human[
            (summary_human.TF==tf) & 
            (summary_human.score >= 0.6)
            ].NoPeak_motif.values
    )

    matches_new = np.unique(
        summary_human[
            (summary_human.score >= 0.7) & 
            np.logical_not(summary_human.NoPeak_motif.isin(matches_expected))
            ].NoPeak_motif.values
    )

    return matches_expected, matches_new

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
    k = get_input().k

    fasta = f"{genomepy_dir}/{assembly}.fa"
    ReferenceGenome.offset = 25
    genome = ReferenceGenome(fasta) # works to here

    # read in NoPeak motifs
    print("Reading in NoPeak motifs...")
    nopeak_motifs = read_motifs_nopeak(".")

    # Pull IC for each motif and write to dataframe
    ic = pull_ic(nopeak_motifs)
    
    # Import asb data and filter for ASBs within peaks
    print("Reading in ASB data...")
    asb_data = get_asb(".",tf_asb)

    if asb_data.shape[0] == 0:
        print(f"No bQTLs found in peaks for {tf_asb}")
        motif_cols = ['ref_seq','alt_seq','motif','motif_pos','motif_strand','ref_score','alt_score','score_diff','tf_motif','Significant_rho','NoPeak_mapping','Concordant']
        colnames = asb_data.columns.values.tolist()
        colnames.extend(motif_cols)
        out = pd.DataFrame(columns = colnames)

        motif_inf =  pd.DataFrame(columns = ['motif','counts','fraction_of_snps','tf_motif','tf_asb','spearmans_rho','pval','information_content'])

    else:
        # Create score set object through gimmemotifs for all facors and motifs within ASB sites for genome hg19
        print("Creating score set object through gimmemotifs for all factors and motifs in ASB sites...")
        score_sets = { factor : ScoreSet(asb_data, motif_list, genome) for factor, motif_list in nopeak_motifs.items()}

        # Get best motif score with FDR < 0.05
        print("Pulling best motif score with FDR < 0.05 for all motifs in this set...")
        scores = { factor : set.get_best_scores(fpr=0.05) for factor, set in score_sets.items() }

        if scores[tf_motifs].empty:
            print(f"No sequences scored against motifs for {tf_motifs} found with FDR < 0.05")
            exit(0)

        print("Computing spearman's rho correlation coefficient for each motif...")
        correlations = compute_correlations(scores)

        # Filter for motifs with PCC pval <= 0.05
        print("Filter for motifs with a significant correlation coefficient rho...")
        out = scores[tf_motifs].copy()
        filtered_motifs, sig = find_significant_motifs(out, correlations, nopeak_motifs, tf_motifs)

        print("Computing the number of SNPs that fall within each motif...")
        snp_counts = compute_snp_counts(out, asb_data, tf_motifs, tf_asb)

        # Write csv with score_diff and AR information
        out["tf_motif"] = tf_motifs
        out["Significant_rho"] = out["motif"].isin(sig["motif"].unique())

        # If no significant motifs, do not match to JASPAR database
        if (correlations["pval"].isnull().all()) | (not (correlations["pval"] <= 0.05).any()):
            print("No significant motifs")
            out["NoPeak_mapping"] = pd.NA

        # Otherwise - match against JASPAR
        else:
            # Compare against JASPAR motifs
            matches_canonical_tf, accessory_motifs = find_motif_matches(filtered_motifs, tf_motifs)

            # Add NoPeak motif information into dataframe
            nopeak_map = []
            for i in out.motif.values:
                # If motif not significant, label as NA
                is_NS = (
                    pd.isnull(correlations[correlations.motif == i].pval.values) | 
                    (correlations[correlations.motif == i].pval.values > 0.05)
                    )[0]
                query = '_'.join(i.split('_')[0:3]) # remove trailing motif characters
                if is_NS:
                    nopeak_map.append(pd.NA)
                elif query in matches_canonical_tf:
                    nopeak_map.append(f"Matches expected JASPAR motif")
                elif query in accessory_motifs:
                    nopeak_map.append("New accessory motif")
                else:
                    nopeak_map.append("No match found in JASPAR database")
            out["NoPeak_mapping"] = nopeak_map
        
        # Filter bQTLs and add information into output file
        total_snps = len(out.ID.unique())
        filtered = filter_bqtls(out)
        filtered_bqtls = len(filtered.ID.unique())
        print(
            f"ASBs factor: {tf_asb}, Motifs factor: {tf_motifs}, Significant SNPs: {total_snps}, Filtered bQTLs: {filtered_bqtls}"
        ) 

        # Add Concordance information into output table
        out["Concordant"] = out["ID"].isin(filtered["ID"].unique())

        # Concatenate all motif dataframe to export as one file
        motif_inf = snp_counts.set_index('motif').join(correlations.set_index('motif')).join(ic.set_index('motif'))
        motif_inf = motif_inf.drop(['factor'],axis=1)

    # Export all CSV files
    print("Saving files to csv...")
    outname = "_".join([tf_motifs, "ASBs", "NoPeak", "Motifs", "AR_score_diff.csv"])
    out.to_csv(f"{outname}", index = False)
    print("Written file:", outname)
    pd.DataFrame(motif_inf).to_csv(f"NoPeak_Motif_metadata_{tf_motifs}.csv", index = True)
    print(f"Written file: NoPeak_Motif_metadata_{tf_motifs}.csv")
