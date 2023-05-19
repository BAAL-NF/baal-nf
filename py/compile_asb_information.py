from argparse import ArgumentParser
from pathlib import Path
import dask.dataframe as dd
import pandas as pd
import numpy as np

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--tf", type = str, help = "Gene symbol for TF with ASB events"
    )

    return parser.parse_args()

def get_asb(input_path, tf):
    data = dd.read_csv(f"{input_path}/*{tf}.withPeaks.csv", include_path_column=True)
    data = data[data.peak]
    data["cell_line"] = data.apply(lambda row: Path(row.path).name.split("_")[0],  axis=1)
    data["tf"] = data.apply(lambda row: Path(row.path).name.split(".")[0].split("_")[1],  axis=1)
    data = data.compute().reset_index(drop=True)

    return data

def get_nopeak(input_path, tf):
    data = pd.read_csv(f"{input_path}/{tf}_ASBs_NoPeak_Motifs_AR_score_diff.csv")
    return data

def get_jaspar(input_path, tf):
    data = pd.read_csv(f"{input_path}/{tf}_ASBs_JASPAR_Motifs_AR_score_diff.csv")
    return data

if __name__ == '__main__':
    tf = get_input().tf

    asb_data = get_asb(".", tf)
    nopeak_data = get_nopeak(".", tf)
    jaspar_data = get_jaspar(".", tf)

    concordant_jaspar = jaspar_data[(jaspar_data.Concordant) & (jaspar_data.Significant_rho)].copy()
    concordant_nopeak = nopeak_data[(nopeak_data.Concordant) & (nopeak_data.Significant_rho)].copy()
    is_concordant = [ [0]*9 for i in range(asb_data.shape[0])]
    for index, snp in enumerate(asb_data.ID.values):
        cell_line = asb_data.iloc[index]['cell_line']
        if snp in concordant_jaspar.ID.values:
            row = concordant_jaspar.index[
                (concordant_jaspar['ID']==snp) & 
                (concordant_jaspar['cell_line'] == cell_line)
                ]
            sub = concordant_jaspar.loc[row].copy()
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[['ID','CHROM','POS','REF','ALT','cell_line','Concordant']].drop_duplicates()
            out['motifs'] = motifs
            out['Motif_group'] = 'JASPAR'
            is_concordant[index] = out

        elif snp in concordant_nopeak.ID.values:
            row = concordant_nopeak.index[
                (concordant_nopeak['ID'] == snp) &
                (concordant_nopeak['cell_line'] == cell_line) &
                (concordant_nopeak['NoPeak_mapping'] == "Matches expected JASPAR motif")
                ]
            group = "NoPeak_JASPAR_match"
            if len(row)==0:
                row = concordant_nopeak.index[
                    (concordant_nopeak['ID'] == snp) &
                    (concordant_nopeak['cell_line'] == cell_line) &
                    (concordant_nopeak['NoPeak_mapping'] == "New accessory motif")
                    ]
                group = "NoPeak_accessory_motif"
            if len(row)==0:
                row = concordant_nopeak.index[
                    (concordant_nopeak['ID'] == snp) &
                    (concordant_nopeak['cell_line'] == cell_line) &
                    (concordant_nopeak['NoPeak_mapping'] == "No match found in JASPAR database")
                    ]
                group = "NoPeak_denovo" 
            sub = concordant_nopeak.loc[row].copy()
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[['ID','CHROM','POS','REF','ALT','cell_line','Concordant']].drop_duplicates()
            out['motifs'] = motifs
            out['Motif_group'] = group
            is_concordant[index] = out
        else:
            sub = asb_data[['ID','CHROM','POS','REF','ALT','cell_line']]
            out = sub.loc[[index]]
            out['Concordant'] = False
            out['motifs'] = "N/A"
            out['Motif_group'] = "N/A"
            is_concordant[index] = out

    compiled = pd.concat(is_concordant)
    final = asb_data.merge(compiled, how = 'inner')

    final.to_csv(f"{tf}.withPeaks.withMotifs.csv", index=False)

