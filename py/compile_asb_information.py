from argparse import ArgumentParser
from pathlib import Path
import dask.dataframe as dd
import pandas as pd
import numpy as np
import sys

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--tf", type = str, help = "Gene symbol for TF with ASB events"
    )
    parser.add_argument(
        "--mode", type = str, nargs = '?', const = "allMotifs", help = "Mode to compile under, depending on which motifs were found for the TF-of-interest"
    )

    return parser.parse_args()

def get_asb(input_path, tf):
    data = dd.read_csv(f"{input_path}/*{tf}.withPeaks.csv", include_path_column=True)
    data["cell_line"] = data.apply(lambda row: Path(row.path).name.split("_")[0],  axis=1)
    data["tf"] = data.apply(lambda row: Path(row.path).name.split(".")[0].split("_")[1],  axis=1)
    peak_data = data[data.peak].copy()
    data = data.compute().reset_index(drop=True)
    peak_data = peak_data.compute().reset_index(drop=True)

    return peak_data, data

def get_nopeak(input_path, tf, motif_mode):
    if ((motif_mode == "noNoPeak") | (motif_mode == "noMotifs")):
        colnames = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'REF.counts', 'ALT.counts', 'Total.counts', 
        'AR', 'RMbias', 'RAF', 'Bayes_lower', 'Bayes_upper', 'Bayes_SD', 'conf_0.99_lower', 'conf_0.99_upper', 
        'Corrected.AR', 'isASB', 'peak', 'path', 'cell_line', 'tf', 'ref_seq', 'alt_seq', 'motif', 'motif_pos', 
        'motif_strand', 'ref_score', 'alt_score', 'score_diff', 'tf_motif', 'Significant_rho', 'NoPeak_mapping', 'Concordant']
        data = pd.DataFrame(columns = colnames)
    elif ((motif_mode == "allMotifs") | (motif_mode == "noJaspar")):
        data = pd.read_csv(f"{input_path}/{tf}_ASBs_NoPeak_Motifs_AR_score_diff.csv")
    else:
        sys.stderr("No valid mode given, Please specify either allMotifs, noMotifs, noNoPeak or noJaspar.")
        exit(1)
    
    return data

def get_jaspar(input_path, tf, motif_mode):
    if ((motif_mode == "noJaspar") | (motif_mode == "noMotifs")):
        colnames = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'REF.counts', 'ALT.counts', 'Total.counts', 
        'AR', 'RMbias', 'RAF', 'Bayes_lower', 'Bayes_upper', 'Bayes_SD', 'conf_0.99_lower', 'conf_0.99_upper', 
        'Corrected.AR', 'isASB', 'peak', 'path', 'cell_line', 'tf', 'ref_seq', 'alt_seq', 'motif', 'motif_pos', 
        'motif_strand', 'ref_score', 'alt_score', 'score_diff', 'tf_motif', 'Significant_rho', 'NoPeak_mapping', 'Concordant']
        data = pd.DataFrame(columns = colnames)
    elif ((motif_mode == "allMotifs") | (motif_mode == "noNoPeak")):
        data = pd.read_csv(f"{input_path}/{tf}_ASBs_JASPAR_Motifs_AR_score_diff.csv")
    else:
        sys.stderr("No valid mode given, Please specify either allMotifs, noMotifs, noNoPeak or noJaspar.")
        exit(1)
    
    return data

if __name__ == '__main__':
    tf = get_input().tf
    mode = get_input().mode

    peak_asb_data, asb_data = get_asb(".", tf)
    nopeak_data = get_nopeak(".", tf, mode)
    jaspar_data = get_jaspar(".", tf, mode)

    concordant_jaspar = jaspar_data[(jaspar_data.Concordant) & (jaspar_data.Significant_rho)].copy()
    concordant_nopeak = nopeak_data[(nopeak_data.Concordant) & (nopeak_data.Significant_rho)].copy()

    # Generate maps_to_jaspar and maps_to_nopeak
    if jaspar_data.shape[0] != 0:
        maps_to_jaspar = jaspar_data[jaspar_data.Significant_rho].copy()
    else:
        maps_to_jaspar = jaspar_data.copy()
    if nopeak_data.shape[0] != 0:
        maps_to_nopeak = nopeak_data[nopeak_data.Significant_rho].copy()
    else:
        maps_to_nopeak = nopeak_data.copy()

    is_concordant = [ [0]*10 for i in range(peak_asb_data.shape[0])]
    for index, snp in enumerate(peak_asb_data.ID.values):
        cell_line = peak_asb_data.iloc[index]['cell_line']
        # First determine whether snp is concordant in JASPAR motifs
        if snp in concordant_jaspar.ID.values:
            row = concordant_jaspar.index[
                (concordant_jaspar['ID']==snp) & 
                (concordant_jaspar['cell_line'] == cell_line)
                ]
            sub = concordant_jaspar.loc[row].copy()
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[['ID','CHROM','POS','REF','ALT','cell_line','Significant_rho','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
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
            out = sub[['ID','CHROM','POS','REF','ALT','cell_line','Significant_rho','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = group
            is_concordant[index] = out

        # Now if snp still maps to a significant rho but not concordant, let's record this
        elif snp in np.append(maps_to_jaspar.ID.values, maps_to_nopeak.ID.values):
            rowjasp = maps_to_jaspar.index[
                (maps_to_jaspar['ID'] == snp) &
                (maps_to_jaspar['cell_line'] == cell_line)
            ]
            rownopeak = maps_to_nopeak.index[
                (maps_to_nopeak['ID'] == snp) &
                (maps_to_nopeak['cell_line'] == cell_line)
            ]
            subjasp = maps_to_jaspar.loc[rowjasp].copy()
            subnopeak = maps_to_nopeak.loc[rownopeak].copy()
            # Define motif group
            if (subjasp.shape[0] != 0) & (subnopeak.shape[0] != 0):
                group = "JASPAR and NoPeak"
            elif (subjasp.shape[0] != 0):
                group = "JASPAR"
            else:
                group = "NoPeak"
            sub = pd.concat([subjasp,subnopeak])
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[['ID','CHROM','POS','REF','ALT','cell_line','Significant_rho','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = group
            is_concordant[index] = out

        else:
            sub = peak_asb_data[['ID','CHROM','POS','REF','ALT','cell_line']]
            out = sub.loc[[index]]
            out['Concordant'] = pd.NA
            out['Motifs'] = pd.NA
            out['Motif_group'] = pd.NA
            # Here setting to NA because it could either be Significant_rho = False OR
            # it doesn't map to a motif, so will choose to be more general
            out['Significant_rho'] = pd.NA 
            is_concordant[index] = out

    # Join with original full asb dataframe
    if len(is_concordant)!=0: #ie. if no SNPs found within peaks
        compiled = pd.concat(is_concordant)
        merged = peak_asb_data.merge(compiled, how = 'inner')
        final = asb_data.merge(merged, how = 'left', on = asb_data.columns.tolist())
    else:
        final = asb_data.copy()
        final['Concordant'] = pd.NA
        final['Motifs'] = pd.NA
        final['Motif_group'] = pd.NA
        final['Significant_rho'] = pd.NA

    # Order column names the same
    colnames = asb_data.columns.tolist() + ['Significant_rho','Concordant','Motifs','Motif_group']
    final = final[colnames]
    
    final.to_csv(f"{tf}.withPeaks.withMotifs.csv", index=False)
