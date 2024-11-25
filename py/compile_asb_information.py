from argparse import ArgumentParser
from pathlib import Path
import dask.dataframe as dd
import pandas as pd
import numpy as np
import sys
import os

def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "--tf", type = str, help = "Gene symbol for TF with ASB events"
    )
    parser.add_argument(
        "--mode", type = str, nargs = '?', default = "allMotifs", help = "Mode to compile under, depending on which motifs were found for the TF-of-interest"
    )
    parser.add_argument(
        "--out_dir", type = str, nargs = '?', default = ".", help = "Output directory"
    ) 

    return parser.parse_args()

def get_asb(input_path, tf):
    data = dd.read_csv(f"{input_path}/*{tf}.withPeaks.csv", include_path_column=True)
    data["cell_line"] = data.apply(lambda row: Path(row.path).name.split("_")[0],  axis=1)
    data["tf"] = data.apply(lambda row: Path(row.path).name.split(".")[0].split("_")[1],  axis=1)
    data = data.drop(['path'], axis = 1)
    peak_data = data[data.peak].copy()
    data = data.compute().reset_index(drop=True)
    peak_data = peak_data.compute().reset_index(drop=True)

    return peak_data, data

def get_nopeak(input_path, tf, motif_mode):
    if ((motif_mode == "noNoPeak") | (motif_mode == "noMotifs")):
        colnames = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'REF.counts', 'ALT.counts', 'Total.counts', 
        'AR', 'RMbias', 'RAF', 'Bayes_lower', 'Bayes_upper', 'Bayes_SD', 'conf_0.99_lower', 'conf_0.99_upper', 
        'Corrected.AR', 'isASB', 'peak', 'path', 'cell_line', 'tf', 'ref_seq', 'alt_seq', 'motif', 'motif_pos', 
        'motif_strand', 'ref_score', 'alt_score', 'score_diff', 'tf_motif', 'High_quality_motif', 'NoPeak_mapping', 'Concordant']
        data = pd.DataFrame(columns = colnames)
    elif ((motif_mode == "allMotifs") | (motif_mode == "noJaspar")):
        data = pd.read_csv(f"{input_path}/{tf}_ASBs_NoPeak_Motifs_AR_score_diff.csv")
    else:
        sys.stderr("No valid mode given, Please specify either allMotifs, noMotifs, noNoPeak or noJaspar.")
        exit(1)
    
    data = data.convert_dtypes()
    data = data[data.tf == tf]
    return data

def get_jaspar(input_path, tf, motif_mode):
    if ((motif_mode == "noJaspar") | (motif_mode == "noMotifs")):
        colnames = ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'REF.counts', 'ALT.counts', 'Total.counts', 
        'AR', 'RMbias', 'RAF', 'Bayes_lower', 'Bayes_upper', 'Bayes_SD', 'conf_0.99_lower', 'conf_0.99_upper', 
        'Corrected.AR', 'isASB', 'peak', 'path', 'cell_line', 'tf', 'ref_seq', 'alt_seq', 'motif', 'motif_pos', 
        'motif_strand', 'ref_score', 'alt_score', 'score_diff', 'tf_motif', 'High_quality_motif', 'NoPeak_mapping', 'Concordant']
        data = pd.DataFrame(columns = colnames)
    elif ((motif_mode == "allMotifs") | (motif_mode == "noNoPeak")):
        data = pd.read_csv(f"{input_path}/{tf}_ASBs_JASPAR_Motifs_AR_score_diff.csv")
    else:
        sys.stderr("No valid mode given, Please specify either allMotifs, noMotifs, noNoPeak or noJaspar.")
        exit(1)
    
    data = data.convert_dtypes()
    data = data[data.tf == tf]
    return data

def subset_bqtls(jaspar, nopeak, tf, concordant = True, isasb = True):
    sub_jaspar = jaspar[
        (jaspar.High_quality_motif) & 
        (jaspar.isASB == isasb) &
        (jaspar.tf == tf)
        ].copy()
    sub_nopeak = nopeak[
        (nopeak.High_quality_motif) & 
        (nopeak.isASB == isasb) &
        (nopeak.tf == tf) &
        (nopeak.NoPeak_mapping != 'Matches expected JASPAR motif')
        ].copy()

    if isasb:
        sub_jaspar = sub_jaspar[sub_jaspar.Concordant == concordant]
        sub_nopeak = sub_nopeak[sub_nopeak.Concordant == concordant]

    return sub_jaspar, sub_nopeak 


if __name__ == '__main__':
    tf = get_input().tf
    mode = get_input().mode
    out_dir = get_input().out_dir

    peak_asb_data, asb_data = get_asb(".", tf)
    nopeak_data = get_nopeak(".", tf, mode)
    jaspar_data = get_jaspar(".", tf, mode)

    # Filter for concordant & discordant bQTLs
    concordant_jaspar, concordant_nopeak = subset_bqtls(jaspar_data, nopeak_data, tf, concordant = True, isasb = True)
    discordant_jaspar, discordant_nopeak = subset_bqtls(jaspar_data, nopeak_data, tf, concordant = False, isasb = True)

    # Identify SNPs that map to high-quality motifs but are not ASB sites
    snps_jaspar, snps_nopeak = subset_bqtls(jaspar_data, nopeak_data, tf, isasb = False)

    is_concordant = [ [0]*25 for i in range(peak_asb_data.shape[0]) ]
    for index, snp in enumerate(peak_asb_data.ID.values):
        cell_line = peak_asb_data.iloc[index]['cell_line']
        # ASB is concordant in JASPAR motifs
        if snp in concordant_jaspar[concordant_jaspar.cell_line == cell_line].ID.values:
            row = concordant_jaspar.index[
                (concordant_jaspar['ID']==snp) & 
                (concordant_jaspar['cell_line'] == cell_line)
                ]
            sub = concordant_jaspar.loc[row].copy()
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[peak_asb_data.columns.tolist() + ['High_quality_motif','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = 'JASPAR'
            if out.shape[0] > 1:
                print(index)
                raise ValueError("More than one transcription factor present in this table")
            is_concordant[index] = out

        # ASB is concordant in NoPeak accessory or de novo motif
        elif snp in concordant_nopeak[concordant_nopeak.cell_line == cell_line].ID.values:
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
                    (concordant_nopeak['NoPeak_mapping'] == "de novo NoPeak motif")
                ]
                group = "NoPeak_denovo" 
            sub = concordant_nopeak.loc[row].copy()
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[peak_asb_data.columns.tolist() + ['High_quality_motif','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = group
            if out.shape[0] > 1:
                print(index)
                raise ValueError("More than one transcription factor present in this table")
            is_concordant[index] = out

        # ASB maps to high quality motif but is discordant 
        elif snp in np.append(discordant_jaspar[discordant_jaspar.cell_line == cell_line].ID.values, 
                                discordant_nopeak[discordant_nopeak.cell_line == cell_line].ID.values):
            
            rowjasp = discordant_jaspar.index[
                (discordant_jaspar['ID'] == snp) &
                (discordant_jaspar['cell_line'] == cell_line)
            ]
            rownopeak = discordant_nopeak.index[
                (discordant_nopeak['ID'] == snp) &
                (discordant_nopeak['cell_line'] == cell_line)
            ]
            subjasp = discordant_jaspar.loc[rowjasp].copy()
            subnopeak = discordant_nopeak.loc[rownopeak].copy()
            # Define motif group
            if (subjasp.shape[0] != 0) & (subnopeak.shape[0] != 0):
                group = "JASPAR and NoPeak"
            elif (subjasp.shape[0] != 0):
                group = "JASPAR"
            else:
                group = "NoPeak"
            sub = pd.concat([subjasp,subnopeak])
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[peak_asb_data.columns.tolist() + ['High_quality_motif','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = group
            if out.shape[0] > 1:
                print(index)
                raise ValueError("More than one transcription factor present in this table")
            is_concordant[index] = out

        # SNP maps to high quality motif but is not ASB site
        elif snp in np.append(snps_jaspar[snps_jaspar.cell_line == cell_line].ID.values, 
                                snps_nopeak[snps_nopeak.cell_line == cell_line].ID.values): 
            rowjasp = snps_jaspar.index[
                (snps_jaspar['ID'] == snp) &
                (snps_jaspar['cell_line'] == cell_line)
            ]
            rownopeak = snps_nopeak.index[
                (snps_nopeak['ID'] == snp) &
                (snps_nopeak['cell_line'] == cell_line)
            ]
            subjasp = snps_jaspar.loc[rowjasp].copy()
            subnopeak = snps_nopeak.loc[rownopeak].copy()
            # Define motif group
            if (subjasp.shape[0] != 0) & (subnopeak.shape[0] != 0):
                group = "JASPAR and NoPeak"
            elif (subjasp.shape[0] != 0):
                group = "JASPAR"
            else:
                group = "NoPeak"
            sub = pd.concat([subjasp,subnopeak])
            motifs = ";".join(np.unique(sub.motif.values))
            out = sub[peak_asb_data.columns.tolist() + ['High_quality_motif','Concordant']].drop_duplicates()
            out['Motifs'] = motifs
            out['Motif_group'] = group
            is_concordant[index] = out 

        # SNP/ASB does not map to high quality motif
        else:
            sub = peak_asb_data.copy()
            out = sub.loc[[index]]
            # Here we set Significant_rho to NA because it could either be Significant_rho = False OR
            # it doesn't map to a motif, so we choose to be more general
            out['High_quality_motif'] = pd.NA
            out['Concordant'] = pd.NA
            out['Motifs'] = pd.NA
            out['Motif_group'] = pd.NA
            is_concordant[index] = out

    # Join with original full asb dataframe
    if len(is_concordant)!=0: 
        compiled = pd.concat(is_concordant, ignore_index = True, sort = False)
        final = asb_data.merge(compiled, how = 'left', on = asb_data.columns.tolist())
        final = final.convert_dtypes() # Convert NAs to be consistent
    else: #ie. if no SNPs found within peaks
        final = asb_data.copy()
        final['High_quality_motif'] = pd.NA
        final['Concordant'] = pd.NA
        final['Motifs'] = pd.NA
        final['Motif_group'] = pd.NA

    # Order column names the same
    colnames = asb_data.columns.tolist() + ['High_quality_motif','Concordant','Motifs','Motif_group']
    final = final[colnames]
    
    os.makedirs(out_dir, exist_ok = True)
    final.to_csv(f"{out_dir}/{tf}.withPeaks.withMotifs.csv", index=False)
