#!/usr/bin/env python3
"""
TCGA-LIHC Somatic Mutation Flanking Peptide Construction Pipeline
Author: [Merve S Uzunoglu / GitHub mervesuzunoglu]
Description: Programmatically extracts missense variant consequence vectors from 
             the NCI GDC, leverages an optimized local transcript sequence cache 
             via the Ensembl REST API, and constructs 31-mer flanking wild-type 
             vs. mutant peptide windows for downstream HLA binding calculations.
"""

import os
import sys
import json
import time
import re
import logging
import requests
import pandas as pd

# Configure structured logging for pipeline traceability
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_somatic_flanking_peptides(input_csv_path, output_directory_path, batch_size=250):
    """
    Translates case-level somatic mutations into 31-mer flanking peptide windows by
    reconciling pre-computed GDC VEP transcript data with Ensembl reference proteomes.
    """
    os.makedirs(output_directory_path, exist_ok=True)
    destination_file = os.path.join(output_directory_path, "constructed_flanking_peptides.csv")
    
    # Ingest baseline coordinate dataset
    mut_df = pd.read_csv(input_csv_path)
    unique_ssm_ids = mut_df['SSM_Core_ID'].dropna().unique().tolist()
    logging.info(f"Loaded source file. Parsing {len(unique_ssm_ids)} unique somatic mutation cores...")

    gdc_ssms_url = "https://api.gdc.cancer.gov/ssms"
    protein_sequence_cache = {}  # Resolves the N+1 network bottleneck via localized memory lookup
    peptide_records = []

    for i in range(0, len(unique_ssm_ids), batch_size):
        batch = unique_ssm_ids[i:i+batch_size]
        
        # Configure specific fields to bypass verbose unneeded JSON payloads
        payload = {
            "filters": json.dumps({"op": "in", "content": {"field": "ssm_id", "value": batch}}),
            "expand": "consequence.transcript",
            "fields": (
                "ssm_id,"
                "consequence.transcript.aa_change,"
                "consequence.transcript.transcript_id,"
                "consequence.transcript.gene_symbol,"
                "consequence.transcript.consequence_type"
            ),
            "format": "JSON",
            "size": str(batch_size)
        }
        
        try:
            response = requests.post(gdc_ssms_url, json=payload, timeout=30)
            response.raise_for_status()
            hits = response.json().get("data", {}).get("hits", [])
            
            for hit in hits:
                ssm_id = hit.get("ssm_id")
                consequences = hit.get("consequence", [])
                
                for cons in consequences:
                    tx = cons.get("transcript", {})
                    
                    if tx.get("consequence_type") == "missense_variant":
                        aa_change = tx.get("aa_change")
                        tx_id = tx.get("transcript_id")
                        gene_symbol = tx.get("gene_symbol")
                        
                        if aa_change and tx_id:
                            # Standardize GDC variant string and extract components via regex
                            clean_aa = str(aa_change).replace("p.", "").strip()
                            match = re.match(r'([a-zA-Z]+)(\d+)([a-zA-Z]+)', clean_aa)
                            
                            if match:
                                wt_aa = match.group(1).upper()
                                aa_start = int(match.group(2))
                                mut_aa = match.group(3).upper()
                                
                                # Fetch from memory cache or query Ensembl REST API defensively
                                if tx_id not in protein_sequence_cache:
                                    seq_url = f"https://rest.ensembl.org/sequence/id/{tx_id}?type=protein"
                                    seq_res = requests.get(seq_url, headers={"Accept": "application/json"}, timeout=15)
                                    
                                    if seq_res.status_code == 200:
                                        protein_sequence_cache[tx_id] = seq_res.json().get("seq", "")
                                    else:
                                        protein_sequence_cache[tx_id] = ""
                                    time.sleep(0.1)  # API pacing compliance
                                
                                full_protein_sequence = protein_sequence_cache[tx_id]
                                
                                if full_protein_sequence:
                                    pos = aa_start - 1  # Standard 0-based indexing shift
                                    
                                    # Construct symmetric context windows (15aa flanking regions)
                                    start_idx = max(0, pos - 15)
                                    end_idx = min(len(full_protein_sequence), pos + 16)
                                    
                                    wt_flank = full_protein_sequence[start_idx:end_idx]
                                    
                                    # Formulate mutated sequence counterpart
                                    mut_protein_list = list(full_protein_sequence)
                                    if pos < len(mut_protein_list):
                                        mut_protein_list[pos] = mut_aa
                                    mut_flank = "".join(mut_protein_list)[start_idx:end_idx]
                                    
                                    # Retain only structural segments valid for downstream tiling
                                    if len(wt_flank) >= 9:
                                        peptide_records.append({
                                            "SSM_Core_ID": ssm_id,
                                            "Hugo_Symbol": gene_symbol,
                                            "Ensembl_Transcript_ID": tx_id,
                                            "Protein_Mutation_Position": aa_start,
                                            "Wild_Type_Residue": wt_aa,
                                            "Mutant_Residue": mut_aa,
                                            "WT_31mer_Flank": wt_flank,
                                            "MUT_31mer_Flank": mut_flank
                                        })
                                        break  # Move to the next unique variant footprint
                                        
            logging.info(f"Processed block iteration {i} to {min(i+batch_size, len(unique_ssm_ids))}. Unique proteins cached: {len(protein_sequence_cache)}")
        except Exception as batch_error:
            logging.error(f"Batch iteration starting at {i} failed: {batch_error}")
            continue

    peptides_df = pd.DataFrame(peptide_records)
    if not peptides_df.empty:
        # Cross-reference back against the source cohort table to re-link patient barcodes
        final_merged_df = pd.merge(mut_df, peptides_df, on="SSM_Core_ID", how="inner")
        final_merged_df.drop_duplicates(subset=['Genomic_DNA_Change', 'WT_31mer_Flank'], inplace=True)
        final_merged_df.to_csv(destination_file, index=False)
        logging.info(f"Pipeline complete. Formatted 31-mer flanking dataset written to: {destination_file}")
        return final_merged_df
    else:
        logging.error("Failed to compile any valid sequence matrices.")
        return pd.DataFrame()

if __name__ == "__main__":
    # Detect environment configuration automatically
    if 'google.colab' in sys.modules:
        from google.colab import drive
        drive.mount('/content/drive')
        input_csv = '/content/drive/MyDrive/HCC_Viral_Immunology_Study/Step2_Host_Neoantigens/filtered_missense_mutations.csv'
        output_workspace = '/content/drive/MyDrive/HCC_Viral_Immunology_Study/Step3_Peptide_Construction'
    else:
        input_csv = "./data_mutation_output/filtered_missense_mutations.csv"
        output_workspace = "./peptide_construction_output"
        
    run_peptide_extraction = process_somatic_flanking_peptides(input_csv_path=input_csv, output_directory_path=output_workspace)
