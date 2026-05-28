#!/usr/bin/env python3
"""
TCGA-LIHC Somatic Mutation Epitope Tiling & Input Chunking Pipeline
Author: [Merve S Uzunoglu / GitHub mervesuzunoglu]
Description: Extracts all overlapping 9-mer, 10-mer (MHC-I), and 15-mer (MHC-II)
             peptide windows containing somatic missense variants from 31-mer flanking 
             sequences.Deduplicates arrays and splits them into clean, newline-separated 
             terminal text files (.pep) for local high-throughput prediction loops.
            
"""

import os
import sys
import logging
import pandas as pd

# Configure structured logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def tile_sliding_window(wt_seq, mut_seq, window_length, target_idx):
    """
    Slides a window of a specified length across the mutated position to capture
    all overlapping fragments containing the altered amino acid.
    """
    fragments = set()
    if pd.isna(wt_seq) or pd.isna(mut_seq) or len(wt_seq) < target_idx + 1:
        return fragments
        
    start_pos = max(0, target_idx - window_length + 1)
    end_pos = min(target_idx, len(wt_seq) - window_length)
    
    for i in range(start_pos, end_pos + 1):
        wt_pep = wt_seq[i : i + window_length]
        mut_pep = mut_seq[i : i + window_length]
        
        # Guard layer to drop truncated edge fragments
        if len(wt_pep) == window_length and len(mut_pep) == window_length:
            fragments.add(wt_pep)
            fragments.add(mut_pep)
    return fragments

def run_epitope_tiling_pipeline(input_csv, output_dir, chunk_size=5000):
    """
    Ingests 31-mer flanking peptide data strings, runs dual-track tiling,
    and partitions outputs into performance-optimized terminal input blocks.
    """
    os.makedirs(output_dir, exist_ok=True)
    c1_dir = os.path.join(output_dir, "Class_I_Chunks")
    c2_dir = os.path.join(output_dir, "Class_II_Chunks")
    os.makedirs(c1_dir, exist_ok=True)
    os.makedirs(c2_dir, exist_ok=True)

    # Ingest Step 3 substrate coordinates
    df = pd.read_csv(input_csv)
    logging.info(f"Loaded {len(df)} flanking records. Beginning overlapping sequence extraction...")

    class_1_peptides = set()
    class_2_peptides = set()
    mutation_index = 15  # Fixed central index of the somatic residue in a 31-mer flank

    for index, row in df.iterrows():
        wt_flank = row['WT_31mer_Flank']
        mut_flank = row['MUT_31mer_Flank']
        
        # Track 1: MHC Class I (9-mer and 10-mer cytotoxic spaces)
        for k in [9, 10]:
            class_1_peptides.update(tile_sliding_window(wt_flank, mut_flank, k, mutation_index))
            
        # Track 2: MHC Class II (15-mer helper spaces)
        class_2_peptides.update(tile_sliding_window(wt_flank, mut_flank, 15, mutation_index))

    logging.info(f"Deduplication complete. Unique MHC-I: {len(class_1_peptides):,} | Unique MHC-II: {len(class_2_peptides):,}")

    def write_pep_files(peptides_set, target_path, prefix):
        peptides_list = list(peptides_set)
        counter = 0
        for i in range(0, len(peptides_list), chunk_size):
            chunk = peptides_list[i : i + chunk_size]
            file_name = os.path.join(target_path, f"{prefix}_chunk_{counter:03d}.pep")
            with open(file_name, 'w') as out_f:
                for peptide in chunk:
                    out_f.write(f"{peptide}\n")
            counter += 1
        logging.info(f"Successfully compiled {counter} raw terminal files (.pep) for {prefix}.")

    # Write newline-delimited strings to local storage directories
    write_pep_files(class_1_peptides, c1_dir, "Class_I")
    write_pep_files(class_2_peptides, c2_dir, "Class_II")
    logging.info("Step 4 Epitope Tiling and Batch Partitioning complete.")

if __name__ == "__main__":
    # Autodetect environmental pipeline execution pathways
    if 'google.colab' in sys.modules:
        from google.colab import drive
        drive.mount('/content/drive')
        input_file = '/content/drive/MyDrive/HCC_Viral_Immunology_Study/Step3_Peptide_Construction/constructed_flanking_peptides.csv'
        output_workspace = '/content/drive/MyDrive/HCC_Viral_Immunology_Study/Step4_Epitope_Tiling'
    else:
        input_file = "./peptide_construction_output/constructed_flanking_peptides.csv"
        output_workspace = "./epitope_tiling_output"
        
    run_epitope_tiling_pipeline(input_csv=input_file, output_dir=output_workspace)
