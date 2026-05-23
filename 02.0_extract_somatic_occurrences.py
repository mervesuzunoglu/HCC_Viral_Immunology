#!/usr/bin/env python3
"""
TCGA-LIHC Somatic Mutation Occurrence Extraction Pipeline
Author: [Your Name/GitHub Handle]
Description: Programmatically streams case-specific somatic missense mutation 
             instances from the live NCI GDC REST API, flattens the nested 
             genomic coordinate architecture, and exports a clean variant matrix.
"""

import os
import json
import logging
import requests
import pandas as pd

# Configure structured logging for production traceability
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def fetch_gdc_mutation_occurrences(project_id="TCGA-LIHC", size_limit=5000):
    """
    Queries the live GDC ssm_occurrences endpoint to map individual variant instances
    directly to active patient case tracking barcodes.
    """
    endpoint_url = "https://api.gdc.cancer.gov/ssm_occurrences"
    
    # Construct precise payload filtering parameters
    query_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "ssm.consequence.transcript.consequence_type", "value": ["missense_variant"]}}
        ]
    }
    
    payload = {
        "filters": json.dumps(query_filters),
        "expand": "case,ssm",
        "format": "JSON",
        "size": str(size_limit)
    }
    
    logging.info(f"Initiating remote handshake with GDC Occurrence Index for project {project_id}...")
    try:
        response = requests.post(endpoint_url, json=payload, timeout=60)
        response.raise_for_status()
    except requests.exceptions.RequestException as error:
        logging.error(f"API connection failed: {error}")
        raise
        
    raw_hits = response.json().get("data", {}).get("hits", [])
    logging.info(f"Successfully retrieved {len(raw_hits)} raw occurrence footprints.")
    return raw_hits

def process_and_flatten_records(raw_hits):
    """
    Parses complex nested JSON nodes, isolates absolute genomic coordinate changes,
    and maps them against clinical submitter identifiers.
    """
    flattened_records = []
    
    for hit in raw_hits:
        case_node = hit.get("case", {})
        patient_barcode = case_node.get("submitter_id")
        
        ssm_node = hit.get("ssm", {})
        genomic_dna_change = ssm_node.get("genomic_dna_change")
        ssm_core_id = ssm_node.get("ssm_id")
        occurrence_id = hit.get("ssm_occurrence_id")
        
        # Enforce strict field validation before downstream parsing inclusion
        if patient_barcode and genomic_dna_change:
            flattened_records.append({
                "patient_barcode": str(patient_barcode).strip().upper(),
                "Genomic_DNA_Change": str(genomic_dna_change).strip(),
                "Variant_Classification": "Missense_Mutation",
                "Occurrence_ID": str(occurrence_id).strip() if occurrence_id else None,
                "SSM_Core_ID": str(ssm_core_id).strip() if ssm_core_id else None
            })
            
    df = pd.DataFrame(flattened_records)
    if not df.empty:
        # Algorithmic pruning of transcript or redundant cohort entry duplications
        before_prune = df.shape[0]
        df.drop_duplicates(subset=['patient_barcode', 'Genomic_DNA_Change'], inplace=True)
        logging.info(f"Matrix flattening complete. Pruned duplicates down from {before_prune} to {df.shape[0]} unique vectors.")
    else:
        logging.warning("Parsing matrix generated an empty dataset matrix profile.")
        
    return df

def main():
    # Setup adaptive environment checking for cloud volume mounts vs local deployment
    import sys
    if 'google.colab' in sys.modules:
        try:
            from google.colab import drive
            logging.info("Google Colab detected. Mounting cloud storage volume...")
            drive.mount('/content/drive', force_remount=True)
            output_dir = '/content/drive/MyDrive/HCC_Viral_Immunology_Study/Step2_Host_Neoantigens'
        except ImportError:
            output_dir = './data_mutation_output'
    else:
        output_dir = './data_mutation_output'
        
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, "filtered_missense_mutations.csv")
    
    try:
        raw_data = fetch_gdc_mutation_occurrences(project_id="TCGA-LIHC", size_limit=5000)
        flattened_df = process_and_flatten_records(raw_data)
        
        if not flattened_df.empty:
            flattened_df.to_csv(output_file_path, index=False)
            logging.info(f"Somatic mutation registry successfully locked at: {output_file_path}")
        else:
            logging.error("Pipeline termination: No valid records compiled.")
    except Exception as pipeline_error:
        logging.error(f"Execution halted due to runtime exception: {pipeline_error}")

if __name__ == "__main__":
    main()
