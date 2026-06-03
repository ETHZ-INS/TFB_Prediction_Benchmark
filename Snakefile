#!/usr/bin/env snakemake -s
##
## Benchmarking TF-Binding Prediction tools

#import os.path as op
#import pandas as pd
#import glob as glob
#from functools import partial

configfile: "config.yaml"

envvars:
"ZENODO_TOKEN"

METHODS = config["Training_Prediction"]["methods"]

# TODO: 
# - sif-files for singularity containers ? 
# download embeddings, cofactormap (zenodo), fpdata (viestra)

import os
ANNOTATION_DIR = os.path.abspath(config.get("global", {}).get("annotation_dir", "data/annotation"))

rule all:
    input:
        "Comparison/output/compute_labels.html",
        "Comparison/output/performance_pr.html",
        "Comparison/output/performance_comp_res.html",
        "Comparison/output/performance_pr_plots.html",

rule download_annotation_data:
  input:
  output:
    add_chip_data=f"{ANNOTATION_DIR}/ChIP_mapped.h5",
    dhs_coords=f"{ANNOTATION_DIR}/dhs.bed",
    mae_object=f"{ANNOTATION_DIR}/mae.rds"
  params:
    out_dir=ANNOTATION_DIR,
    record_id="18198234",
    chip_file="ChIP_mapped.h5",
    dhs_file="dhs.bed.gz",
    mae_file="mae.rds"
  shell:
      """
        mkdir -p {params.out_dir}
        # TODO: Just for testing purposes --------------------------------------
        URL_PRIVATE_CHIP="https://zenodo.org/api/records/18198234/draft/files/ChIP_mapped.h5/content"
        URL_PRIVATE_DHS="https://zenodo.org/api/records/18198234/draft/files/PREs.bed.gz/content"
        URL_PRIVATE_MAE="https://zenodo.org/api/records/18198234/draft/files/mae_object.rds/content"
        
        wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir}/{params.chip_file}" "$URL_PRIVATE_CHIP"
        wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir}/{params.dhs_file}" "$URL_PRIVATE_DHS"
        wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir}/{params.mae_file}" "$URL_PRIVATE_MAE"
      
        gunzip "{params.out_dir}/{params.dhs_file}"

        # TODO: Download footprints, profiles (?), cofactor map, and embeddings.rds (last two put to zenodo)
        # blacklisted files / resized chromosomes ?
      """
      
rule create_dhs_margins:
  input:
    dhs_coords=f"{ANNOTATION_DIR}/dhs.bed"
  output:
    dhs_margin_coords=f"{ANNOTATION_DIR}/dhs_margin.bed"
  params:
    margin=201
  script: "src/create_margins.R"
    
module Training_Prediction:
    snakefile: "Training_Prediction/Snakefile"
    config: dict(config["Training_Prediction"], **config["global"], **{"module_prefix": "Training_Prediction"})
    prefix: "Training_Prediction"

module Comparison:
    snakefile: "Comparison/Snakefile"
    config: dict(config["Comparison"], **config["global"], **{"module_prefix": "Comparison"})
    prefix: "Comparison"

use rule * from Training_Prediction
use rule common_prereg_data from Training_Prediction with:
    input:
      add_chip_data=rules.download_annotation_data.output.add_chip_data,
      dhs_coords=rules.download_annotation_data.output.dhs_coords

if "TFBlearner" in METHODS:
  use rule TFBlearner_download_gam from Training_Prediction with:
    input:
      dhs_coords=rules.download_annotation_data.output.dhs_coords
        
  use rule TFBlearner_subset_atac_data from Training_Prediction with:
    input:
      all_raw_atac="Training_Prediction/data/common/04_all_raw_atac.tsv",
      atac_merged_bam="Training_Prediction/data/common/03_merged_cleaned_atac/merged_{CellularContext}.bam",
      dhs_path=rules.create_dhs_margins.output.dhs_margin_coords
  
  use rule TFBlearner_addSupport_download_data from Training_Prediction with:
      input:
        dhs_coords=rules.download_annotation_data.output.dhs_coords,
  
  use rule TFBlearner_addSupport_construct_object from Training_Prediction with:
       input:
        add_motif_data=rules.TFBlearner_addSupport_download_data.output.add_motif_data,
        motif_add_base_dir=rules.TFBlearner_addSupport_download_data.output.add_motif_match_data,
        chip_add_path=rules.download_annotation_data.output.add_chip_data,
        mae_add_path=rules.download_annotation_data.output.mae_object,
        add_meta_path=rules.TFBlearner_addSupport_download_data.output.add_meta_data,
        combinations="Training_Prediction/data/common/04_train_test_combinations.tsv",
        mae_basic=rules.TFBlearner_construct_object.output.mae,
      



use rule * from Comparison
use rule comparison_compute_labels from Comparison with:
    input:
      prediction_list="Training_Prediction/all_prediction_list.txt",
      add_chip_data=rules.download_annotation_data.output.add_chip_data,
      dhs_coords=rules.download_annotation_data.output.dhs_coords,
      mae_object=rules.download_annotation_data.output.mae_object,  
      combinations="Training_Prediction/data/common/04_train_test_combinations.tsv",
      merged_peaks_dir="Training_Prediction/data/common/02_merged_chip_peaks",
      dhs_motifs_report=rules.Common_get_matches_dhs.output.html,
      full_chr_motifs_report=rules.Common_get_matches_genome.output.html,

use rule comparison_performance_metrics from Comparison with:
    input:
      dhs_coords=rules.download_annotation_data.output.dhs_coords,
      combinations="Training_Prediction/data/common/04_train_test_combinations.tsv",
      labels_dhs=rules.comparison_compute_labels.output.labels_dhs,
      labels_full_chr=rules.comparison_compute_labels.output.labels_full_chr,
      mae_object=rules.download_annotation_data.output.mae_object,  
      add_chip_data=rules.download_annotation_data.output.add_chip_data,
      full_chr_bin_ranges=rules.comparison_compute_labels.output.full_chr_bin_ranges,
      full_chr_overlaps=rules.comparison_compute_labels.output.full_chr_overlaps,

use rule comparison_computational_resources from Comparison with:
    input:
      prediction_list="Training_Prediction/all_prediction_list.txt",
      combinations="Training_Prediction/data/common/04_train_test_combinations.tsv",
