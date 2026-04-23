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

# TODO: 
# - sif-files for singularity containers ? 
# make modules dependend on Training_Prediction & Comparison
# create dhs margin

rule download_annotation_data:
  input:
  output:
    add_chip_data="data/annotation/ChIP_mapped.h5",
    dhs_coords="data/annotation/dhs.bed"
  params:
    out_dir="data/annotation",
    record_id="18198234",
    chip_file="ChIP_mapped.h5",
    dhs_file="dhs.bed.gz"
  shell:
      """
        mkdir -p {params.out_dir}
        # TODO: Just for testing purposes --------------------------------------
        URL_PRIVATE_CHIP="https://zenodo.org/api/records/18198234/draft/files/ChIP_mapped.h5/content"
        URL_PRIVATE_DHS="https://zenodo.org/api/records/18198234/draft/files/PREs.bed.gz/content"
        
        wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir}/{params.chip_file}" "$URL_PRIVATE_CHIP"
        wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir}/{params.dhs_file}" "$URL_PRIVATE_DHS"
        gunzip "{params.out_dir}/{params.dhs_file}"
      """
      
rule create_dhs_margins:
  input:
    dhs_coords="data/annotation/dhs.bed"
  output:
    dhs_margin_coords="data/annotation/dhs_margin.bed"
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
      
use rule TFBlearner_download_gam from Training_Prediction with:
  input:
    dhs_coords=rules.download_annotation_data.output.dhs_coords
      
use rule TFBlearner_subset_atac_data from Training_Prediction with:
  input:
    all_raw_atac="Training_Prediction/data/common/04_all_raw_atac.tsv",
    atac_merged_bam="Training_Prediction/data/common/03_merged_cleaned_atac/merged_{CellularContext}.bam",
    dhs_path=rules.create_dhs_margins.output.dhs_margin_coords
    
# use rule TFBlearner_addSupport_download_data from Training_Prediction with:
#     input:
#       add_chip_data=rules.download_annotation_data.output.add_chip_data,
#       dhs_coords=rules.download_annotation_data.output.dhs_coords

use rule * from Comparison
use rule compute_metrics from Comparison with:
    input:
        prediction_list="Training_Prediction/all_prediction_list.txt"
        
rule all:
     default_target: True
     input:
       "data/annotation/ChIP_mapped.h5",
       "data/annotation/dhs.bed",
       "Training_Prediction/all_prediction_list.txt",
       "Comparison/output/performance_pr.html"
