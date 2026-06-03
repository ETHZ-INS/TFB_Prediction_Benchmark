import os
ANNOTATION_DIR = os.path.abspath(config.get("global", {}).get("annotation_dir", "data/annotation"))

os.makedirs('Training_Prediction/data/TFBlearner', exist_ok=True)

envvars:
"ZENODO_TOKEN"

NSUBREADS=250000000

def tfblearner_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    combs = combs_list()
    tfs = df["TF"].unique().tolist()
    
    proc_preds = [
       f"output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_prediction.html"
      for tf in tfs
    ]
    TFBlearner_outputs=proc_preds

    if "TFBlearner.addSupport" in METHODS:
      proc_addSupport_preds = [
           f"output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_addSupport_prediction.html"
          for tf in tfs
      ]
      TFBlearner_outputs=(TFBlearner_outputs+proc_addSupport_preds)

    return TFBlearner_outputs
  
def read_samples_df():
    ckpt = checkpoints.common_get_combinations.get()
    return pd.read_csv(ckpt.output[0], sep="\t")  
  
def contexts_list():
    return read_samples_df()["CellularContext"].unique().tolist()
  
def combs_list():
    df = read_samples_df()
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    combs = df["comb"].unique().tolist()

    return combs
  
def motifs_list_tfblearner(wildcards):
    ckpt = checkpoints.TFBlearner_motif_placeholder.get()
    df = pd.read_csv(ckpt.output["placeholder_motif_files"], sep="\t")
    motifs = df["motif"].unique().tolist()
    
    return motifs
  
  
ALL_TARGET_BUILDERS.append(tfblearner_targets)

#TODO: download other common annotation data
rule TFBlearner_download_gam:
  input:
    dhs_coords=f"{ANNOTATION_DIR}/dhs.bed",
  output:
    gam="data/TFBlearner/gam_rep.rds",
    embeddings=f"{ANNOTATION_DIR}/embeddings.rds",
    cofactor_map=f"{ANNOTATION_DIR}/topInteractors.rds",
    footprints=f"{ANNOTATION_DIR}/consensus_footprints_and_collapsed_motifs_hg38.bed",
  params:
    out_dir=prefix_path("data/TFBlearner")
  benchmark: "benchmarks/TFBlearner/download_gam.tsv"
  log: "logs/TFBlearner/download_gam.log"
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec > {log} 2>&1
      mkdir -p {params.out_dir}
      URL_PRIVATE_GAM="https://zenodo.org/api/records/18198234/draft/files/reproducibilityGAM.rds/content"
      wget --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.gam}" "$URL_PRIVATE_GAM"

      URL_PRIVATE_EMBD="https://zenodo.org/api/records/18198234/draft/files/sequenceEmbeddings.rds/content"
      wget --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.embeddings}" "$URL_PRIVATE_EMBD"

      URL_PRIVATE_CFCTR="https://zenodo.org/api/records/18198234/draft/files/topInteractors.rds/content"
      wget --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.cofactor_map}" "$URL_PRIVATE_CFCTR"
      
      URL_FP_DATA="https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_collapsed_motifs_hg38.bed.gz"
      wget --output-document="{output.footprints}" "$URL_FP_DATA"    
    """
  
rule TFBlearner_derive_chip_labels:
  input:
    gam="data/TFBlearner/gam_rep.rds",
    combinations="data/common/04_train_test_combinations.tsv",
    merged_chip_peaks="data/common/02_merged_chip_peaks/filtered_merged_peaks_signal_{comb}.bed",
    noncons_chip_peaks="data/common/02_merged_chip_peaks/filtered_nonCons_peaks_{comb}.bed",
    dhs_coords=f"{ANNOTATION_DIR}/dhs.bed",
  output:
    chip_labels="data/TFBlearner/chip_labels/labels_{comb}.tsv",
  benchmark: "benchmarks/TFBlearner/prep_drive_chip_labels_{comb}.tsv"
  log: "logs/TFBlearner/derive_chip_labels_{comb}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  script: "../src/TFBlearner/TFBlearner_ChIP_rep_prob.R"

rule TFBlearner_subset_atac_data:
  input:
    all_raw_atac="data/common/04_all_raw_atac.tsv",
    atac_merged_bam="data/common/03_merged_cleaned_atac/merged_{CellularContext}.bam",
    dhs_path=f"{ANNOTATION_DIR}/dhs_margin.bed",
  output:
    atac_subset_bed="data/TFBlearner/atac_subsetted/subset_{CellularContext}.bed",
  params:
    out_dir=prefix_path("data/TFBlearner/atac_subsetted"),
    n_sub_reads=NSUBREADS
  benchmark: "benchmarks/TFBlearner/prep_subset_atac_data_{CellularContext}.tsv"
  log: "logs/TFBlearner/subset_atac_data_{CellularContext}.log"
  threads: 4
  container: "tfbench-preproc.sif"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      samtools sort -n {input.atac_merged_bam} -o {params.out_dir}/sorted_merged_{wildcards.CellularContext}.bam
      bedtools bamtobed -i {params.out_dir}/sorted_merged_{wildcards.CellularContext}.bam -bedpe > {params.out_dir}/merged_{wildcards.CellularContext}.bedpe
      rm {params.out_dir}/sorted_merged_{wildcards.CellularContext}.bam
      awk '$1 != "." && $9 != "."' {params.out_dir}/merged_{wildcards.CellularContext}.bedpe > {params.out_dir}/merged2_{wildcards.CellularContext}.bedpe
      rm {params.out_dir}/merged_{wildcards.CellularContext}.bedpe
      awk '{{print $1, $2, $6}}'  OFS="\t" {params.out_dir}/merged2_{wildcards.CellularContext}.bedpe | sort -k1,1 -k2,2n > {params.out_dir}/merged_fragments_{wildcards.CellularContext}.bed
      rm {params.out_dir}/merged2_{wildcards.CellularContext}.bedpe
      
      # filter for min-width of 10
      awk '($3 - $2) >= 10 {{ print $1 "\t" $2 "\t" $3 }}' {params.out_dir}/merged_fragments_{wildcards.CellularContext}.bed > {params.out_dir}/filtered_merged_fragments_{wildcards.CellularContext}.bed
      rm {params.out_dir}/merged_fragments_{wildcards.CellularContext}.bed
      
      bedtools intersect -a {params.out_dir}/filtered_merged_fragments_{wildcards.CellularContext}.bed -b {input.dhs_path} -u > {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed
      rm {params.out_dir}/filtered_merged_fragments_{wildcards.CellularContext}.bed
    
      shuf -n {params.n_sub_reads} {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed > {params.out_dir}/subset_{wildcards.CellularContext}.bed
      rm {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed
      #touch {output.atac_subset_bed}
    """
    
checkpoint TFBlearner_motif_placeholder:
   input:
     combinations="data/common/04_train_test_combinations.tsv",
   output:
     placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
   params:
     out_dir=prefix_path("data/TFBlearner/motifs"),
     motif_models_path=f"{ANNOTATION_DIR}/motif_models.rds",
   benchmark: "benchmarks/TFBlearner/prep_process_motif_placeholders.tsv"
   log:  "logs/TFBlearner/process_motif_placeholders.log"
   threads: 4
   container: "tfbench-tfblearner.sif"
   script: "../src/TFBlearner/TFBlearner_01_motif_placeholder.R"
   
ruleorder: TFBlearner_motif_placeholder > TFBlearner_process_motifs > TFBlearner_subset_atac_data

rule TFBlearner_process_motifs:
  input:
     placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
  output:
    processed_motifs="data/TFBlearner/motifs/{motif}.rds",
  params:
    out_dir=prefix_path("data/TFBlearner/motifs"),
    motif_models_path=f"{ANNOTATION_DIR}/motif_models.rds",
    genome_path=f"{ANNOTATION_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa",
    dhs_path=f"{ANNOTATION_DIR}/dhs.bed",
  benchmark: "benchmarks/TFBlearner/prep_process_motifs_{motif}.tsv"
  log: "logs/TFBlearner/process_motifs_{motif}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  script: "../src/TFBlearner/TFBlearner_01_motif_processing.R"

rule TFBlearner_construct_object:
  input:
    placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
    combinations="data/common/04_train_test_combinations.tsv",
    motifs=lambda wildcards: expand("data/TFBlearner/motifs/{motif}.rds", 
                                    motif=motifs_list_tfblearner(wildcards)),
    atac_subset_bed=lambda wildcards: expand("data/TFBlearner/atac_subsetted/subset_{CellularContext}.bed", 
                                             CellularContext=contexts_list()),
    chip_labels=lambda wildcards: expand("data/TFBlearner/chip_labels/labels_{comb}.tsv", 
                                         comb=combs_list())
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_construct_object.html", 
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/01_mae.rds",
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_02_construct_object.Rmd"),
    out_dir_report=prefix_path("output/TFBlearner/{HOLDOUT}"),
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/basic_object"),
    chip_dir=prefix_path("data/TFBlearner/chip_labels"),
    atac_dir=prefix_path("data/TFBlearner/atac_subsetted"),
    motifs_dir=prefix_path("data/TFBlearner/motifs")
  benchmark: "benchmarks/TFBlearner/prep_construct_object_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/construct_object_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      mkdir -p "{params.out_dir_report}"
      mkdir -p "{params.out_dir}"
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outDir="{params.out_dir}",
                                                 combinationsPath="{input.combinations}",
                                                 chipDir="{params.chip_dir}",
                                                 atacDir="{params.atac_dir}",
                                                 motifsDir="{params.motifs_dir}",
                                                 dhsPath="{ANNOTATION_DIR}/dhs.bed",
                                                 threads="{threads}"))'
    """
    
rule TFBlearner_site_features:
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_construct_object.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/01_mae.rds",
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_03_site_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/02_mae_site.rds",
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_03_site_features.Rmd"),
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/basic_object"),
    footprint_path=f"{ANNOTATION_DIR}/consensus_footprints_and_collapsed_motifs_hg38.bed",
    emb_path=f"{ANNOTATION_DIR}/embedding.rds",
  benchmark: "benchmarks/TFBlearner/prep_site_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/site_features_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outPath="{output.mae}",
                                                 outDir="{params.out_dir}",
                                                 maePath="{input.mae}",
                                                 footprintPath="{params.footprint_path}",
                                                 embPath="{params.emb_path}",
                                                 threads="{threads}"))'
    """
    
rule TFBlearner_context_features: 
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_03_site_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/02_mae_site.rds",
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_context_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/03_mae_context.rds",
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_04_context_features.Rmd"),
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/basic_object"),
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_context_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/context_features_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outPath="{output.mae}",
                                                 outDir="{params.out_dir}",
                                                 maePath="{input.mae}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """
    
rule TFBlearner_feature_matrix: 
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_context_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/03_mae_context.rds",
    combinations="data/common/04_train_test_combinations.tsv",
  output:
    fm_train_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTrainSe.rds"),
    fm_test_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTestSe.rds"),
    fm_train_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/train_feature_matrix_{tf}.h5"),
    #fm_test_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/test_feature_matrix_{tf}.h5"),
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_feature_matrix.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_05_feature_matrix.Rmd"),
    out_dir=prefix_path("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object"),
    profile_dir="",
    cofactor_map_path=f"{ANNOTATION_DIR}/topInteractors.rds",
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_feature_matrix_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/feature_matrix_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outDir="{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 combinationsPath="{input.combinations}",
                                                 cofactorMapPath="{params.cofactor_map_path}",
                                                 maePath="{input.mae}",
                                                 computeProfile=TRUE,
                                                 profilesDir="{params.profile_dir}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """

rule TFBlearner_training:
  input:
    fm_train="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTrainSe.rds",
    fm_train_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/train_feature_matrix_{tf}.h5",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_feature_matrix.html"
  output:
    model_basic="models/TFBlearner/{HOLDOUT}/{tf}_basic_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_training.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_06_training.Rmd"),
    out_dir=prefix_path("models/TFBlearner/{HOLDOUT}"),
    seed=SEED
  benchmark: "benchmarks/TFBlearner/train_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/train_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outPath="{output.model_basic}",
                                                 outDir="{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 fmTrainPath="{input.fm_train}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """  

rule TFBlearner_prediction:  
  input:
    fm_test_basic="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTestSe.rds",
    #fm_test_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/test_feature_matrix_{tf}.h5",
    model_basic="models/TFBlearner/{HOLDOUT}/{tf}_basic_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_training.html"
  output:
    pred_basic="predictions/TFBlearner/{HOLDOUT}/{tf}_pred_basic_model.rds",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_prediction.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_07_prediction.Rmd"),
    out_dir=prefix_path("predictions/TFBlearner/{HOLDOUT}"),
    val_chrs=VAL_CHRS
  benchmark: "benchmarks/TFBlearner/pred_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/pred_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript -e 'rmarkdown::render("{params.script}", 
                                    "html_document",
                                     knit_root_dir=getwd(), 
                                     output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                     params=list(outPath="{output.pred_basic}",
                                                 outDir="{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 fmTestPath="{input.fm_test_basic}",
                                                 modelPath="{input.model_basic}",
                                                 valChrs="{params.val_chrs}",
                                                 threads="{threads}"))'
    """

rule TFBlearner_addSupport_download_data:
  input:
    dhs_coords=f"{ANNOTATION_DIR}/dhs.bed",
  output:
    profiles=directory("data/TFBlearner/addSupport/precomputed_profiles"),
    add_motif_data="data/TFBlearner/addSupport/Motif_mapped.h5",
    add_motif_match_data=directory("data/TFBlearner/addSupport/motif_matches"),
    add_meta_data="data/TFBlearner/addSupport/meta.tsv",
    add_assoc_data="data/TFBlearner/addSupport/Activity.Association_mapped.h5",
  params:
    out_dir_parent=prefix_path("data/TFBlearner/addSupport"),
    out_dir_name_profiles="precomputed_profiles",
    out_dir_names_motif_maches="motif_matches"
  benchmark: "benchmarks/TFBlearner/download_profiles.tsv"
  log: "logs/TFBlearner/download_addSupport_data.log"
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec > {log} 2>&1
      mkdir -p {params.out_dir_parent}
      mkdir -p {params.out_dir_parent}/{params.out_dir_name_profiles}
      mkdir -p {params.out_dir_parent}/{params.out_dir_names_motif_maches}

      URL_PRIVATE_PROFILES="https://zenodo.org/api/records/18198234/draft/files/insertions.tar/content"
      wget --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir_parent}/insertions.tar" "$URL_PRIVATE_PROFILES"
      tar -xf "{params.out_dir_parent}/insertions.tar" -C "{params.out_dir_parent}/{params.out_dir_name_profiles}" --strip-components=1
      rm "{params.out_dir_parent}/insertions.tar"
      
      URL_PRIVATE_ASSOC="https://zenodo.org/api/records/18198234/draft/files/Activity.Association_mapped.h5/content"
      wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.add_assoc_data}" "$URL_PRIVATE_ASSOC"

      URL_PRIVATE_MOTIFS="https://zenodo.org/api/records/18198234/draft/files/Motif_mapped.h5/content"
      wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.add_motif_data}" "$URL_PRIVATE_MOTIFS"
        
      URL_PRIVATE_MOTIFS_Matches="https://zenodo.org/api/records/18198234/draft/files/motifs.tar/content"
      wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{params.out_dir_parent}/motifs.tar" "$URL_PRIVATE_MOTIFS_Matches"
      tar -xf "{params.out_dir_parent}/motifs.tar" -C "{params.out_dir_parent}/{params.out_dir_names_motif_maches}" --strip-components=1
      rm "{params.out_dir_parent}/motifs.tar"

      URL_PRIVATE_META="https://zenodo.org/api/records/18198234/draft/files/Supplementary_table2_aggregated.tsv/content"
      wget --quiet --header="Authorization: Bearer $ZENODO_TOKEN" --output-document="{output.add_meta_data}" "$URL_PRIVATE_META"
    """

  '''
  # https://zenodo.org/uploads/18198234
  #rule TFBlearner_addSupport_download_data:
    input:
    output:
      add_chip_data="data/TFBlearner/addSupport/ChIP_mapped.h5",
      add_motif_data="data/TFBlearner/addSupport/Motif_mapped.h5"
    params:
      out_dir="data/TFBlearner/addSupport",
      record_id="18198234",
      chip_file="ChIP_mapped.h5",
      motifs_file="Motif_mapped.h5",
    container: "tfbench-tfblearner.sif"
    shell:
      """
        ZENODO_BASE_URL="https://zenodo.org"
        mkdir -p {params.out_dir}
        
        #URL_PUBLIC_CHIP="${ZENODO_BASE_URL}"/records/{params.record_id}/files/{params.chip_file}?download=1"
        #URL_PUBLIC_MOTIFS="${ZENODO_BASE_URL}"/records/{params.record_id}/files/{params.motifs_file}?download=1"
        #URL_IPS
        #wget --quiet --output-document="${params.out_dir}/{params.chip_file}" "${URL_PUBLIC_CHIP}"; then
        #wget --quiet --output-document="${params.out_dir}/{params.motifs_file}" "${URL_PUBLIC_MOTIFS}"
        
        # PRIVATE_DOWNLOAD_URL=$(wget --quiet --header="Authorization: Bearer ${ZENODO_TOKEN}" -O - "${API_URL}" | jq -r ".files[] | select(.filename==\"${FILENAME}\") | .links.download")
        
        # TODO: Just for testing purposes --------------------------------------
        #URL_PRIVATE_CHIP="${ZENODO_BASE_URL}"/api/deposit/depositions/{params.record_id}/files/{params.chip_file}?download=1"
        #URL_PRIVATE_MOTIFS="${ZENODO_BASE_URL}"/api/deposit/depositions/{params.record_id}/files/{params.motifs_file}?download=1"
        
        URL_PRIVATE_CHIP="https://zenodo.org/api/records/18198234/draft/files/ChIP_mapped.h5/content"
        URL_PRIVATE_MOTIFS="https://zenodo.org/api/records/18198234/draft/files/Motif_mapped.h5/content"
        
        wget --quiet --header="Authorization: Bearer ${ZENODO_TOKEN}" --output-document="${params.out_dir}/{params.chip_file}" "${URL_PRIVATE_CHIP}"
        wget --quiet --header="Authorization: Bearer ${ZENODO_TOKEN}" --output-document="${params.out_dir}/{params.motifs_file}" "${URL_PUBLIC_MOTIFS}"
      """
  '''

rule TFBlearner_addSupport_construct_object:
  input:
    add_motif_data="data/TFBlearner/addSupport/Motif_mapped.h5",
    motif_add_base_dir="data/TFBlearner/addSupport/motif_matches",
    combinations="data/common/04_train_test_combinations.tsv",
    add_meta_path="data/TFBlearner/addSupport/meta.tsv",
    chip_add_path=f"{ANNOTATION_DIR}/ChIP_mapped.h5",
    mae_add_path=f"{ANNOTATION_DIR}/mae.rds",
    mae_basic="data/TFBlearner/objects/{HOLDOUT}/basic_object/01_mae.rds",
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_addSupport_construct_object.html", 
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/01_maeAddSupport.rds"
  params: 
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/additional_support_object"),
    script=prefix_path("src/TFBlearner/TFBlearner_02_construct_addSupport_object.Rmd"),
  benchmark: "benchmarks/TFBlearner/prep_construct_addSupport_object_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_construct_object_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
      exec 2>>{log} 
      mkdir -p "{params.out_dir}"
      Rscript -e 'rmarkdown::render("{params.script}",
                                    "html_document",
                                    knit_root_dir=getwd(), 
                                    output_file=basename("{output.html}"),
                                    output_dir=dirname("{output.html}"),
                                    params=list(maePath="{input.mae_basic}",
                                                combinationsPath="{input.combinations}",
                                                chipAddPath="{input.chip_add_path}",
                                                maeAddSuportPath="{input.mae_add_path}",
                                                metaAddPath="{input.add_meta_path}",
                                                motifAddPath="{input.add_motif_data}",
                                                motifAddBaseDir="{input.motif_add_base_dir}",
                                                outDir="{params.out_dir}",
                                                threads="{threads}"))'
      """ 

rule TFBlearner_addSupport_site_features:
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_addSupport_construct_object.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/01_maeAddSupport.rds"
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_03_addSupport_site_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/02_maeAddSupport_site.rds",
  params:
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/additional_support_object"),
    script=prefix_path("src/TFBlearner/TFBlearner_03_site_features.Rmd"),
    footprint_path=f"{ANNOTATION_DIR}/consensus_footprints_and_collapsed_motifs_hg38.bed",
    emb_path=f"{ANNOTATION_DIR}/embedding.rds",
  benchmark: "benchmarks/TFBlearner/prep_addSupport_site_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_site_features_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("{params.script}", 
                                      "html_document",
                                      knit_root_dir=getwd(), 
                                      output_file=basename("{output.html}"),
                                      output_dir=dirname("{output.html}"),
                                      params=list(outPath="{output.mae}",
                                                  outDir="{params.out_dir}",
                                                  maePath="{input.mae}",
                                                  footprintPath="{params.footprint_path}",
                                                  embPath="{params.emb_path}",
                                                  threads="{threads}"))'
    """
  
rule TFBlearner_addSupport_context_features: 
  input:
    add_assoc_path="data/TFBlearner/addSupport/Activity.Association_mapped.h5",
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_03_addSupport_site_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/02_maeAddSupport_site.rds",
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_addSupport_context_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/03_maeAddSupport_context.rds",
  params:
    mae_add_path=f"{ANNOTATION_DIR}/mae.rds",
    out_dir=prefix_path("data/TFBlearner/objects/{HOLDOUT}/additional_support_object"),
    script=prefix_path("src/TFBlearner/TFBlearner_04_context_features_addSupport.Rmd"),
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_addSupport_context_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_context_features_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("{params.script}", 
                                      "html_document",
                                      knit_root_dir=getwd(), 
                                      output_file=basename("{output.html}"),
                                      output_dir=dirname("{output.html}"),
                                      params=list(outPath="{output.mae}",
                                                  outDir="{params.out_dir}",
                                                  maePath="{input.mae}",
                                                  addAssocPath="{input.add_assoc_path}",
                                                  maeAddSuportPath="{params.mae_add_path}",
                                                  seed="{params.seed}",
                                                  threads="{threads}"))'
    """
  
rule TFBlearner_addSupport_feature_matrix: 
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_addSupport_context_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/03_maeAddSupport_context.rds",
    combinations="data/common/04_train_test_combinations.tsv",
    profile_dir="data/TFBlearner/addSupport/precomputed_profiles",
  output:
    fm_train_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTrainSe.rds"),
    fm_test_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTestSe.rds"),
    fm_train_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/train_feature_matrix_{tf}.h5"),
    #fm_test_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/test_feature_matrix_{tf}.h5"),
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_addSupport_feature_matrix.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_05_feature_matrix.Rmd"),
    out_dir=prefix_path("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object"),
    cofactor_map_path=f"{ANNOTATION_DIR}/topInteractors.rds",
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_addSupport_feature_matrix_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_feature_matrix_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("{params.script}", 
                                      "html_document",
                                      knit_root_dir=getwd(), 
                                      output_file=basename("{output.html}"),
                                      output_dir=dirname("{output.html}"),
                                      params=list(outDir="{params.out_dir}",
                                                  tfName="{wildcards.tf}",
                                                  cofactorMapPath="{params.cofactor_map_path}",
                                                  combinationsPath="{input.combinations}",
                                                  maePath="{input.mae}",
                                                  computeProfile=FALSE,
                                                  profilesDir="{input.profile_dir}",
                                                  seed="{params.seed}",
                                                  threads="{threads}"))'
    """
    
rule TFBlearner_addSupport_training:
  input:
    fm_train="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTrainSe.rds",
    fm_train_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/train_feature_matrix_{tf}.h5",
    #html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_addSupport_feature_matrix.html"
  output:
    model_addSupport="models/TFBlearner/{HOLDOUT}/{tf}_addSupport_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_addSupport_training.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_06_training.Rmd"),
    out_dir=prefix_path("models/TFBlearner/{HOLDOUT}"),
    seed=SEED
  benchmark: "benchmarks/TFBlearner/train_addSupport_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_train_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
        exec 2>>{log} 
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("{params.script}", 
                                      "html_document",
                                     knit_root_dir=getwd(), 
                                       output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                       params=list(outPath="{output.model_addSupport}",
                                                   outDir="{params.out_dir}",
                                                   tfName="{wildcards.tf}",
                                                   fmTrainPath="{input.fm_train}",
                                                   seed="{params.seed}",
                                                   threads="{threads}"))'
    """  
  
rule TFBlearner_addSupport_prediction:  
  input:
    fm_test_addSupport="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTestSe.rds",
    #fm_test_addSupport_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/test_feature_matrix_{tf}.h5",
    model_addSupport="models/TFBlearner/{HOLDOUT}/{tf}_addSupport_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_addSupport_training.html"
  output:
    pred_addSupport="predictions/TFBlearner/{HOLDOUT}/{tf}_pred_addSupport_model.rds",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_addSupport_prediction.html"
  params:
    script=prefix_path("src/TFBlearner/TFBlearner_07_prediction.Rmd"),
    out_dir=prefix_path("predictions/TFBlearner/{HOLDOUT}"),
    val_chrs=VAL_CHRS
  benchmark: "benchmarks/TFBlearner/pred_addSupport_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/addSupport_pred_{tf}_{HOLDOUT}.log"
  threads: 4
  container: "tfbench-tfblearner.sif"
  shell:
    """
        exec 2>>{log} 
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("{params.script}", 
                                      "html_document",
                                     knit_root_dir=getwd(), 
                                       output_file=basename("{output.html}"),
                                     output_dir=dirname("{output.html}"),
                                       params=list(outPath="{output.pred_addSupport}",
                                                   outDir="{params.out_dir}",
                                                   tfName="{wildcards.tf}",
                                                   fmTestPath="{input.fm_test_addSupport}",
                                                   modelPath="{input.model_addSupport}",
                                                   valChrs="{params.val_chrs}",
                                                   threads="{threads}"))'
    """
