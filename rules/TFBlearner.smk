os.makedirs('data/TFBlearner', exist_ok=True)


NSUBREADS=250000000

def tfblearner_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    combs = combs_list()
    tfs = df["TF"].unique().tolist()
    
    proc_preds = [
       f"output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_prediction.html"
      for tf in tfs
    ]
    
    proc_addSupport_preds = [
         f"output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_addSupport_prediction.html"
        for tf in tfs
    ]

    TFBlearner_outputs=(proc_preds+proc_addSupport_preds)

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
  
def motifs_list_tfblearner():
    ckpt = checkpoints.TFBlearner_motif_placeholder.get()
    df = pd.read_csv(ckpt.output["placeholder_motif_files"], sep="\t")
    motifs = df["motif"].unique().tolist()

    return motifs  
  
  
ALL_TARGET_BUILDERS.append(tfblearner_targets)

rule TFBlearner_derive_chip_labels:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    merged_chip_peaks="data/common/02_merged_chip_peaks/filtered_merged_peaks_signal_{comb}.bed",
    noncons_chip_peaks="data/common/02_merged_chip_peaks/filtered_nonCons_peaks_{comb}.bed"
  output:
    chip_labels="data/TFBlearner/chip_labels/labels_{comb}.tsv",
  params:
    dhs_path="data/annotation/dhs.bed",
    gam_path="data/TFBlearner/gam_rep.rds",
    out_dir="data/TFBlearner/chip_labels",
  benchmark: "benchmarks/TFBlearner/prep_drive_chip_labels_{comb}.tsv"
  log: "logs/TFBlearner/derive_chip_labels_{comb}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript --vanilla src/TFBlearner/TFBlearner_ChIP_rep_prob.R {input.merged_chip_peaks} {params.gam_path} {params.dhs_path} {params.out_dir} {input.noncons_chip_peaks} 
    """

rule TFBlearner_subset_atac_data:
  input:
    all_raw_atac="data/common/04_all_raw_atac.tsv",
    atac_merged_bam="data/common/03_merged_cleaned_atac/merged_{CellularContext}.bam",
  output:
    atac_subset_bed="data/TFBlearner/atac_subsetted/subset_{CellularContext}.bed",
  params:
    out_dir="data/TFBlearner/atac_subsetted",
    dhs_path="data/annotation/dhs_margin.bed",
    n_sub_reads=NSUBREADS
  benchmark: "benchmarks/TFBlearner/prep_subset_atac_data_{CellularContext}.tsv"
  log: "logs/TFBlearner/subset_atac_data_{CellularContext}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
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
      
      bedtools intersect -a {params.out_dir}/filtered_merged_fragments_{wildcards.CellularContext}.bed -b {params.dhs_path} -u > {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed
      rm {params.out_dir}/filtered_merged_fragments_{wildcards.CellularContext}.bed
    
      shuf -n {params.n_sub_reads} {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed > {params.out_dir}/subset_{wildcards.CellularContext}.bed
      rm {params.out_dir}/subset_merged_fragments_{wildcards.CellularContext}.bed
    """
    
checkpoint TFBlearner_motif_placeholder:
   input:
     combinations="data/common/04_train_test_combinations.tsv",
   output:
     placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
   params:
     out_dir="data/TFBlearner/motifs",
     motif_models_path="data/motifs/motif_models.rds",
   benchmark: "benchmarks/TFBlearner/prep_process_motif_placeholders.tsv"
   log:  "logs/TFBlearner/process_motif_placeholders.log"
   threads: 4
   # conda: "../envs/tfblearner.yml"
   script: "../src/TFBlearner/TFBlearner_01_motif_placeholder.R"
   
rule TFBlearner_process_motifs:
  input:
     placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
     placeholder_motif="data/TFBlearner/motifs/{motif}.tsv",
  output:
    processed_motifs="data/TFBlearner/motifs/{motif}.rds",
  params:
    out_dir="data/TFBlearner/motifs",
    motif_models_path="data/motifs/motif_models.rds",
    genome_path="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa",
    dhs_path="data/annotation/dhs.bed",
  benchmark: "benchmarks/TFBlearner/prep_process_motifs_{motif}.tsv"
  log: "logs/TFBlearner/process_motifs_{motif}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  script: "../src/TFBlearner/TFBlearner_01_motif_processing.R"

rule TFBlearner_construct_object:
  input:
    placeholder_motif_files="data/TFBlearner/01_placeholder_motif_files.tsv",
    combinations="data/common/04_train_test_combinations.tsv",
    motifs=lambda wildcards: expand("data/TFBlearner/motifs/{motif}.rds", 
                                    motif=motifs_list_tfblearner()),
    atac_subset_bed=lambda wildcards: expand("data/TFBlearner/atac_subsetted/subset_{CellularContext}.bed", 
                                             CellularContext=contexts_list()),
    chip_labels=lambda wildcards: expand("data/TFBlearner/chip_labels/labels_{comb}.tsv", 
                                         comb=combs_list())
  output:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_construct_object.html", 
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/01_mae.rds",
  params:
    out_dir="data/TFBlearner/objects/{HOLDOUT}/basic_object"
  benchmark: "benchmarks/TFBlearner/prep_construct_object_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/construct_object_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      mkdir -p "output/TFBlearner"
      mkdir -p "{params.out_dir}"
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_02_construct_object.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outDir="{params.out_dir}",
                                                 combinationsPath="{input.combinations}",
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
    out_dir="data/TFBlearner/objects/{HOLDOUT}/basic_object",
    footprint_path="data/TFBlearner/consensus_footprints_and_collapsed_motifs_hg38.bed",
    emb_path="data/TFBlearner/embedding.rds",
  benchmark: "benchmarks/TFBlearner/prep_site_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/site_features_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_03_site_features.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outPath="../{output.mae}",
                                                 outDir="../{params.out_dir}",
                                                 maePath="../../{input.mae}",
                                                 footprintPath="../../{params.footprint_path}",
                                                 embPath="../../{params.emb_path}",
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
    out_dir="data/TFBlearner/objects/{HOLDOUT}/basic_object",
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_context_features_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/context_features_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_04_context_features.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outPath="../{output.mae}",
                                                 outDir="../{params.out_dir}",
                                                 maePath="../../{input.mae}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """
    
rule TFBlearner_feature_matrix: 
  input:
    html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_context_features.html",
    mae="data/TFBlearner/objects/{HOLDOUT}/basic_object/03_mae_context.rds"
  output:
    fm_train_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTrainSe.rds"),
    fm_test_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTestSe.rds"),
    fm_train_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/train_feature_matrix_{tf}.h5"),
    fm_test_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/test_feature_matrix_{tf}.h5"),
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_feature_matrix.html"
  params:
    out_dir="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object",
    profile_dir="data/TFBlearner/precomputed_profiles",
    cofactor_map_path="data/TFBlearner/topInteractors.rds",
    seed=SEED
  benchmark: "benchmarks/TFBlearner/prep_feature_matrix_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/feature_matrix_{tf}_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_05_feature_matrix.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outDir="../../{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 cofactorMapPath="../../{params.cofactor_map_path}",
                                                 maePath="../../{input.mae}",
                                                 computeProfile=TRUE,
                                                 profilesDir="{params.profile_dir}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """
# TODO: switch back from output to params
rule TFBlearner_training:
  input:
    fm_train="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTrainSe.rds",
    fm_train_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/train_feature_matrix_{tf}.h5",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_feature_matrix.html"
  output:
    model_basic="models/TFBlearner/{HOLDOUT}/{tf}_basic_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_training.html"
  params:
    out_dir="models/TFBlearner/{HOLDOUT}",
    seed=SEED
  benchmark: "benchmarks/TFBlearner/train_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/train_{tf}_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_06_training.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outPath="../{output.model_basic}",
                                                 outDir="../{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 fmTrainPath="../../{input.fm_train}",
                                                 seed="{params.seed}",
                                                 threads="{threads}"))'
    """  

rule TFBlearner_prediction:  
  input:
    fm_test_basic="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/{tf}_fmTestSe.rds",
    fm_test_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/basic_object/test_feature_matrix_{tf}.h5",
    model_basic="models/TFBlearner/{HOLDOUT}/{tf}_basic_model.txt",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_training.html"
  output:
    pred_basic="predictions/TFBlearner/{HOLDOUT}/{tf}_pred_basic_model.rds",
    html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_prediction.html"
  params:
    out_dir="predictions/TFBlearner/{HOLDOUT}",
    val_chrs=VAL_CHRS
  benchmark: "benchmarks/TFBlearner/pred_{tf}_{HOLDOUT}.tsv"
  log: "logs/TFBlearner/pred_{tf}_{HOLDOUT}.log"
  threads: 4
  # conda: "../envs/tfblearner.yml"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_07_prediction.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.html}",
                                     params=list(outPath="../../{output.pred_basic}",
                                                 outDir="../{params.out_dir}",
                                                 tfName="{wildcards.tf}",
                                                 fmTestPath="../../{input.fm_test_basic}",
                                                 modelPath="../../{input.model_basic}",
                                                 valChrs="{params.val_chrs}",
                                                 threads="{threads}"))'
    """
    
if "TFBlearner.addSupport" in METHODS:
  
  rule TFBlearner_addSupport_construct_object:
    input:
      combinations="data/common/04_train_test_combinations.tsv",
      mae_basic="data/TFBlearner/objects/{HOLDOUT}/basic_object/01_mae.rds",
    output:
      html="output/TFBlearner/{HOLDOUT}/TFBlearner_02_addSupport_construct_object.html", 
      mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/01_maeAddSupport.rds"
    params:
      out_dir="data/TFBlearner/objects/{HOLDOUT}/additional_support_object"
    benchmark: "benchmarks/TFBlearner/prep_construct_addSupport_object_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_construct_object_{HOLDOUT}.log"
    threads: 4
    # conda: "../envs/tfblearner.yml"
    shell:
      """
        exec 2>>{log} 
        mkdir -p "output/TFBlearner"
        mkdir -p "{params.out_dir}"
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_02_construct_addSupport_object.Rmd",
                                      "html_document", 
                                       output_file="../../{output.html}",
                                       params=list(outDir="{params.out_dir}",
                                                   maePath="{input.mae_basic}",
                                                   combinationsPath="{input.combinations}",
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
      out_dir="data/TFBlearner/objects/{HOLDOUT}/additional_support_object",
      footprint_path="data/TFBlearner/consensus_footprints_and_collapsed_motifs_hg38.bed",
      emb_path="data/TFBlearner/embedding.rds",
    benchmark: "benchmarks/TFBlearner/prep_addSupport_site_features_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_site_features_{HOLDOUT}.log"
    threads: 4
    # conda: "../envs/tfblearner.yml"
    shell:
     """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_03_site_features.Rmd", 
                                      "html_document", 
                                       output_file="../../{output.html}",
                                       params=list(outPath="../{output.mae}",
                                                   outDir="../{params.out_dir}",
                                                   maePath="../../{input.mae}",
                                                   footprintPath="../../{params.footprint_path}",
                                                   embPath="../../{params.emb_path}",
                                                   threads="{threads}"))'
      """
  
  rule TFBlearner_addSupport_context_features: 
    input:
      html="output/TFBlearner/{HOLDOUT}/TFBlearner_03_addSupport_site_features.html",
      mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/02_maeAddSupport_site.rds",
    output:
      html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_addSupport_context_features.html",
      mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/03_maeAddSupport_context.rds",
    params:
      out_dir="data/TFBlearner/objects/{HOLDOUT}/additional_support_object",
      seed=SEED
    benchmark: "benchmarks/TFBlearner/prep_addSupport_context_features_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_context_features_{HOLDOUT}.log"
    threads: 4
  # conda: "../envs/tfblearner.yml"
    shell:
     """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_04_context_features_addSupport.Rmd", 
                                      "html_document", 
                                       output_file="../../{output.html}",
                                       params=list(outPath="../../{output.mae}",
                                                   outDir="../../{params.out_dir}",
                                                   maePath="../../{input.mae}",
                                                   seed="{params.seed}",
                                                   threads="{threads}"))'
      """
  
  rule TFBlearner_addSupport_feature_matrix: 
    input:
      html="output/TFBlearner/{HOLDOUT}/TFBlearner_04_addSupport_context_features.html",
      mae="data/TFBlearner/objects/{HOLDOUT}/additional_support_object/03_maeAddSupport_context.rds"
    output:
      fm_train_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTrainSe.rds"),
      fm_test_basic=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTestSe.rds"),
      fm_train_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/train_feature_matrix_{tf}.h5"),
      fm_test_h5=temp("data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/test_feature_matrix_{tf}.h5"),
      #html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_addSupport_feature_matrix.html"
    params:
      html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_05_addSupport_feature_matrix.html",
      out_dir="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object",
      profile_dir="data/TFBlearner/addSupport/precomputed_profiles",
      cofactor_map_path="data/TFBlearner/topInteractors.rds",
      seed=SEED
    benchmark: "benchmarks/TFBlearner/prep_addSupport_feature_matrix_{tf}_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_feature_matrix_{tf}_{HOLDOUT}.log"
    threads: 4
    # conda: "../envs/tfblearner.yml"
    shell:
      """
        exec 2>>{log} 
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_05_feature_matrix.Rmd", 
                                      "html_document", 
                                       output_file="../../{params.html}",
                                       params=list(outDir="../../{params.out_dir}",
                                                   tfName="{wildcards.tf}",
                                                   cofactorMapPath="../../{params.cofactor_map_path}",
                                                   maePath="../../{input.mae}",
                                                   computeProfile=FALSE,
                                                   profilesDir="{params.profile_dir}",
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
      out_dir="models/TFBlearner/{HOLDOUT}",
      seed=SEED
    benchmark: "benchmarks/TFBlearner/train_addSupport_{tf}_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_train_{tf}_{HOLDOUT}.log"
    threads: 4
    # conda: "../envs/tfblearner.yml"
    shell:
      """
        exec 2>>{log} 
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_06_training.Rmd", 
                                      "html_document", 
                                       output_file="../../{output.html}",
                                       params=list(outPath="../{output.model_addSupport}",
                                                   outDir="../{params.out_dir}",
                                                   tfName="{wildcards.tf}",
                                                   fmTrainPath="../../{input.fm_train}",
                                                   seed="{params.seed}",
                                                   threads="{threads}"))'
      """  
  
  rule TFBlearner_addSupport_prediction:  
    input:
      fm_test_addSupport="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/{tf}_fmTestSe.rds",
      fm_test_addSupport_h5="data/TFBlearner/Feature_matrix/{HOLDOUT}/additional_support_object/test_feature_matrix_{tf}.h5",
      model_addSupport="models/TFBlearner/{HOLDOUT}/{tf}_addSupport_model.txt",
      html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_06_addSupport_training.html"
    output:
      pred_addSupport="predictions/TFBlearner/{HOLDOUT}/{tf}_pred_addSupport_model.rds",
      html="output/TFBlearner/{HOLDOUT}/{tf}_TFBlearner_07_addSupport_prediction.html"
    params:
      out_dir="predictions/TFBlearner/{HOLDOUT}",
      val_chrs=VAL_CHRS
    benchmark: "benchmarks/TFBlearner/pred_addSupport_{tf}_{HOLDOUT}.tsv"
    log: "logs/TFBlearner/addSupport_pred_{tf}_{HOLDOUT}.log"
    threads: 4
    # conda: "../envs/tfblearner.yml"
    shell:
      """
        exec 2>>{log} 
        mkdir -p {params.out_dir}
        Rscript -e 'rmarkdown::render("src/TFBlearner/TFBlearner_07_prediction.Rmd", 
                                      "html_document", 
                                       output_file="../../{output.html}",
                                       params=list(outPath="../../{output.pred_addSupport}",
                                                   outDir="../{params.out_dir}",
                                                   tfName="{wildcards.tf}",
                                                   fmTestPath="../../{input.fm_test_addSupport}",
                                                   modelPath="../../{input.model_addSupport}",
                                                   valChrs="{params.val_chrs}",
                                                   threads="{threads}"))'
      """
