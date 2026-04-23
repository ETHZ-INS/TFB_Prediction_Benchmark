import pandas as pd

os.makedirs('data/maxATAC', exist_ok=True)

EPOCHS = config["maxatac_epochs"]

def maxatac_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    pairs = df[['TF', 'CellularContext']].drop_duplicates()
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
  
    combsTest = df[df["set"] == "testData" ]["comb"].unique().tolist()
    tfs, contexts = zip(*(s.rsplit("_", 1) for s in combsTest)) if combsTest else ([], [])
    
    pred_all_combs = [
      f"predictions/maxATAC/{HOLDOUT}/{EPOCHS}/all/{tf}/{context}_pred.bw"
      for tf, context in zip(tfs, contexts)
    ]
    pred_dhs_combs = [
      f"predictions/maxATAC/{HOLDOUT}/{EPOCHS}/dhs/{tf}/{context}_pred.bw"
      for tf, context in zip(tfs, contexts)
    ]
    max_atac_outputs = (pred_all_combs+pred_dhs_combs)                    
    
    return max_atac_outputs
  
  
ALL_TARGET_BUILDERS.append(maxatac_targets)

def get_samples_df(which="all"):
    ckpt1 = checkpoints.common_get_combinations.get()   # waits for checkpoint
    ckpt2 = checkpoints.common_get_encode_chip_coverage_ids.get()
    
    if which=="atac":
      df=pd.read_csv(ckpt1.output[1], sep="\t")
    elif which=="chip":
      df=pd.read_csv(ckpt2.output[0], sep="\t")
    else:
      df=pd.read_csv(ckpt1.output[0], sep="\t")
      
    return df

def contexts_list():
    return get_samples_df()["CellularContext"].unique().tolist()

def combs_train_list():
    df = get_samples_df()
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    combs = df[df["set"]== "trainData"]["comb"].unique().tolist()

    return combs
  
# Input function for merge rule (per group)
def group_accession_atac_context_inputs(wildcards):
    df = get_samples_df(which="atac")
    accession_subset = df.loc[df.CellularContext == wildcards.CellularContext, "accession"].tolist()

    return expand("data/maxATAC/01_prepared_atac/{accession}_IS_slop20_RP20M.bw", accession=accession_subset)
  
def group_accession_chip_context_inputs(wildcards):
    df = get_samples_df(which="chip")
    df = df[df["set"]== "trainData"]
    
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")

    accession_subset = df.loc[df.comb == wildcards.comb, "Accession"].tolist()

    return expand("data/maxATAC/01_prepared_chip/{accession}.bw", accession=accession_subset)
  
rule maxATAC_download_annotation_data:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
  output:
    downloaded_data="data/maxATAC/download_annot.txt",
  container: "tfbench-maxatac.sif",
  benchmark: "benchmarks/maxATAC/prep_download_annotation_data.tsv"
  log: "logs/maxATAC/download_annotation_data.log"
  params:
    data_dir="data/maxATAC",
  threads: 1
  shell:
    r"""
      export HOME="$PWD"
      mkdir -p "$HOME/opt"
      maxatac data --output "$HOME/opt" &> "{log}"
      touch "{output.downloaded_data}"
    """
    
rule maxATAC_prepare_atac:
  input:
    downloaded_data="data/maxATAC/download_annot.txt",
    raw_bam="data/common/01_raw_atac/{accession}.bam"
  output:
    bw="data/maxATAC/01_prepared_atac/{accession}_IS_slop20_RP20M.bw"
  params:
    out_dir="data/maxATAC/01_prepared_atac",
    blacklist_bed="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bed",
    blacklist_bw="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bw",
    chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    train_chrs=AUT_CHRS
  threads: 4
  benchmark: "benchmarks/maxATAC/prep_prepare_atac_{accession}.tsv"
  log: "logs/maxATAC/prepare_atac_{accession}.log"
  container: "tfbench-maxatac.sif",
  shell:
    """
      mkdir -p {params.out_dir}
      maxatac prepare -i $(pwd)/{input.raw_bam} \
      -o {params.out_dir} \
      --prefix {wildcards.accession} \
      --blacklist {params.blacklist_bed} \
      --blacklist_bw {params.blacklist_bw} \
      --chrom_sizes {params.chrom_size} \
      --threads {threads} \
      --chromosomes {params.train_chrs} &> {log}
    """

rule maxATAC_average_atac:
    input:
      group_accession_atac_context_inputs
    output:
      averaged_bw="data/maxATAC/02_averaged_atac/{CellularContext}.bw"
    params:
      merged_dir="data/maxATAC/02_averaged_atac",
      chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
      train_chrs=AUT_CHRS
    benchmark: "benchmarks/maxATAC/prep_average_atac_{CellularContext}.tsv"
    log: "logs/maxATAC/average_atac_{CellularContext}.log"
    container: "tfbench-maxatac.sif",
    threads: 4
    shell:
      """
        mkdir -p {params.merged_dir}
        maxatac average -i {input} \
        -n {wildcards.CellularContext} \
        -o {params.merged_dir} \
        --chrom_sizes {params.chrom_size} \
        --chromosomes {params.train_chrs} \
        --threads {threads} &> {log}
      """
        
rule maxATAC_normalize_atac:
    input:
      averaged_bw="data/maxATAC/02_averaged_atac/{CellularContext}.bw"
    output:
      norm_bw="data/maxATAC/03_normalized_atac/{CellularContext}_minmax.bw"
    params:
      out_dir="data/maxATAC/03_normalized_atac",
      chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
      blacklist_bw="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bw",
      train_chrs=AUT_CHRS
    benchmark: "benchmarks/maxATAC/prep_normalize_atac_{CellularContext}.tsv"
    log: "logs/maxATAC/normalize_atac_{CellularContext}.log"
    container: "tfbench-maxatac.sif",
    threads: 4
    shell:
      """
        mkdir -p {params.out_dir}
        maxatac normalize -i {input.averaged_bw} \
        --name {wildcards.CellularContext}_minmax \
        -o {params.out_dir} \
        --method min-max \
        --max_percentile 99 \
        --blacklist_bw {params.blacklist_bw} \
        --chrom_sizes {params.chrom_size} \
        --chromosomes {params.train_chrs} \
        --threads {threads} &> {log}
      """

rule maxATAC_convert_chip_coverage_tracks:
  input:
    downloaded_data="data/maxATAC/download_annot.txt",
    bam_file="data/common/01_raw_chip_coverage/{accession}.bam",
  output:
    bw_file="data/maxATAC/01_prepared_chip/{accession}.bw",
  params:
    out_dir="data/maxATAC/01_prepared_chip",
    chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    autosomes_path="data/annotation/autosomes.bed",
  benchmark: "benchmarks/maxATAC/prep_convert_chip_coverage_tracks_{accession}.tsv"
  log: "logs/maxATAC/convert_chip_coverage_tracks_{accession}.log"
  container: "tfbench-maxatac.sif",
  threads: 4
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      bedtools genomecov -bg -ibam {input.bam_file} > {params.out_dir}/{wildcards.accession}.bedGraph
      sort -k1,1 -k2,2n -o {params.out_dir}/{wildcards.accession}.bedGraph {params.out_dir}/{wildcards.accession}.bedGraph
        
      bedGraphToBigWig {params.out_dir}/{wildcards.accession}.bedGraph {params.chrom_size} {output.bw_file}
      rm {params.out_dir}/{wildcards.accession}.bedGraph
    """

rule maxATAC_average_chip_coverage_tracks:
  input:
    group_accession_chip_context_inputs
  output:
    averaged_bw="data/maxATAC/02_averaged_chip/{comb}.bw"
  params:
    out_dir="data/maxATAC/02_averaged_chip",
    chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    train_chrs=AUT_CHRS
  benchmark: "benchmarks/maxATAC/prep_average_chip_coverage_tracks_{comb}.tsv"
  log: "logs/maxATAC/average_chip_coverage_tracks_{comb}.log"
  container: "tfbench-maxatac.sif",
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}
      maxatac average -i {input} \
      -n {wildcards.comb} \
      -o {params.out_dir} \
      --chrom_sizes {params.chrom_size} \
      --chromosomes {params.train_chrs} \
      --threads {threads} &> {log}
    """

rule maxATAC_normalize_chip_coverage_tracks:
  input:
    averaged_bw="data/maxATAC/02_averaged_chip/{comb}.bw"
  output:
    norm_bw="data/maxATAC/03_normalized_chip/{comb}_minmax.bw"
  params:
    out_dir="data/maxATAC/03_normalized_chip",
    chrom_size="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    blacklist_bw="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bw",
    train_chrs=AUT_CHRS
  benchmark: "benchmarks/maxATAC/prep_normalize_chip_coverage_tracks_{comb}.tsv"
  log: "logs/maxATAC/normalize_chip_coverage_tracks_{comb}.log"
  container: "tfbench-maxatac.sif",
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}
      maxatac normalize -i {input.averaged_bw} \
      --name {wildcards.comb}_minmax \
      -o {params.out_dir} \
      --method min-max \
      --max_percentile 99 \
      --blacklist_bw {params.blacklist_bw} \
      --chrom_sizes {params.chrom_size} \
      --chromosomes {params.train_chrs} \
      --threads {threads} &> {log}
    """

rule maxATAC_prepare_metadata:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
  params:
    holdout=HOLDOUT,
    out_dir="data/maxATAC/metadata",
    out_dir_report="output/maxATAC",
  benchmark: "benchmarks/maxATAC/prepare_{HOLDOUT}_metadata.tsv"
  output:
    meta_report="output/maxATAC/01_metaData_{HOLDOUT}_maxATAC.html"
  benchmark: "benchmarks/maxATAC/prep_prepare_metadata_{HOLDOUT}.tsv"
  log: "logs/maxATAC/prepare_metadata_{HOLDOUT}.log"
  threads: 1
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      mkdir -p {params.out_dir_report}
      Rscript -e 'rmarkdown::render("src/utils/maxATAC/01_metaData_maxATAC.Rmd", 
                                      "html_document", 
                                       output_file="../../../{output.meta_report}",
                                       params=list(outDir="{params.out_dir}",
                                                   combDir="{input.combinations}",
                                                   holdOut="{params.holdout}"))'
    """

rule maxATAC_train:
  input:
    meta_report="output/maxATAC/01_metaData_{HOLDOUT}_maxATAC.html",
    atac_bw=lambda wildcards: expand("data/maxATAC/03_normalized_atac/{cellularContext}_minmax.bw", 
                                      cellularContext=contexts_list()),
    chip_bam=lambda wildcards: expand("data/maxATAC/03_normalized_chip/{comb}_minmax.bw", 
                                      comb=combs_train_list())
  output:
    models="models/maxATAC/{HOLDOUT}/{EPOCHS}/all/{tf}/best_epoch.txt",
  container: 
    lambda wc: "tfbench-maxatac_gpu.sif" if int(wc.EPOCHS) >= 100 else "tfbench-maxatac.sif"
  resources:
    gpu=lambda wc, attempt: 1 if config["maxatac_epochs"] >= 100 else 0, #TODO: do also for the image used
  #container: "tfbench-maxatac_gpu.sif",
  params:
    genome="hg38",
    holdout=HOLDOUT,
    model_dir="models/maxATAC/{HOLDOUT}/{EPOCHS}/all",
    meta_dir="data/maxATAC/metadata",
    sequence="data/maxATAC/maxatac/data/hg38/hg38.2bit",
    blacklist="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bed",
    chrom_sizes="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    epochs=EPOCHS,
    train_chrs=AUT_CHRS,
    seed=SEED,
  benchmark: "benchmarks/maxATAC/train_{HOLDOUT}_{tf}_{EPOCHS}.tsv"
  log: "logs/maxATAC/train_{HOLDOUT}_{tf}_{EPOCHS}.log"
  threads: 4,
  shell:
    """
      mkdir -p {params.model_dir}/{wildcards.tf}
      maxatac train --genome {params.genome} \
                    --sequence {params.sequence} \
                    --meta_file {params.meta_dir}/maxATAC_meta_{params.holdout}_{wildcards.tf}.tsv \
                    --output {params.model_dir}/{wildcards.tf} \
                    --prefix {wildcards.tf} \
                    --blacklist {params.blacklist} \
                    --chrom_sizes {params.chrom_sizes} \
                    --chromosomes {params.train_chrs} \
                    --epochs {params.epochs} \
                    --seed {params.seed} \
                    --plot False \
                    --threads {threads} \
                    --multiprocessing True &> {log}
    """
    
rule maxATAC_predict_all:
  input:
    models="models/maxATAC/{HOLDOUT}/{EPOCHS}/all/{tf}/best_epoch.txt",
  output:
    predictions="predictions/maxATAC/{HOLDOUT}/{EPOCHS}/all/{tf}/{context}_pred.bw",
  container: "tfbench-maxatac.sif",
  params:
    chrom_sizes="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    blacklist="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bed",
    sequence="data/maxATAC/maxatac/data/hg38/hg38.2bit",
    out_dir="predictions/maxATAC/{HOLDOUT}/{EPOCHS}/all",
    model_dir="models/maxATAC/{HOLDOUT}/{EPOCHS}/all",
    val_chrs=VAL_CHRS
  resources:
    gpu=lambda wc, attempt: 1 if config["maxatac_epochs"] >= 100 else 0,
  benchmark: "benchmarks/maxATAC/pred_all_{HOLDOUT}_{tf}_{context}_{EPOCHS}.tsv"
  log: "logs/maxATAC/predict_all_{HOLDOUT}_{tf}_{context}_{EPOCHS}.log"
  threads: 4
  shell:
    """
      # Read the best model
      model_summary_path="{params.model_dir}/{wildcards.tf}/best_epoch.txt"
      model_name=$(head -n 1 "$model_summary_path")

      mkdir -p {params.out_dir}
      export CUDA_VISIBLE_DEVICES=""; export OMP_NUM_THREADS={threads}; \
      maxatac predict --model "$model_name" \
                      --signal data/maxATAC/03_normalized_atac/{wildcards.context}_minmax.bw \
                      --name {wildcards.context}_pred \
                      --output {params.out_dir}/{wildcards.tf} \
                      --chromosomes {params.val_chrs} \
                      --sequence {params.sequence} \
                      --blacklist {params.blacklist} \
                      --chrom_sizes {params.chrom_sizes} \
                      --batch_size 3000 \
                      --threads {threads} \
                      --skip_call_peaks &> {log}
    """

rule maxATAC_predict_dhs:
  input:
    models="models/maxATAC/{HOLDOUT}/{EPOCHS}/all/{tf}/best_epoch.txt",
  output:
    predictions="predictions/maxATAC/{HOLDOUT}/{EPOCHS}/dhs/{tf}/{context}_pred.bw",
  container: "tfbench-maxatac.sif",
  benchmark: "benchmarks/maxATAC/pred_dhs_{HOLDOUT}_{tf}_{context}_{EPOCHS}.tsv"
  log: "logs/maxATAC/predict_dhs_{HOLDOUT}_{tf}_{context}_{EPOCHS}.log"
  params:
    dhs_regions="data/annotation/dhs.bed",
    chrom_sizes="data/maxATAC/maxatac/data/hg38/hg38.chrom.sizes",
    blacklist="data/maxATAC/maxatac/data/hg38/hg38_maxatac_blacklist.bed",
    sequence="data/maxATAC/maxatac/data/hg38/hg38.2bit",
    out_dir="predictions/maxATAC/{HOLDOUT}/{EPOCHS}/dhs",
    model_dir="models/maxATAC/{HOLDOUT}/{EPOCHS}/all",
    val_chrs=VAL_CHRS,
  resources:
    gpu=lambda wc, attempt: 1 if config["maxatac_epochs"] >= 100 else 0,
  threads: 4
  shell:
    """
      # Read the best model
      model_summary_path="{params.model_dir}/{wildcards.tf}/best_epoch.txt"
      model_name=$(head -n 1 "$model_summary_path")

      mkdir -p {params.out_dir}
      export CUDA_VISIBLE_DEVICES=""; export OMP_NUM_THREADS={threads}; \
      maxatac predict --model "$model_name" \
                      --signal data/maxATAC/03_normalized_atac/{wildcards.context}_minmax.bw \
                      --name {wildcards.context}_pred \
                      --output {params.out_dir}/{wildcards.tf} \
                      --bed {params.dhs_regions} \
                      --chromosomes {params.val_chrs} \
                      --sequence {params.sequence} \
                      --blacklist {params.blacklist} \
                      --chrom_sizes {params.chrom_sizes} \
                      --batch_size 3000 \
                      --threads {threads} \
                      --skip_call_peaks &> {log}
    """
