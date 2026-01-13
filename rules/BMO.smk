os.makedirs('data/BMO', exist_ok=True)

def bmo_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    
    combsTest = df[df["set"] == "testData" ]["comb"].unique().tolist()
    tfs, contexts = zip(*(s.rsplit("_", 1) for s in combsTest)) if combsTest else ([], [])
    
    pred_combs = [
      f"predictions/BMO/{tf}_{cellularContext}.all.bed"
      for tf, cellularContext in zip(tfs, contexts)
    ]
    
    return pred_combs
  
ALL_TARGET_BUILDERS.append(bmo_targets)  

rule BMO_prepare_motif_txt:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_models="data/motifs/motif_models.rds"
  params:
    out_dir="data/BMO"
  log: "logs/BMO/prepare_motif_{cellularContext}_txt.log"
  benchmark: "benchmarks/BMO/prep_motif_{cellularContext}_txt.tsv"
  threads: 1
  output:
    motif_list="data/BMO/motif_{cellularContext}.txt"
  script: "../src/BMO/BMO_prepare_motif_list.R"
    
rule BMO_prepare_config:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_matching_report="output/FIMO/01_motif_matches_all.html",
    config_yaml="rules/BMO/config/config.yaml",
    motif_list="data/BMO/motif_{cellularContext}.txt",
    atac_bam_dir="data/common/03_merged_cleaned_atac/merged_{cellularContext}.bam",
    atac_peak_dir="data/common/03_peaks_atac/filtered_merged_{cellularContext}_peaks.narrowPeak",
  params:
    motif_dir="data/motifs/motif_matches_genome",
    pred_dir="predictions/BMO",
    out_dir="data/BMO",
    bmo_dir="rules/BMO",
  threads: 4
  log: "logs/BMO/prepare_config_{cellularContext}.log"
  benchmark: "benchmarks/BMO/prep_config_{cellularContext}.tsv"
  output: 
    config_yaml="data/BMO/input_{cellularContext}_config.yaml"
  conda: "../envs/bmo.yml"
  script: "../src/BMO/BMO_prepare_config.R"

rule BMO_run_analysis:
  input:
    config_yaml="data/BMO/input_{cellularContext}_config.yaml",
  output:
    predictions=directory("predictions/BMO/{cellularContext}")
  params:
    snakefile="rules/BMO/Snakefile",
    out_dir="predictions/BMO",
    inner_wd="rules/BMO"
  conda: "../envs/bmo.yml"
  threads: 4
  log: "logs/BMO/run_analysis_{cellularContext}.log"
  benchmark: "benchmarks/BMO/train_pred_{cellularContext}.tsv"
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.cellularContext}
      snakemake -s {params.snakefile} -d {params.inner_wd}  -j {threads} --resources io_limit=2 --forceall --configfile {input.config_yaml} &> {log}
    """
    
rule BMO_move_results:
  input:
    predictions="predictions/BMO/{cellularContext}"
  output:
    prediction="predictions/BMO/{tf}_{cellularContext}.all.bed"
  benchmark: "benchmarks/BMO/move_results_{tf}_{cellularContext}.tsv"
  log: "logs/BMO/move_results_{tf}_{cellularContext}.log"
  shell:
    """
      mv {input.predictions}/bmo/{wildcards.cellularContext}/all/{wildcards.tf}.all.bed {output.prediction} &> {log}
    """
