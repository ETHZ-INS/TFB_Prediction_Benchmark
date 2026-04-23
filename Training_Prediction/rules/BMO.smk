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

'''
rule BMO_download:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
  output:
    config_yaml="src/BMO/config/config.yaml"
  params:
    out_dir="src"
  log: "logs/BMO/download_bmo.log"
  benchmark: "benchmarks/BMO/download_bmo.tsv"
  container: "tfbench-bmo.sif"
  threads: 1
  shell: 
    """
      exec 2>>{log}
      git clone https://github.com/ParkerLab/BMO.git {params.out_dir}/BMO
    """
    
          #mv BMO ./{params.out_dir}
    # clone specifically the release, however that release has some version mismatches
    # wget -O {params.out_dir}/v1.0.zip https://github.com/ParkerLab/BMO/archive/refs/tags/v1.0.zip
    # unzip -o {params.out_dir}/v1.0.zip -d {params.out_dir}
    # rm {params.out_dir}/v1.0.zip
''' 

rule BMO_prepare_motif_txt:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_models="data/annotation/motif_models.rds"
  params:
    out_dir="data/BMO"
  log: "logs/BMO/prepare_motif_{cellularContext}_txt.log"
  benchmark: "benchmarks/BMO/prep_motif_{cellularContext}_txt.tsv"
  container: "tfbench-bmo.sif"
  threads: 1
  output:
    motif_list="data/BMO/motif_{cellularContext}.txt"
  script: "../src/utils/BMO/BMO_prepare_motif_list.R"
  
rule BMO_prepare_config:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_matching_report="output/FIMO/01_motif_matches_all.html",
    #config_yaml=rules.BMO_download.output.config_yaml,
    motif_list="data/BMO/motif_{cellularContext}.txt",
    atac_bam_dir="data/common/03_merged_cleaned_atac/merged_{cellularContext}.bam",
    atac_peak_dir="data/common/03_peaks_atac/filtered_merged_{cellularContext}_peaks.narrowPeak",
  params:
    config_yaml="/opt/BMO/config/config.yaml",
    motif_dir="data/motifs/motif_matches_genome",
    pred_dir="predictions/BMO",
    out_dir="data/BMO",
    bmo_dir="/opt/BMO", # src/BMO
  threads: 4
  log: "logs/BMO/prepare_config_{cellularContext}.log"
  benchmark: "benchmarks/BMO/prep_config_{cellularContext}.tsv"
  output: 
    config_yaml="data/BMO/input_{cellularContext}_config.yaml"
  container: "tfbench-bmo.sif"
  script: "../src/utils/BMO/BMO_prepare_config.R"

rule BMO_run_analysis:
  input:
    config_yaml="data/BMO/input_{cellularContext}_config.yaml",
  output:
    predictions=directory("predictions/BMO/{cellularContext}")
  wildcard_constraints:
    cellularContext="[^./]+"
  params:
    snakefile="/opt/BMO/Snakefile",
    out_dir="predictions/BMO",
  threads: 4
  log: "logs/BMO/run_analysis_{cellularContext}.log"
  benchmark: "benchmarks/BMO/train_pred_{cellularContext}.tsv"
  container: "tfbench-bmo.sif"
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.cellularContext}
      snakemake -s {params.snakefile} -d "$PWD" -j {threads} --resources io_limit=2 --forceall --rerun-incomplete --configfile {input.config_yaml} &> {log}
    """
    # -d {params.inner_wd}
    
rule BMO_move_results:
  input:
    predictions="predictions/BMO/{cellularContext}"
  output:
    prediction="predictions/BMO/{tf}_{cellularContext}.all.bed"
  benchmark: "benchmarks/BMO/move_results_{tf}_{cellularContext}.tsv"
  container: "tfbench-bmo.sif"
  log: "logs/BMO/move_results_{tf}_{cellularContext}.log"
  shell:
    """
      mv {input.predictions}/bmo/{wildcards.cellularContext}/all/{wildcards.tf}.all.bed {output.prediction} &> {log}
    """
