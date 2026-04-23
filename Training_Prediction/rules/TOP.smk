os.makedirs('data/TOP', exist_ok=True)

# change here path to your bwtool install, in case of not using the container
######################################################################################
BWTOOLPATH="/usr/local/bin/bwtool"
######################################################################################

def combs_train_list():
    df = read_samples_df()
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    combs = df[df["set"]== "trainData"]["comb"].unique().tolist()

    return combs

def top_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    
    combs = df["comb"].unique().tolist()
    tfs = df["TF"].unique().tolist()
    
    combsTest = df[df["set"] == "testData" ]["comb"].unique().tolist()
    tfs, contexts = zip(*(s.rsplit("_", 1) for s in combsTest)) if combsTest else ([], [])
    
    pred_logistic_files = [
      f"predictions/TOP/{HOLDOUT}/logistic/pred_{tf}_{context}.rds"
      for tf, context in zip(tfs, contexts)
    ]
    
    return pred_logistic_files
 
ALL_TARGET_BUILDERS.append(top_targets) 
  
def read_samples_df(which=None):
    ckpt1 = checkpoints.common_get_combinations.get()   # waits for checkpoint
    ckpt2 = checkpoints.common_get_encode_chip_coverage_ids.get()
    
    if which=="atac":
      df=pd.read_csv(ckpt1.output[1], sep="\t")
    elif which=="chip":
      df=pd.read_csv(ckpt2.output[0], sep="\t")
    else:
      df=pd.read_csv(ckpt1.output[0], sep="\t")
      
    return df

def group_accession_chip_context_inputs(wildcards):
    df = read_samples_df(which="chip")
   
    accession_subset = df.loc[df.TF == wildcards.tf and df.CellularContext==wildcards.cellularContext, "Accession"].tolist()

    return expand("data/common/01_raw_chip_coverage/{accession}.bw", accession=accession_subset)
  
ALL_TARGET_BUILDERS.append(top_targets)  

rule TOP_prepare_motif_matches:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_matching_report="output/FIMO/01_motif_matches_all.html", # TODO: Switch to all
  output:
    cand_sites="data/TOP/motif_matches/candidate_sites_{tf}.tsv",
  params:
    black_list_path="data/annotation/blacklist.bed",
    motif_match_dir="data/motifs/motif_matches_genome"
  threads: 4
  benchmark: "benchmarks/TOP/prep_prepare_motif_matches_{tf}.tsv"
  log: "logs/TOP/prepare_motif_matches_{tf}.log"
  container: "tfbench-top.sif",
  script: "../src/utils/TOP/TOP_process_candidate_sites.R"
  
rule TOP_count_genome_cuts:
  input:
    atac_bam_file="data/common/03_merged_cleaned_atac/merged_{cellularContext}.bam"
  output:
    atac_fwd_bw_file="data/TOP/atac_cut_counts/cuts_{cellularContext}.fwd.genomecounts.bw",
    atac_rev_bw_file="data/TOP/atac_cut_counts/cuts_{cellularContext}.rev.genomecounts.bw",
  params:
    chrom_sizes="data/annotation/reChr.sizes",
    out_dir="data/TOP/atac_cut_counts"
  benchmark: "benchmarks/TOP/prep_count_genome_cuts_{cellularContext}.tsv"
  log: "logs/TOP/count_genome_cuts_{cellularContext}.log"
  threads: 4
  container: "tfbench-top.sif",
  script: "../src/utils/TOP/TOP_count_genome_cuts.R"

rule TOP_normalize_count_matrix:
  input:
    cand_sites="data/TOP/motif_matches/candidate_sites_{tf}.tsv",
    atac_fwd_bw_file="data/TOP/atac_cut_counts/cuts_{cellularContext}.fwd.genomecounts.bw",
    atac_rev_bw_file="data/TOP/atac_cut_counts/cuts_{cellularContext}.rev.genomecounts.bw",
    idxstats_file="data/common/03_merged_cleaned_atac/merged_{cellularContext}.idxstats.txt",
  output:
    norm_count_matrix="data/TOP/atac_cut_counts/count_matrix_{tf}_{cellularContext}.rds"
  params:
    bwtool_path=BWTOOLPATH
  benchmark: "benchmarks/TOP/prep_normalize_count_matrix_{tf}_{cellularContext}.tsv"
  log: "logs/TOP/normalize_count_matrix_{tf}_{cellularContext}.log"
  threads: 4
  container: "tfbench-top.sif"
  script: "../src/utils/TOP/TOP_normalize_count_matrix.R"

rule TOP_assemble_matrix_logistic:
  input:
    cand_sites="data/TOP/motif_matches/candidate_sites_{tf}.tsv",
    chip_peak_file="data/common/02_merged_chip_peaks/filtered_merged_peaks_{tf}_{cellularContext}.bed",
    binned_atac_mat="data/TOP/atac_cut_counts/count_matrix_{tf}_{cellularContext}.rds",
  output:
    feature_matrix=temp("data/TOP/{HOLDOUT}/feature_matrix_logistic/feature_matrix_{tf}_{cellularContext}.rds")
  benchmark: "benchmarks/TOP/prep_assemble_matrix_logistic_{HOLDOUT}_{tf}_{cellularContext}.tsv"
  log: "logs/TOP/assemble_matrix_logistic_{HOLDOUT}_{tf}_{cellularContext}.log"
  threads: 4
  container: "tfbench-top.sif",
  script: "../src/utils/TOP/TOP_assemble_matrix.R"

rule TOP_train_logistic:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    feat_mats=lambda wildcards: expand("data/TOP/{holdout}/feature_matrix_logistic/feature_matrix_{comb}.rds", 
                                       comb=combs_train_list(), holdout=HOLDOUT),
  output:
    posterior_samples="models/TOP/{HOLDOUT}/TOP_samples_logistic.rds",
    posterior_mean="models/TOP/{HOLDOUT}/TOP_mean_logistic.rds"
  params:
    seed=SEED,
    out_dir="models/TOP/{HOLDOUT}/logistic",
    feat_mat_dir="data/TOP/{HOLDOUT}/feature_matrix_logistic",
    train_chrs=AUT_CHRS
  benchmark: "benchmarks/TOP/train_logistic_{HOLDOUT}.tsv"
  log: "logs/TOP/train_logistic_{HOLDOUT}.tsv"
  threads: 4
  container: "tfbench-top.sif",
  script: "../src/utils/TOP/TOP_train_logistic.R"

rule TOP_predict_logistic:
  input:
    posterior_samples="models/TOP/{HOLDOUT}/TOP_samples_logistic.rds",
    posterior_mean="models/TOP/{HOLDOUT}/TOP_mean_logistic.rds",
    feature_matrix="data/TOP/{HOLDOUT}/feature_matrix_logistic/feature_matrix_{tf}_{cellularContext}.rds"
  output:
    pred_file="predictions/TOP/{HOLDOUT}/logistic/pred_{tf}_{cellularContext}.rds"
  params:
    out_dir="predictions/TOP/{HOLDOUT}/logistic",
  benchmark: "benchmarks/TOP/pred_logistic_{HOLDOUT}_{tf}_{cellularContext}.tsv"
  log: "logs/TOP/predict_logistic_{HOLDOUT}_{tf}_{cellularContext}.tsv"
  threads: 4
  container: "tfbench-top.sif",
  script: "../src/utils/TOP/TOP_pred_logistic.R"
