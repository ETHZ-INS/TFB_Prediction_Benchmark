os.makedirs('data/Catchitt', exist_ok=True)

BINSIZE=50

def catchitt_targets(ckpt):
    df = pd.read_csv(ckpt.output.combinations, sep="\t")
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    
    combs = df["comb"].unique().tolist()
    tfs = df["TF"].unique().tolist()
    
    # combsTrain = df[df["set"] == "trainData" ]["comb"].unique().tolist()
    # proc_combs = [
    #     f"data/Catchitt/chip/{comb}/Labels.tsv.gz"
    #     for comb in combsTrain
    # ]
    # 
    # proc_contexts = [
    #   f"data/Catchitt/atac/{context}/Chromatin_accessibility.tsv.gz"
    #   for context in contexts_list()
    # ]
    # 
    # proc_tfs = [
    #   f"data/Catchitt/motifs/{tf}/Motif_scores.tsv.gz"
    #   for tf in tfs
    # ]
    # 
    # html_reports = [ "output/Catchitt/Catchitt_convert_motifs.html" ]
    # 
    # tfs = [x for x in tfs if x not in to_drop]  
    # models = [
    #   f"models/Catchitt/{HOLDOUT}/{tf}/Classifiers.xml"
    #   for tf in tfs
    # ]
    

    testSub =  df[df["set"] == "testData" ]
    combsTest = testSub["comb"].unique().tolist()
    tfs, contexts = zip(*(s.rsplit("_", 1) for s in combsTest)) if combsTest else ([], [])

    pred_combs = [
       f"predictions/Catchitt/{HOLDOUT}/all/{tf}/{context}_pred.bw"
       for tf, context in zip(tfs, contexts)
    ]
    catchitt_outputs = pred_combs

    return catchitt_outputs
  
ALL_TARGET_BUILDERS.append(catchitt_targets)  
  
# These can be defined in the main Snakefile at the end
def read_samples_df():
    ckpt = checkpoints.common_get_combinations.get()
    return pd.read_csv(ckpt.output[0], sep="\t")  
  
def combs_train_list():
    df = read_samples_df()
    df["comb"] = df["TF"].astype("string").fillna("NA") + "_" + df["CellularContext"].astype("string").fillna("NA")
    combs = df[df["set"]== "trainData"]["comb"].unique().tolist()

    return combs
  
def contexts_list():
    return read_samples_df()["CellularContext"].unique().tolist()
  
def contexts_train_list():
    df= read_samples_df()
    contexts = df[df["set"] == "trainData" ]["CellularContext"].unique().tolist()
    
    return contexts

def tfs_list():
    return read_samples_df()["TF"].unique().tolist()
  
def tfs_list2():
    ckpt = checkpoints.Catchitt_convert_motifs_jaspar.get()
    print(read_samples_df()["TF"].unique().tolist())
    return read_samples_df()["TF"].unique().tolist()

rule Catchitt_process_atac:
  input:
    combs="data/common/04_train_test_combinations.tsv",
    bam_file="data/common/03_merged_cleaned_atac/merged_{cellularContext}.bam",
  output:
    catchitt_atac="data/Catchitt/atac/{cellularContext}/Chromatin_accessibility.tsv.gz",
  params:
    out_dir="data/Catchitt/atac",
    bin_size=BINSIZE,
    java_opts="-Xms48g -Xmx48g"
  benchmark: "benchmarks/Catchitt/prep_process_atac_{cellularContext}.tsv"
  log: "logs/Catchitt/process_atac_{cellularContext}.log"
  #conda: "../envs/Catchitt.yml"
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.cellularContext}
      
      export _JAVA_OPTIONS=
      export JAVA_TOOL_OPTIONS=
      export JDK_JAVA_OPTIONS=
      java {params.java_opts} -jar Catchitt.jar access d="BAM/SAM" i={input.bam_file} b={params.bin_size} outdir={params.out_dir}/{wildcards.cellularContext} &> {log}
    """
    
rule Catchitt_derive_chip_labels:
  input:
    cons_peaks="data/common/02_merged_chip_peaks/filtered_merged_peaks_{comb}.bed",
    nonCons_peaks="data/common/02_merged_chip_peaks/filtered_nonCons_peaks_{comb}.bed"
  output:
     catchitt_labels="data/Catchitt/chip/{comb}/Labels.tsv.gz", 
  params:
    out_dir="data/Catchitt/chip",
    fai="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa.fai",
    java_opts="-Xms48g -Xmx48g",
    bin_size=BINSIZE
  benchmark: "benchmarks/Catchitt/prep_derive_chip_labels_{comb}.tsv"
  log: "logs/Catchitt/derive_chip_labels_{comb}.log"
  #conda: "../envs/Catchitt.yml"
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.comb}
      # $(dirname {output.catchitt_labels})
      export _JAVA_OPTIONS=
      export JAVA_TOOL_OPTIONS=
      export JDK_JAVA_OPTIONS=

      java {params.java_opts} -jar Catchitt.jar labels b={params.bin_size} \
      c={input.cons_peaks} \
      r={input.nonCons_peaks} \
      f={params.fai} \
      outdir={params.out_dir}/{wildcards.comb} &> {log}
    """
    
checkpoint Catchitt_convert_motifs_jaspar:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    motif_models="data/motifs/motif_models.rds",
  output:
    motif_conversion_report="output/Catchitt/Catchitt_convert_motifs.html",
  params:
    out_dir="output/Catchitt",
    motif_jaspar_dir="data/motifs/jaspar"
  benchmark: "benchmarks/Catchitt/prep_convert_motifs_jaspar.tsv"
  log: "logs/Catchitt/convert_motifs_jaspar.log"
  #conda: "../envs/Catchitt.yml"
  threads: 1
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      mkdir -p {params.motif_jaspar_dir}
      Rscript -e 'rmarkdown::render("src/Catchitt/Catchitt_convert_motifs.Rmd", 
                                    "html_document", 
                                     output_file="../../{output.motif_conversion_report}",
                                     params=list(outDir="{params.motif_jaspar_dir}",
                                                 motifModelsPath="{input.motif_models}"))'
    """    

rule Catchitt_get_motifs_scores:
  input:
    motif_conversion_report="output/Catchitt/Catchitt_convert_motifs.html",
    motif_jaspar="data/motifs/jaspar/{tf}.jaspar",
    #motif_jaspar=lambda wildcards: expand("data/motifs/jaspar/{tf}.jaspar", 
    #                                      tf=tfs_list())
  output:
    catchitt_motifscores="data/Catchitt/motifs/{tf}/Motif_scores.tsv.gz"
  params:
    out_dir="data/Catchitt/motifs",
    bin_size=BINSIZE,
    java_opts="-Xms48g -Xmx48g",
    fai="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa.fai",
    fa="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa",
  benchmark: "benchmarks/Catchitt/prep_get_motifs_scores_{tf}.tsv"
  log: "logs/Catchitt/get_motifs_scores_{tf}.log"
  #conda: "../envs/Catchitt.yml"
  threads: 4
  shell:
    """
      echo "{input.motif_jaspar}"
      mkdir -p {params.out_dir}/{wildcards.tf}
      
      export _JAVA_OPTIONS=
      export JAVA_TOOL_OPTIONS=
      export JDK_JAVA_OPTIONS=
      
      java {params.java_opts} -jar Catchitt.jar motif m=Jaspar \
      j="data/motifs/jaspar/{wildcards.tf}.jaspar" \
      g={params.fa} \
      f={params.fai} \
      b={params.bin_size} \
      outdir={params.out_dir}/{wildcards.tf} \
      threads={threads} &> {log}
    """

rule Catchitt_merge_train:
  input:
    combinations="data/common/04_train_test_combinations.tsv",
    catchitt_atac=lambda wildcards: expand("data/Catchitt/atac/{cellularContext}/Chromatin_accessibility.tsv.gz", 
                                           cellularContext=contexts_train_list()),
    catchitt_labels=lambda wildcards: expand("data/Catchitt/chip/{comb}/Labels.tsv.gz", 
                                             comb=combs_train_list()),
    motif_scores=lambda wildcards: expand("data/Catchitt/motifs/{tf}/Motif_scores.tsv.gz", 
                                          tf=tfs_list())
  output:
    merged_labels=temp("data/Catchitt/train/{HOLDOUT}/{tf}/Chromatin_accessibility.tsv.gz"),
    merged_motif_scores=temp("data/Catchitt/train/{HOLDOUT}/{tf}/Motif_scores.tsv.gz"),
    merged_atac=temp("data/Catchitt/train/{HOLDOUT}/{tf}/Labels.tsv.gz")
  params:
    temp_train_dir="data/Catchitt/train/{HOLDOUT}",
  benchmark: "benchmarks/Catchitt/prep_merge_train_{HOLDOUT}_{tf}.tsv"
  log: "logs/Catchitt/merge_train_{HOLDOUT}_{tf}.log"
  threads: 4
  script: "../src/Catchitt/Catchitt_merge_train_data.R"
  
rule Catchitt_train:
  input:
    merged_atac="data/Catchitt/train/{HOLDOUT}/{tf}/Chromatin_accessibility.tsv.gz",
    merged_motif_scores="data/Catchitt/train/{HOLDOUT}/{tf}/Motif_scores.tsv.gz",
    merged_labels="data/Catchitt/train/{HOLDOUT}/{tf}/Labels.tsv.gz",
  output:
    models="models/Catchitt/{HOLDOUT}/{tf}/Classifiers.xml"
  params:
    out_dir="models/Catchitt/{HOLDOUT}",
    fa_path="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa.fai",
    bin_size=BINSIZE,
    java_opts="-Xms100g -Xmx100g",
    train_chrs=",".join(AUT_CHRS)
  benchmark: "benchmarks/Catchitt/train_{HOLDOUT}_{tf}.tsv"
  log: "logs/Catchitt/train_{HOLDOUT}_{tf}.log"
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.tf}
      export _JAVA_OPTIONS=
      export JAVA_TOOL_OPTIONS=
      export JDK_JAVA_OPTIONS=
      java {params.java_opts} -jar Catchitt.jar itrain a={input.merged_atac} \
      m={input.merged_motif_scores} \
      l={input.merged_labels} \
      f={params.fa_path} \
      t="{params.train_chrs}" \
      b={params.bin_size} \
      threads={threads} \
      outdir={params.out_dir}/{wildcards.tf} &> {log}
    """

rule Catchitt_predict_all:
  input:
    model="models/Catchitt/{HOLDOUT}/{tf}/Classifiers.xml",
    atac_file="data/Catchitt/atac/{context}/Chromatin_accessibility.tsv.gz",
    motif_scores="data/Catchitt/motifs/{tf}/Motif_scores.tsv.gz",
  output:
    predictions="predictions/Catchitt/{HOLDOUT}/all/{tf}/{context}_pred.bw"
  params:
    out_dir="predictions/Catchitt/{HOLDOUT}/all",
    fai="data/annotation/Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa.fai",
    java_opts="-Xms48g -Xmx48g",
    val_chrs=",".join(VAL_CHRS)
  benchmark: "benchmarks/Catchitt/pred_{HOLDOUT}_{tf}_{context}.tsv"
  log: "logs/Catchitt/predict_all_{HOLDOUT}_{tf}_{context}.log"
  threads: 4
  shell:
    """
      mkdir -p {params.out_dir}/{wildcards.tf}
      export _JAVA_OPTIONS=
      export JAVA_TOOL_OPTIONS=
      export JDK_JAVA_OPTIONS=
      
      java {params.java_opts} -jar Catchitt.jar predict \
      c={input.model} \
      a={input.atac_file} \
      m={input.motif_scores} \
      p="{params.val_chrs}" \
      f={params.fai} \
      outdir={params.out_dir}/{wildcards.tf} &> {log}
    """
