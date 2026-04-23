library(data.table)
library(yaml)

outFile <- snakemake@output[["config_yaml"]]
combinations <- fread(snakemake@input[["combinations"]])
testCombinations <- subset(combinations, set=="testData")

configPath <- snakemake@params[["config_yaml"]]
config <- read_yaml(configPath)
cellularContext <- unlist(tstrsplit(outFile, split="_", keep=2))

if(!dir.exists(dirname(outFile))) dir.create(dirname(outFile))

config$results <- file.path(snakemake@params[["pred_dir"]], cellularContext) #file.path("../../..", snakemake@params[["pred_dir"]], cellularContext)
config$motif_dir <- snakemake@params[["motif_dir"]] #file.path("../../..", snakemake@params[["motif_dir"]])
config$motif_file <- snakemake@input[["motif_list"]] #file.path("..", gsub("data/BMO/", "", snakemake@input[["motif_list"]]))
config$bmo_dir <- snakemake@params[["bmo_dir"]] #file.path("../..", snakemake@params[["bmo_dir"]]) # ""
config$motif_ext <- "bed"
  
config$samples <- NULL
config$samples[[cellularContext]]$bamfile <- snakemake@input[["atac_bam_dir"]] #file.path("../../..", snakemake@input[["atac_bam_dir"]])
config$samples[[cellularContext]]$peakfile <- snakemake@input[["atac_peak_dir"]] #file.path("../../..", snakemake@input[["atac_peak_dir"]])
write_yaml(config, outFile)
