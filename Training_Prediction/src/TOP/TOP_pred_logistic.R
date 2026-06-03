log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)


library(TOP)
library(data.table)

outDir <- snakemake@params[["out_dir"]]
tf <- unlist(tstrsplit(basename(snakemake@output[["pred_file"]]), split="_", 
                       keep=2))
context <- unlist(tstrsplit(basename(snakemake@output[["pred_file"]]), split="_", 
                            keep=3))
posteriorMean <- readRDS(snakemake@input[["posterior_mean"]])

featMat <- readRDS(snakemake@input[["feature_matrix"]])
# Q: Does chip_label column need to be removed explicitly?
TOPResult <- predict_TOP(featMat, 
                         TOP_coef=posteriorMean,
                         tf_name=tf,
                         cell_type=context,
                         use_model='ATAC',
                         level='best',
                         logistic_model = TRUE)
if(!dir.exists(outDir)) dir.create(outDir)
saveRDS(TOPResult, snakemake@output[["pred_file"]])


