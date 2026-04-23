# cp <- Sys.getenv("CONDA_PREFIX")
# .libPaths(c(file.path(cp, "lib", "R", "library"), .libPaths()))
# print(cp)
# cp <- Sys.getenv("CONDA_PREFIX")
# library(withr)
# with_envvar(
#   c(
#     JAGS_MODULE_PATH = file.path(cp, "lib", "JAGS", "modules-4"),
#     LD_LIBRARY_PATH = paste(file.path(cp,"lib"), file.path(cp,"lib","R","lib"),
#                             Sys.getenv("LD_LIBRARY_PATH"), sep=":")
#   ),
#   { library(rjags) }
# )
#install.packages("rjags")

#load.module("basemod", path = file.path(Sys.getenv("CONDA_PREFIX"), "lib", "JAGS", "modules-4"))
library(TOP)
library(data.table)

tf <- tstrsplit(basename(snakemake@output[["posterior_samples"]]),
                split="_", keep=2)
combinations <- fread(snakemake@input[["combinations"]])
trainData <- subset(combinations, set=="trainData")
trainData <- unique(trainData, by=c("TF", "CellularContext"))

trainContexts <-  unique(trainData$CellularContext)
featMatDir <- snakemake@params[["feat_mat_dir"]]
trainChrs <- snakemake@params[["train_chrs"]]
seed <- snakemake@params[["seed"]]
outDir <- snakemake@params[["out_dir"]]

# assemble the training data
tfCellTable <- data.table(tf_name=trainData$TF, 
                          cell_type=trainData$CellularContext)
tfCellTable[,data_file:=file.path(featMatDir,
                                  paste("feature_matrix",
                                        tf_name, cell_type, sep="_"))]
tfCellTable[,data_file:=paste0(data_file, ".rds")]


tfCellTable <- as.data.frame(tfCellTable)
trainData <- assemble_training_data(tfCellTable,
                                    logistic_model = TRUE,
                                    chip_col = 'chip_label',
                                    training_chrs = trainChrs,
                                    seed=seed)

if(!dir.exists(outDir)) dir.create(outDir)
allTOPSamples <- fit_TOP_M5_model(all_training_data=trainData,
                                  logistic_model = TRUE,
                                  n_iter = 5000,
                                  n_burnin = 2000,
                                  n_chains = 3,
                                  n_thin = 2,
                                  n_cores=snakemake@threads,
                                  outdir = outDir,
                                  return_type="samplefiles")
#saveRDS(allTOPSamples, "/mnt/germain/esonder/benchmark_finalized/benchmark_final_2/models/TOP/logistic/test_output.rds")
#allTOPSamples <- list.files(dirname(snakemake@output[["posterior_samples"]]),
#                            pattern="posterior", full.names=TRUE)
TOPSamples <- combine_TOP_samples(allTOPSamples)
TOPMeanCoef <- extract_TOP_mean_coef(TOPSamples,
                                     assembled_training_data=trainData)
saveRDS(TOPSamples, snakemake@output[["posterior_samples"]])
saveRDS(TOPMeanCoef, snakemake@output[["posterior_mean"]])
