library(TOP)
library(data.table)
library(GenomicRanges)

sites <- fread(snakemake@input[["cand_sites"]])
atacMat <- readRDS(snakemake@input[["binned_atac_mat"]])
chipPeaksPath <- snakemake@input[["chip_peak_file"]]
siteRanges <- makeGRangesFromDataFrame(as.data.frame(sites))
peakRanges <- makeGRangesFromDataFrame(as.data.frame(fread(chipPeaksPath,
                                                           col.names=c("chr",
                                                                       "start",
                                                                       "end"))))
chipLabels <- data.frame(chip_label=as.integer(overlapsAny(siteRanges, peakRanges)))

# this fails if not all columns of a bed file are present, which we do not have after the merging.
# however what it does is equivalent to the above overlapsAny
#chipLabels <- add_chip_peak_labels_to_sites(sites, 
#                                            chip_peak_file=chipPeaksPath)
#colnames(chipLabels) <- "chip_label"
colnames(atacMat) <-  paste0('bin', 1:ncol(atacMat))
featMat <- cbind(sites, atacMat)
featMat <- cbind(featMat, chipLabels)

saveRDS(featMat, snakemake@output[["feature_matrix"]])