library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)
library(TFBlearner)
library(data.table)
library(GenomicRanges)
library(BiocParallel)
library(universalmotif)
library(TFBSTools)

log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")

setDTthreads(as.integer(snakemake@threads))
dhsCoords <- fread(snakemake@params[["dhs_path"]], col.names=c("chr", "start", "end"))
dhsCoords <- dhsCoords[chr %in% paste0("chr", 1:22)]
dhsCoords <- makeGRangesFromDataFrame(as.data.frame(dhsCoords))

motifModels <- readRDS(snakemake@params[["motif_models_path"]])
outDir <- snakemake@params[["out_dir"]]
fa <- FaFile(snakemake@params[["genome_path"]])

tfName <- snakemake@wildcards[["motif"]]
tfName <- gsub(".", "-", tfName, fixed=TRUE)
motifModel <- motifModels[tfName]
names(motifModel) <- gsub("-", ".", names(motifModel), fixed=TRUE)
mm <- convert_motifs(motifModel, "TFBSTools-PFMatrix")
mms <- lapply(mm, toPWM, type="prob")
mms <- do.call(PWMatrixList, mms)

mm <- TFBlearner::prepMotifs(dhsCoords, mms, outputFolder=outDir, genome=fa, BPPARAM=SerialParam())
#saveRDS(mm, snakemake@output[["processed_motifs"]])
