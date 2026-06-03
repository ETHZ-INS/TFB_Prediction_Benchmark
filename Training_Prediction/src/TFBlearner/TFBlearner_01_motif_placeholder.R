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
motifModels <- readRDS(snakemake@params[["motif_models_path"]])
outDir <- snakemake@params[["out_dir"]]
if(!dir.exists(outDir)) dir.create(outDir)

motifNames <- names(motifModels)
motifNames <- gsub("-", ".", motifNames)
for(tfName in motifNames){
  fwrite(data.table(tfName=tfName), 
         file=file.path(outDir, paste0(tfName, ".tsv")),
         sep="\t")
}
fwrite(data.table(motif=motifNames), 
       file=snakemake@output[["placeholder_motif_files"]],
       sep="\t")