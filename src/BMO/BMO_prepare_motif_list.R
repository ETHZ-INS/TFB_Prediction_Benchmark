library(data.table)

motifModels <- readRDS(snakemake@input[["motif_models"]])
outPath <- snakemake@output[["motif_list"]]
combs <- fread(snakemake@input[["combinations"]])
testComb <- subset(combs, set=="testData")

motifDt <- data.table(motif=unique(names(motifModels)))
motifDt <- subset(motifDt, motif %in% testComb$TF)

if(grepl("Jurkat", outPath)){
  motifDt <- subset(motifDt, motif!="RXRA")
}

write.table(motifDt, file=outPath, quote=FALSE, 
            row.names=FALSE, col.names=FALSE)