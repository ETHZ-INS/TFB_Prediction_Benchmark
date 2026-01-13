library(data.table)
library(universalmotif)

combinations <- fread(snakemake@input[["combinations"]])
tf <- basename(dirname(snakemake@output[["merged_motif_scores"]]))
tempDir <- file.path(snakemake@params[["temp_train_dir"]], tf)
trainData <- subset(combinations, TF==tf & set=="trainData")

# TODO: Sort within the loops!

# 1. Merge ATAC-Data
chromDts <- list()
for(context in unique(trainData$CellularContext)){
  chromDt <- fread(file.path("data/Catchitt/atac", context, 
                              "Chromatin_accessibility.tsv.gz"))
  setorder(chromDt, V1, V2)
  chromDts[[context]] <- chromDt
}
mergedChromDt <- rbindlist(chromDts)
mergedChromDt <- subset(mergedChromDt, V1 %in% paste0("chr", 1:22))

if(!dir.exists(tempDir)) dir.create(tempDir)
fwrite(mergedChromDt, file.path(tempDir,
                                "Chromatin_accessibility.tsv.gz"), 
       sep="\t",
       compress="gzip",
       row.names=FALSE, col.names=FALSE, quote=FALSE,
       encoding="UTF-8", na="NaN")
rm(mergedChromDt)

# 2. Merge Motif-Data
motifDts <- list()
for(i in 1:length(trainData$CellularContext)){
  motifDt <- fread(file.path("data/Catchitt/motifs", tf, "Motif_scores.tsv.gz"))
  setorder(motifDt, V1, V2)
  motifDts[[i]] <- motifDt
}
mergedMotifDt <- rbindlist(motifDts)
#mergedMotifDt[,V1:=paste0("chr", V1)]
mergedMotifDt <- subset(mergedMotifDt, V1 %in% paste0("chr", 1:22))
#mergedMotifDt[,V3:=fifelse(is.na(V3), 0, V3)]
#mergedMotifDt[,V4:=fifelse(is.na(V4), 0, V4)]
fwrite(mergedMotifDt, file.path(tempDir,
                                "Motif_scores.tsv.gz"), 
       sep="\t",
       compress="gzip", row.names=FALSE, col.names=FALSE, quote=FALSE,
       encoding="UTF-8", na="NaN")
rm(mergedMotifDt)

# 3. Merge ChIP-Labels
labelDts <- list()
for(context in unique(trainData$CellularContext)){
  labelDt <- fread(file.path("data/Catchitt/chip", 
                             paste(tf, context, sep="_"), 
                             "Labels.tsv.gz"))
  setorder(labelDt, V1, V2)
  labelDts[[context]] <- labelDt
}
labelDt <- rbindlist(labelDts)
labelDt <- subset(labelDt, V1 %in% paste0("chr", 1:22))
fwrite(labelDt, file.path(tempDir, "Labels.tsv.gz"), 
       sep="\t",
       compress="gzip", row.names=FALSE, col.names=FALSE, quote=FALSE, 
       encoding="UTF-8", na="NaN")
rm(labelDt)

if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
  logf <- snakemake@log[[1]]
  con <- file(logf, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message", append = TRUE)
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(con), silent = TRUE)
  }, add = TRUE)
  options(error = function() {
    traceback()
    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(con), silent = TRUE)
    quit(status = 1)
  })
}