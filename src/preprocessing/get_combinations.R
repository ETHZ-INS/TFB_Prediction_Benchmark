library(data.table)

outCombPath <- snakemake@output[["combinations"]]
holdOut <- snakemake@params[["holdout"]]
combData <- fread(snakemake@input[["downloaded_chip_peak_ids"]])
combData <- unique(combData, by=c("TF", "CellularContext", "set"))[,c("TF",
                                                                      "CellularContext",
                                                                      "set"), 
                                                                   with=FALSE]
if(holdOut=="train"){
  combData <- subset(combData, set!="testData")
  combData[,set:=fifelse(CellularContext %in% c("HepG2", "MCF-7"), 
                         "testData", "trainData")]
}
write.table(combData, file=outCombPath, quote=FALSE, row.names=FALSE, sep="\t")

outAtacPath <- snakemake@output[["all_raw_atac"]]
atacData <- fread(snakemake@input[["downloaded_atac_ids"]])
atacData <- atacData[,c("accession", "CellularContext"), with=FALSE]
atacData[,accession:=fifelse(CellularContext=="HeLa-S3", "merged_HeLa-S3", accession)]
atacData[,accession:=fifelse(CellularContext=="H1", "merged_H1", accession)]
atacData[,accession:=fifelse(CellularContext=="Jurkat", "merged_Jurkat", accession)]

write.table(atacData, file=outAtacPath, quote=FALSE, row.names=FALSE, sep="\t")