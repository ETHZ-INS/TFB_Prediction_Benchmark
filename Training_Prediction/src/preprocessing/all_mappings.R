library(data.table)

args <- commandArgs(trailingOnly=TRUE)
allCombs <- fread(args[1])

rawATACEncode <- fread(args[2])
rawATACEncode <- subset(rawATACEncode, !is.na(accession))
rawATACEncode <- rawATACEncode[,c("accession", "CellularContext"), with=FALSE]

# TODO: Refactor this and read directly from the sraRunTables.R !!!
rawATACNonEncode <- data.table(accession=c("merged_H1", 
                                           "merged_Jurkat", 
                                           "merged_HeLa-S3"),
                               CellularContext=c("H1", 
                                                 "Jurkat",
                                                 "HeLa-S3"))
rawATACAll <- rbind(rawATACEncode, 
                    rawATACNonEncode)
rawATACAll <- subset(rawATACAll, 
                     CellularContext %in% unique(allCombs$CellularContext))

outDir <- args[3]
rawEncodeATACDir <- args[4]
rawNonEncodeATACDir <- args[5]
outFile <- args[6]

if(!dir.exists(outDir)) dir.create(outDir)

for(file in list.files(rawEncodeATACDir, full.names=TRUE)){
  file.copy(file, file.path(outDir, basename(file)))
}

for(file in list.files(rawNonEncodeATACDir, full.names=TRUE)){
  file.copy(file, file.path(outDir, basename(file)))
}

write.table(rawATACAll, file=outFile, 
            quote=FALSE, row.names=FALSE, sep="\t")