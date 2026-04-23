library(httr)
library(jsonlite)
library(data.table)

#args = commandArgs(trailingOnly=TRUE)
# trainData <- args[1]
# testData <- args[2]
# tfsSub <- args[3]
# outPath <- args[4]

trainData <- snakemake@input[["train_data_identifiers"]]
tfsSub <- snakemake@params[["tfs_sub"]]

outPath <- snakemake@output[["downloaded_chip_coverage_ids"]]

trainData <- fread(trainData)
trainData$set <- "trainData"

idDt <- trainData

if(!is.null(tfsSub)){
  idDt <- subset(idDt, TF %in% tfsSub)  
}

metaCovAll <- data.table()
for(id in idDt$Accession){
  set <- subset(idDt, Accession==id)$set
  cellularContext <- subset(idDt, Accession==id)$CellularContext
  statusAllowed <- "released"
  assemblyAllowed <- "GRCh38"
  fileType <- "bam"
  outputType <- "alignments"
  
  tryCatch({
    # get ENCODE meta-data
    url <- paste0("https://www.encodeproject.org/experiments/", 
                  id, "/?format=json")
    response <- GET(url)
    content <- httr::content(response, "text")
    data <- fromJSON(content)
    filesDt <- as.data.table(data$files)
  
    # get the most recent data
    filesDt[,date_created_short:=tstrsplit(date_created, split="T", keep=1)]
    filesDt[,date_created_short:=as.Date(filesDt$date_created_short)]
    mostRecentDate <- max(filesDt$date_created_short)

    metaCov <- subset(filesDt, date_created_short==mostRecentDate &
                               file_type==fileType & 
                               output_type==outputType & 
                               !grepl("unfiltered", output_type) &
                               !grepl("unfiltered", file_type) &
                               status==statusAllowed & 
                               assembly==assemblyAllowed)
    if(nrow(metaCov)==0){
      metaCov <- subset(filesDt, 
                          file_type==fileType & 
                          output_type==outputType & 
                          !grepl("unfiltered", output_type) &
                          !grepl("unfiltered", file_type) &
                          status==statusAllowed & 
                          assembly==assemblyAllowed)
    }
    
    expSpecDt <- subset(idDt, Accession==id)

    metaCov$TF <- unique(expSpecDt$TF)
    metaCov$CellularContext <- unique(expSpecDt$CellularContext)
    metaCov$Exp_Accession <- id
    metaCov$set <- set
    setnames(metaCov, "accession", "Accession")
    metaCov <- metaCov[,c("Accession", "href", "TF", "CellularContext", 
                          "Exp_Accession", "set",
                          "output_type",	"file_type"), with=FALSE]
    metaCov <- subset(metaCov, !grepl("unfiltered", output_type))

    metaCovAll <- rbind(metaCovAll, metaCov, fill=TRUE, use.names=TRUE)
    write.table(metaCovAll, file=outPath, quote=FALSE, 
                sep="\t", row.names=FALSE)},
    error=function(cond){
      NULL
    })
}
