library(httr)
library(jsonlite)
library(data.table)

trainData <- snakemake@input[["train_data_identifiers"]]
testData <- snakemake@input[["test_data_identifiers"]]
tfsSub <- snakemake@params[["tfs_sub"]]

outPath <- snakemake@output[["downloaded_chip_peak_ids"]]

trainData <- fread(trainData)
trainData$set <- "trainData"
testData <- fread(testData)
testData$set <- "testData"
testData <- subset(testData, Type=="ChIP")
testData <- testData[,colnames(trainData), with=FALSE]

idDt <- rbind(trainData, testData, use.names=TRUE)

if(!is.null(tfsSub)){
  idDt <- subset(idDt, TF %in% tfsSub)  
}

metaPeaksAll <- data.table()
for(id in idDt$Accession){
  set <- subset(idDt, Accession==id)$set
  cellularContext <- subset(idDt, Accession==id)$CellularContext
  statusAllowed <- fifelse(cellularContext=="Jurkat", "archived","released")
  fileType <- fifelse(cellularContext=="Jurkat", "bed bed3+","bed narrowPeak")
  assemblyAllowed <- fifelse(cellularContext=="Jurkat", "hg19", "GRCh38")
  outputType <- fifelse(cellularContext=="Jurkat", 
                        "optimal IDR thresholded peaks",
                        "IDR thresholded peaks")
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

    subFilesDt <- subset(filesDt, date_created_short==mostRecentDate &
                                  file_type==fileType & 
                                  status==statusAllowed & 
                                  assembly==assemblyAllowed)

    consPeaks <- subset(subFilesDt, 
                        output_type=="conservative IDR thresholded peaks")
    consPeaks <- consPeaks[,c("accession", "output_type", 
                              "file_type", "href"), with=FALSE]
    consPeaks$peak_type <- "conservative"
    libPeaks <- subset(subFilesDt, output_type==outputType)
    if(nrow(libPeaks)>1 & 
       "Rep 1, Rep 2" %in% unique(libPeaks$biological_replicates_formatted)){
        libPeaks <- subset(libPeaks, 
                           biological_replicates_formatted=="Rep 1, Rep 2")
    }
    
    #TODO: Add Optimal IDR thresholded peaks (if these exist)
    libPeaks <- libPeaks[,c("accession", "output_type", 
                            "file_type", "href"), with=FALSE]
    libPeaks$peak_type <- "liberal"

    metaPeaks <- rbind(consPeaks, libPeaks)
    expSpecDt <- subset(idDt, Accession==id)

    metaPeaks$TF <- unique(expSpecDt$TF)
    metaPeaks$CellularContext <- unique(expSpecDt$CellularContext)
    metaPeaks$Exp_Accession <- id
    metaPeaks$set <- set

    metaPeaksAll <- rbind(metaPeaksAll, metaPeaks, fill=TRUE, use.names=TRUE)
    write.table(metaPeaksAll, file=outPath, quote=FALSE, 
                sep="\t", row.names=FALSE)},
    error=function(cond){
      NULL
    })
}
