library(httr)
library(jsonlite)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
idSetPath <- snakemake@input[["downloaded_chip_peak_ids"]]
outPath <- snakemake@output[["downloaded_atac_ids"]]

idSets <- fread(idSetPath, header=TRUE)
idSets <- split(idSets, by="set")

# get all ATAC-seq dataset identifiers on ENCODE
allExpDt <- fread("https://www.encodeproject.org/report.tsv?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=ATAC-seq&files.run_type=paired-ended",
                  skip=1)
allExpDt <- subset(allExpDt, Organism=="Homo sapiens")

metaAtacAll <- data.table()
for(idDt in idSets){
  set <- unique(idDt$set)
  cellularContexts <- unique(idDt$CellularContext)
  
  for(cellularContext in cellularContexts){
    if(!(cellularContext %in% c("H1", "HeLa-S3", "Jurkat"))){
      expIdDt <- subset(allExpDt, `Biosample term name`==cellularContext)
      for(id in unique(expIdDt$Accession)){
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
                                      file_type=="bam" & 
                                      status=="released" & 
                                      output_type=="alignments" &
                                      assembly=="GRCh38")
  
        subFilesDt <- subFilesDt[,c("accession", "output_type", 
                                    "file_type", "href"), with=FALSE]
      subFilesDt$CellularContext <- cellularContext
      subFilesDt$Exp_Accession <- id
      subFilesDt$set <- set
      
      metaAtacAll <- rbind(metaAtacAll, subFilesDt, fill=TRUE, use.names=TRUE)
     }
    }
    else{
      subFilesDt <- data.table()
      subFilesDt$CellularContext <- cellularContext
      if(cellularContext=="H1"){
        subFilesDt$Exp_Accession <- "GSM2584741"}
      else if(cellularContext=="Jurkat"){
        subFilesDt$Exp_Accession <- "GSM3693103"
      }
      else if(cellularContext=="HeLa-S3"){
        subFilesDt$Exp_Accession <- "GSM4764093"
      }
      
      subFilesDt$set <- set
      metaAtacAll <- rbind(metaAtacAll, subFilesDt, fill=TRUE, use.names=TRUE)
    }
  }
}

metaAtacAll <- metaAtacAll[, c("accession",
                               "output_type",
                               "file_type",
                               "href",	
                               "CellularContext",	
                               "Exp_Accession",	
                               "set"), with=FALSE]
write.table(metaAtacAll, file=outPath, quote=FALSE, 
            sep="\t", row.names=FALSE)
