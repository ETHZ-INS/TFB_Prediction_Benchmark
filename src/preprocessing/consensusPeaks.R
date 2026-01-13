.mergeConservative <- function(bedFilePaths, minOverlap=1){
  
  beds <- lapply(bedFilePaths, function(bedFile){
    bed <- fread(bedFile, select=c(1:3,9,10), 
                 col.names = c("chr", "start", "end", "qValue", "dist_center"))
    bed[,center:=start+dist_center]
  })
  names(beds) <- names(bedFilePaths)
  bed <- rbindlist(beds, idcol="experiment")
  bed$peak_id <- 1:nrow(bed)
  bed$type <- "union"
  bed[,center:=as.numeric(center)]

  bedRanges <- makeGRangesFromDataFrame(as.data.frame(bed), keep.extra.columns=TRUE)
  mergedRanges <- reduce(bedRanges)
  qValue <- TFBlearner::genomicRangesMapping(mergedRanges, bed, 
                                             scoreCol="qValue", 
                                             byCols="type",
                                             aggregationFun=mean)
  center <- TFBlearner::genomicRangesMapping(mergedRanges, bed, 
                                             scoreCol="center", 
                                             byCols="type",
                                             aggregationFun=median)
  
  mergedRanges$qValue <- qValue[,1,drop=TRUE]
  mergedRanges$center <- as.integer(center[,1,drop=TRUE])
  mergedRanges
}

pasteOrder <- function(x1, x2){
  x1 <- unlist(x1)
  x2 <- unlist(x2)
  if(grepl("_", x1) & grepl("_", x2)){
    x1 <- unlist(tstrsplit(x1,split="_")) 
    x2 <- unlist(tstrsplit(x2, split="_"))
  }
  x <- c(x1, x2)
  x <- as.integer(x)
  label <- paste(x[order(x)], collapse="_")
  return(label)
}

.mergeBoundaries <- function(bed, minOverlap, rounds){
  
  if(rounds>0)
  {
    bed2 <- copy(bed)
    setkey(bed2, chr, start, end)
    if(!("experiment" %in% colnames(bed))){
        bed$experiment <- 1
        bed2$experiment <- 2
    }
    
    ovBed <- foverlaps(bed, bed2, 
                       by.x=c("chr", "start", "end"),
                       by.y=c("chr", "start", "end"), 
                       minoverlap=minOverlap, 
                       nomatch=NULL)
    
    ovBed <- subset(ovBed, peak_id!=i.peak_id & 
                           experiment!=i.experiment & 
                           peak_label!=i.peak_label)
    if(nrow(ovBed)>0){
  
    ovBed[,ov_dist:=pmax(end-i.start,
                         i.end-start)]#add a column ov_dist to ovBed 

    ovBed <- subset(ovBed, ov_dist>=minOverlap) #get rid of overlaps that are smaller than the minOverlap, if min=1 then all overlaps are kept
    ovBed$ov_id <- 1:nrow(ovBed) #add overlap id
    ovBed[,peak_label_new:=pasteOrder(peak_label, i.peak_label), by=ov_id] #add a column peak label that is like (peak1_peak2) that overlap
    ovBed$peak_label <- NULL #get rid of old peak label column
    ovBed[,peak_label:=peak_label_new] #update old peak label with new one

    # merge the peak boundaries
    pairOvs <- ovBed[,.(chr=chr,
                        qValue=mean(c(qValue, i.qValue)),
                        start=median(c(start, i.start)),
                        end=median(c(end, i.end)),
                        center=median(c(center,i.center)),
                        peak_label=unique(peak_label)),
                     by=.(ov_id)] #merging boundaries

    pairOvs <- unique(pairOvs, by=c("chr", "start", "end")) #get rid of peaks with the same boundaries (duplicates)
    pairOvs <- pairOvs[,.(chr=chr,
                          start=median(start),
                          qValue=mean(qValue),
                          end=median(end),
                          center=median(center)),
                       by=.(peak_label)] #comes into play after 1st iteration and merges the rest of overlap combinations
    
    pairOvs$peak_id <- 1:nrow(pairOvs)
    
    return(.mergeBoundaries(pairOvs, minOverlap, rounds-1))
    }
    else{
      bed$is_rep <- FALSE
      bed
    }
  }
  else{
    bed$is_rep <- TRUE
    return(bed)
  }
}

.mergeOptimal <- function(bedFilePaths, minOverlap=1, mergeRounds=2){
  
  beds <- lapply(unique(names(bedFilePaths)), function(id){
    bed <- lapply(bedFilePaths[names(bedFilePaths)==id], fread, 
                  select=c(1:3,9,10), 
                  col.names=c("chr", "start", "end", "qValue", "dist_center"))
    bed <- rbindlist(bed)
    bed[,center:=start+dist_center]
    bed <- unique(bed)
    bed 
  })
  names(beds) <- unique(names(bedFilePaths))
  bed <- rbindlist(beds, idcol="experiment") 
  bed <- unique(bed)
  
  bed$peak_id <- 1:nrow(bed) 
  bed$peak_label <- 1:nrow(bed) 
  ovBed <- .mergeBoundaries(bed, minOverlap, mergeRounds)
  
  ovBed <- subset(ovBed, is_rep)
  ovBed <- makeGRangesFromDataFrame(as.data.frame(ovBed), keep.extra.columns=TRUE)
}

consensusPeaks <- function(consBedFilePaths, 
                           optBedFilePaths,
                           minOverlap=1){
  
  if(length(consBedFilePaths)>0){
    consMergedPeaks <- .mergeConservative(consBedFilePaths, 
                                          minOverlap=minOverlap)
    consMergedPeaks$type <- "conservative"}
  else{
    consMergedPeaks <- GRanges()
  }
  
  if(length(optBedFilePaths)>0){
    optMergedPeaks <- .mergeOptimal(optBedFilePaths, minOverlap=1, 
                                    mergeRounds=min(length(unique(names(optBedFilePaths))), 2))
    optMergedPeaks <- unique(optMergedPeaks)
    if(length(optMergedPeaks)>0){
      optMergedPeaks <- optMergedPeaks[!overlapsAny(optMergedPeaks, consMergedPeaks)]
      if(length(optMergedPeaks)>0){
        optMergedPeaks$type <- "merged_optimal"
      }
    }
    
  }
  else{
    optMergedPeaks <-  GRanges()
  }
  
  mergedPeaks <- c(consMergedPeaks, optMergedPeaks)
  mergedPeaks
}
