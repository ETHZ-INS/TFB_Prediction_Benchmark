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
  qValue <- genomicRangesMapping(mergedRanges, bed, 
                                             scoreCol="qValue", 
                                             byCols="type",
                                             aggregationFun=mean)
  center <- genomicRangesMapping(mergedRanges, bed, 
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

# from TFBlearner: 
# copied here to save all dependencies for the singularity/docker image
genomicRangesMapping <- function(refRanges,
                                 assayTable,
                                 byCols=c("tf_name",
                                          "cellular_context"),
                                 scoreCol=NULL,
                                 aggregationFun=NULL,
                                 minoverlap=1,
                                 type=c("any", "start", "end",
                                        "within", "equal"),
                                 shift=FALSE,
                                 chunk=NULL,
                                 BPPARAM=BiocParallel::SerialParam()){
  
  type <- match.arg(type, choices=c("any","start", "end",
                                    "within", "equal"))
  
  # TODO: - add warning for integer overflows - data.table size
  assayTable <- .processData(assayTable, readAll=TRUE, shift=shift,
                             seqLevelStyle=seqlevelsStyle(refRanges))
  colsDepth <- unique(assayTable[[byCols[1]]])
  
  if(is.null(chunk)){
    chunk <- fifelse(nrow(assayTable)>1e7, TRUE, FALSE)
  }
  
  if(chunk){
    chrLevelsRef <- levels(seqnames(refRanges))
    assayTable <- subset(assayTable, chr %in% chrLevelsRef)
    assayTable[,chr:=factor(chr, levels=chrLevelsRef)]
    assayTable[[byCols[1]]] <- factor(assayTable[[byCols[1]]])
    if(length(byCols)>1){
      assayTable[[byCols[2]]] <- factor(assayTable[[byCols[2]]])}
    
    assayTable <- split(assayTable, by="chr")
    refRangesList <- S4Vectors::split(refRanges, seqnames(refRanges))
    
    assayTable <- assayTable[names(refRangesList)]
    overlapTable <- mapply(genomicRangesMapping,
                           refRangesList, assayTable,
                           MoreArgs=list(byCols=byCols, scoreCol=scoreCol,
                                         aggregationFun=aggregationFun,
                                         chunk=FALSE, shift=shift,
                                         type=type,
                                         BPPARAM=BPPARAM),
                           SIMPLIFY=FALSE)
    
    # retrieve original order
    refRangesList <- Reduce("c", refRangesList[-1], refRangesList[[1]])
    ind <- GenomicRanges::findOverlaps(refRangesList, refRanges,
                                       select="first", type="equal")
    
    rbindFill  <- function(mat1, mat2){
      
      if(is.null(mat1)) mat1 <- Matrix::Matrix(nrow=0, ncol=0)
      if(is.null(mat2)) mat2 <- Matrix::Matrix(nrow=0, ncol=0)
      
      allCols <- union(colnames(mat1), colnames(mat2))
      diffCols1 <- setdiff(allCols, colnames(mat1))
      diffCols2 <- setdiff(allCols, colnames(mat2))
      
      # get missing columns
      mat1Missing <- Matrix::Matrix(0, nrow=nrow(mat1), ncol=length(diffCols1),
                                    dimnames=list(NULL, diffCols1))
      mat1 <- cbind(mat1, mat1Missing)
      
      mat2Missing <- Matrix::Matrix(0, nrow=nrow(mat2), ncol=length(diffCols2),
                                    dimnames=list(NULL, diffCols2))
      mat2 <- cbind(mat2, mat2Missing)
      
      mat1 <- mat1[,allCols, drop=FALSE]
      mat2 <- mat2[,allCols, drop=FALSE]
      rbind(mat1, mat2)
    }
    
    if(length(byCols)>1){
      overlapTable <- lapply(colsDepth, function(col){
        tablesChr <- lapply(overlapTable,
                            function(tables) tables[[col]])
        tablesChr <- Reduce("rbindFill", tablesChr[-1],
                            tablesChr[[1]])
        tablesChr <- tablesChr[ind,,drop=FALSE]
        tablesChr})
      names(overlapTable) <- colsDepth
    }
    else{
      overlapTable <- Reduce("rbindFill", overlapTable[-1], overlapTable[[1]])
      overlapTable <- overlapTable[ind,,drop=FALSE]
    }
    
    return(overlapTable)
  }
  
  seqNamesCol <- "chr"
  startCol <- "start"
  endCol <- "end"
  
  if(sum(!(byCols %in% colnames(assayTable)))>0 | is.null(byCols)){
    stop("byCols needed to be provided and column names of assayTable.")
  }
  
  # attribute generic names to dimensionalities
  if(length(byCols)==2)
  {
    setnames(assayTable, byCols, c("col_depth", "col_width"))
    byCols <- c("col_depth", "col_width")
    multiTf <- TRUE
  }
  else
  {
    setnames(assayTable, byCols, c("col_width"))
    byCols <- c("col_width")
    multiTf <- FALSE
  }
  
  # get dimensions of tables
  nRefs <- length(refRanges)
  
  if(is.factor(assayTable$col_width)){
    colsWidth <- levels(assayTable$col_width)
  }
  else{
    colsWidth <- unique(assayTable$col_width)
  }
  nColsWidth <- length(colsWidth)
  # convert to integer for speed-up
  assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=colsWidth))]
  
  if(is.factor(assayTable$col_depth)){
    colsDepth <- levels(assayTable$col_depth)
  }
  else{
    colsDepth <- unique(assayTable$col_depth)
  }
  
  # convert to GRanges for faster overlap finding
  if("width" %in% colnames(assayTable)) assayTable$width <- NULL
  if("strand" %in% colnames(assayTable)) assayTable$strand <- NULL
  
  assayRanges <- .dtToGr(assayTable, seqCol="chr", addMetaCols=TRUE)
  
  # find overlaps with ref. coordinates
  overlapTable <- as.data.table(GenomicRanges::findOverlaps(refRanges,
                                                            assayRanges,
                                                            type=type,
                                                            minoverlap=minoverlap,
                                                            ignore.strand=TRUE))
  rm(refRanges, assayRanges)
  
  # retrieve tf and cell type ids
  overlapTable <- cbind(overlapTable$queryHits,
                        assayTable[overlapTable$subjectHits,
                                   c(byCols, scoreCol),
                                   with=FALSE])
  rm(assayTable)
  
  threads <- floor(getDTthreads())/BPPARAM$workers
  
  if(multiTf)
  {
    setkey(overlapTable, V1, col_width)
    if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- split(overlapTable, by=c("col_depth"))
    
    overlapTable <- BiocParallel::bplapply(overlapTable, function(table,
                                                                  scoreCol,
                                                                  aggregationFun,
                                                                  nRefs,
                                                                  colsWidth,
                                                                  threads){
      
      data.table::setDTthreads(threads)
      
      if(is.null(scoreCol) | is.null(aggregationFun)){
        table <- table[,.(value=.N),
                       by=c("V1", "col_width")]}
      else{
        table <- table[,.(value=aggregationFun(scoreCol)),
                       by=c("V1", "col_width")]}
      
      nColsWidth <- length(colsWidth)
      
      # convert to sparse matrix
      table <- sparseMatrix(i=table$V1,
                            j=as.integer(table$col_width),
                            dims=c(nRefs, nColsWidth),
                            x=table$value)
      colnames(table) <- colsWidth
      return(table)}, scoreCol=scoreCol, aggregationFun=aggregationFun,
      nRefs=nRefs, colsWidth=colsWidth, threads=threads,
      BPPARAM=BPPARAM)
  }
  else
  {
    # setkeys for speed-up
    overlapTable[,V1:=as.integer(V1)]
    setkey(overlapTable, col_width, V1)
    
    # overlap with ref. coordinates
    if(is.null(scoreCol) | is.null(aggregationFun)){
      overlapTable <- overlapTable[,.(scoreCol=.N),
                                   by=c("col_width", "V1")]}
    else{
      setnames(overlapTable, scoreCol, "scoreCol")
      overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
                                   by=c("col_width", "V1")]}
    
    # convert to sparse matrix
    overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1,
                                         j=overlapTable$col_width,
                                         dims=c(nRefs, nColsWidth),
                                         x=overlapTable$scoreCol)
    overlapTable <- Matrix::Matrix(overlapTable)
    
    colnames(overlapTable) <- colsWidth
  }
  
  # add combinations with zero overlaps
  missingDepthCols <- setdiff(colsDepth, names(overlapTable))
  if(length(missingDepthCols)>0){
    missingMat <- Matrix(0,nrow=nRefs, ncol=nColsWidth, doDiag=FALSE)
    colnames(missingMat) <- colsWidth
    missingTables <- replicate(length(missingDepthCols),
                               missingMat)
    names(missingTables) <- missingDepthCols
    overlapTable <- c(overlapTable, missingTables)
  }
  
  gc()
  return(overlapTable)
}

.processData <- function(data, readAll=FALSE, shift=FALSE,
                         subSample=NULL, seqLevelStyle="UCSC"){
  if(is.character(data)){
    if(grepl(".bam", basename(data), fixed=TRUE))
    {
      param <- Rsamtools::ScanBamParam(what=c('pos', 'qwidth', 'isize'))
      readPairs <- GenomicAlignments::readGAlignmentPairs(data, param=param)
      
      # get fragment coordinates from read pairs
      seqDat <- GRanges(seqnames(readPairs@first),
                        IRanges(start=pmin(GenomicAlignments::start(readPairs@first),
                                           GenomicAlignments::start(readPairs@last)),
                                end=pmax(GenomicAlignments::end(readPairs@first),
                                         GenomicAlignments::end(readPairs@last))),
                        strand=GenomicAlignments::strand(readPairs))
      seqDat <- granges(seqDat, use.mcols=TRUE)
      seqDat <- as.data.table(seqDat)
      setnames(seqDat, c("seqnames"), c("chr"))
    }
    else if(grepl(".bed", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data, stringsAsFactors=TRUE)
      else{
        
        readBed <- function(data){
          tryCatch(
            {
              seqDat <- fread(data, select=c(1:3,6),
                              col.names=c("chr", "start", "end", "strand"),
                              stringsAsFactors=TRUE)
              return(seqDat)},
            error = function(cond){
              seqDat <- fread(data, select=c(1:3),
                              col.names=c("chr", "start", "end"),
                              stringsAsFactors=TRUE)
              return(seqDat)
            })}
        seqDat <- readBed(data)
      }
    }
    else if(grepl(".tsv", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data, stringsAsFactors=TRUE)
      else{
        seqDat <- fread(data, select=c(1:3),
                        col.names=c("chr", "start", "end"),
                        stringsAsFactors=TRUE)}
      if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    }
    else if(grepl(".rds", basename(data), fixed=TRUE)){
      seqDat <- as.data.table(readRDS(data))
      if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    }
  }
  else{
    seqDat <- as.data.table(data)
    if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    seqDat$chr <- factor(seqDat$chr)
  }
  
  if(!is.null(subSample) & is.numeric(subSample)){
    message("Subsampling file")
    subSample <- as.integer(subSample)
    seqDat <- seqDat[sample(1:nrow(seqDat), min(nrow(seqDat), subSample)),]
  }
  
  # Match seqlevelstyle to reference
  if((sum(grepl("chr", levels(seqDat$chr)))==0 & seqLevelStyle=="UCSC") |
     (sum(grepl("chr", levels(seqDat$chr)))>0 & seqLevelStyle=="NCBI")){
    tmpgr <- GRanges(levels(seqDat[["chr"]]),
                     IRanges(seq_along(levels(seqDat[["chr"]])), width=2L))
    seqlevelsStyle(tmpgr) <- seqLevelStyle
    levels(seqDat$chr) <- seqlevels(tmpgr)
  }
  
  # Insert ATAC shift
  if(shift){
    seqDat[, start:=start+4L]
    seqDat[, end:=end-4L]
  }
  else if(shift){
    warning("Did not shift as no column named strand was not found")
  }
  
  seqDat[, start:=as.integer(start)]
  seqDat[, end:=as.integer(end)]
  if("width" %in% colnames(seqDat)) seqDat$width <- NULL
  
  return(seqDat)
}

.dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end",
                    strandCol="strand", stranded=FALSE, addMetaCols=FALSE){
  dt <- copy(dt)
  setnames(dt, seqCol, "seqnames", skip_absent = TRUE)
  
  if(stranded) strand <- dt[[strandCol]] else strand <- NULL
  
  gr <- GRanges(seqnames=dt[["seqnames"]],
                strand=strand,
                ranges=IRanges(start=dt[[startCol]], end=dt[[endCol]]))
  
  if(startCol==endCol)
  {
    gr <- GPos(seqnames=dt[["seqnames"]],
               strand=strand,
               pos=dt[[startCol]])
  }
  
  if(addMetaCols){
    metaCols <- dt[,setdiff(colnames(dt),
                            c(seqCol, startCol, endCol, strandCol,
                              "seqnames", "chr")),with=FALSE]
    mcols(gr) <- metaCols
  }
  
  return(gr)
}