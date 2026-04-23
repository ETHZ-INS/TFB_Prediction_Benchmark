library(data.table)
library(gam)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
chIPPeaksPath <- args[1]
gamPath <- args[2]
dhsPath <- args[3]
outDir <- args[4]
nonconsChIPPeaksPath <- args[5]

modRep <- readRDS(gamPath)
refCoords <- fread(dhsPath, col.names=c("chr", "start", "end", "width", "strand"))
refCoords <- makeGRangesFromDataFrame(as.data.frame(refCoords))
labelDt <- fread(chIPPeaksPath, col.names=c("chr", "start", "end", "qValue"))

fileName <- basename(chIPPeaksPath)
strs <- unlist(tstrsplit(fileName, split="_", keep=5:6))
tf <- strs[[1]]
context <- gsub(".bed", "", strs[[2]])
labelDt$rep_prob <- predict(modRep, 
                            newdata=data.frame(qValue=labelDt$qValue),
                            type="response")
labelDt[,mid:=floor((end-start)/2)+start]
labelMidDt <- labelDt[,c("chr", "mid"), with=FALSE]
labelMidDt[,start:=mid-15]
labelMidDt[,end:=mid+15]
labelDt$overlaps_mid <- overlapsAny(makeGRangesFromDataFrame(as.data.frame(labelMidDt)),
                                    refCoords)
labelDt$overlaps_flank <- overlapsAny(makeGRangesFromDataFrame(as.data.frame(labelDt)),
                                      refCoords)
labelDt[,is_uncertain:=fifelse(overlaps_flank & !overlaps_mid, TRUE, FALSE)]

# label non-conservative peaks as uncertain
nonConsPeaks <- fread(nonconsChIPPeaksPath, 
                      col.names=c("chr", "start", "end"), 
                      select=1:3)
nonConsPeaks <- makeGRangesFromDataFrame(as.data.frame(nonConsPeaks))
refCoordsNoPeaks <- subsetByOverlaps(refCoords,
                                     makeGRangesFromDataFrame(as.data.frame(labelDt)),
                                     invert=TRUE)
refCoordsPeaks <- subsetByOverlaps(refCoords,
                                   makeGRangesFromDataFrame(as.data.frame(labelDt)))

# get ambigous peaks
uncertainPeaks <- subsetByOverlaps(nonConsPeaks, 
                                   refCoordsNoPeaks)
# remove those that overlap of the define training peaks
uncertainPeaks <- subsetByOverlaps(uncertainPeaks, 
                                   refCoordsPeaks, invert=TRUE)
uncertainPeaksDt <- as.data.table(uncertainPeaks)
setnames(uncertainPeaksDt, "seqnames", "chr")
uncertainPeaksDt$is_uncertain <- TRUE
# does not really matter as its marked as uncertain and removed from training
uncertainPeaksDt$rep_prob <- 1 
labelDt <- rbind(labelDt, uncertainPeaksDt, use.names=TRUE, fill=TRUE)

write.table(labelDt, file.path(outDir, 
                               paste0(paste("labels", tf, 
                                            context, sep="_"), ".tsv")), 
            row.names=FALSE, quote=FALSE, sep="\t")