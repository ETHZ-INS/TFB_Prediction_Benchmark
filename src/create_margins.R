library(data.table)

margin <- as.integer(snakemake@params[["margin"]])
dhsCoords <- fread(snakemake@input[["dhs_coords"]])
dhsCoords[,start:=V2-margin]
dhsCoords[,end:=V3+margin]
fwrite(dhsCoords, snakemake@output[["dhs_margin_coords"]], 
       col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")


