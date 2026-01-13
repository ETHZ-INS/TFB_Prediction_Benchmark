library(TOP)
library(data.table)

candSites <- as.data.frame(fread(snakemake@input[["cand_sites"]]))
idxstatsFile <- snakemake@input[["idxstats_file"]]

genomeCountDir <- dirname(snakemake@input[["atac_fwd_bw_file"]])
genomeCountName <- unlist(tstrsplit(basename(snakemake@input[["atac_fwd_bw_file"]]), 
                                    split=".", fixed=TRUE, keep=1))

countMatrix <- get_sites_counts(candSites,
                                genomecount_dir=genomeCountDir,
                                genomecount_name=genomeCountName,
                                tmpdir=genomeCountDir,
                                bwtool_path=snakemake@params[["bwtool_path"]])
binnedMat <- normalize_bin_transform_counts(countMatrix, 
                                            idxstats_file=idxstatsFile, 
                                            ref_size=5e7,
                                            transform='asinh')
outDir <- dirname(snakemake@output[["norm_count_matrix"]])
if(!dir.exists(outDir)) dir.create(outDir)
saveRDS(binnedMat, snakemake@output[["norm_count_matrix"]])