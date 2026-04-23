library(TOP)

chipBamFiles <- snakemake@input[["bam_files"]]
candSites <- snakemake@input[["cand_sites"]]
# could be provided as a seperate input
chipIdxStatFiles <- gsub(".bam", ".idxstats.txt",snakemake@input[["bam_files"]])

chromSizesPath <- snakemake@params[["chrom_sizes_path"]]
sitesChip <- count_normalize_chip(candSites,
                                  chip_bam_files=chipBamFiles,
                                  chip_idxstats_files=chipIdxStatFiles,
                                  chrom_size_file=chromSizesPath)

outDir <- dirname(snakemake@output[["norm_count_matrix"]])
if(!dir.exists(outDir)) dir.create(outDir)
saveRDS(sitesChip, snakemake@output[["norm_count_matrix"]])