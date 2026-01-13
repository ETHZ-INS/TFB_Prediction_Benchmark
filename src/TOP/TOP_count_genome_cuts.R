library(TOP)
library(data.table)

outDir <- snakemake@params[["out_dir"]]
if(!dir.exists(outDir)) dir.create(outDir)
print(basename(snakemake@output[["atac_fwd_bw_file"]]))
outName <- unlist(tstrsplit(basename(snakemake@output[["atac_fwd_bw_file"]]), 
                            split=".", fixed=TRUE, keep=1))
count_genome_cuts(bam_file=snakemake@input[["atac_bam_file"]], 
                  chrom_size_file=snakemake@params[["chrom_sizes"]], 
                  data_type='ATAC',
                  shift_ATAC=TRUE,
                  shift_ATAC_bases = c(4L, -4L),
                  outdir=outDir,
                  outname=outName)
