log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink(type = "output")
  sink(type = "error")
  close(log_con)
}, add = TRUE)

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
