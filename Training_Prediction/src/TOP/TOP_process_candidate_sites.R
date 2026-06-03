log_file <- snakemake@log[[1]]
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)


library(TOP)
library(data.table)

outPath <- snakemake@output[["cand_sites"]]
comps <- tstrsplit(snakemake@output[["cand_sites"]], split="_")
tf <- gsub(".tsv", "", comps[[length(comps)]])

fimoFilePath <- file.path(snakemake@params[["motif_match_dir"]], 
                          paste0(tf, "_fimo.tsv"))
blackListPath <- snakemake@params[["black_list_path"]]

blackListPath
file.exists(blackListPath)

sites <- process_candidate_sites(fimo_file=fimoFilePath,
                                 thresh_pValue=1e-04,
                                 chr_order=paste0("chr", 1:22),
                                 flank=100,
                                 blacklist_file=blackListPath)
write.table(sites, file=outPath, 
            sep="\t", quote=FALSE, row.names=FALSE)