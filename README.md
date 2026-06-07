# Transcription Factor Binding Benchmark 

The protocol for this benchmark and the datasets have been preregistered under: https://osf.io/kv2tw/overview.    
ATAC- & ChIP-seq datasets are downloaded based on the specified TFs & cellular contexts.

The benchmark can be run with the following command:

```bash
snakemake --singularity-args "--nv --bind $PWD:$PWD --pwd $PWD" --use-singularity --rerun-triggers mtime --cores 16 --resources gpu=2 mem_gb=120 --configfile config.yaml
```

**Note:**    
* Adapt cores, gpus, and memory according to available resources.   
* Files of several 100s of GBs will be downloaded.   
* Alternatively, it can be run with the specified config ([`config.yaml`](./config.yaml)). This allows for subsetting methods and TF-cellular context combinations (untested!).   

The workflow uses centralized directories for downloading and processing the preregistered benchmark ATAC- & ChIP-seq datasets. All sub-workflows of the benchmarked methods start based on these commonly preprocessed files.

## Required Files 

### 1. Annotation Files

All Snakemake sub-workflows expect annotation files to exist centrally in the `data/annotation/` directory relative to top level directory. 
Further following files need to be available in `data/annotation/` before running the full pipeline:

- Reference genome assembly: **`Homo_sapiens.GRCh38.dna.primary_assembly.autosomes.chr.fa`** (and `.fai`).    
  *Downloaded from: [Ensembl (HTTP)](http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) and subsetted to autosomes.*
- Bowtie2 index: **`BOWTIE2Index_chr/genome`**     
  *Generated via `bowtie2-build` from the primary assembly above)*
- Chromosome sizes: **`reChr.sizes`**   
- Minimal blacklisted regions set: **`blacklist.bed`**  
  *Downloaded from: https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz*
- Extended blacklisted regions set:  **`blacklist_extended.bed`**    
  *Downloaded from: https://github.com/MiraldiLab/maxATAC_data/blob/main/hg38/hg38_maxatac_blacklist.bed. 
   These were used for the final evaluation of all methods, by removing all overlapping regions from computing performance metrics.*
- TruSeq adapters: **`truseqPE.fa`**      
  *Downloaded from: https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq2-PE.fa*

### 2. Singularity Images

- All singularity images of the benchmarked methods and the preprocessing are deposited under: [https://doi.org/10.5281/zenodo.20582127](https://zenodo.org/records/20582127)

Alternatively, the singularity images can also be generated from the dockerfiles specified in [`./Dockerfiles/`](./Dockerfiles/). 


## Benchmark Workflow Graph

![Benchmark Rulegraph](rulegraph.png)
