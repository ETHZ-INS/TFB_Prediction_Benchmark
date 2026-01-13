Configurations (Methods, TFs) of the benchmark can be defined in the [config.yaml](config.yaml).

The workflow starts with downloading the ENCODE/GEO IDs of the preregistered datasets from the [preregistered protocol](https://osf.io/kv2tw/files/osfstorage) on OSF.
Subsequently, the specified data is downloaded from ENCODE/GEO and commonly preprocessed for all methods.

![Snakemake workflow for benchmarking](rulegraph.png)


Currently methods are run based on local installs / conda environments. 
TODO:        
1. Separate docker containers for each method      
2. Push extra data required by methods to Zenodo to allow for downloads during Snakemake execution     
