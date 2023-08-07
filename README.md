# rlooper-plasmids

Basically just a refined version of [this repo](https://github.com/EthanHolleman/rlooperVariableTopology). Main changes are 

- Runs on new HPC cluster; Franklin
- Uses Snakemake parameter exploration functions to allow for easily setting large number of run parameters (Do this by changing
  values in `resources/rlooper_params.tsv`
- Output filepaths include parameters used in each run
- Additional python scripts to process `wig` bp predictions into tsv files with the parameters used in each run in long format
  for easy plotting with `ggplot2`
- Runs on any sequence (in theory) with a correct R-loop fasta header that is located in the `resources/sequences` directory

