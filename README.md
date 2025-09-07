# inversionSimulation
### This is a snakemake pipeline for simulating recurrent inversion events using structured coalescent process and performing the same evolutionary approach proposed in Porubsky et al. (2022) to detect and infer recurrent inversion loci.

## Required software
#### Python (v3.10)
#### Snakemake (v7.2)
#### R


## Input data
#### All required input data are listed in the config.yaml file.

## Usage
```
snakemake -j 10 -w 60 -kp
```
or using the run_slurm.sh to submit jobs over a SLURM system.
