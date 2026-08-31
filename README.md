# inversionSimulation
### This is a snakemake pipeline for simulating recurrent inversion events using structured coalescent process and performing the same evolutionary approach proposed in Porubsky et al. (2022) to detect and infer recurrent inversion loci.

## Required software
#### Python (v3.10)
#### Snakemake (v7.2)
#### R
#### msprime (v2.1.2)
#### IQ-TREE (v3.1.3)

## Input data
#### All required input data are listed in the config.yaml file.

## Usage
```
# Use one of the following models. Need to modify *model_manifest* in the config.yaml accordingly before running.

# For running single event simulations

snakemake -j 10 -w 60 -kp -s Snakefile.single.snake

# For running recurrent event simulations
snakemake -j 10 -w 60 -kp -s Snakefile.recurrent.snake
```
or using the run_slurm.sh to submit jobs over a SLURM system.
