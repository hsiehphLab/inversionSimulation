#!/usr/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hsiehph@umn.edu
#SBATCH -o msmc2snake.out
#SBATCH -e msmc2snake.err
#SBATCH -p sioux,msismall
set -euo pipefail

snakemake --use-conda --profile profile/ "$@"

