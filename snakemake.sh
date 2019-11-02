#!/bin/bash
#SBATCH --mem=2G
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=4-0
#SBATCH --mail-type=END

mkdir "output"
mkdir "output/slurm_out"
source ~/.bashrc
conda activate snakemake
echo "Executing paired-end workflow"
snakemake -s Snakefile -p -j 9999 --latency-wait 120 --use-conda --cluster-config config/cluster.json --cluster "sbatch -c {cluster.cpu} --mem={cluster.mem} -t {cluster.time} -p {cluster.partition} -o {cluster.output} -e {cluster.error} --mail-type=FAIL"
