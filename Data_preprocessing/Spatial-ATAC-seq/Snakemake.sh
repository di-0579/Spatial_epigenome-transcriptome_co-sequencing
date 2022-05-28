#!/bin/bash
#SBATCH --partition=general 
#SBATCH --job-name=Snakemake
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem=64g
#SBATCH --time=120:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=your email
#SBATCH -o Snakemake.%j.out
#SBATCH -e Snakemake.%j.err

SLURM_ARGS="-p {cluster.partition} -J {cluster.job-name} -n {cluster.ntasks} -c {cluster.cpus-per-task} \
--mem={cluster.mem} -t {cluster.time} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user} \
-o {cluster.output} -e {cluster.error}"

snakemake -j 20 --cluster-config cluster.json --cluster "sbatch $SLURM_ARGS"

