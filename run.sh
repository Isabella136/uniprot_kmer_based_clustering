#!/bin/bash

#SBATCH --job-name=uniprot_clustering
#SBATCH --output=log/uniprot_clustering.out%j
#SBATCH --error=log/uniprot_clustering.err%j
#SBATCH --time=12:00:00
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=512gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

cargo run --release -- uniprot_arg.fasta 32