#!/bin/bash
#SBATCH -t 24:00:00 --mem=150G --nodes=1 --cpus-per-task=32 --mail-type=ALL --mail-user=rvalerischacin@cvm.tamu.edu --account=132769842952


module load Anaconda3
source activate mashtree_myco

module load GCC/10.3.0
module load Mash/2.3
module load BioPerl/1.7.8


mashtree_bootstrap.pl --numcpus 32 --reps 100 /scratch/user/robertv/allgenomic_fna/*.fna -- --min-depth 0 > mycobovis_mashtree.dnd
