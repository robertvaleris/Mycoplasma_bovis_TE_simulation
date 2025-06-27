#!/bin/bash
#SBATCH -t 30:00:00 --mem=150G --nodes=1 --cpus-per-task=32 --mail-type=ALL --mail-user=rvalerischacin@cvm.tamu.edu --account=132769842952


module load GCC/11.2.0
module load OpenMPI/4.1.1
module load Kraken2/2.1.2


#add PSV Mbovis genomes to database
for file in /scratch/user/robertv/Relabeled_PSV_Mbovis_genomes/*ado	echo $file	kraken2-build --add-to-library $file --db /scratch/user/robertv/databases/kraken2/Mycoplasmoidales_Mbovis_dbdone# Add Mycoplasmoidales genomesfor file in /scratch/user/robertv/Labeled_Mycoplasmoidales_genbank_genomes/*ado	echo $file	kraken2-build --add-to-library $file --db /scratch/user/robertv/databases/kraken2/Mycoplasmoidales_Mbovis_dbdone# Build databasekraken2-build --build --db /scratch/user/robertv/databases/kraken2/Mycoplasmoidales_Mbovis_db
