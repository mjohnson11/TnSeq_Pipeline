#!/bin/bash
#SBATCH -J annotate  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-00:30              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=3500               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/annotate.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/annotate.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc01
source activate milo_simple_conda5
module load bowtie2
python annotate.py
rm tmp.fasta
rm tmp_bowtie_output.tsv
