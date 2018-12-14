#!/bin/bash
#SBATCH -J GC_BFA_parse  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-05:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=20000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/GC_BFA_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/GC_BFA_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc01
source activate milo_simple_conda5

python BT_parse.py ../accessory_files/GC_Exp_Demult.csv ../../BT_Bioinformatic_Work/GC_BFA_output/ ../raw_sequencing_data/GC_BFA/ "${SLURM_ARRAY_TASK_ID}" GC_BFA
