#!/bin/bash
#SBATCH -J TP_BFA_cluster  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-05:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=30000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/TP_BFA_cluster_%A_%a.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/TP_BFA_cluster_%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc01
source activate milo_simple_conda5

python BT_combine_and_cluster.py ../accessory_files/TP_BFA_assay_file.csv ../../BT_Bioinformatic_Work/TP_output/TP_BFA/ ../../BT_Bioinformatic_Work/TP_output/TP_BC_Assoc/TP_CS_BC_Assoc_filtered_clusters.csv "${SLURM_ARRAY_TASK_ID}"
