#!/bin/bash
#SBATCH -J BT_bc_parse_cluster_071918  #job name for array
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-10:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=60000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../shell_outputs/BT_bc_parse_cluster_071918.out      # File to which STDOUT will be written
#SBATCH -e ../../shell_outputs/BT_bc_parse_cluster_071918.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=milo.s.johnson.13@gmail.com  # Email to which notifications will be sent

module load Anaconda3/5.0.1-fasrc01
source activate milo_simple_conda5

python3 BT_BC_association_parse.py ../../raw_sequencing_data/BT1_BC_Assoc/ ../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_ ../accessory_files/BT_BC_Assoc_Demult.csv -genomic_bp_cutoff 30 -paired_end

python3 BT_BC_association_cluster.py ../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_counts.csv ../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_unfiltered_clusters.csv ../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_filtered_clusters.csv ../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_excluded_bcs.txt

sbatch --array=1-10 BT1_BFA_parsing_jabbaray.sh
