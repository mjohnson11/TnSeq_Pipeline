# Tnseq pipeline scripts
## Scripts:
1. BC association:
  * BT_BC_association_parse.py
  * BT_BC_association_cluster.py
2. BC reading:
  * BT_parse.py
  * BT_combine_and_cluster.py
3. Annotation
  * annotate.py
4. S Estimation
  * simple_s_estimation.py
  * joint_s_inference.py
  * TP_s_measure.py
  * BT_s_measure.py
  * GC_s_measure.py
5. Analysis
  * calculate_edge_stats.py
  * TP_analysis.py
  * BT_analysis.py


Steps necessary to recreate data:
### make directories (these go in the parent directory of the github project) 
bash making_directory_structure.sh

### download raw_sequencing_data into the parent directory of the github project from XXX

### this does the BC assocation for the first experiment and then submits sbatch --array=1-10 BT1_BFA_parsing_jabbaray.sh to get the BFA parsing started
sbatch BT_BC_assoc_runner.sh

### this does the BC assocation for the second experiment and then submits sbatch --array=1-20 TP_BFA_parsing_jabbaray.sh to get the BFA parsing started
sbatch TP_BC_assoc_runner.sh

### BFA barcode clustering and file combining
sbatch --array=1-2 BT1_BFA_clustering_jabbaray.sh
sbatch --array=1-4 TP_BFA_clustering_jabbaray.sh

### Mutation annotation
sbatch run_annotation.sh

### Measuring fitness effects
sbatch --array=0-4 run_s_measure_BT.sh
sbatch --array=0-35 run_s_measure_TP.sh

### Run analysis scripts
sbatch analysis_runner_TP.sh
sbatch analysis_runner_BT.sh
sbatch DFE_analysis_runner.sh
sbatch run_QTL_GO_analysis.sh

### GC experiment calls
sbatch --array=1-20 GC_BFA_parsing_jabbaray.sh
sbatch --array=1-2 GC_BFA_clustering_jabbaray.sh
sbatch run_s_measure_GC.sh

### Plotting figures and making tables is all done in the jupyter notebooks TNSEQ_PLOTTING.ipynb and GC_plotting.ipynb
