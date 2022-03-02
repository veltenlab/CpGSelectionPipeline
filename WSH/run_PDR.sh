#$ -q long-sl7
#$ -N PDR_GSM1274425
#$ -e /users/lvelten/project/Methylome/analysis/selection_pipeline/WSH/GSM1274425/PDR_2/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/selection_pipeline/WSH/GSM1274425/PDR_2/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=92G #Amount of memory
#$ -l h_rt=720:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 6

source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/rnbeads/
Rscript /users/lvelten/project/Methylome/src/selection_pipeline/WSHScripts/scores/PDR/compute_PDR.R PDR_GSM1274425 /users/lvelten/project/Methylome/analysis/selection_pipeline/WSH/GSM1274425/PDR_2/ /users/lvelten/project/Methylome/data/external/Cabezas/raw/GSM1274425/merged/merged_sorted.bam /users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20211004_reduced/cluster_run/preprocessing_RnBSet/ 6