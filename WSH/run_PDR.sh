#$ -q long-sl7
#$ -N PDR_HSC_1
#$ -e WSH/HSC_1/PDR/stderr.txt #The file that stderr will get written to
#$ -o WSH/HSC_1/PDR/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=92G #Amount of memory
#$ -l h_rt=720:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 12

conda activate selection_pipeline
Rscript ../WSHScripts/scores/PDR/compute_PDR.R PDR_HSC_1 WSH/HSC_1/PDR/ HSC_1/HSC_1.bam rnb_report/cluster_run/preprocessing_RnBSet/ 12
