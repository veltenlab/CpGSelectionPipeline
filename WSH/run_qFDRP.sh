#$ -q long-sl7
#$ -N qFDRP_GSM1274425
#$ -e WSH/GSM1274425/stderr.txt #The file that stderr will get written to
#$ -o WSH/GSM1274425/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=92G #Amount of memory
#$ -l h_rt=720:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 10

conda activate selection_pipeline
Rscript ../WSHScripts/scores/qFDRP/compute_qFDRP.R qFDRP_GSM1274425 WSH/GSM1274425/ GSM1274425/GSM1274425.bam RnBeads/rnb_report/cluster_run/preprocessing_RnBSet/ 10
