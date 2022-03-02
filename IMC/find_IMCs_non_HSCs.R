##################################### find_IMS.R #####################################
#' This script is used to deteremine intermediately methylated sites from the RnBeads
#' set, focusing only on the HSCs

.libPaths(c('/users/mscherer/conda/envs/rnbeads/lib/R/library/', .libPaths()))
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

rnb_set_path <- '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20211004/cluster_run/preprocessing_RnBSet/'
out_folder <- '/users/lvelten/project/Methylome/analysis/selection_pipeline/IMS/'
all_dmrs <- c('/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_HSC.csv',
              '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_MPP.csv',
              '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_MPP1.csv',
              '/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_MPP2.csv')
config_file <- '/users/lvelten/project/Methylome/src/selection_pipeline/config.yaml'
config <- yaml.load_file(config_file)

rnb_set <- load.rnb.set(rnb_set_path)
anno_fr <- annotation(rnb_set)
anno_gr <- makeGRangesFromDataFrame(anno_fr)
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
dmrs_gr <- resize(dmrs_gr, width = 500, fix = 'center')
op <- findOverlaps(anno_gr, dmrs_gr)
anno_gr <- anno_gr[-queryHits(op)]
meth_data <- meth(rnb_set)[, c("GSM1274427",
                               "GSM1274428",
                               "GSM1274429",
                               "GSM1274430",
                               "GSM1274431",
                               "GSM1274432",
                               "GSM1274433",
                               "GSM1274434",
                               "GSM1274435")]
meth_data <- meth_data[-queryHits(op), ]
is_intermediate <- apply(meth_data, 1, function(x){
  all(x>0.25&x<0.75)
})
meth_data <- meth_data[is_intermediate, ]
anno_gr <- anno_gr[is_intermediate]
covg_data <- covg(rnb_set)[, c("GSM1274427",
                               "GSM1274428",
                               "GSM1274429",
                               "GSM1274430",
                               "GSM1274431",
                               "GSM1274432",
                               "GSM1274433",
                               "GSM1274434",
                               "GSM1274435")]
covg_data <- covg_data[-queryHits(op), ][is_intermediate, ]
mean_covg <- rowMeans(covg_data)
too_high <- mean_covg>quantile(mean_covg, .95)
mean_covg <- mean_covg[!too_high]
meth_data <- meth_data[!too_high, ]
anno_gr <- anno_gr[!too_high]
meth_data_fr <- data.frame(Chromosome=seqnames(anno_gr),
                    Start=start(anno_gr),
                    End=end(anno_gr),
                    MeanMeth=rowMeans(meth_data),
                    MeanCovg=mean_covg)
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
write.csv(res, paste0(out_folder, '/IMS_annotated_non_HSCs.csv'))
