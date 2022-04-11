##################################### find_IMS.R #####################################
#' This script is used to deteremine intermediately methylated sites from the RnBeads
#' set, focusing only on the HSCs

library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

rnb_set_path <- 'RnBeads/rnb_report_20211004_reduced/cluster_run/preprocessing_RnBSet/'
out_folder <- 'IMC/'
all_dmrs <- c('DMRs/high_filtered_HSC.csv',
              'DMRs/high_filtered_preB.csv',
              'DMRs/high_filtered_naiveB.csv',
              'DMRs/high_filtered_gcB.csv',
              'DMRs/high_filtered_memB.csv')
pdrs <- c('WSH/HSC_1/PDR/PDR_HSC_1',
          'WSH/HSC_2/PDR/PDR_HSC_2')
pdr_annotations <- c('WSH/HSC_1/PDR/annotation.RData',
                     'WSH/HSC_2/PDR/annotation.RData')
config_file <- '../config.yaml'
config <- yaml.load_file(config_file)
load(pdr_annotations[1])
for(i in 2:length(pdr_annotations)){
  last_anno <- annotation
  load(pdr_annotations[i])
  op <- findOverlaps(last_anno, annotation)
  last_anno <- last_anno[queryHits(op)]
}
all_pdr <- matrix(nrow=length(last_anno), ncol=length(pdrs))
for(i in 1:length(pdr_annotations)){
  load(pdr_annotations[i])
  op <- findOverlaps(last_anno, annotation)
  pdr <- read.csv(pdrs[i])
  all_pdr[queryHits(op), i] <- pdr[subjectHits(op), ]
}
annotation <- last_anno

source('../scripts/checkForCutSite.R')

rnb_set <- load.rnb.set(rnb_set_path)
anno_fr <- annotation(rnb_set)
anno_gr <- makeGRangesFromDataFrame(anno_fr)
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
dmrs_gr <- resize(dmrs_gr, width = 500, fix = 'center')
op <- findOverlaps(anno_gr, dmrs_gr)
anno_gr <- anno_gr[-queryHits(op)]
meth_data <- meth(rnb_set)[, c("HSC_1",
                               "HSC_2")]
meth_data <- meth_data[-queryHits(op), ]
is_intermediate <- apply(meth_data, 1, function(x){
  all(x>0.25&x<0.75)
})
meth_data <- meth_data[is_intermediate, ]
anno_gr <- anno_gr[is_intermediate]
covg_data <- covg(rnb_set)[, c("HSC_1",
                               "HSC_2")]
covg_data <- covg_data[-queryHits(op), ][is_intermediate, ]
mean_covg <- rowMeans(covg_data)
too_high <- mean_covg>quantile(mean_covg, .95)
mean_covg <- mean_covg[!too_high]
meth_data <- meth_data[!too_high, ]
anno_gr <- anno_gr[!too_high]
op <- findOverlaps(anno_gr, annotation)
pdrs <- rowMeans(all_pdr)
pdrs <- pdrs[subjectHits(op)]
meth_data_fr <- data.frame(Chromosome=seqnames(anno_gr),
                    Start=start(anno_gr),
                    End=end(anno_gr),
                    MeanMeth=rowMeans(meth_data),
                    MeanCovg=mean_covg)
meth_data_fr[queryHits(op), 'PDR'] <- pdrs
res <- checkForCutSite(na.omit(meth_data_fr),
                       number=750,
                       config=config_file, 
                       sort.col='PDR',
                       decreasing=FALSE)
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
write.csv(res, paste0(out_folder, '/IMC_annotated_all.csv'))
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/IMC_annotated.csv'))
