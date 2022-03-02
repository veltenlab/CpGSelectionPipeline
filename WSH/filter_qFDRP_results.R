##################################### annotate_WSH.R #####################################
#' With this script, we filter for potentially interesting variable regions and annotate them
#' according to further epigenomic annotations such as TFBSs or ChIP-seq data

library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

res_folder <- 'WSH/GSM1274424/'
sample_name <- 'qFDRP_GSM1274424'
out_folder <- 'WSH/'
all_dmrs <- c('DMRs/non_pairwise/high_HSCs.csv',
              'DMRs/non_pairwise/high_MPP.csv',
              'DMRs/non_pairwise/high_MPP1.csv',
              'DMRs/non_pairwise/high_MPP2.csv')
imcs <- c('IMC/IMC_annotated_all.csv')
pdrs <- c('WSH/GSM1274424/PDR/PDR_GSM1274424.csv',
          'WSH/GSM1274425/PDR/PDR_GSM1274425.csv',
          'WSH/GSM1274426/PDR/PDR_GSM1274426.csv')
pdr_annotations <- c('WSH/GSM1274424/PDR/annotation.RData',
                     'WSH/GSM1274425/PDR/annotation.RData',
                     'WSH/GSM1274426/PDR/annotation.RData')
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
pdr_annotation <- last_anno
# Include high PDR here
config_file <- '../config.yaml'
config <- yaml.load_file(config_file)

source('../scripts/checkForCutSite.R')

qfdrp <- read.csv(paste0(res_folder, sample_name, '.csv'))
load(paste0(res_folder, 'annotation.RData'))
nas <- is.na(qfdrp$x)
qfdrp <- qfdrp[!nas, ]
annotation <- annotation[!nas]
is_1 <- qfdrp==1
annotation <- annotation[!is_1]
qfdrp <- qfdrp[!is_1]
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
dmrs_gr <- resize(dmrs_gr, width = 500, fix = 'center')
op <- findOverlaps(annotation,dmrs_gr)
qfdrp <- qfdrp[-queryHits(op)]
annotation <- annotation[-queryHits(op)]
imcs <- read.csv(imcs)
imcs_gr <- makeGRangesFromDataFrame(imcs, end='Start')
imcs_gr <- resize(imcs_gr, width = 500, fix = 'center')
op <- findOverlaps(annotation, imcs_gr)
qfdrp <- qfdrp[-queryHits(op)]
annotation <- annotation[-queryHits(op)]
op <- findOverlaps(annotation, pdr_annotation)
qfdrp <- qfdrp[queryHits(op)]
pdrs <- rowMeans(all_pdr[subjectHits(op), ])
qfdrp <- data.frame(Chromosome=seqnames(annotation),
                    Start=start(annotation),
                    End=end(annotation),
                    qFDRP=qfdrp,
                    PDR=pdrs)
res <- checkForCutSite(na.omit(qfdrp),
                       number=230,
                       config=config_file, 
                       sort.col=c('qFDRP', 'PDR'))
out_folder <- file.path(out_folder, paste0('filtered_', sample_name))
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/filtered_', sample_name, '.csv'))
