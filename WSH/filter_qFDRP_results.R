##################################### annotate_WSH.R #####################################
#' With this script, we filter for potentially interesting variable regions and annotate them
#' according to further epigenomic annotations such as TFBSs or ChIP-seq data

#.libPaths(c('/users/mscherer/R/R-4.1.2/', .libPaths()))
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

sample_names <- c('HSC_1',
                  'HSC_2')
res_folder <- 'WSH/'
out_folder <- 'WSH/'
all_dmrs <- c('RnBeads/DMRs/high_HSCs.csv',
              'RnBeads/DMRs/high_preB.csv',
              'RnBeads/DMRs/high_naiveB.csv', 
              'RnBeads/DMRs/high_gcB.csv',
              'RnBeads/DMRs/high_memB.csv')
imcs <- c('IMS/IMS_annotated_all.csv')
pdrs <- c('WSH/HSC_1/PDR/PDR_HSC_1.csv',
          'WSH/HSC_2/PDR/PDR_HSC_2.csv')
pdr_annotations <- c('WSH/HSC_1/PDR/annotation.RData',
                     'WSH/HSC_2/PDR/annotation.RData')
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
config_file <- '..//config.yaml'
config <- yaml.load_file(config_file)

source('../scripts/checkForCutSite.R')

qfdrps <- lapply(sample_names, function(x){
  read.csv(paste0(res_folder, x, '/qFDRP_' , x, '.csv'))$x
 })
qfdrps <- data.frame(do.call(cbind, qfdrps))
colnames(qfdrps) <- sample_names
qfdrp <- rowMeans(qfdrps)
load(paste0(paste0(res_folder, sample_names[1], '/'), 'annotation.RData'))
nas <- is.na(qfdrp)
qfdrp <- qfdrp[!nas]
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
                       number=300,
                       config=config_file, 
                       sort.col=c('qFDRP', 'PDR'),
                       use.extended = TRUE)
out_folder <- file.path(out_folder, paste0('filtered_qFDRP'))
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/filtered_qFDRP.csv'))
