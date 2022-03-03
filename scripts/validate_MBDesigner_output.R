validateMBDesigner <- function(csv,
                            config='../config.yaml'){
  
  require(yaml)
  require(data.table)
  require(RnBeads)
  require(RnBeads.mm10)
  
  suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
  genome <- BSgenome.Mmusculus.UCSC.mm10
  dat <- read.csv(csv)
  hox_cluster <- c('HoxA'=GRanges('chr6:52150000-52275000'),
                   'HoxB'=GRanges('chr11:96280000-96375000'),
                   'HoxC'=GRanges('chr15:102910000-103040000'),
                   'HoxD'=GRanges('chr2:74665000-74765000'))
  config <- yaml.load_file(config)
  cut.seq <- DNAString(config[['general']][['cut_seq']])
  max.cpgs <- as.numeric(config[['general']][['max_cpgs']])
  hema.motifs <- read.csv(config[['annotations']][['hema_tf_motifs']])
  pu1.motif <- DNAString(config[['general']][['pu1_motif']])
  pu1.chip <- read.table(config[['annotations']][['pu1_chip']])
  pu1.chip$V2 <- as.numeric(pu1.chip$V2)
  pu1.chip$V3 <- as.numeric(pu1.chip$V3)
  pu1.chip <- makeGRangesFromDataFrame(na.omit(pu1.chip), seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V6')
  gen.version <- config[['general']][['genome']]
  variable_genes <- as.character(read.csv(config[['annotations']][['variable_genes']])[, 2])
  if(gen.version=='mm10'){
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    genome <- BSgenome.Mmusculus.UCSC.mm10
  }else if(gen.version=='hg19'){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    genome <- BSgenome.Hsapiens.UCSC.hg19
  }else{
    stop("Genome version '", gen.version, "' not supported")
  }
  i <- 1
  num <- nrow(dat)
  dat$cutsiteInRegion <- NA
  dat$gene <- NA
  dat$promoter <- NA
  dat$enhancer_annotation <- NA
  dat$enhancer_annotation_gene_name <- NA
  dat$pu1_motif <- NA
  dat$pu1_chip <- NA
  dat$Hox <- NA
  dat$GCContent <- NA
  dat$CpGCount <- NA
  dat$HhaCutsites <- NA
  dat$ClosestVariableGene <- NA
  dat$ClosestVariableGeneDistance <- NA
  dat$AgingDMC <- NA
  dat[hema.motifs$TF] <- NA
  tfs.files.all <- list.files(config[['annotations']][['tf_chip']],full.names=TRUE)
  tfs.all <- gsub('.bed', '', list.files(config[['annotations']][['tf_chip']]))
  for(tf in tfs.all){
    dat[,tf] <- NA
  }
  cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG", assembly = gen.version)))
  genes <- unlist(rnb.get.annotation('genes',assembly = gen.version))
  names(genes) <- gsub('chr[[:alnum:]][[:punct:]]', '', names(genes))
  names(genes) <- gsub('chr[[:alnum:]][[:alnum:]][[:punct:]]', '', names(genes))
  promoters <- unlist(rnb.get.annotation('promoters',assembly = gen.version))
  reg.elements <- read.csv(config[['annotations']][['enhancer_catalog']])
  row.names(reg.elements) <- paste0(reg.elements$Chr,'_',reg.elements$Start)
  reg.elements.tab <- read.table(config[['annotations']][['enhancer_catalog_bed']],sep='\t')
  reg.elements.gr <- makeGRangesFromDataFrame(reg.elements.tab,seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V4')                    
  values(reg.elements.gr) <- reg.elements.tab$V5
  aging_dmcs <- read.csv(config[['annotations']][['aging_dmcs']])
  aging_dmcs <- makeGRangesFromDataFrame(aging_dmcs, seqnames.field = 'Chromosome',
                                         start.field='Start',
                                         end.field='Start',
                                         keep.extra.columns=TRUE)
  while(i < num){
    chr <- as.character(dat[i,"chr"])
    seq <- genome[[chr]]
    start <- as.numeric(as.character(dat[i,"amplicon_start"]))
    end <- as.numeric(as.character(dat[i,"amplicon_end"]))
    cpg_count <- length(findOverlaps(GRanges(paste0(chr, ':',  start, '-', end)),
                                     cpgs))/2
    dat[i, 'CpGCount'] <- cpg_count
    region.seq <- seq[start:end]
    res.extended <- matchPattern(cut.seq, region.seq)
    dat[i, 'HhaCutsites'] <- length(res.extended)
    dat$cutsiteInRegion[i] <- paste0(start(res.extended), '-', end(res.extended))
    freq <- alphabetFrequency(region.seq)
    dat$GCContent[i] <- (freq['C']+freq['G'])/length(region.seq)
    region.gr <- GRanges(paste0(chr,':', start,'-', end))
    op <- findOverlaps(region.gr,genes)
    if(length(op)>0){
      dat$gene[i] <- paste(values(genes)$symbol[subjectHits(op)],collapse=',')
    }else{
      dat$gene[i] <- NA
    }
    closest_variable_gene <- distance(region.gr, genes[variable_genes])
    dat$ClosestVariableGene[i] <- paste(values(genes[variable_genes])$symbol[which.min(closest_variable_gene)],collapse=',')
    dat$ClosestVariableGeneDistance[i] <- min(closest_variable_gene, na.rm = TRUE)
    op <- findOverlaps(region.gr, aging_dmcs)
    if(length(op)>0){
      dat$AgingDMC[i] <- values(aging_dmcs)$mean.diff[subjectHits(op)]
    }
    op <- findOverlaps(region.gr,promoters)
    if(length(op)>0){
      dat$promoter[i] <- paste(values(promoters)$symbol[subjectHits(op)],collapse=',')
    }else{
      dat$promoter[i] <- NA
    }
    op <- which(unlist(sapply(hox_cluster, function(x,y){
      length(queryHits(findOverlaps(x,y)))>0
    }, y=region.gr)))
    if(length(op)>0){
      dat$Hox[i] <- names(hox_cluster)[op]
    }
    op <- findOverlaps(region.gr,reg.elements.gr)
    if(length(op)>0){
      reg.element <- unname(unlist(values(reg.elements.gr))[subjectHits(op)])
      dat$enhancer_annotation[i] <- paste(reg.elements[reg.element,'Annotation'],collapse=',')
      dat$enhancer_annotation_gene_name[i] <- paste(reg.elements[reg.element,'Gene.Name'],collapse=',')
    }else{
      dat$enhancer_annotation[i] <- NA
      dat$enhancer_annotation_gene_name[i] <- NA
    }
    for(j in 1:length(tfs.all)){
      tf.gr <- makeGRangesFromDataFrame(read.table(tfs.files.all[[j]],sep = '\t'),seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
      op <- findOverlaps(region.gr,tf.gr)
      if(length(op)>0){
        dat[i,tfs.all[j]] <- tfs.all[j]
      }else{
        dat[i,tfs.all[j]] <- NA
      }
    }
    pattern.match <- matchPattern(pu1.motif, region.seq, max.mismatch = 1)
    if(length(pattern.match)>0){
      dat[i, 'pu1_motif'] <- 'PU1'
    }
    op <- findOverlaps(region.gr, pu1.chip)
    if(length(op)>0){
      dat[i, 'pu1_chip'] <- 'PU1'
    }
    for(j in 1:nrow(hema.motifs)){
      tf.name <- hema.motifs[j , 'TF']
      pattern.match <- matchPattern(hema.motifs[j, 'Motif'], region.seq, max.mismatch = 1)
      if(length(pattern.match)>0){
        dat[i, tf.name] <- tf.name
      }
    }
    i <- i+1
  }
  return(dat)
}

mb_id <- '1234'
res <- validateMBDesigner(paste0('MB_output/Designer/', mb_id, '-design-summary.csv'))
input <- read.csv('MB_input/frame_cleartext.csv')
output.gr <- makeGRangesFromDataFrame(res, seqnames.field = 'chr', start.field = 'amplicon_start', end.field = 'amplicon_end')
min_dists <- c()
for(i in 1:length(output.gr)){
  min_dists <- c(min_dists, min(distance(output.gr[i], output.gr[-i]), na.rm = TRUE))
}
if(sum(min_dists<1000)>50){
  stop('Many closeby CpGs')
}
input.gr <- makeGRangesFromDataFrame(input)
op <- findOverlaps(input.gr, output.gr, ignore.strand=T)
if(any(!(1:nrow(res)%in%subjectHits(op)))){
  print('Something is wrong')
}
res <- data.frame(res[subjectHits(op), ], input[queryHits(op), ])
write.csv(res, paste0('MB_output/Designer/', mb_id, '-design-summary_annotated.csv'))

library(ggsci)
plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_blank(),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))
plot_base <- data.frame(Type=gsub('[[:digit:]]', '', res$Type),
                      AgingDMC=ifelse(is.na(res$AgingDMC), 'NoDMC', 'DMC'),
                      Enhancer=ifelse(is.na(res$enhancer_annotation), 'NoEnhancer', 'Enhancers'))
to_plot <- c()
for(ty in unique(plot_base$Type)){
  to_plot <- rbind(to_plot, data.frame(plyr::count(plot_base[plot_base$Type%in%ty, 'AgingDMC']), Type=ty))
}
plot <- ggplot(to_plot, aes(x="", y=freq, fill=x))+geom_bar(stat="identity", width=1, color="white")+plot_theme+scale_fill_tron()+ylab("")+xlab("")+facet_wrap(Type~., ncol=3)
ggsave('amplicons_aging.pdf',
       height=100,
       width=100,
       unit='mm')

to_plot <- c()
for(ty in unique(plot_base$Type)){
  to_plot <- rbind(to_plot, data.frame(plyr::count(plot_base[plot_base$Type%in%ty, 'Enhancer']), Type=ty))
}
plot <- ggplot(to_plot, aes(x="", y=freq, fill=x))+geom_bar(stat="identity", width=1, color="white")+plot_theme+scale_fill_tron()+ylab("")+xlab("")+facet_wrap(Type~., ncol=3)
ggsave('amplicons_enhancer.pdf',
       height=100,
       width=100,
       unit='mm')

dat <- res
types <- gsub('[0-9]', '', dat$Type)
types[grep('high', types)] <- 'DMC'
types[types%in%'IMR'] <- 'IMC'
to_plot <- plyr::count(types)
colnames(to_plot) <- c('Type','Count')
plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_blank(),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))
to_plot$Type <- factor(to_plot$Type, levels=c('Non_cut',
                                              'Always_methylated',
                                              'Always_unmethylated',
                                              'DMC',
                                              'IMC',
                                              'WSH'))
to_plot <- to_plot[order(to_plot$Type), ]
to_plot$prop <- to_plot$Count/sum(to_plot$Count) *100
to_plot$ypos <- cumsum(to_plot$prop)-0.5*to_plot$prop
plot <- ggplot(to_plot, aes(x="", y=Count, fill=Type))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+scale_fill_tron()+ylab("")+xlab("")+geom_text(aes(label = Count),
                                                                                     position = position_stack(vjust = 0.5))
ggsave('amplicon_type_distribution.png',
       height=150,
       width=150,
       unit='mm')

plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text.x = element_text(angle=45, hjust=1, color='black', size=8),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))

dat <- read.csv('4135-design-summary_annotated.csv')
types <- gsub('[0-9]', '', dat$Type)
types[grep('high', types)] <- 'DMC'
types[types%in%'IMR'] <- 'IMC'
dat$Type <- types
tf_names <- colnames(dat)
start_tf <- which(tf_names%in%'CpGCount')+1
end_tf <- which(tf_names%in%"Scl")
tf_names <- tf_names[start_tf:end_tf]
counts_tfs <- unlist(lapply(tf_names, function(x){
  sum(!is.na(dat[, x]))
}))
to_plot <- data.frame(Count=counts_tfs, TF=tf_names)
to_plot$TF <- factor(to_plot$TF, levels=to_plot$TF[order(to_plot$Count, decreasing=T)])
plot <- ggplot(to_plot, aes(x=TF, y=Count))+geom_bar(stat="identity", width=1, color="white")+
  plot_theme
ggsave('amplicon_TF_all_distribution.png',
       height=75,
       width=150,
       unit='mm')

to_plot <- c()
for(type in unique(dat$Type)){
  sel_dat <- subset(dat, subset=Type==type)
  counts_tfs <- unlist(lapply(tf_names, function(x){
    sum(!is.na(sel_dat[, x]))
  }))
  to_plot <- rbind(to_plot, data.frame(Type=type, TF=factor(tf_names, levels=tf_names[order(counts_tfs, decreasing=TRUE)]), Count=counts_tfs))
}
plot <- ggplot(to_plot, aes(x=TF, y=Count))+geom_bar(stat="identity", width=1, color="white")+
  plot_theme+facet_wrap(Type~., ncol=3)
ggsave('amplicon_TF_stratified.png',
       height=150,
       width=300,
       unit='mm')
