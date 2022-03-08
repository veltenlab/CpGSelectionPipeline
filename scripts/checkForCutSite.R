#' checkForCutSite
#'
#' This function checks for the presence of a cutsite in an amplicon using the configuration provided in the config file
#' @param input The input data frame
#' @param number The number of sites/amplicons to return
#' @param config The parsed configuration
#' @param sort.col The sorting column for the input table
#' @param decreasing Decreasing order?
#' @param use.extended Should the cutsite be padded by two Cs to check for no additional CpGs
#' @return A modified version of \code{input} with only those lines containing a cutsite for the specified enzyme
checkForCutSite <- function(dat,
                            number,
                            config='../config.yaml',
                            sort.col=c('mean.diff', 'mean.diff.1', 'mean.diff.2'),
                            decreasing=TRUE,
                            use.extended=FALSE){
  require(yaml)
  require(data.table)
  require(RnBeads)
  require(RnBeads.mm10)
  
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
  }else if(gen.version=='mm10'){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    genome <- BSgenome.Hsapiens.UCSC.hg19
  }else{
    stop("Genome version '", gen.version, "' not supported")
  }
  if(sort.col=='random'){
    dat <- dat[sample(1:nrow(dat), nrow(dat)),]
  }else{
    dat <- dat[order(rowMeans(dat[, sort.col, drop=FALSE]), decreasing=decreasing),]
  }
  i <- 1
  num <- 0
  vec.all <- c()
  dat$regionStart <- NA
  dat$regionEnd <- NA
  dat$cutsiteInRegion <- NA
  dat$gene <- NA
  dat$promoter <- NA
  dat$enhancer_annotation <- NA
  dat$enhancer_annotation_gene_name <- NA
  dat$pu1_motif <- NA
  dat$pu1_chip <- NA
  dat$Hox <- NA
  dat$GCContent <- NA
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
  while(num<number && i<nrow(dat)){
    chr <- as.character(dat[i,"Chromosome"])
    start <- as.numeric(as.character(dat[i,"Start"]))
    seq <- genome[[chr]]
    region.start <- start-2
    region.end <- start+2
    sel.seq <- seq[region.start:region.end]
    res <- matchPattern(cut.seq,sel.seq)
    extended.region.start <- region.start-200
    extended.region.end <- region.end+200
    cpg_count <- length(findOverlaps(GRanges(paste0(chr, ':',  extended.region.start, '-', extended.region.end)),
                                     cpgs))/2
    if(cpg_count<max.cpgs){
      extended.region.seq <- seq[extended.region.start:extended.region.end]
      res.extended <- matchPattern(cut.seq,extended.region.seq)
      if(length(res)==1 & length(res.extended)==1){
        if(use.extended){
          extended.cut1 <- matchPattern(DNAString(paste0('C', config[['general']][['cut_seq']])), sel.seq)
          extended.cut2 <- matchPattern(DNAString(paste0(config[['general']][['cut_seq']], 'G')), sel.seq)
          if(length(extended.cut1)>0 | length(extended.cut2)>0){
            i <- i+1
            next
          }
        }
        dat$regionStart[i] <- region.start
        dat$regionEnd[i] <- region.end
        dat$cutsiteInRegion[i] <- paste0(start(res), '-', end(res))
        freq <- alphabetFrequency(extended.region.seq)
        dat$GCContent[i] <- (freq['C']+freq['G'])/length(extended.region.seq)
        region.gr <- GRanges(paste0(chr,':',region.start,'-',region.end))
        op <- findOverlaps(region.gr,genes)
        closest_variable_gene <- distance(region.gr, genes[variable_genes])
        dat$ClosestVariableGene[i] <- paste(values(genes[variable_genes])$symbol[which.min(closest_variable_gene)],collapse=',')
        dat$ClosestVariableGeneDistance[i] <- min(closest_variable_gene, na.rm = TRUE)
        if(length(op)>0){
          dat$gene[i] <- paste(values(genes)$symbol[subjectHits(op)],collapse=',')
        }else{
          dat$gene[i] <- NA
        }
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
        tf.region.start <- region.start-10
        tf.region.end <- region.end+10
        tf.region.seq <- seq[tf.region.start:tf.region.end]
        pattern.match <- matchPattern(pu1.motif, tf.region.seq, max.mismatch = 1)
        if(length(pattern.match)>0){
          dat[i, 'pu1_motif'] <- 'PU1'
        }
        op <- findOverlaps(region.gr, pu1.chip)
        if(length(op)>0){
          dat[i, 'pu1_chip'] <- 'PU1'
        }
        for(j in 1:nrow(hema.motifs)){
          tf.name <- hema.motifs[j , 'TF']
          pattern.match <- matchPattern(hema.motifs[j, 'Motif'], tf.region.seq, max.mismatch = 1)
          if(length(pattern.match)>0){
            dat[i, tf.name] <- tf.name
          }
        }
        vec.all <- c(vec.all,i)
        num <- num+1
        print(num)
      }
    }
    i <- i+1
  }
  dat <- dat[vec.all,]
  return(dat)
}
