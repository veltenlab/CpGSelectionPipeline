suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
pars <- ArgumentParser()
pars$add_argument("-o","--output",type="character",
    help="The output file")
pars$add_argument("-c","--config",type="character",
    help="Path to the configuration file")
args <- pars$parse_args()
config <- read_yaml(args$config)
cut.seq <- DNAString(config[['general']][['cut_seq']])
gen.version <- config[['general']][['genome']]
if('non_cut'%in%names(config)){
    config <- config[['non_cut']]
}else{
    stop("Missing required key 'non-cut' in the configuration file")
}
max.amps <- as.numeric(config[['n_amplicons']])
aci1.site <- DNAString('CCGC')
allowed.cpg.count <- 1:as.numeric(config[['max_cpgs']])
ampli.mean <- as.numeric(config[['length']])
cgis <- as.numeric(config[['n_cgis']])
proms <- as.numeric(config[['n_promoters']])
genes <- as.numeric(config[['n_genes']])
none <- as.numeric(config[['none']])
if(gen.version=='mm10'){
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    genome <- BSgenome.Mmusculus.UCSC.mm10
    all.chroms <- paste0("chr",1:19)
    all.types <- c(rep("cpgislands",cgis),
        rep("promoters",proms),
        rep("genes",genes),
        rep("none",none))
}else if(geno.version=='hg19'){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    suppressPackageStartupMessages(library(RnBeads.hg19))
    genome <- BSgenome.Hsapiens.UCSC.hg19
    rnb.load.annotation.from.db('ensembleRegBuildBPall',assembly=gen.version)
    ctcf <- as.numeric(config[['n_ctcf']])
    proximal <- as.numeric(config[['n_proximal']])
    tss <- as.numeric(config[['n_tss']])
    all.chroms <- paste0("chr",1:22)
    all.types <- c(rep("cpgislands",cgis),
        rep("promoters",proms),
        rep("genes",genes),
        rep("ensembleRegBuildBPctcf",ctcf),
        rep("ensembleRegBuildBPproximal",tss),
        rep("ensembleRegBuildBPtss",proximal),
        rep("none",none))
}else{
    stop("Unsupported genome, currently only ' hg19'and 'mm10'")
}
cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG",assembly=gen.version)))
if('db_snp'%in%names(config)){
    db.snp <- fread(config[['db_snp']])
    db.snp <- makeGRangesFromDataFrame(db.snp,start.field='POS',end.field='POS',seqnames.field='#CHROM')
    seqlevelsStyle(db.snp) <- 'UCSC'
}else{
    db.snp <- GRanges()
}

n.amps <- 0
ampli.pos <- c()
while(n.amps < max.amps){
    chr <- sample(all.chroms,1)
    ampli.size <- ampli.mean    
    #ampli.size <- round(rnorm(1,ampli.mean,5))
    seq <- genome[[chr]]
    type <- sample(all.types,1)
    if(type=="none"){
        pos <- sample(1:length(seq),1)
        if((pos+ampli.size) > length(seq)){
            next
        }
        sel.seq <- seq[pos:(pos+ampli.size)]
        sel.seq.hha <- seq[(pos-50):(pos+ampli.size+50)]
        freq <- alphabetFrequency(sel.seq)
        gc.cont <- (freq['C']+freq['G'])/length(sel.seq)
        if(gc.cont < 0.45 || gc.cont > 0.65){
            next
        }
        ampli.code <- paste0(chr,":",pos,"-",pos+ampli.size)
        if(length(matchPattern(cut.seq,sel.seq.hha))>0){
            next
        }else if(gen.version=='hg19' && length(findOverlaps(GRanges(ampli.code),rnb.get.annotation('ensembleRegBuildBPall',assembly=gen.version)))>0){
            next
        }else if(gen.version=='hg19' && length(findOverlaps(GRanges(ampli.code),db.snp))>0){
            next
        }else{
            cpg.count <- findOverlaps(GRanges(ampli.code),cpgs,ignore.strand=T)
            cpg.count <- length(cpg.count)/2
            if(!(cpg.count%in%allowed.cpg.count)){
                next
            }
            ampli.pos <- rbind(ampli.pos,c(ampli.code,chr,pos,pos+ampli.size,cpg.count,
                length(matchPattern(aci1.site,sel.seq))==0,type,gc.cont))
        }
    }else{
        if(!(type %in% rnb.region.types(assembly=gen.version))){
            rnb.load.annotation.from.db(type,assembly=gen.version)
        }
        an.data <- rnb.get.annotation(type,assembly=gen.version)[[chr]]
        sel.pos <- sample(an.data,1)
        sel.pos <- resize(sel.pos,ampli.size+1)
        sel.seq <- seq[start(sel.pos):end(sel.pos)]
        sel.seq.hha <- seq[(start(sel.pos)-50):(end(sel.pos)+50)]
        freq <- alphabetFrequency(sel.seq)
        gc.cont <- (freq['C']+freq['G'])/length(sel.seq)
        if(gc.cont < config[['min_gc']] || gc.cont > config[['max_gc']]){
            next
        }
        if(length(matchPattern(cut.seq,sel.seq.hha))>0){
            next
        }else if(length(findOverlaps(sel.pos,db.snp))>0){
            next
        }else{
            ampli.code <- paste0(chr,":",start(sel.pos),"-",end(sel.pos))
            cpg.count <- findOverlaps(GRanges(ampli.code),cpgs,ignore.strand=T)
            cpg.count <- length(cpg.count)/2
            if(!(cpg.count%in%allowed.cpg.count)){
                next
            }
            ampli.pos <- rbind(ampli.pos,c(ampli.code,chr,start(sel.pos),end(sel.pos),cpg.count,
                length(matchPattern(aci1.site,sel.seq))==0,type,gc.cont))
        }
    }
    all.types <- all.types[-(which(all.types %in% type)[1])]
    n.amps <- n.amps+1
}

ampli.pos <- as.data.frame(ampli.pos)
colnames(ampli.pos) <- c("Location","Chromosome","Start","End","CpGCount","AciSite","Type","GCContent")
write.csv(ampli.pos,file.path(args$output))

