##################################### select_always_meth_unmeth.R #####################################
#' This file computes sites that are always unmethylated/methylated in all the comparisons that
#' were conducted with RnBeads.
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

report <- 'RnBeads/rnb_report/'
output <- 'meth_unmeth/'
config_file <- '../config.yaml'
config <- yaml.load_file(config_file)


all.comparisons <- list.files(file.path(report,'differential_methylation_data'), full.names=TRUE, pattern = 'diffMethTable_site')
system(paste0('rm -rf ', output,'/always_meth.csv'))
system(paste0('rm -rf ', output,'/always_unmeth.csv'))
dmrs <- lapply(all.comparisons,function(comp){
    print(comp)
    diff.table <- as.data.frame(fread(comp))
    group.names <- colnames(diff.table)[grepl('sd.',colnames(diff.table))]
    group.names <- gsub('sd.','', group.names)
    print(paste('Processing diffMethTable for ', group.names[1], 'and', group.names[2]))
    sds.first <- diff.table[[paste0('sd.',group.names[1])]]
    sds.second <- diff.table[[paste0('sd.',group.names[2])]]
    se.sd.first <- sd(sds.first)#/sqrt(length(sds.first))
    se.sd.second <- sd(sds.second)#/sqrt(length(sds.second))
    low.sds <- (sds.first<(mean(sds.first)+2*se.sd.first)) & (sds.second<(mean(sds.second)+2*se.sd.second))
    diff.table <- diff.table[low.sds,]
    row.names(diff.table) <- paste0(diff.table$Chromosome, '_', diff.table$Start)
    diff.positive <- diff.table[diff.table[[paste0('mean.',group.names[1])]]>config[['dmcs']][['methylated']]&
                                diff.table[[paste0('mean.',group.names[2])]]>config[['dmcs']][['methylated']],]
    diff.positive$type <- 'always_meth'
    row.names(diff.positive) <- paste0(diff.positive$Chromosome, '_', diff.positive$Start)
    file.positive <- file.path(output,'always_meth.csv')
    if(file.exists(file.positive)){
        all.positive <- read.csv(file.positive,row.names = 1)
        if(nrow(all.positive)<4){
            stop('Please delete the always_meth file in the output directory')
        }
        both <- intersect(row.names(diff.positive),row.names(all.positive))
        all.positive <- data.frame(all.positive[both,],diff.positive[both,])
        write.csv(all.positive,file.positive)
    }else{
        write.csv(diff.positive,file.positive)
    }
    diff.negative <- diff.table[diff.table[[paste0('mean.',group.names[1])]]<config[['dmcs']][['unmethylated']]&
                                    diff.table[[paste0('mean.',group.names[2])]]<config[['dmcs']][['unmethylated']],]
    diff.negative$type <- 'always_unmeth'
    row.names(diff.negative) <- paste0(diff.negative$Chromosome, '_', diff.negative$Start)
    file.negative <- file.path(output,'always_unmeth.csv')
    if(file.exists(file.negative)){
        all.negative <- read.csv(file.negative, row.names=1)
        if(nrow(all.negative)<4){
            stop('Please delete the always_meth file in the output directory')
        }
        both <- intersect(row.names(diff.negative),row.names(all.negative))
        all.negative <- data.frame(all.negative[both,],diff.negative[both,])
        write.csv(all.negative,file.negative)
    }else{
        write.csv(diff.negative,file.negative)
    }
    return(list(diff.positive,diff.negative))
})

source('../checkForCutSite.R')

always_meth <- list.files(output,pattern='always_meth.csv', full.names=TRUE)
always_meth <- checkForCutSite(read.csv(always_meth), config=config_file, number=75, sort.col='random', decreasing = FALSE)
write.csv(always_meth, file.path(output, 'always_meth_filtered.csv'))
always_unmeth <- list.files(output,pattern='always_unmeth.csv', full.names=TRUE)
always_unmeth <- checkForCutSite(read.csv(always_unmeth), config=config_file, number=75, sort.col='random', decreasing = FALSE)
write.csv(always_unmeth, file.path(output, 'always_unmeth_filtered.csv'))

