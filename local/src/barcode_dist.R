#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
args <- commandArgs(trailingOnly = T)
input <- args[1]
outplot <- args[2]
outnum <- args[3]

save.image('pippo.Rdata')
samples <- unlist(strsplit(input, ','))
samples_names <- unlist(lapply(samples, function(x) {y <- basename(x)}))

data <- read.table(samples[1], sep="\t", header=FALSE)
names(data) <- c('id','sequence','count')
data$sample <- samples_names[1]

for (i in seq(2, length(samples))) {
    d <- read.table(samples[i], sep="\t", header=FALSE)
    names(d) <- c('id','sequence','count')
    d$sample <- samples_names[i]
    data <- rbind(data, d)
}

infobarcodes <- function(sample, data) {
  d <- data[data$sample == sample, ]
  return(c(length(unique(d$sequence)),c(length(d$sequence)), sum(d$count), nrow(d)))
}

numbers <- as.data.frame(t(sapply(samples, infobarcodes, data)))
colnames(numbers) <- c('a','n_seen_barcodes','reads_with_barcodes','a1')
numbers$a <- NULL
numbers$a1 <- NULL


ggplot(data, aes(x=count)) + geom_histogram(bins=50)+facet_wrap(~sample)+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+theme_bw()+xlab('barcodes counts')+ylab('')
ggsave(outplot)

write.table(numbers, file=outnum, sep="\t", quote=FALSE)


#d <- read.table(infile, sep="\t", header=FALSE, comment.char='')
#names(d) <- c('id','sequence','count')
#infobarcodes <- function(data) {
#  return(c(length(unique(d$sequence)),c(length(d$sequence)), sum(d$count))
#}

#numbers <- infobarcodes(d)
#colnames(numbers) <- c('n_barcodes','unique_n_seen_barcodes','total_read_barcodes')
#numbers$a <- NULL
#numbers$a1 <- NULL


#plotbarcodes <- function(sample, data) {
#  d <- data[data$sample == sample, ]
#  
#}

#ggplot(d, aes(x=count)) + geom_histogram(bins=50)+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+theme_bw()+xlab('barcodes counts')+ylab('')+ggtitle(sample)