library(corrplot)
library(reshape)
library(ggplot2)
library(scales)

jsplot_f <- snakemake@output[['jsplot']]
freqs_f <- snakemake@output[['freqplot']]
counts_f <- snakemake@output[['countplot']]
n_samples <- as.numeric(snakemake@params[['nsamples']])

image <- snakemake@input[['Rimage']]
load(image)

seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==n_samples])
data_common <- data[data$sequence %in% keepbc,]
get_fraction <- function(sample, counts) {
  res <- counts[counts$sample == sample,]
  tot <- sum(res$count)
  res$fraction <- res$count / tot
  return(res)
}
ggplot(data_common, aes(x=count)) +
  geom_histogram(bins=50)+facet_wrap(~sample, ncol=3)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+xlab('n.reads')+ylab('n.barcodes')

ggsave(counts_f)
freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]

ggplot(freqs, aes(x=fraction)) + 
  geom_histogram(bins=50)+facet_wrap(~sample, ncol=3)+
  theme_bw()+xlab('fraction tot reads')+ylab('n.barcodes')+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

ggsave(freqs_f)


cast <- cast(fr, id ~ sample, value="fraction")

js <- function(p,q) {
  n <- 0.5 * (p + q)
  JS <- 0.5 * (sum(p * log(p / n)) + sum(q * log(q / n)))
  JS
}
dc <- as.data.frame(cast)
rownames(dc) <- dc$id
dc$id <- NULL
allc <- colnames(dc)
pairs <- expand.grid(unique(allc),unique(allc))

jsd <- function(pair,data) {
  pair <- unlist(pair)
  js(data[, colnames(data)==pair[1]], data[,colnames(data)==pair[2]])
}

jsdiv <- apply(pairs, 1, jsd, dc)
jm <- matrix(jsdiv, nrow=length(allc))


colnames(jm) <- unique(allc)
rownames(jm) <- unique(allc)
pdf(jsplot_f)
corrplot(jm, is.corr=FALSE)
graphics.off()
#summary(jm[upper.tri(jm)])