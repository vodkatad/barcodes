library(corrplot)
library(reshape)
library(ggplot2)
library(scales)

jsplot_f <- snakemake@output[['jsplot']]

image <- snakemake@input[['Rimage']]
load(image)

all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL

all_counts_wide <- all_counts_wide + 1

freqs_cast <- apply(all_counts_wide, 2, function(x) x/sum(x))

#cast <- cast(fr, id ~ sample, value="fraction")

js <- function(p,q) {
  n <- 0.5 * (p + q)
  JS <- 0.5 * (sum(p * log(p / n)) + sum(q * log(q / n)))
  JS
}
dc <- as.data.frame(freqs_cast)

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