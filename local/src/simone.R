library(reshape)
library(ggplot2)

scatter_f <- snakemake@output[['scatter']]
histo_f <- snakemake@output[['histo']]
#counts_f <- snakemake@output[['countplot']]
#n_samples <- as.numeric(snakemake@params[['nsamples']])

image <- snakemake@input[['Rimage']]
load(image)


all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL

all_counts_wide <- all_counts_wide + 1

#ccast <- cast(data, sequence ~ sample, value="count", fill=0)
#rownames(ccast) <- ccast$sequence
#ccast$sequence <- NULL
freqs_cast <- apply(all_counts_wide, 2, function(x) x/sum(x))
long_fr <- melt(freqs_cast)
colnames(long_fr) <- c('seq','sample','freq')
long_fr$sample <- as.character(long_fr$sample)

#ref <- 'T0'
#set <- 'CTX1' 'NT1' set as wildcards TODO

long_fr$treat <- sapply(long_fr$sample, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
long_fr$rep <- sapply(long_fr$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})

# very bad code following
all_rep <- unique(long_fr$rep)
pdata_all <- NULL
for (i in seq(1, length(all_rep)-1)) {
  w <- all_rep[i]
  subset <- long_fr[long_fr$rep == w,]
  subset <- subset[order(subset$seq),]
  subset_y1 <- subset[subset$treat == 'CTX1','freq']
  subset_y2 <- subset[subset$treat == 'NT1','freq']
  subset_x <- subset[subset$treat=='T0','freq']
  ssubset_y1 <- subset[subset$treat == 'CTX1','seq']
  ssubset_y2 <- subset[subset$treat == 'NT1','seq']
  ssubset_x <- subset[subset$treat=='T0','seq']
  stopifnot(all(ssubset_y1 == ssubset_y2))
  stopifnot(all(ssubset_x == ssubset_y2))
  pdata <- data.frame(x=rep(subset_x, 2), y=c(subset_y1, subset_y2), fill=c(rep('CTX', length(subset_y1)),rep('NT', length(subset_y2))), rep=rep(w, length(subset_x)*2))
  if (i == 1) {
    pdata_all <- pdata
  } else {
    pdata_all <- rbind(pdata_all, pdata)
  }
}
ggplot(data=pdata_all, aes(x=x, y=y, color=fill))+geom_point()+facet_wrap(~rep)+scale_color_manual(values=c('red','grey'))+theme_bw()+xlab('f0')+ylab('ft')+ggtitle('w1 vs t0')
ggsave(scatter_f)


# keep only > 0.001 and produce histograms
#thr <- 0.00025
thr <- 0.0005

pdata_all_1 <- NULL
pdata_all_2 <- NULL
for (i in seq(1, length(all_rep)-1)) {
  w <- all_rep[i]
  subset <- long_fr[long_fr$rep == w,]
  subset0 <- subset[(subset$treat =='T0' & subset$freq > thr),]
  subset1 <- subset[(subset$treat == "CTX1" & subset$seq %in% subset0$seq),]
  subset2 <- subset[(subset$treat == "NT1" & subset$seq %in% subset0$seq),]
  subset <- rbind(subset1, subset2, subset0)
  subset <- subset[order(subset$seq),]
  subset_y1 <- subset[subset$treat == 'CTX1','freq']
  subset_y2 <- subset[subset$treat == 'NT1','freq']
  subset_x <- subset[subset$treat=='T0','freq']
  seqs <- as.character(subset[subset$treat=='T0','seq'])
  logfr_1 <- log2(subset_y1/subset_x)
  logfr_2 <- log2(subset_y2/subset_x)
  pdata1 <- data.frame(seq=seqs, logfr=logfr_1, fill=rep('CTX', length(logfr_1)),rep=rep(w, length(logfr_1)), stringsAsFactors = FALSE)
  pdata2 <- data.frame(seq=seqs, logfr=logfr_2, fill=rep('NT', length(logfr_1)), rep=rep(w, length(logfr_2)), stringsAsFactors = FALSE)
  if (i == 1) {
    pdata_all_1 <- pdata1
    pdata_all_2 <- pdata2
  } else {
    pdata_all_1 <- rbind(pdata_all_1, pdata1)
    pdata_all_2 <- rbind(pdata_all_2, pdata2)
  }
}

ggplot(data=rbind(pdata_all_1,pdata_all_2), aes(x=logfr, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=15, color='black', data=pdata_all_1)+
  geom_histogram(aes(y=..count..), bins=15, color='black', data=pdata_all_2)+
  facet_wrap(~rep)+
  theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()

ggsave(histo_f)

# then look for enriched barcodes in all replicates
# logFC > 1 in all replicates

selected_1 <- pdata_all_1[pdata_all_1$logfr > 1,]
selected_2 <- pdata_all_2[pdata_all_2$logfr > 1,]

sel <- as.data.frame(table(selected_1$seq))
sel[sel$Freq > 1,]
seq1 <- sel[sel$Freq > 1, 'Var1']

sel2 <- as.data.frame(table(selected_2$seq))
sel2[sel2$Freq > 1, ]
seq2 <- sel2[sel2$Freq > 1,'Var1']

intersect(seq1, seq2)

pp <- long_fr[long_fr$seq=="ATGCCAGAACATAT-CAACCAACGTCACATGCAGTTGTGGTACTG",]
pp$treat <- as.character(pp$treat)
pp[pp$rep=="overall",'treat'] <- 'overall0'
pp$time <- substr(pp$treat, nchar(pp$treat), nchar(pp$treat))

pp$treat <- factor(pp$treat, levels=c('overall', 'preT0','T0','NT0','CTX0','NT1','CTX1','CTX2'))
ggplot(data=pp, aes(y=freq, x=treat, color=rep, group=rep))+geom_line()+geom_point()+theme_bw()