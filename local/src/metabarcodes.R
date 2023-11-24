load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/counts.Rdata')


all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL

all_counts_wide <- all_counts_wide + 1

#ccast <- cast(data, sequence ~ sample, value="count", fill=0)
#rownames(ccast) <- ccast$sequence
#ccast$sequence <- NULL
freqs_cast <- apply(all_counts_wide, 2, function(x) x/sum(x))


# select 500 group of barcodes at random
total_bc <- nrow(freqs_cast)
group <- floor(total_bc/500)
new_order <- sample(1:total_bc)
n_groups <- floor(total_bc/group)
new_freqs <- data.frame(matrix(rep(0, ncol(freqs_cast)*n_groups), ncol=ncol(freqs_cast)))
colnames(new_freqs) <- colnames(freqs_cast)
j <- 1
for (i in seq(1, n_groups)) {
  for (k in seq(1, group)) {
    new_freqs[i,] <- new_freqs[i,] + freqs_cast[new_order[j],]
    j <- j + 1
  }
}
rownames(new_freqs) <- paste0('mbc_', seq(1, nrow(new_freqs)))

long_fr <- melt(as.matrix(new_freqs))
colnames(long_fr) <- c('seq', 'sample','freq')
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
  subset_y1 <- subset[subset$treat == 'CTX','freq'] # CTX for CRC0322, cetuxi for CRC0327
  subset_y2 <- subset[subset$treat == 'NT','freq']
  subset_x <- subset[subset$treat=='T0','freq']
  ssubset_y1 <- subset[subset$treat == 'CTX','seq']
  ssubset_y2 <- subset[subset$treat == 'NT','seq']
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


# keep only > 0.001 and produce histograms
#thr <- 0.00025
thr <- 0.0005

ggplot(data=pdata_all, aes(x=x, y=y, color=fill))+geom_point()+facet_wrap(~rep)+
  scale_color_manual(values=c('red','grey'))+theme_bw()+xlab('Initial frequency')+
  ylab('Final Frequency')+ggtitle('meta')+geom_vline(xintercept=thr)


thr <- 0
pdata_all_1 <- NULL
pdata_all_2 <- NULL
for (i in seq(1, length(all_rep)-1)) {
  w <- all_rep[i]
  subset <- long_fr[long_fr$rep == w,]
  subset0 <- subset[(subset$treat =='T0' & subset$freq > thr),]
  subset1 <- subset[(subset$treat == "CTX" & subset$seq %in% subset0$seq),]
  subset2 <- subset[(subset$treat == "NT" & subset$seq %in% subset0$seq),]
  subset <- rbind(subset1, subset2, subset0)
  subset <- subset[order(subset$seq),]
  subset_y1 <- subset[subset$treat == 'CTX','freq']
  subset_y2 <- subset[subset$treat == 'NT','freq']
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