load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library(UpSetR)

#####
all_counts_wide <- dcast(data, sequence ~ sample, value.var="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL
write.table(all_counts_wide, file='barcodes_counts_wide.tsv.gz', sep="\t", quote=F)
barcode_appeareance <- lapply(unique(data$sample), function(x) {data[data$sample==x,'sequence']})
names(barcode_appeareance)  <- unique(data$sample)

seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==10])
data_common <- data[data$sequence %in% keepbc,]


get_fraction <- function(sample, counts) {
  res <- counts[counts$sample == sample,]
  tot <- sum(res$count)
  res$fraction <- res$count / tot
  return(res)
}

freqs <- lapply(unique(data$sample), get_fraction, data)
freqs <- do.call(rbind, freqs)



cast <- dcast(freqs, sequence ~ sample, value.var="fraction")
rownames(cast) <- cast$sequence
cast$sequence <- NULL

reps <- c('T0','NT','CTX')

averaged_freqs <- sapply(reps, function(x) { d <- cast[,grepl(x ,colnames(cast), fixed=T)]; apply(d, 1, function(x) {mean(x[!is.na(x)])})  })
averaged_freqs[rownames(averaged_freqs)=="AGCAGGCGAAGTTA-ACGTTGCAGTGTTGACGTCAACTGACTGCA",] #? not == to Simone what did he do?
averaged_freqs <- as.data.frame(averaged_freqs)
fcs <- data.frame(row.names=rownames(averaged_freqs), ctx_t0 = averaged_freqs$CTX/averaged_freqs$T0, nt_t0=averaged_freqs$NT/averaged_freqs$T0)
all <- names(cast[complete.cases(cast),])
repu <- fcs[!is.na(fcs$ctx_t0) & !is.na(fcs$nt_t0),]
head(repu[repu$nt_t0> 3 & repu$ctx_t0  < 3 & rownames(repu) %in% all,])
# remove low frequency TODO f0 < 0.005

#Shared selected barcodes (Log Ratio > 1.02)

###

fcgen <- function(replicate, model, counts, cond1, cond2) {
  nt <- counts[counts$sample == paste0(model, "_", replicate, cond1),]
  cetuxi <- counts[counts$sample == paste0(model, "_", replicate, cond2),]
  nt <- nt[order(nt$sequence),]
  cetuxi <- cetuxi[order(cetuxi$sequence),]
  if (! all(cetuxi$sequence == cetuxi$NT)) {
    return(NA)
  }
  res <- log(cetuxi$fraction/ nt$fraction)
  names(res) <- cetuxi$id
  return(res)
}


fcs_t0_ctx <- sapply(c('A','B','C'), fcgen, 'CRC0322', freqs, "_T0", "_CTX")
colnames(fcs_t0_ctx) <- c('repA', 'repB', 'repC')

m2 <- melt(fcs_t0_ctx)
m2$class <- ifelse(m2$X1 %in% rev$id, 'atleast','other')

ggplot(data=m2, aes(value, y=..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)


fcs_t0_nt <- sapply(c('A','B','C'), fcgen, 'CRC0322', freqs, "_T0", "_NT")
colnames(fcs_t0_nt) <- c('repA', 'repB', 'repC')

m3 <- melt(fcs_t0_nt)
m3$class <- ifelse(m3$X1 %in% rev$id, 'atleast','other')

ggplot(data=m3, aes(value, ..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)

#### pseudocounts
ccast <- dcast(data, sequence ~ sample, value.var="count", fill=1)
rownames(ccast) <- ccast$sequence
ccast$sequence <- NULL

reps <- c('T0','NT','CTX','overall')
freqs_cast <- apply(ccast, 2, function(x) x/sum(x))
averaged_freqs <- sapply(reps, function(x) { d <- freqs_cast[,grepl(x ,colnames(freqs_cast), fixed=T), drop=FALSE]; apply(d, 1, function(x) {mean(x)})  })
averaged_freqs[rownames(averaged_freqs)=="AGCAGGCGAAGTTA-ACGTTGCAGTGTTGACGTCAACTGACTGCA",] #? not == to Simone what did he do?
averaged_freqs <- as.data.frame(averaged_freqs)
nodrift <- averaged_freqs[averaged_freqs$T0 > 0.0005,]

fcs <- data.frame(row.names=rownames(nodrift), ctx_t0 = nodrift$CTX/nodrift$T0, nt_t0=nodrift$NT/nodrift$T0, t0_overall=nodrift$T0/nodrift$overall)
m <- melt(fcs)
ggplot(data=m, aes(x=log(value), fill=variable))+geom_density(alpha=0.7)+theme_bw()


fcs <- data.frame(row.names=rownames(averaged_freqs), ctx_t0 = averaged_freqs$CTX/averaged_freqs$T0, nt_t0=averaged_freqs$NT/averaged_freqs$T0, t0_overall=averaged_freqs$T0/averaged_freqs$overall )
m <- melt(fcs)
ggplot(data=m, aes(x=log(value), fill=variable))+geom_density(alpha=0.7)+theme_bw()


fcs_all <- data.frame(row.names=rownames(freqs_cast), ctx_A = freqs_cast[,1]/freqs_cast[,3],nt_A = freqs_cast[,2]/freqs_cast[,3],
                      ctx_B = freqs_cast[,4]/freqs_cast[,6],nt_B = freqs_cast[,5]/freqs_cast[,6],
                      ctx_C = freqs_cast[,7]/freqs_cast[,9],nt_C = freqs_cast[,8]/freqs_cast[,9])
lfcs_all <- log(fcs_all)

lfcs_all[rownames(lfcs_all) %in% rownames(nodrift) & lfcs_all$ctx_A > 1.02 & lfcs_all$ctx_B > 1.02 & lfcs_all$ctx_C > 1.02,]

#
wanted <- c('AGCAGGCGAAGTTA-ACGTTGCAGTGTTGACGTCAACTGACTGCA','AGTTTCCTGCGTGT-GTGTACACACACACGTCATGACGTGTACGT','AGTTTCCTGCGTGT-GTTGCACATGTGACTGCACAGTCAGTACCA','ATGCCAGAACATAT-CACACAACGTACACGTCATGCACATGTGCA','ATGCCAGAACATAT-CAGTTGTGCAACACGTTGCAACTGTGCATG')

## pseudocounts like Simone and plots
data$count <- data$count + 1
ccast <- dcast(data, sequence ~ sample, value.var="count", fill=0)
rownames(ccast) <- ccast$sequence
ccast$sequence <- NULL
freqs_cast <- apply(ccast, 2, function(x) x/sum(x))
long_fr <- melt(freqs_cast)
colnames(long_fr) <- c('seq','sample','freq')
long_fr$sample <- as.character(long_fr$sample)
long_fr[grepl('overall', long_fr$sample),"sample"] <- "CRC0322_overall_overall"
long_fr$treat <- sapply(long_fr$sample, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
long_fr$rep <- sapply(long_fr$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
long_fr$treat <- factor(long_fr$treat, levels=c('overall','T0','NT', "CTX"))
ggplot(data=long_fr, aes(x=log(freq), fill=rep))+geom_density(alpha=0.5)+facet_wrap(~treat)+theme_bw()
ggplot(data=long_fr, aes(x=log(freq), fill=treat))+geom_density(alpha=0.5)+theme_bw()+facet_wrap(~rep)
seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==10])
head(keepbc)
long_fr_common <- long_fr[long_fr$seq %in% keepbc,]
ggplot(data=long_fr_common, aes(x=log(freq), fill=treat))+geom_density(alpha=0.5)+theme_bw()+facet_wrap(~rep)
ggplot(data=long_fr_common, aes(x=log(freq), fill=rep))+geom_density(alpha=0.5)+facet_wrap(~treat)+theme_bw()


