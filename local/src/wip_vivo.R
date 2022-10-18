load('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/invivo_feb2022/counts.Rdata')

library(ggpubr)
library(ggplot2)
library(vegan)
library(reshape)

dd <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(dd) <- dd$sequence
dd$sequence <- NULL
Hall <- diversity(t(dd))

plotd <- data.frame(entropy=Hall, sample=as.character(names(Hall)), stringsAsFactors = FALSE)
plotd$condition <- sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][1]; return(y[1])})
ggplot(plotd, aes(x=condition, y=entropy, fill=condition))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c('white',"darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))

all_counts_wide <- dd
all_counts_wide <- all_counts_wide + 1

freqs_cast <- as.data.frame(apply(all_counts_wide, 2, function(x) x/sum(x)))



library(MASS)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#https://slowkow.com/notes/ggplot2-color-by-density/
compare <- function(x, y, log, nx, ny) {
  zx <- length(x[x==0.000000001])
  zy <- length(y[y==0.000000001])
  zz <- x==0.000000001 & y==0.000000001
  zz <- sum(zz)
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  print(ggplot(d, aes(x=x, y=y, color=density)) +geom_point()+theme_bw()+theme(text = element_text(size=20))+xlab(nx)+ylab(ny)+labs(caption=paste0(pe$estimate, ', pval=', pe$p.value))+scale_color_viridis_c())
  return(c(zx, zy, zz))
}

compare(freqs_cast$basale_1, freqs_cast$basale_2, TRUE, 'basale_1', 'basale_2')
compare(freqs_cast$basale_1, freqs_cast$vitro_1, TRUE, 'basale_1', 'vitro_1')
compare(freqs_cast$basale_1, freqs_cast$topo_1, TRUE, 'basale_1', 'topo_1')
compare(freqs_cast$vitro_1, freqs_cast$topo_1, TRUE, 'vitro_1', 'topo_1')
compare(freqs_cast$topo_3, freqs_cast$topo_1, TRUE, 'topo_2', 'topo_1')


##########3 baloon

# comment block to work on common only
seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==10])
data_common <- data[data$sequence %in% keepbc,]
get_fraction <- function(sample, counts) {
  res <- counts[counts$sample == sample,]
  tot <- sum(res$count)
  res$fraction <- res$count / tot
  return(res)
}

freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
head(freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]

cast <- cast(fr, id ~ sample, value="fraction")
rownames(cast) <- cast$id
cast$id <- NULL
### common only block

### all block 
all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL

all_counts_wide <- all_counts_wide + 1

cast <- as.data.frame(apply(all_counts_wide, 2, function(x) x/sum(x)))

###

ll <- apply(cast, 1, function(x) {any(x > 0.001)})

sel <- cast[ll,]

sel$id <- rownames(sel)
rev <- melt(sel)

#freqs$atleast <- ifelse(freqs$id %in% rev$id, 'atleast','other')
#ggplot(freqs, aes(x=reorder(id,-fraction),fraction,fill=atleast)) + geom_bar(stat="identity")+facet_wrap(~sample, ncol=3)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank, axis.ticks.x=element_blank())+theme_bw()
#ggplot(freqs, aes(x=reorder(id,-count),log10(count),fill=atleast)) + geom_bar(stat="identity")+ theme(axis.title.x=element_blank(), axis.text.x=element_blank, axis.ticks.x=element_blank())+facet_wrap(~sample, ncol=3)+theme_bw()


x <- 'variable'
y <- 'id'
size <- 'value'
fill <- 'id'
barcode <- 'barcode'
#rev$variable <- factor(rev$variable, levels=c('CRC0322_overall', 'CRC0322_A_T0','CRC0322_B_T0','CRC0322_C_T0','CRC0322_A_NT','CRC0322_B_NT', 'CRC0322_C_NT','CRC0322_A_CTX','CRC0322_B_CTX', 'CRC0322_C_CTX'))
rev$variable <- as.factor(rev$sample)

#rev$barcode <- as.factor(sample(nrow(rev)))
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                         axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                         axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# order on freq in basale 1

o <- rev[rev$variable=="basale_1",]
o <- o[order(o$value),]
rev$id <- factor(rev$id, levels=o$id)


ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### freqs with zero to study shared vivo/etc
all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL
freqs_cast <- as.data.frame(apply(all_counts_wide, 2, function(x) x/sum(x) ))

topi <- rownames(freqs_cast)[apply(freqs_cast[, grepl('topo', colnames(freqs_cast))], 1, function(x) {all(x!=0)} )]
notopi <- rownames(freqs_cast)[apply(freqs_cast[, !grepl('topo', colnames(freqs_cast))], 1, function(x) {any(x!=0)} )]
solotopi <- setdiff(topi, notopi)
data <- data.frame(id=rownames(freqs_cast), freqs=freqs_cast[, grepl('topo', colnames(freqs_cast))])
data$solotopi <- ifelse(data$id %in% solotopi, 'yes', 'no')
table(data$solotopi)

ggplot(data=data, aes(y=freqs.topo_1, x=solotopi))+geom_boxplot(varwidth = TRUE)+scale_y_log10()




topi <- rownames(freqs_cast)[apply(freqs_cast[, !grepl('topo', colnames(freqs_cast))], 1, function(x) {all(x!=0)} )]
notopi <- rownames(freqs_cast)[apply(freqs_cast[, grepl('topo', colnames(freqs_cast))], 1, function(x) {any(x!=0)} )]
solotopi <- setdiff(topi, notopi)
data <- data.frame(id=rownames(freqs_cast), freqs=freqs_cast[, !grepl('topo', colnames(freqs_cast))])
data$solotopi <- ifelse(data$id %in% solotopi, 'yes', 'no')
table(data$solotopi)

ggplot(data=data, aes(y=freqs.basale_1, x=solotopi))+geom_boxplot(varwidth = TRUE)+scale_y_log10()

#########

library(reshape)
library(ggplot2)

load('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/invivo_feb2022/counts.Rdata')


dd <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(dd) <- dd$sequence
dd$sequence <- NULL
all_counts_wide <- dd
all_counts_wide <- all_counts_wide + 1

freqs_cast <- as.data.frame(apply(all_counts_wide, 2, function(x) x/sum(x)))


get_ave <- function(name, data) {
  sub <- data[,grepl(name, colnames(data), fixed=TRUE)]
  return(apply(sub, 1, mean))
}

ns <- c('basale', 'topo', 'vitro')

freq_ave <- as.data.frame(sapply(ns, get_ave, freqs_cast))

freq_ave_s <- freq_ave[freq_ave$basale > 0.0005,]

logfc_1 <- data.frame(logfc=log(freq_ave_s$topo/freq_ave_s$basale), fill=rep('vivo', nrow(freq_ave_s)))
logfc_2 <- data.frame(logfc=log(freq_ave_s$vitro/freq_ave_s$basale), fill=rep('vitro', nrow(freq_ave_s)))


ggplot(data=rbind(logfc_1, logfc_2), aes(x=logfc, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=15, color='black', data=logfc_1)+
  geom_histogram(aes(y=..count..), bins=15, color='black', data=logfc_2)+
  theme_bw()+scale_fill_manual(values=c('grey','red'))+coord_flip()+ylab('n barcodes')

logfc_1 <- data.frame(logfc=log(freq_ave$topo/freq_ave$basale), fill=rep('vivo', nrow(freq_ave)))
logfc_2 <- data.frame(logfc=log(freq_ave$vitro/freq_ave$basale), fill=rep('vitro', nrow(freq_ave)))


ggplot(data=rbind(logfc_1, logfc_2), aes(x=logfc, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=15, color='black', data=logfc_1)+
  geom_histogram(aes(y=..count..), bins=15, color='black', data=logfc_2)+
  theme_bw()+scale_fill_manual(values=c('grey','red'))+coord_flip()


ref <- "basale_1"
logfc <- function(col, reference, data) {
  return(log(data[,col]/data[, reference]))
}

topos <- grep('topo_', colnames(freqs_cast), fixed=TRUE)

singlelfc <- as.data.frame(sapply(topos, logfc, ref, freqs_cast))
colnames(singlelfc) <- colnames(freqs_cast)[topos]
rownames(singlelfc) <- rownames(freqs_cast)
singlelfc[, ref] <- freqs_cast[, ref]

ss <- singlelfc[singlelfc[,ref] > 0.0005,]

cor.test(ss$topo_1, ss$basale_1)
ss$id <- rownames(ss)
longdf <- melt(ss, measure.vars = colnames(freqs_cast)[topos])

ggplot(data=longdf, aes(x=basale_1,y=value))+geom_point()+facet_wrap(~variable) + theme_bw()



topos <- grep('vitro_', colnames(freqs_cast), fixed=TRUE)

singlelfc <- as.data.frame(sapply(topos, logfc, ref, freqs_cast))
colnames(singlelfc) <- colnames(freqs_cast)[topos]
rownames(singlelfc) <- rownames(freqs_cast)
singlelfc[, ref] <- freqs_cast[, ref]

ss2 <- singlelfc[singlelfc[,ref] > 0.0005,]

ss2$id <- rownames(ss2)
cor.test(ss2$vitro_1, ss2$basale_1)
longdf2 <- melt(ss2, measure.vars = colnames(freqs_cast)[topos])

ggplot(data=longdf2, aes(x=basale_1,y=value))+geom_point()+facet_wrap(~variable) + theme_bw()


cc <- cor(freqs_cast)
library(corrplot)
corrplot(cc)


all(ss2$id==ss$id)
sss <- cbind(ss, ss2)


cor.test(sss$topo_1, sss$vitro_1)

cor.test(sss$topo_1, sss$topo_3)

cor.test(sss$vitro_1, sss$vitro_2)

compare(sss$topo_1, sss$topo_3, FALSE, 'topo_1', 'topo_3')
compare(sss$topo_1, sss$vitro_1, FALSE, 'topo_1', 'vitro_1')
compare(sss$vitro_2, sss$vitro_1, FALSE, 'vitro_2', 'vitro_1')




