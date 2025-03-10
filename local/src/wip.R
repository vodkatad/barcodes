library(scales)

library(ggplot2)
s <- c("bulk_1","bulk_2","bulk_3", "linea_S", "linea_M", "linea_L", "Basale_0X002", "Basale_0X006", "Basale_0X007", "Phy_0X001", "Phy_0X008", "Phy_0X011", "Cetux_0x003", "Cetux_0x004", "Cetux_0x005", "Cetux_0x009", "Cetux_0x0012", "Cetux_0x0013")
o <- LETTERS[1:18]
map <- data.frame(o=o, s=s, stringsAsFactors = F)

setwd('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/sept2023')
load('counts.Rdata')
tot <- read.table('../sept2023_phix/info.tsv', sep="\t", header=F, stringsAsFactors = F)
tot$o <- sapply(strsplit(tot$V1, "_"), function(x){x[[5]][1]})

tot <- merge(map, tot, by="o")


m <- merge(tot, numbers, by.x="s", by.y="row.names")
m <- m[c(4,5,6,1,2,3,16,17,18,7,8,9,10,11,12,13,14,15),]
m$frac_barcodes <- m$reads_with_barcodes/m$V3
m$group <- sapply(strsplit(m$s, "_"), function(x){x[[1]][1]})
m$s <- factor(m$s, levels = m$s)
ggplot(data=m)+geom_col(aes(x=s, y=n_seen_barcodes, fill=group))+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text=element_text(size=15))

ggplot(data=m)+geom_col(aes(x=s, y=frac_barcodes, fill=group))+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(data, aes(x=count)) + geom_histogram(bins=50)+facet_wrap(~factor(sample, levels=m$s), ncol=3)+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                                                                                    labels = trans_format("log10", math_format(10^.x)))+theme_bw()+xlab('n.reads')+ylab('n.barcodes')
                                      
                                                                                                                   
n_samples <- 18                                                                                                                   
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
  geom_histogram(bins=50)+facet_wrap(~factor(sample, levels=m$s), ncol=3)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+xlab('n.reads')+ylab('n.barcodes')


library(reshape)
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


base <-  'Basale_0X006'
phy <-  'Phy_0X008'
cetux <-  'Cetux_0x0013'

get_freqs <- function(b, p, c, data) {
  subset <- data[data$sample %in% c(b, p, c),]
  subset <- subset[order(subset$seq),]
  subset_y1 <- subset[subset$sample == c,'freq'] 
  subset_y2 <- subset[subset$sample == p,'freq']
  subset_x <- subset[subset$sample==b,'freq']
  ssubset_y1 <- subset[subset$sample == c,'seq']
  ssubset_y2 <- subset[subset$sample == p,'seq']
  ssubset_x <- subset[subset$sample==b,'seq']
  stopifnot(all(ssubset_y1 == ssubset_y2))
  stopifnot(all(ssubset_x == ssubset_y2))
  data.frame(x=rep(subset_x, 2), y=c(subset_y1, subset_y2), 
             fill=c(rep('Cetux', length(subset_y1)),rep('Phy', length(subset_y2))), 
                    rep=rep(paste(b, p, c, sep='-'), length(subset_x)*2))
}



pd <- get_freqs(base, phy, cetux, long_fr)


base <-  'Basale_0X007'
phy <-  'Phy_0X001'
cetux <-  'Cetux_0x003'
pd2 <- get_freqs(base, phy, cetux, long_fr)

pd <- rbind(pd, pd2)
thr <- 0.0005
ggplot(data=pd, aes(x=x, y=y, color=fill))+geom_point()+facet_wrap(~rep)+
  scale_color_manual(values=c('red','grey'))+theme_bw()+xlab('Initial frequency')+
  ylab('Final Frequency')+ggtitle('CRC0322')+geom_vline(xintercept=thr)


get_log <- function(b, p, c, data) {
  subset <- data[data$sample %in% c(b, c, p),]
  subset0 <- subset[(subset$sample == b & subset$freq > thr),]
  subset1 <- subset[(subset$sample == c & subset$seq %in% subset0$seq),]
  subset2 <- subset[(subset$sample == p  & subset$seq %in% subset0$seq),]
  subset <- rbind(subset1, subset2, subset0)
  subset <- subset[order(subset$seq),]
  subset_y1 <- subset[subset$sample == c,'freq']
  subset_y2 <- subset[subset$sample == p,'freq']
  subset_x <- subset[subset$sample == b,'freq']
  seqs <- as.character(subset[subset$sample == b,'seq'])
  logfr_1 <- log2(subset_y1/subset_x)
  logfr_2 <- log2(subset_y2/subset_x)
  pdata1 <- data.frame(seq=seqs, logfr=logfr_1, fill=rep('Cetux', length(logfr_1)),rep=rep(paste(b, p, c, sep='-'), length(logfr_1)), stringsAsFactors = FALSE)
  pdata2 <- data.frame(seq=seqs, logfr=logfr_2, fill=rep('Phy', length(logfr_1)), rep=rep(paste(b, p, c, sep='-'), length(logfr_2)), stringsAsFactors = FALSE)
  return(list(pdata2, pdata1))
}

base <-  'Basale_0X006'
phy <-  'Phy_0X008'
cetux <-  'Cetux_0x0013'
pd <- get_log(base, phy, cetux, long_fr)

ggplot(data=rbind(pd[[1]],pd[[2]]), aes(x=logfr, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=20, color='black', data=pd[[1]])+
  geom_histogram(aes(y=..count..), bins=20, color='black', data=pd[[2]])+
  facet_wrap(~rep)+
  theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()

base <-  'Basale_0X007'
phy <-  'Phy_0X001'
cetux <-  'Cetux_0x0013'
pd <- get_log(base, phy, cetux, long_fr)

ggplot(data=rbind(pd[[1]],pd[[2]]), aes(x=logfr, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=20, color='black', data=pd[[1]])+
  geom_histogram(aes(y=..count..), bins=20, color='black', data=pd[[2]])+
  facet_wrap(~rep)+
  theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()


base <-  'Basale_0X002'
phy <-  'Phy_0X001'
cetux <-  'Cetux_0x0013'
pd <- get_log(base, phy, cetux, long_fr)

ggplot(data=rbind(pd[[1]],pd[[2]]), aes(x=logfr, fill=fill))+
  geom_histogram(aes(y=-1*..count..), bins=20, color='black', data=pd[[1]])+
  geom_histogram(aes(y=..count..), bins=20, color='black', data=pd[[2]])+
  facet_wrap(~rep)+
  theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()
library(ggpubr)
wanted_c <- c( 'Cetux_0x0012',  'Cetux_0x0013',  'Cetux_0x003',  'Cetux_0x004', 'Cetux_0x005', 'Cetux_0x009')
base <-  'Basale_0X006'
phy <-  'Phy_0X008'

plotlist <- list()
logsc <- list()
logsp <- list()
i <- 1
for (c in wanted_c) {
  pd <- get_log(base, phy, c, long_fr)
  p <- ggplot(data=rbind(pd[[1]],pd[[2]]), aes(x=logfr, fill=fill))+
    geom_histogram(aes(y=-1*..count..), bins=20, color='black', data=pd[[1]])+
    geom_histogram(aes(y=..count..), bins=20, color='black', data=pd[[2]])+
    facet_wrap(~rep)+
    theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()
  plotlist[[i]] <- p
  logsc[[i]] <- pd[[2]]
  logsp[[i]] <- pd[[1]]
  i <- i +1
}

ggarrange(plotlist=plotlist)


indc <- lapply(logsc, function(x){x[x$logfr>1, 'seq']})
names(indc) <- wanted_c
indp <- lapply(logsp, function(x){x[x$logfr>1, 'seq']})

names(indc) <- wanted_c

upset(fromList(indc))
upset(fromList(indc), text.scale=2)

indc2 <- indc
indc2$Cetux_0x009 <- NULL
indc2$Cetux_0x004 <- NULL
indc2$Cetux_0x005 <- NULL
common <- Reduce(intersect, indc2) 
common[common %in% indp[[1]] ]


wanted_c <- c( 'Cetux_0x0012',  'Cetux_0x0013',  'Cetux_0x003',  'Cetux_0x004', 'Cetux_0x005', 'Cetux_0x009')
base <-  'Basale_0X007'
phy <-  'Phy_0X008'

plotlist <- list()
logsc <- list()
logsp <- list()
i <- 1
for (c in wanted_c) {
  pd <- get_log(base, phy, c, long_fr)
  p <- ggplot(data=rbind(pd[[1]],pd[[2]]), aes(x=logfr, fill=fill))+
    geom_histogram(aes(y=-1*..count..), bins=20, color='black', data=pd[[1]])+
    geom_histogram(aes(y=..count..), bins=20, color='black', data=pd[[2]])+
    facet_wrap(~rep)+
    theme_bw()+scale_fill_manual(values=c('red','grey'))+coord_flip()
  plotlist[[i]] <- p
  logsc[[i]] <- pd[[2]]
  logsp[[i]] <- pd[[1]]
  i <- i +1
}

ggarrange(plotlist=plotlist)


indc <- lapply(logsc, function(x){x[x$logfr>1, 'seq']})
names(indc) <- wanted_c
indp <- lapply(logsp, function(x){x[x$logfr>1, 'seq']})

names(indc) <- wanted_c

upset(fromList(indc), text.scale=2)

indc2 <- indc
indc2$Cetux_0x009 <- NULL
indc2$Cetux_0x004 <- NULL
indc2$Cetux_0x005 <- NULL
common2 <- Reduce(intersect, indc2) 
common2[common2 %in% indp[[1]] ]

### baloons
freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
head(freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]

cast <- cast(fr, id ~ sample, value="fraction")
rownames(cast) <- cast$id
cast$id <- NULL

ll <- apply(cast, 1, function(x) {any(x > 0.005)})

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
rev$variable <- factor(rev$sample, levels=m$s)

#rev$barcode <- as.factor(sample(nrow(rev)))
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# order on freq in basale 1

o <- rev[rev$sample=="Basale_0X002",]
o <- o[order(o$value),]
rev$id <- factor(rev$id, levels=o$id)
rev$variable <- rev$sample
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##################### average bulk e basali
#cast <- as.data.frame(freqs_cast) # all pc

#cast$Basale_0X002 <- NULL
avgBasali <- rowMeans(cast[, grepl('Basale', colnames(cast))])

avgBulk <- rowMeans(cast[, grepl('bulk', colnames(cast))])

keep <- row.names(cast)[avgBulk > 0.0005]
#keep <- row.names(cast)[avgBulk > 0.00005]

sel <- cast[rownames(cast) %in% keep,]
sel$basale <- rowMeans(sel[, grepl('Basale', colnames(sel))])
sel$bulk <- rowMeans(sel[, grepl('bulk', colnames(sel))])
logfc <- function(col, ref, data) {
  log2(data[,col] / data[,ref])
}

cet_lfc <- as.data.frame(sapply(colnames(sel)[grepl('Cetux', colnames(sel))], logfc, 'bulk', sel))
phy_lfc <- as.data.frame(sapply(colnames(sel)[grepl('Phy', colnames(sel))], logfc, 'bulk', sel))
rownames(cet_lfc) <- row.names(sel)
rownames(phy_lfc) <- row.names(sel)
### plot lfc
plot_dist <- function(df) {
  df$id <- row.names(df)
  pd <- melt(df)
  print(ggplot(data=pd, aes(x=value, color=variable))+geom_density(aes(y= ..count.. ))+theme_bw()+
    theme(text=element_text(size=15)))
}

plot_dist(cet_lfc)
plot_dist(phy_lfc)

cetux <- c("Cetux_0x0012", "Cetux_0x0013", "Cetux_0x003",  "Cetux_0x004")
cet_lfc_b <- cet_lfc

cet_lfc <- cet_lfc[, cetux]
phy <- 'Phy_0X008'
phy_lfc <- phy_lfc[, phy, drop=F]
nc <- apply(cet_lfc, 1, function(x) {sum(x > 1)})
np <- apply(phy_lfc, 1, function(x) {sum(x > 1)})

table(nc)
table(np)


upc <- setdiff(names(nc)[nc==4],names(np)[np==1])

dfc <- as.data.frame(sel[rownames(sel) %in% upc,c("bulk_1","bulk_2","bulk_3", "linea_S", "linea_M", "linea_L", "Basale_0X002", "Basale_0X006", "Basale_0X007", "Phy_0X001", "Phy_0X008", "Phy_0X011", "Cetux_0x003", "Cetux_0x004", "Cetux_0x005", "Cetux_0x009", "Cetux_0x0012", "Cetux_0x0013")])

dfc$id <- rownames(dfc)
m <- melt(dfc)

m$s <- factor(m$variable, levels=c("bulk_1","bulk_2","bulk_3", "linea_S", "linea_M", "linea_L", "Basale_0X002", "Basale_0X006", "Basale_0X007", "Phy_0X001", "Phy_0X008", "Phy_0X011", "Cetux_0x003", "Cetux_0x004", "Cetux_0x005", "Cetux_0x009", "Cetux_0x0012", "Cetux_0x0013"))
ggplot(data=m, aes(x=s, y=value, colour=id, group=id))+geom_point()+geom_line()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text=element_text(size=15))

##

### do we see same n. of up if we randomize?
set.seed(42)
n <- c()
for (i in seq(1, 1000)) {
  rsel <- sel
  cets <- colnames(sel)[grepl('Cetux', colnames(sel))]
  for (c in cets) {
    rsel[,c] <- sample(sel[,c], size=nrow(sel))
  }
  
  cet_lfc <- as.data.frame(sapply(colnames(rsel)[grepl('Cetux', colnames(rsel))], logfc, 'bulk', rsel))
  phy_lfc <- as.data.frame(sapply(colnames(rsel)[grepl('Phy', colnames(rsel))], logfc, 'bulk', rsel))
  rownames(cet_lfc) <- row.names(sel)
  rownames(phy_lfc) <- row.names(sel)
  
  colnames(cet_lfc) <- paste0('rand_', colnames(cet_lfc))
  
  cetux <- paste0("rand_", c("Cetux_0x0012", "Cetux_0x0013", "Cetux_0x003",  "Cetux_0x004"))
  cet_lfc <- cet_lfc[, cetux]
  phy <- 'Phy_0X008'
  phy_lfc <- phy_lfc[, phy, drop=F]
  nc <- apply(cet_lfc, 1, function(x) {sum(x > 1)})
  np <- apply(phy_lfc, 1, function(x) {sum(x > 1)})
  
  upc <- setdiff(names(nc)[nc==4],names(np)[np==1])
  n <- c(n, length(upc))
}
table(n)
basale <- colnames(sel)[grepl('Basale', colnames(sel))]
#phy <- colnames(sel)[grepl('Phy', colnames(sel))]
allc <- colnames(sel)[grepl('Cetux', colnames(sel))]
dfc <- as.data.frame(sel[rownames(sel) %in% upc, c(cetux, basale, phy)])
dfc$id <- rownames(dfc)
m <- melt(dfc)
m$s <- factor(m$variable, levels = c(basale, phy, cetux))

ggplot(data=m, aes(x=s, y=value, colour=id, group=id))+geom_point()+geom_line()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text=element_text(size=15))

##

upc_s <- setdiff(names(nc)[nc==4],names(np)[np>=1])

dfc <- as.data.frame(sel[rownames(sel) %in% upc_s, c(cetux, basale, phy)])
dfc$id <- rownames(dfc)
m <- melt(dfc)
m$s <- factor(m$variable, levels = c(basale, phy, cetux))

ggplot(data=m, aes(x=s, y=value, colour=id, group=id))+geom_point()+geom_line()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text=element_text(size=15))




