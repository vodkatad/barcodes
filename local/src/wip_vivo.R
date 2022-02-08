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
