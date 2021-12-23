load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(scales)
library(reshape)
library(pheatmap)
library(UpSetR)

#####
all_counts_wide <- dcast(data, sequence ~ sample, value.var="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL
write.table(all_counts_wide, file='barcodes_counts_wide.tsv.gz', sep="\t", quote=F)
barcode_appeareance <- lapply(unique(data$sample), function(x) {data[data$sample==x,'sequence']})
names(barcode_appeareance)  <- unique(data$sample)
upset(fromList(barcode_appeareance), order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), sets=unique(data$sample), keep.order=TRUE)
upset(fromList(barcode_appeareance), order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), sets=unique(data$sample)[c(2,3,4)], keep.order=TRUE)
seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==10])
data_common <- data[data$sequence %in% keepbc,]
get_fraction <- function(sample, counts) {
  res <- counts[counts$sample == sample,]
  tot <- sum(res$count)
  res$fraction <- res$count / tot
  return(res)
}
ggplot(data_common, aes(x=count)) + geom_histogram(bins=50)+facet_wrap(~sample, ncol=3)+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
labels = trans_format("log10", math_format(10^.x)))+theme_bw()+xlab('n.reads')+ylab('n.barcodes')
freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
head(freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]

ggplot(freqs, aes(x=fraction)) + geom_histogram(bins=50)+facet_wrap(~sample, ncol=3)+theme_bw()+xlab('fraction tot reads')+ylab('n.barcodes')+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                                                                                                                             labels = trans_format("log10", math_format(10^.x)))
#
#

library(reshape2)
cast <- dcast(fr, id ~ sample, value.var="fraction")
rownames(cast) <- cast$id
cast$id <- NULL

ll <- apply(cast, 1, function(x) {any(x > 0.001)})
# noone 0.01
# 239 0.001
# 4402 0.00001
sel <- cast[ll,]
#sel <- cast
sel$id <- rownames(sel)
#sel <- head(sel, n=100)
rev <- melt(sel)

freqs$atleast <- ifelse(freqs$id %in% rev$id, 'atleast','other')
ggplot(freqs, aes(x=reorder(id,-fraction),fraction,fill=atleast)) + geom_bar(stat="identity")+facet_wrap(~sample, ncol=3)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank, axis.ticks.x=element_blank())+theme_bw()
ggplot(freqs, aes(x=reorder(id,-count),log10(count),fill=atleast)) + geom_bar(stat="identity")+ theme(axis.title.x=element_blank(), axis.text.x=element_blank, axis.ticks.x=element_blank())+facet_wrap(~sample, ncol=3)+theme_bw()

#pheatmap(sel)

#library(RColorBrewer)
#n <- 100
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


#palette <- brewer.pal(n=12, 'Paired')
#all <- palette
#for (i in seq(1,10)) {
#  all <- c(all, palette)
#}


x <- 'variable'
y <- 'id'
size <- 'value'
fill <- 'id'
barcode <- 'barcode'
rev$variable <- factor(rev$variable, levels=c('CRC0322_overall', 'CRC0322_A_T0','CRC0322_B_T0','CRC0322_C_T0','CRC0322_A_NT','CRC0322_B_NT', 'CRC0322_C_NT','CRC0322_A_CTX','CRC0322_B_CTX', 'CRC0322_C_CTX'))

#rev <- rev[order(rev$variable),]
#o <- rev[rev$variable=="CRC0322_overall",]
#o <- o[order(o$value),]
#rev$id <- factor(rev$id, levels=o$id)

#"CRC0322_overall","CRC0322_A_T0","CRC0322_A_NT","CRC0322_A_CTX", "CRC0322_B_T0","CRC0322_B_NT","CRC0322_B_CTX", "CRC0322_C_T0","CRC0322_C_NT","CRC0322_C_CTX"

#ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
#                                                data = rev, size = size, fill = fill, 
#                                                color=fill, alpha=0.7) + scale_size(range=c(0,15))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
#                  axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
#                  axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rev$barcode <- as.factor(sample(nrow(rev)))
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = barcode, 
                                          color=barcode, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                          axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                          axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

o <- rev[rev$variable=="CRC0322_overall",]
o <- o[order(o$value),]
rev$id <- factor(rev$id, levels=o$id)
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = barcode, 
                                          color=barcode, alpha=0.7) + scale_size(range=c(-1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                         axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                         axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# diversity

library(vegan)

dd <- dcast(data, sequence ~ sample, value.var="count", fill=0)
dd$sequence <- NULL
Hall <- diversity(t(dd))
plotd <- data.frame(entropy=Hall, sample=as.character(names(Hall)), stringsAsFactors = FALSE)
plotd$condition <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
plotd$rep <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
plotd[is.na(plotd$condition), ]$condition = 'overall'
plotd$condition <- factor(plotd$condition, levels=c('overall','T0','NT','CTX'))
ggplot(plotd, aes(x=condition, y=entropy, fill=condition))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c('white',"darkgoldenrod","darkgreen","red"))+theme(text=element_text(size=15))


dd <- dcast(data_common, sequence ~ sample, value.var="count", fill=0)
dd$sequence <- NULL
H <- diversity(t(dd))
plotd <- data.frame(entropy=H, sample=as.character(names(H)), stringsAsFactors = FALSE)
plotd$condition <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
plotd$rep <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
plotd[is.na(plotd$condition), ]$condition = 'overall'
plotd$condition <- factor(plotd$condition, levels=c('overall','T0','NT','CTX'))
ggplot(plotd, aes(x=condition, y=entropy, fill=condition))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c('white',"darkgoldenrod","darkgreen","red"))+theme(text=element_text(size=15))


data_small <- data[data$id %in% rev$id,]
dd <- dcast(data_small, sequence ~ sample, value.var="count", fill=0)
dd$sequence <- NULL
H <- diversity(t(dd))
plotd <- data.frame(entropy=H, sample=as.character(names(H)), stringsAsFactors = FALSE)
plotd$condition <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][3]; return(y[1])})
plotd$rep <-  sapply(plotd$sample, function(x) {y<-strsplit(x, '_')[[1]][2]; return(y[1])})
plotd[is.na(plotd$condition), ]$condition = 'overall'
plotd$condition <- factor(plotd$condition, levels=c('overall','T0','NT','CTX'))
ggplot(plotd, aes(x=condition, y=entropy, fill=condition))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c('white',"darkgoldenrod","darkgreen","red"))+theme(text=element_text(size=15))

## fcs

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

fcs_t0_ctx <- sapply(c('A','B','C'), fcgen, 'CRC0322', freqs, "_NT", "_CTX")
colnames(fcs_t0_ctx) <- c('repA', 'repB', 'repC')

m2 <- melt(fcs_t0_ctx)
m2$class <- ifelse(m2$X1 %in% rev$id, 'atleast','other')

ggplot(data=m2, aes(value, ..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)


fcs_t0_NT <- sapply(c('A','B','C'), fcgen, 'CRC0322', freqs, "_T0", "_NT")
colnames(fcs_t0_NT) <- c('repA', 'repB', 'repC')

m2 <- melt(fcs_t0_NT)
m2$class <- ifelse(m2$X1 %in% rev$id, 'atleast','other')

ggplot(data=m2, aes(value, ..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)

t0_CTX_NT <- cbind(fcs_t0_ctx, fcs_t0_NT)
corrplot(cor(t0_CTX_NT))
t0_CTX_NT$atleast <- ifelse(rownames(t0_CTX_NT) %in% sel$id, 'yes','no')
ggplot(t0_CTX_NT, aes(x=A_CTX, y=B_CTX, color=atleast)) +geom_point()+scale_color_manual(vaues=c('black','red'))+theme_bw()
ggplot(t0_CTX_NT, aes(x=A_CTX, y=A_NT, color=atleast)) +geom_point()+scale_color_manual(values=c('black','red'))+theme_bw()

n <- t0_CTX_NT[,c(1,2,3,4,5,6)]
ann <- t0_CTX_NT[,7, drop=F]
pheatmap(n, annotation_row = ann)

#########
castcorr <- cast
castcorr <- castcorr/cast[,10]
castcorr <- castcorr[,-10]
castcorr$id <- rownames(castcorr)
rev <- melt(castcorr)


x <- 'variable'
y <- 'id'
size <- 'value'
fill <- 'id'
barcode <- 'barcode'
rev$variable <- factor(rev$variable, levels=c('CRC0322_A_T0','CRC0322_B_T0','CRC0322_C_T0','CRC0322_A_NT','CRC0322_B_NT', 'CRC0322_C_NT','CRC0322_A_CTX','CRC0322_B_CTX', 'CRC0322_C_CTX'))
rev$barcode <- as.factor(sample(nrow(rev)))

o <- rev[rev$variable=="CRC0322_A_T0",]
o <- o[order(o$value),]
rev$id <- factor(rev$id, levels=o$id)
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = barcode, 
                                          color=barcode, alpha=0.7) + scale_size(range=c(1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                         axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                         axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(1,7))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                     axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                     axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### correct counts by overall

#cast <- dcast(data_common, id ~ sample, value.var="count") # mmm, no, will get < 1 wtf?

rev_nt_ctx <- sapply(c('A','B','C'), fcgen, 'CRC0322', rev, "_NT", "_CTX")
freqs_nt_ctx <- sapply(c('A','B','C'), fcgen, 'CRC0322', freqs, "_NT", "_CTX")

posrev <- apply(rev_nt_ctx, 1, function(x) all(x>0))
posfreq <- apply(freqs_nt_ctx, 1, function(x) all(x>0))

# stupida ovviamente i relativi rimangono uguali

# FC vs overall?
fcoverall <- function(replicate, model, counts, cond2) {
  nt <- counts[counts$sample == paste0(model, "_overall"),]
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

freqs_ctx_overall <- sapply(c('A','B','C'), fcoverall, 'CRC0322', freqs, "_CTX")
m2 <- melt(freqs_ctx_overall)
m2$class <- ifelse(m2$X1 %in% rev$id, 'atleast','other')
ggplot(data=m2, aes(value, y=..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)

### 
js <- function(p,q) {
  n <- 0.5 * (p + q)
  JS <- 0.5 * (sum(p * log(p / n)) + sum(q * log(q / n)))
  JS
}
dc <- as.data.frame(cast)
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
corrplot(jm, is.corr=F)
summary(jm[upper.tri(jm)])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00802 0.02718 0.04682 0.05186 0.07226 0.12890 
# entropy without selected atleast

#mah

###

freqs_ctx_overall <- sapply(c('A','B','C'), fcoverall, 'CRC0322', freqs, "_T0")
m2 <- melt(freqs_ctx_overall)
m2$class <- ifelse(m2$X1 %in% rev$id, 'atleast','other')
ggplot(data=m2, aes(value, y=..count.., color=class))+geom_density()+theme_bw()+facet_wrap(~X2)



###
pieplot <- function(sample, data, thr) {
  f <- data[data$sample==sample,'fraction']
  f <- f[f>thr]
  #pie(f, labels=NA)
  df <- data.frame(f=f, fill=as.factor(seq(1, length(f))))
  ggplot(df, aes(x="", y=f, fill=fill)) + geom_bar(stat="identity", width=1)+
   theme_classic() + 
  theme(axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(), legend.position="none")+scale_fill_discrete() # coord_polar("y", start=0)+
}

pieplot('CRC0322_overall',freqs, 0.0005)
pieplot('CRC0322_overall',freqs, 0.001)
ggplot(fr, aes(x=log10(fraction))) + geom_density()+facet_wrap(~sample, ncol=3)

###########
# 1/x^2
curve(1/x^2, from=1, to=50, , xlab="x", ylab="y")

# cumulative: how many clones larger than this one do we have?
cumfreq <- function(sample, data) {
  fr1 <- data[data$sample==sample,]; 
  excum <- sapply(1:nrow(fr1),function(i) sum(fr1[i,'fraction'] <= fr1$fraction))
  return(data.frame(row.names=fr1$id, cum=excum, fraction=fr1$fraction))
}

allcum <- lapply(unique(fr$sample), cumfreq, fr)

points <- function(data) {
  pd <- as.data.frame(data)
  ggplot(pd, aes(x = fraction, y = cum)) + geom_point() + theme_bw() + ylab('N.clones larger')+xlab('Clone size')
}

all_p <- lapply(allcum, points)

#p+stat_smooth(method = "lm", formula = y ~ 1/x, size = 1, color="red")
cowplot::plot_grid(plotlist = all_p,nrow = 3, labels=unique(fr$sample))

#### sottoriva's plots
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

compare(ccast$CRC0322_A_CTX, ccast$CRC0322_A_NT, TRUE, 'A_CTX', 'A_NT')
compare(ccast$CRC0322_A_CTX, ccast$CRC0322_A_T0, TRUE, 'A_CTX', 'A_T0')
compare(ccast$CRC0322_A_CTX, ccast$CRC0322_overall, TRUE, 'A_CTX', 'overall')

compare(cast$CRC0322_A_CTX, cast$CRC0322_A_NT, TRUE, 'A_CTX', 'A_NT')
compare(cast$CRC0322_A_CTX, cast$CRC0322_A_T0, TRUE, 'A_CTX', 'A_T0')
compare(cast$CRC0322_A_CTX, cast$CRC0322_overall, TRUE, 'A_CTX', 'overall')

comparepc <- function(x, y, log, nx, ny, pc) {
  zx <- length(x[x==0])
  zy <- length(y[y==0])
  zz <- x==0 & y==0
  zz <- sum(zz)
  x[x==0] <- pc
  y[y==0] <- pc
  if (log) {
    x <- log10(x)
    y <- log10(y)
  }
  pe <- cor.test(x, y)
  d <- data.frame(x=x, y=y, density=get_density(x,y, n=100))
  print(ggplot(d, aes(x=x, y=y, color=density)) +geom_point()+theme_bw()+theme(text = element_text(size=20))+xlab(nx)+ylab(ny)+labs(caption=paste0(pe$estimate, ', pval=', pe$p.value))+scale_color_viridis_c())
  return(c(zx, zy, zz))
}



all_freqs <- lapply(unique(data$sample), get_fraction, data)
all_freqs <- do.call(rbind, all_freqs)
all_frac_wide <- dcast(all_freqs, sequence ~ sample, value.var="fraction")
rownames(all_frac_wide) <- all_frac_wide$sequence
all_frac_wide$sequence <- NULL
all_frac_wide[is.na(all_frac_wide)] <- 0.000000001


compare(all_frac_wide$CRC0322_A_T0, all_frac_wide$CRC0322_overall, TRUE, 'CRC0322_T0', 'overall')

comparepc(all_counts_wide$CRC0322_A_T0, all_counts_wide$CRC0322_overall, TRUE, 'CRC0322_T0', 'overall')
####