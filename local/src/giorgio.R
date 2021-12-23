load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library(UpSetR)

#####
# all_counts_wide <- dcast(data, sequence ~ sample, value.var="count", fill=0)
# rownames(all_counts_wide) <- all_counts_wide$sequence
# all_counts_wide$sequence <- NULL
# write.table(all_counts_wide, file='barcodes_counts_wide.tsv.gz', sep="\t", quote=F)
# barcode_appeareance <- lapply(unique(data$sample), function(x) {data[data$sample==x,'sequence']})
# names(barcode_appeareance)  <- unique(data$sample)


giorgio <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/clonetracker_barcodes_count.csv',sep="\t", header=FALSE, stringsAsFactors = FALSE)
map <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/all_seq_id.tsv',sep="\t", header=FALSE, stringsAsFactors = FALSE)

colnames(giorgio) <- c('bc', 'CTXp_original_pool1','CTXp_original_pool2')
colnames(map) <- c('bc', 'seq')

ov <- data[data$sample =="CRC0322_overall",]


giorgio2 <- giorgio[grepl('bc14-00', giorgio$bc, fixed=T),]
giorgio3 <- giorgio[grepl('bc14-010', giorgio$bc, fixed=T),]
giorgio_pool1 <- rbind(giorgio2, giorgio3)


giorgio2 <- giorgio[grepl('bc14-01', giorgio$bc, fixed=T),]
giorgio3 <- giorgio2[!grepl('bc14-010', giorgio2$bc, fixed=T),]
giorgio_pool2 <- giorgio3

dim(giorgio_pool2)
dim(giorgio_pool1)
summary(giorgio_pool2)
summary(giorgio_pool1)


pool1 <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/pool1.tsv',sep="\t", header=FALSE, stringsAsFactors = FALSE)
pool2 <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/pool2.tsv',sep="\t", header=FALSE, stringsAsFactors = FALSE)
# all_counts_wide <- dcast(data, sequence ~ sample, value.var="count", fill=0)
# pdf <- data.frame(row.names=all_counts_wide$sequence, reads=all_counts_wide$CRC0322_overall)
# pdf <- pdf[order(pdf$reads),, drop=F]
# pdf$x <- seq(1, nrow(pdf))
# pdf$cumsum <- cumsum(pdf$reads)
# 
pdf <- ov[order(ov$count),]
colnames(pdf) <- c('id','seq','reads','sample')
pdf$cumsum <- cumsum(pdf$reads)
pdf$x <- seq(1, nrow(pdf))
#xt <- as.data.frame(table(pdf$reads))
#xt$Var1 <- as.numeric(as.character(xt$Var1))
#ggplot(data=xt, aes(x=Var1, y=Freq))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('log10(reads)')+ylab('nBarcodes')+scale_x_log10()
ggplot(data=pdf, aes(x=x, y=reads))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()

pdf1 <- giorgio[order(giorgio$CTXp_original_pool1),c(1,2)]
pdf1$x <- seq(1, nrow(pdf1))
colnames(pdf1)[2] <- 'reads'
ggplot(data=pdf1, aes(x=x, y=reads))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()


pdf1 <- giorgio[order(giorgio$CTXp_original_pool2),c(1,3)]
pdf1$x <- seq(1, nrow(pdf1))
colnames(pdf1)[2] <- 'reads'
ggplot(data=pdf1, aes(x=x, y=reads))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('log10(reads)')+ylab('nBarcodes')+scale_y_log10()
# 
# 

bc <- strsplit(ov$sequence, '-')
ov$bc1 <- sapply(bc, function(x){x[[1]][1]})
ov$bc2 <- sapply(bc, function(x){x[[2]][1]})
library(Biostrings)
ov$rbc1 <- sapply(ov$bc1, function(x) { as.character(reverseComplement(DNAString(x))) })
ov$rbc2 <- sapply(ov$bc2, function(x) { as.character(reverseComplement(DNAString(x))) })
ov$revbc <- paste0(ov$rbc1, '-', ov$rbc2)

m <- merge(ov, map, by.x="revbc", by.y="seq")
m2 <- merge(m, giorgio_pool1, by.x="bc", by.y="bc")
m2$lc0 <- log10(m2$count+1)
m2$lc1 <- log10(m2$CTXp_original_pool1+1)
m2$lc2 <- log10(m2$CTXp_original_pool2+1)
ggplot(data=m2, aes(x=lc0, y=lc1))+geom_point()+current_theme
ggplot(data=m2, aes(x=lc0, y=lc2))+geom_point()+current_theme


library(reshape)
#m2 <- m2[order(m2$count),]
m2 <- m2[order(m2$CTXp_original_pool1),]
or <- m2$bc
m3 <- m2[,c('bc','count','CTXp_original_pool1')]
colnames(m3) <- c('id','overall_CRC0322','original_pool1')
m3$id <- as.character(m3$id)
m3$original_pool1 <- as.character(m3$original_pool1)
m3$overall_CRC0322 <- as.character(m3$overall_CRC0322)
long <- melt(m3, id="id")
long1 <- long[long$variable=="original_pool1",]
long1 <- long1[match(or, long1$id),]
long1$x <- seq(1, nrow(long1))
long2 <- long[long$variable=="overall_CRC0322",]
long2 <- long2[match(or, long2$id),]
long2$x <- seq(1, nrow(long2))
long <- rbind(long1, long2)
long$value <- as.numeric(as.character(long$value))
ggplot(data=long, aes(x=x, y=value, color=variable))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('log10(reads)')+ylab('nBarcodes')+scale_y_log10()+theme(axis.text.x=element_blank())


####
library(Biostrings)
prepare_pool <- function(ov) {
  colnames(ov) <- c('id', 'sequence', 'count')
  bc <- strsplit(ov$sequence, '-')
  ov$bc1 <- sapply(bc, function(x){x[[1]][1]})
  ov$bc2 <- sapply(bc, function(x){x[[2]][1]})
  ov$rbc1 <- sapply(ov$bc1, function(x) { as.character(reverseComplement(DNAString(x))) })
  ov$rbc2 <- sapply(ov$bc2, function(x) { as.character(reverseComplement(DNAString(x))) })
  ov$revbc <- paste0(ov$rbc1, '-', ov$rbc2)
  m <- merge(ov, map, by.x="revbc", by.y="seq")
  m
}

merge_with_g <- function(m, giorgio_df, giorgio_col) {
  m2 <- merge(m, giorgio_df, by.x="bc", by.y="bc")
  m2$lc0 <- log10(m2$count+1)
  m2$lc1 <- log10(m2[, giorgio_col]+1)
  print(cor.test(m2[,'count'], m2[,giorgio_col]))
  print(ggplot(data=m2, aes_string(x='lc0', y='lc1'))+geom_point())
  m2
}

mpool1 <- prepare_pool(pool1)
mp1 <- merge_with_g(mpool1, giorgio_pool1, 'CTXp_original_pool1')

mpool2 <- prepare_pool(pool2)

mp2 <- merge_with_g(mpool2, giorgio_pool2, 'CTXp_original_pool2')
mp3 <- merge_with_g(mpool1, giorgio_pool2, 'CTXp_original_pool2')




#####
merge_with_self <- function(m, o) {
  m2 <- merge(m, o, by="id")
  m2$lc0 <- log10(m2[,'count.x']+1)
  m2$lc1 <- log10(m2[,'count.y']+1)
  print(cor.test(m2[,'count.x'], m2[,'count.y']))
  print(ggplot(data=m2, aes_string(x='lc0', y='lc1'))+geom_point())
  m2
}


mp4 <- merge_with_self(ov, mpool1)
mp5 <- merge_with_self(ov, mpool2)

#################### new slides

load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library(UpSetR)
pool1 <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/pool1.tsv',sep="\t", header=FALSE, stringsAsFactors = FALSE)
pool2 <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/giorgio_library/pool2.tsv',sep="\t", header=FALSE, stringsAsFactors = FALSE)
ov <- data[data$sample =="CRC0322_overall",]

colnames(pool1) <- c('id_lib', 'sequence_lib', 'count_lib')
m <- merge(pool1, ov, by.x="id_lib", by.y="id")

ggplot(data=m, aes(x=log10(count_lib), y=log10(count)))+geom_point()+current_theme
cor.test(log10(m$count_lib), log10(m$count))


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

compare(m$count_lib, m$count, TRUE, 'log10 reads library', 'log10 reads infection')

pool1 <- pool1[order(pool1$count_lib),]
pool1$x <- seq(1, nrow(pool1))
#pool1$common <- ifelse(pool1$id_lib %in% m$id_lib, 'yes', 'no')
#ggplot(data=pool1, aes(x=x, y=count_lib, fill=common))+geom_point(pch=21, colour="black")+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()+scale_colour_manual(values=c('white','red'))

ggplot(data=pool1, aes(x=x, y=count_lib))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()




ov <- ov[order(ov$count),]
ov$x <- seq(1, nrow(ov))
#pool1$common <- ifelse(pool1$id_lib %in% m$id_lib, 'yes', 'no')
#ggplot(data=pool1, aes(x=x, y=count_lib, fill=common))+geom_point(pch=21, colour="black")+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()+scale_colour_manual(values=c('white','red'))

ggplot(data=ov, aes(x=x, y=count))+geom_point()+theme_bw()+theme(text = element_text(size = 15))+xlab('barcodes')+ylab('log10(reads)')+scale_y_log10()


ggplot(data=ov, aes(x=count))+geom_histogram()+current_theme+scale_x_log10()+xlab('reads')

pool1$common <- NULL
ov$sample <- NULL
colnames(pool1) <- colnames(ov)
ov$kind <- 'infection'
pool1$kind <- 'library'
d <- rbind(ov, pool1)

bin_inf = diff(range(ov$count))/30
bin_lib = diff(range(pool1$count))/30

ggplot(data=pool1, aes(x=count))+geom_histogram()+current_theme+scale_x_log10()+xlab('log10 reads')+ylab("n.barcodes")
ggplot(data=ov, aes(x=count))+geom_histogram()+current_theme+scale_x_log10()+xlab('log10 reads')+ylab("n.barcodes")


ggplot(d, aes(count, fill=kind)) +  geom_histogram(data=subset(d, kind=="infection"), aes( y= - ..count..)) +  geom_histogram(data=subset(d, kind=="library"), aes( y= ..count.. ))+theme_bw()+scale_x_log10()