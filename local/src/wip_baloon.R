load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(scales)
library(reshape)
library(pheatmap)
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
fr <- freqs(,c('id', 'sample', 'fraction'))
fr <- freqs[,c('id', 'sample', 'fraction')]

library(reshape2)
cast <- dcast(fr, id ~ sample, value.var="fraction")
rownames(cast) <- cast$id
cast$id <- NULL

ll <- apply(cast, 1, function(x) {any(x > 0.0001)})

sel <- cast[ll,]
sel$id <- rownames(sel)
#sel <- head(sel, n=100)
rev <- melt(sel)



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


o <- rev[rev$variable=="CRC0322_overall",]
o <- o[order(o$value),]
rev$id <- factor(rev$id, levels=o$id)

#"CRC0322_overall","CRC0322_A_T0","CRC0322_A_NT","CRC0322_A_CTX", "CRC0322_B_T0","CRC0322_B_NT","CRC0322_B_CTX", "CRC0322_C_T0","CRC0322_C_NT","CRC0322_C_CTX"

#ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
#                                                data = rev, size = size, fill = fill, 
#                                                color=fill, alpha=0.7) + scale_size(range=c(0,15))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
#                  axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
#                  axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rev$barcode <- as.factor(sample(nrow(rev)))
ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = barcode, 
                                          color=barcode, alpha=0.7) + scale_size(range=c(-2,10))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                      axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                      axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(rev, aes_string(x, y)) + geom_exec(geom_point, 
                                          data = rev, size = size, fill = fill, 
                                          color=fill, alpha=0.7) + scale_size(range=c(-2,10))+ theme_minimal() + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                          axis.title.y = ggplot2::element_blank(), axis.text.y=element_blank(),
                                                                                                                          axis.ticks.y=element_blank(), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ok we do need enrichments
