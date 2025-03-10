library(reshape)
library(pheatmap)
library(ggplot2)
setwd('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/invivo_feb2022')
load('counts.Rdata')

n_samples <- 10
seen_bc <- table(data$sequence)
keepbc <- names(seen_bc[seen_bc==n_samples])
data_common <- data[data$sequence %in% keepbc,]
get_fraction <- function(sample, counts) {
  res <- counts[counts$sample == sample,]
  tot <- sum(res$count)
  res$fraction <- res$count / tot
  return(res)
}

freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]
cast <- cast(fr, id ~ sample, value="fraction")
rownames(cast) <- cast$id
cast$id <- NULL

avg <- rowMeans(cast)

histogram(-log10(avg))

length(avg)
table(avg > 2.23*10**-4)
table(avg > 5.543*10**-4)

common_thr <- as.data.frame(cast[avg > mean(avg),]) # do we have a bias here in the selection for the basale or..?

pheatmap(common_thr, show_rownames = FALSE)

avg_basale <- rowMeans(common_thr[, grepl('basale', colnames(common_thr))])
avg_vivo <- rowMeans(common_thr[, grepl('topo', colnames(common_thr))])
avg_vitro <- rowMeans(common_thr[, grepl('vitro', colnames(common_thr))])

histogram(-log10(avg_basale))
histogram(-log10(avg_vivo))
histogram(-log10(avg_vitro))

lfc_vivo <- avg_vivo / avg_basale

vivo_up <- common_thr[lfc_vivo > 5,]

pheatmap(vivo_up, show_rownames = FALSE)

lfc_vitro <- avg_vitro / avg_basale

lfc <- data.frame(lfc =c(lfc_vivo, lfc_vitro), sample=c(rep('vivo', length(avg_vitro)), rep('vitro', length(avg_vivo))))

ggplot(data=lfc, aes(y=lfc, x=sample))+geom_violin()+theme_bw(base_size=20)+geom_boxplot(width=0.1)

lfc_up <- lfc[lfc_vivo > 5,]
ggplot(data=lfc_up, aes(y=lfc, x=sample))+geom_violin()+theme_bw(base_size=20)+geom_boxplot(width=0.1)

vivo_up <- common_thr[lfc_vivo > 5,]

pheatmap(vivo_up, show_rownames = FALSE)

mat_upup <- common_thr[lfc_vivo > 15,]

avg_basale <- rowMeans(mat_upup[, grepl('basale', colnames(mat_upup))])
avg_vivo <- rowMeans(mat_upup[, grepl('topo', colnames(mat_upup))])
avg_vitro <- rowMeans(mat_upup[, grepl('vitro', colnames(mat_upup))])

pd <- data.frame(id=rownames(mat_upup), basale=avg_basale, vivo=avg_vivo, vitro=avg_vitro)

mee <- melt(pd)
mee$variable=factor(mee$variable, levels=c('basale', 'vitro', 'vivo'))
ggplot(data=mee, aes(x=variable, y=value, group=id))+geom_point()+geom_line()+theme_bw(base_size=20)

mat_upup <- common_thr[lfc_vivo > 5,]

avg_basale <- rowMeans(mat_upup[, grepl('basale', colnames(mat_upup))])
avg_vivo <- rowMeans(mat_upup[, grepl('topo', colnames(mat_upup))])
avg_vitro <- rowMeans(mat_upup[, grepl('vitro', colnames(mat_upup))])

pd <- data.frame(id=rownames(mat_upup), basale=avg_basale, vivo=avg_vivo, vitro=avg_vitro)

mee <- melt(pd)
mee$variable=factor(mee$variable, levels=c('basale', 'vitro', 'vivo'))
ggplot(data=mee, aes(x=variable, y=value, group=id))+geom_point()+geom_line()+theme_bw(base_size=20)

