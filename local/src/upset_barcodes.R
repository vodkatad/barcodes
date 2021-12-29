library(UpSetR)
library(reshape)

upsetplot_f <- snakemake@output[['upset']]
wide_f <- snakemake@output[['widecounts']]

image <- snakemake@input[['Rimage']]
load(image)


all_counts_wide <- cast(data, sequence ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$sequence
all_counts_wide$sequence <- NULL
write.table(all_counts_wide, file=wide_f, sep="\t", quote=F)
barcode_appeareance <- lapply(unique(data$sample), function(x) {data[data$sample==x,'sequence']})
names(barcode_appeareance)  <- unique(data$sample)
pdf(upsetplot_f)
upset(fromList(barcode_appeareance), order.by = "freq", text.scale = c(0.5, 0.5, 0.5, 0.5, 1, 0.5), sets=unique(data$sample), keep.order=TRUE, set_size.show=FALSE)
graphics.off()
#upset(fromList(barcode_appeareance), order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 1), sets=unique(data$sample)[c(2,3,4)], keep.order=TRUE)

#seen_bc <- table(data$sequence)
#keepbc <- names(seen_bc[seen_bc==TODO])
#data_common <- data[data$sequence %in% keepbc,]