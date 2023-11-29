# setup common ###########################
library(scales)
library(ggplot2)
s <- c("bulk_1","bulk_2","bulk_3", "linea_S", "linea_M", "linea_L", "Basale_0X002", "Basale_0X006", "Basale_0X007", "Phy_0X001", "Phy_0X008", "Phy_0X011", "Cetux_0x003", "Cetux_0x004", "Cetux_0x005", "Cetux_0x009", "Cetux_0x0012", "Cetux_0x0013")
o <- LETTERS[1:18]
map <- data.frame(o=o, s=s, stringsAsFactors = F)
setwd('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/sept2023')
load('counts.Rdata')
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
freqs <- lapply(unique(data_common$sample), get_fraction, data_common)
freqs <- do.call(rbind, freqs)
head(freqs)
fr <- freqs[,c('id', 'sample', 'fraction')]
cast <- cast(fr, id ~ sample, value="fraction")
rownames(cast) <- cast$id
cast$id <- NULL

# setup all ###########################
library(scales)
library(ggplot2)
s <- c("bulk_1","bulk_2","bulk_3", "linea_S", "linea_M", "linea_L", "Basale_0X002", "Basale_0X006", "Basale_0X007", "Phy_0X001", "Phy_0X008", "Phy_0X011", "Cetux_0x003", "Cetux_0x004", "Cetux_0x005", "Cetux_0x009", "Cetux_0x0012", "Cetux_0x0013")
o <- LETTERS[1:18]
map <- data.frame(o=o, s=s, stringsAsFactors = F)
setwd('/mnt/cold1/snaketree/prj/cellecta_barcode/dataset/sept2023')
load('counts.Rdata')
library(reshape)
all_counts_wide <- cast(data, id ~ sample, value="count", fill=0)
rownames(all_counts_wide) <- all_counts_wide$id
all_counts_wide$id <- NULL

all_counts_wide <- all_counts_wide + 1

freqs_cast <- apply(all_counts_wide, 2, function(x) x/sum(x))

cast <- as.data.frame(freqs_cast)
