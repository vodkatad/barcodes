load('/mnt/trcanmed/snaketree/prj/cellecta_barcode/dataset/CRC0322_cetuxi/pippo.Rdata')
library(ggplot2)
library(reshape2)

ov <- data[data$sample =="CRC0322_overall",]

nreads <- read.table('/mnt/trcanmed/snaketree/prj/cellecta_barcode/local/share/data/n_reads_barcodes_cetuxi.txt', sep="\t", header=T)

ov <- data[data$sample =="CRC0322_overall",]
ov <- rbind(ov, data.frame(id="unaligned", sequence="unalignedreads", 
                           count=nreads[nreads$sample=="CRC0322_overall","unaligned"],
                           sample="CRC0322_overall"))
ov <- ov[order(ov$count),]
ov$index <- ov$count / sum(ov$count)

tot <-  nreads[nreads$sample=="CRC0322_overall","tot"]
subs <- seq(1000000,5000000,by=500000)
toremove <- tot - subs
toremove <- rev(toremove)

subsample <- function(nremove, data) {
  while (nremove > 0) {
    i <- runif(1)
    chosen <- tail(which(data$index <= i), n=1)
    if (length(chosen) == 1) { # cooorner case of chosing 0?
      dc <- data[chosen,'count']
      if (dc >= 1) {
        data[chosen,'count'] <- dc - 1
        ov$index <- ov$count / sum(ov$count)
        nremove <- nremove -1
      }
    }
  }
  print(paste0('done ', nremove))
  return(length(unique(data[data$id != "unaligned" & data$count != 0,'id'])))
}

ov$upto <- cumsum(ov$count)
subsample2 <- function(nremove, data) {
  chosen_removal <- sample.int(sum(data$count), size=nremove, replace=FALSE)
  chosen_removal <- chosen_removal[order(chosen_removal)]
  for (i in seq(1, length(chosen_removal))) {
    chosen <- tail(which(data$upto <= chosen_removal[i]), n=1)
    data[chosen,'count'] <- data[chosen,'count'] - 1
  }
  print(paste0('done2 ', nremove))
  return(length(unique(data[data$id != "unaligned" & data$count != 0,'id'])))
}

n_seen_bc <- sapply(toremove, subsample2, ov)

n_tot <- length(unique(ov[ov$id != "unaligned" & ov$count != 0,'id']))
pd <- data.frame(n_barcode=c(n_seen_bc, n_tot), n_reads=c(rev(subs), tot))
ggplot(data=pd, aes(n_reads, n_barcode))+geom_point()+theme_bw()

model <- lm(data=pd, formula=n_barcode ~ poly(n_reads,3))

xx <- c(pd$n_reads, seq(6000000, 40000000, by=1000000))
fit.ggplot=data.frame(y=predict(model, newdata=data.frame(n_reads=xx)), x=xx)

ggplot(data=pd, aes(n_reads, n_barcode))+geom_point()+theme_bw()+
  geom_line(data=fit.ggplot,aes(x=x,y=y))

xx <- c(pd$n_reads, seq(6000000, 10000000, by=1000000))
fit.ggplot=data.frame(y=predict(model, newdata=data.frame(n_reads=xx)), x=xx)

ggplot(data=pd, aes(n_reads, n_barcode))+geom_point()+theme_bw()+
  geom_line(data=fit.ggplot,aes(x=x,y=y))

# 100000
xx <- c(pd$n_reads, seq(6000000, 20000000, by=1000000))
fit.ggplot=data.frame(y=predict(model, newdata=data.frame(n_reads=xx)), x=xx)

ggplot(data=pd, aes(n_reads, n_barcode))+geom_point()+theme_bw()+
  geom_line(data=fit.ggplot,aes(x=x,y=y))
