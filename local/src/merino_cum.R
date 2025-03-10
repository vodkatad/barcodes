# PDX-110: 16 mice at T1, 17 mice at T2

t1 <- runif(3000, 0.0001, 100)
t2 <- runif(2400, 0.00001, 100)

pdata <- data.frame(x=c(t1, t2), class=c(rep('t1', 3000), rep('t2', 2400)))

library(ggplot2)
ggplot(data=pdata, aes(x=x, color=class))+geom_density()+scale_x_log10()


#represented on square-root scale to highlight dominant clones),

t1 <- runif(5000)*100

ggplot(data=data.frame(x=t1), aes(x=x))+geom_histogram(color='blue', fill='white')+theme_bw()+theme(text=element_text(size = 20))
ggplot(data=data.frame(x=t1), aes(x=x))+geom_histogram(color='blue', fill='white')+theme_bw()+theme(text=element_text(size = 20))+scale_x_log10()
ggplot(data=data.frame(x=t1), aes(x=x))+geom_histogram(color='blue', fill='white')+theme_bw()+theme(text=element_text(size = 20))+scale_x_sqrt(breaks=c(0, 25, 50, 75, 100))

x <- seq(1, 100, by=0.01)
d <- data.frame(x=x, y=x)
ggplot(d, aes(x=x, y=y))+geom_line()+theme_bw()+theme(text=element_text(size = 20))
d <- data.frame(x=x, y=log(x))
ggplot(d, aes(x=x, y=y))+geom_line()+theme_bw()+theme(text=element_text(size = 20))
d <- data.frame(x=x, y=sqrt(x))
ggplot(d, aes(x=x, y=y))+geom_line()+theme_bw()+theme(text=element_text(size = 20))