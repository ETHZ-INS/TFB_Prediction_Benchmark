library(ggplot2)
library(data.table)

dt <- data.table(x=runif(10), y=runif(10))
ggplot(dt, aes(x=x, y=y))+geom_point()+theme_bw()
dir.create("plots/dummy_fig.png")
ggsave("plots/dummy_fig.png")