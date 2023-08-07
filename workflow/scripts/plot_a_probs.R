library(ggplot2)
library(ggpubr)
library(RColorBrewer)



plot.bpprobs <- function(df){


    ggplot(df, aes(x=position, y=prob, color=as.factor(a))) +
        geom_line() + theme_pubr()



}

df <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE)

# subset to just collection of values

a.vals <- c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377)
df <- subset(df, df$a %in% a.vals)

plt <- plot.bpprobs(df)
ggsave(snakemake@output[[1]], plt, width=10, height=10, dpi=300, units='in')