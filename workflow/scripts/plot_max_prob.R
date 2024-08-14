
library(ggplot2)
library(ggpubr)



# plot function for each plasmid and a value
plot.max.prob <- function(sub.df, a_val){

    title <- c('a =', a_val)
    title<- paste(title, sep=' ')

     p <- ggplot(sub.df, aes(x=sigma, y=rloop_prob, color=as.factor(plot_name))) + 
        geom_line(size=1) +
        theme_pubr() + 
        facet_wrap(~orrientation, ncol=1) + 
        labs(x='sigma', y='R-loop prob', title=title) +
        xlim(-0.12, 0)
    
    return(p)

}

df <- read.table(snakemake@input[[1]], sep='\t', header=TRUE)
plasmids <- unique(df$plasmid)
a_vals <- unique(df$a)

pdf(snakemake@output[[1]], width=10, height=10)

print(plasmids)
print(a_vals)

for (each_a in a_vals){

        plot.df <- subset(df, a==each_a)
        p <- plot.max.prob(plot.df, each_a)
        print(p)

}

dev.off()


