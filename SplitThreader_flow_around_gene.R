library(ggplot2)


genes <- c("ERBB2","MYC","TPD52","BCAS1","BRCA1","BRCA2","MET","EGFR","TP53")

for (i in seq(1,length(genes))) {
    filename <- paste("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/SplitThreader.", genes[i], ".flow.csv", sep="")

    # args<-commandArgs(TRUE)
    # output_prefix <- args[1]
    
    theme_set(theme_gray(base_size = 24))
    
    data <- read.csv(filename)
    
    Mb = function(x) {
        tmp = x/1000000
        paste(tmp, "Mb")
    }
    
    chrom.colors = c('#E41A1C', '#A73C52', '#6B5F88', '#3780B3', '#3F918C', '#47A266','#53A651', '#6D8470', '#87638F', '#A5548D', '#C96555', '#ED761C','#FF9508', '#FFC11A', '#FFEE2C', '#EBDA30', '#CC9F2C', '#AD6428','#BB614F', '#D77083', '#F37FB8', '#DA88B3', '#B990A6', '#999999',"black")
    chromosomes <- c(seq(1,22),"X","Y","maxed out")
    
    
    data$chromosome <- factor(data$chromosome,levels=chromosomes)
    head(data)
    
    png(paste(filename,".flow_plot.png", sep=""),1000,1000)
    print(
        ggplot(data,aes(x=distance,fill=chromosome)) + geom_bar() + scale_x_continuous(labels=Mb) + scale_fill_manual(values=chrom.colors,limits=chromosomes) + labs(y="thread count",title=genes[i])
    )
    dev.off()

}




