library(ggplot2)


filename <- "/Users/mnattest/Desktop/SplitThreader_testcases/Lumpy_SplitThreader.karyotype_saved"

data <- read.csv(filename, header=FALSE,sep="\t")
names(data) <- c("path","chrom","length")
head(data)

ggplot(data,aes(fill=chrom,y=length,x=1)) + geom_bar(position="stack",stat="identity") + facet_wrap(~ path)


