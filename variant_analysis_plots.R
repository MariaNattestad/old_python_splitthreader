library(ggplot2)

##################
dataset_files <- c("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7.SRVs.csv","/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7_10_reads.SRVs.csv")
dataset_names <- c("5 reads","10 reads")
output_filename <- "/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles.PLOT.SRVs.Feb7.5_vs_10_reads.png"
##################

##################
# dataset_files <- c("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Oct28.SRVs.csv","/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7_10_reads.SRVs.csv")
# dataset_names <- c("Oct28","Feb7")
# output_filename <- "/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles.PLOT.SRVs.Oct28_vs_Feb7_10_reads.png"
##################


combined_data <- data.frame()
for (i in c(1:length(dataset_files))) {
    filename <- dataset_files[i]
    dataset <- read.csv(filename)
    dataset <- dataset[,c("flow_category","description")]
    dataset$Dataset <- dataset_names[i]
    combined_data <- rbind(combined_data,dataset)
}

nrow(combined_data)

combined_data$flow_category <- factor(combined_data$flow_category, levels = c("Perfect","Great","Poor","Bad"))
combined_data$description <- factor(combined_data$description, levels =  c("exactly 1 CNV on each side", "1 or more explanatory CNVs","explanatory CNVs only on one side","CNVs missing or in wrong direction"))
combined_data$Dataset <- factor(combined_data$Dataset, levels = dataset_names)


theme_set(theme_gray(base_size = 20)) 

png(file=output_filename,width=1000,height=1000)
ggplot(combined_data, aes(x=Dataset,fill=description)) + geom_bar() + facet_grid(. ~ flow_category)
dev.off()

###########################################################################################################################################



##############     CNVs    #################


library(ggplot2)

##################
dataset_files <- c("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7.CNVs.tab","/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7_10_reads.CNVs.tab")
dataset_names <- c("5 reads","10 reads")
output_filename <- "/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles.PLOT.CNVs.Feb7.5_vs_10_reads.png"
##################

##################
# dataset_files <- c("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Oct28.CNVs.tab","/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles_Feb7_10_reads.CNVs.tab")
# dataset_names <- c("Oct28","Feb7")
# output_filename <- "/Applications/XAMPP/htdocs/splitthreader/user_data/example2/CNV_SRV_concordance_Sniffles.PLOT.CNVs.Oct28_vs_Feb7_10_reads.png"
##################


combined_data <- data.frame()
for (i in c(1:length(dataset_files))) {
    filename <- dataset_files[i]
    dataset <- read.csv(filename,sep="\t")
    dataset$Dataset <- dataset_names[i]
    combined_data <- rbind(combined_data,dataset)
}

nrow(combined_data)
head(combined_data)

# combined_data$flow_category <- factor(combined_data$flow_category, levels = c("Perfect","Great","Poor","Bad"))
# combined_data$description <- factor(combined_data$description, levels =  c("exactly 1 CNV on each side", "1 or more explanatory CNVs","explanatory CNVs only on one side","CNVs missing or in wrong direction"))

combined_data$Dataset <- factor(combined_data$Dataset, levels = dataset_names)
combined_data$copy_number_change <- abs(combined_data$left_copy_number - combined_data$right_copy_number)

theme_set(theme_gray(base_size = 20)) 


png(file=output_filename,width=1000,height=1000)
ggplot(combined_data, aes(x=Dataset,fill=SRV_evidence)) + geom_bar() + facet_grid(CN_change_categorical ~ segment_lengths)
dev.off()




# ggplot(combined_data, aes(x=Dataset, fill=SRV_evidence)) + geom_bar() + facet_grid(. ~ SRV_evidence)

# ggplot(combined_data,aes(x=left_copy_number, y= right_copy_number,color=SRV_evidence)) + geom_point() + xlim(0,1000) + ylim(0,1000) #+ facet_grid(. ~ Dataset)


###########################################################################################################################################

