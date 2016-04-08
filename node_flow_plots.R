library(ggplot2)

##################

theme_set(theme_gray(base_size = 24))
# 
# filenames <- c("/Applications/XAMPP/htdocs/splitthreader/user_data/example2/Oct28_Evolution_unfiltered.node_flows.csv","/Applications/XAMPP/htdocs/splitthreader/user_data/example2/Oct28_Evolution_SRV_filtered.node_flows.csv")
# 
# 
# for (i in c(1:length(filenames))) {
#     filename <- filenames[i]


    filename <- "/Applications/XAMPP/htdocs/splitthreader/user_data/example2/SplitThreader.after_balancing.node_flows.csv"
    mydata <- read.csv(filename)
    

    # Calculating differences
    mydata$diff <- mydata$flow_in - mydata$flow_out
    mydata$in_vs_weight <- mydata$flow_in - mydata$node_weight
    mydata$out_vs_weight <- mydata$flow_out - mydata$node_weight



    # Pretty ordering for plot
    mydata$chromosome <- factor(mydata$chromosome, levels = c(seq(1,22),"X","Y"))
    mydata$in_matches_weight <- factor(mydata$in_matches_weight, levels = c("in matches CN: True","in matches CN: False"))
    mydata$out_matches_weight <- factor(mydata$out_matches_weight, levels = c("out matches CN: False","out matches CN: True"))


    message(paste(filename,".PLOT.png", sep=""))
    
    png(paste(filename,".PLOT.png", sep=""),width=1000,height=1000)
    print(ggplot(mydata,aes(x=flow_in,y=flow_out,color=flow_balance)) + geom_point() + xlim(0,1000) + ylim(0,1000) + facet_grid(in_matches_weight ~ out_matches_weight))
    
    dev.off()

    print(ggplot(mydata,aes(x=flow_in,y=flow_out,color=flow_balance)) + geom_point() + xlim(0,1000) + ylim(0,1000) + facet_grid(in_matches_weight ~ out_matches_weight))
    
    # By chromosome
    # ggplot(mydata,aes(x=flow_in,y=flow_out)) + geom_point() + xlim(0,1000) + ylim(0,1000) + facet_wrap(~ chromosome)
# }


