library("ggpubr")

args=commandArgs(T)
setwd(args[3])

Data<-read.table(args[1],header=T,sep='\t')
Data$MutType <- factor(Data$MutType,levels=c("Wildtype","Mut"))

pdf(args[2],width=20,height=8)
p <- ggboxplot(Data, x = "Celltype", y = "CellFrac",color = "MutType", palette = "lancet", bxp.errorbar=TRUE,bxp.errorbar.width=0.25)+stat_compare_means(aes(group=MutType),label = "p.format")+ylim(0,0.8)+ ylab("Cell composition")+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
p
dev.off()
