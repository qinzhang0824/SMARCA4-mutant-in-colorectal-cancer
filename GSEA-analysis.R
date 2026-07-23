library(ggplot2)
setwd('/SMARCA4_Data/7.TCGA_data')
Data<-read.table('gsea_C2.txt',sep="\t",header=T)
Data$Description = factor(Data$Description,levels = Data$Description,ordered = T)

pdf('TCGA.COAD_SMARCA4_Hallmark_GSEA.pdf',width=8.5,height=6)
p <- ggplot(Data,aes(x = NES,y = Description))+geom_point(aes(color = pvalue,size = setSize))+scale_color_gradient(low = "red", high = "blue")+ xlab("Normalized Enrichment Score(NES)")+theme_bw()+guides(color = guide_colorbar(reverse = TRUE))
p
dev.off()
