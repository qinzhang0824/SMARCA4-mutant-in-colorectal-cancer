library(maftools)
library(RColorBrewer)
library(dplyr)
library(scitb)
library(ComplexHeatmap)
library("ggpubr")

args=commandArgs(T)
help_usage <- function(){
        cat ("Usage:\n")
        cat ("Rscript coxph.R <maf> <clinical> <output file>\n")
}
if(args[1] == '-help' || args[1] == '-h' || length(args) != 4){
        help_usage()
        q()
}

swi_gene_list <- c("ARID1A","ARID1B","ARID2","PBRM1","SMARCA2","SMARCA4","SMARCB1","SMARCD1","SMARCE1")
setwd('/SMARCA4_Data/1.simceredx_data')

#############################################################################################
all <-read.table('simceredx_3273s_clinical_info_addSMARCA4.xls',sep='\t',header=T)

all_TMB <-read.table('simceredx_3273s_clinical_info_addSMARCA4_addTMB.xls',sep='\t',header=T)
all_TMB$Tumor_Sample_Barcode <- all_TMB$order_number
all_tmb <- all_TMB[,c(5,9)]
all_addtmb <- left_join(all,all_tmb,by="Tumor_Sample_Barcode")
write.table(all_addtmb,file='simceredx_3273s_clinical_info_addSMARCA4_addTMB.xls',sep='\t',quote=FALSE,row.names = FALSE)


all_MSI <-read.table('simceredx_MSI.txt',sep='\t',header=T)
all_MSI$Tumor_Sample_Barcode <- all_MSI$order_number
all_msi <-all_MSI[,c(8,9)]
all_addtmb_add_msi <- left_join(all_TMB,all_msi,by="Tumor_Sample_Barcode")
write.table(all_addtmb_add_msi,file='simceredx_3273s_clinical_info_addSMARCA4_addTMB_addMSI.xls',sep='\t',quote=FALSE,row.names = FALSE)

mut<-read.table('simceredx_clinical_info_223s.SMARCA4.Mut.xls',sep='\t',header=T)
all <-as.data.frame(all)
all$Group <-ifelse(all$Tumor_Sample_Barcode %in% mut$Tumor_Sample_Barcode,'Mut','Wild')
write.table(all,file='simceredx_3273s_clinical_info_addSMARCA4.xls',sep='\t',quote=FALSE)

colon = read.maf(maf = 'simceredx_somatic_mutation_Func_3273s.maf', cnTable='simceredx_3273s_somatic_CNV_for_landscape.txt',clinicalData = 'simceredx_3273s_clinical_info_addSMARCA4_addTMB_addMSI.xls')
colon = read.maf(maf = 'simceredx_somatic_mutation_Func_3273s.maf',clinicalData = 'simceredx_3273s_clinical_info_addSMARCA4_addTMB_addMSI.xls')

write.mafSummary(maf = colon, basename = 'Coclon_3273s')

wildmaf <- read.maf(maf = 'simceredx_3050s.SMARCA4.Wild_somatic_mutation_Func.maf',clinicalData = 'simceredx_clinical_info_3050s.SMARCA4.Wild.xls',cnTable='simceredx_3050s.SMARCA4.Wild_somatic_CNV_for_landscape.txt' )
write.mafSummary(maf = wildmaf, basename = 'Coclon.SMARCA4.Wild_3050s')

smarca4 = read.maf(maf = 'simceredx_223s.SMARCA4.Mut_somatic_mutation_Func.maf',clinicalData = 'simceredx_SMARCA4.Mut.223s_clinical_info_addSMARCA4_addTMB_addMSI.xls',cnTable='simceredx_223s.SMARCA4.Mut_somatic_CNV_for_landscape.txt')

smarca4 = read.maf(maf = 'simceredx_223s.SMARCA4.Mut_somatic_mutation_Func.maf',clinicalData = 'simceredx_SMARCA4.Mut.223s_clinical_info_addSMARCA4_addTMB_addMSI.xls')
write.mafSummary(maf = smarca4, basename = 'Coclon.SMARCA4.Mut_223s')
########################################################################## SWI
pdf('simceredx_3273s_SWI_SNF_landscape.pdf',height=5,width=30)
oncoplot(maf = colon, altered =TRUE,colors = col,genes=swi_gene_list,bgCol = "#F8F8F7",legend_height=1.5,anno_height=1.0,barcode_mar=10,logColBar=TRUE,fontSize = 1.5,SampleNamefontSize = 1.8,titleFontSize = 1.8,legendFontSize = 1.8,annotationFontSize = 1.8,removeNonMutated=TRUE,fill=TRUE,gene_mar=20,drawRowBar=TRUE,drawColBar=F,showTumorSampleBarcodes = FALSE)
dev.off()


pdf('MSK_5543s_SWI_SNF_landscape.pdf',height=5,width=30)
oncoplot(maf = msk, altered =TRUE,colors = col,genes=swi_gene_list,bgCol = "#F8F8F7",legend_height=1.5,anno_height=1.0,barcode_mar=10,logColBar=TRUE,fontSize = 1.5,SampleNamefontSize = 1.8,titleFontSize = 1.8,legendFontSize = 1.8,annotationFontSize = 1.8,removeNonMutated=TRUE,fill=TRUE,gene_mar=20,drawRowBar=TRUE,drawColBar=F,showTumorSampleBarcodes = FALSE)
dev.off()

sim_swi <-subsetMaf(maf=colon,genes=swi_gene_list)
msk_swi <-subsetMaf(maf=msk,genes=swi_gene_list)

comp <-mafCompare(m1 =colon , m2 =msk, m1Name = 'Our Cohort', m2Name = 'MSK cohort', minMut = 1)
write.table(comp$results,file='simceredx.3273s_vs_MSK.5543s_somatic_SWI.SNF_GeneFreq.xls',sep='\t',quote=FALSE,row.names = FALSE)

pdf('simceredx.3273s_vs_MSK.5543s_somatic_SWI.SNF_Barplot.pdf',height=5,width=15)
coBarplot(m1 =colon , m2 =msk , m1Name = 'Our Cohort', m2Name = 'MSK cohort',genes=rev(swi_gene_list),yLims =c(20,20))
dev.off()
########################################################################## MSK cohort
setwd('/home/data/t020412/Qinzhang_data/SMARCA4_Data/2.msk_chord_2024/')

msk = read.maf(maf = 'msk_chord_2024_mutation_Func_colon_5543s.maf',clinicalData = 'msk_chord_2024_clinical_data_colon_5543s.txt')

msk_smarca4 = read.maf(maf = 'msk_chord_2024_mutation_Func_colon_SMARCA4_Mut.324s.maf',clinicalData = 'msk_chord_2024_clinical_data_colon_SMARCA4.mut.324s.txt')
write.mafSummary(maf = msk_smarca4, basename = 'MSK.SMARCA4.Mut_324s')

pdf('MSK.CHORD_SMARCA4.Mut.324s_lollipopPlot.pdf',height=3,width=12)
lollipopPlot(
  maf = msk_smarca4,
  gene = 'SMARCA4',
  AACol = 'HGVSp',
  showMutationRate = TRUE,
  showDomainLabel=TRUE,labelPos = c(271,910,859)
)
dev.off()

##########################################################################

clin <- read.table('/home/data/t020412/Qinzhang_data/SMARCA4_Data/1.simceredx_data/simceredx_3273s_clinical_info_addSMARCA4_addTMB_addMSI.xls', sep='\t',header=T)

clin$Agegroup <- ifelse(clin$Age>=65,">=65","<65")
write.table(clin,file = "simceredx_3273s_clinical_info_addSMARCA4_addTMB_addMSI_addAgeGroup.xls", col.names = T, row.names = F, quote = F, sep = "\t")


allvars <- c("Age","Agegroup","Gender","CancerType","msi_status","Stage")
fvars <- c("Agegroup","Gender","CancerType","msi_status","Stage")

strata<-"Group"

table1<-scitb1(vars=allvars,fvars=fvars,strata=strata,data=clin)

write.table(table1,file = "Table1_all.3273s_SMARCA4.Mut.223_vs_wild.3050.xls", col.names = T, row.names = F, quote = F, sep = "\t")
####################################################################################### TMB 
clin_sim_TMB <- subset(clin,clin$tmb!="NA")
data <-as.data.frame(clin_sim_TMB)
data$Group <-factor(data$Group,levels=c("Mut","Wild"))
#data$tmb <- log(data$tmb+1,2)

pdf('simceredx_3273s_SMARCA4_mut.vs.wild_TMB_boxplot.pdf',width=3,height=5)
p <- ggboxplot(data, x = "Group", y = "tmb",color = "Group",palette = c("#ED0000","#A58AFF"),add = "jitter",ylab = "TMB")+stat_compare_means()
p
dev.off()

all_msk <-read.table('msk_chord_2024_clinical_data_colon_5543s.txt',sep='\t',header=T)
msk_mut<-read.table('msk_chord_2024_clinical_data_colon_SMARCA4.mut.324s.txt',sep='\t',header=T)
all_msk <-as.data.frame(all_msk)
all_msk$Group <-ifelse(all_msk$Tumor_Sample_Barcode %in% msk_mut$Tumor_Sample_Barcode,'Mut','Wild')
write.table(all_msk,file = "msk_chord_2024_clinical_info_addSMARCA4_addTMB.xls", col.names = T, row.names = F, quote = F, sep = "\t")

all_msk$Group <-factor(all_msk$Group,levels=c("Mut","Wild"))
pdf('MSK_5543s_SMARCA4_mut.vs.wild_TMB_boxplot.pdf',width=3,height=5)
p <- ggboxplot(all_msk, x = "Group", y = "TMB.nonsynonymous",color = "Group",palette = c("#ED0000","#A58AFF"),add = "jitter",ylab = "TMB")+stat_compare_means()
p
dev.off()


#################################################  Sim vs MSK
setwd('/home/data/t020412/Qinzhang_data/SMARCA4_Data/3.chinese_vs_western/')
clin <- read.table('chinese_vs_western_SMARCA4.Mut_clinical.xls', sep='\t',header=T)
clin$Agegroup <- ifelse(clin$Age>=65,">=65","<65")

allvars <- c("Age","Agegroup","Gender","CancerType","msi_status","Stage")
fvars <- c("Agegroup","Gender","CancerType","msi_status","Stage")

strata<-"Population"

table1<-scitb1(vars=allvars,fvars=fvars,strata=strata,data=clin)

write.table(table1,file = "Table2_SMARCA4.Mut_Chinese.223_vs_Western.324s.xls", col.names = T, row.names = F, quote = F, sep = "\t")

##################################################### 
data <- read.table('chinese_Muttype.xls', sep='\t',header=T)

pdf('simceredx_SMARCA4.Mut.223s_MutType_pie.pdf',height=3,width=5)
ggplot(data, aes(x = "", y = summary, fill = ID)) + 
  geom_col(color = "black") + 
  geom_text(aes(label = paste0(Freq, "%")), position = position_stack(vjust = 0.5), size = 4) + 
  scale_fill_brewer(palette = "Paired") + 
  coord_polar("y", start = 0)  + 
  theme_void()
dev.off()

data <- read.table('msk_Muttype.xls', sep='\t',header=T)

pdf('MSK_SMARCA4.Mut.324s_MutType_pie.pdf',height=3,width=5)
ggplot(data, aes(x = "", y = summary, fill = ID)) + 
  geom_col(color = "black") + 
  geom_text(aes(label = paste0(Freq, "%")), position = position_stack(vjust = 0.5), size = 4) + 
  scale_fill_brewer(palette = "Paired") + 
  coord_polar("y", start = 0)  + 
  theme_void()
dev.off()

###########################################################
col = brewer.pal(11,'Paired')[c(1:8,10,11,9)]
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins', 'Amplification','In_Frame_Ins', 'Splice_Site','In_Frame_Del','Deletion','Translation_Start_Site')

Sexcol = brewer.pal(11,'RdBu')[c(4,8)]
names(Sexcol) = c('Female','Male')

msi_col= brewer.pal(8,'Accent')[c(7,5,1)]
names(msi_col)= c('MSS','MSI-H','MSI-L')

age_col=brewer.pal(8,'BrBG')[c(1,4)]
names(age_col)=c('<65','>=65')

annocolors = list(msi_status=msi_col,Agegroup=age_col,Gender=Sexcol)

"#F8F8F7"
"#CCCCCC"
#########################################################################draw picture
pdf('simceredx_SMARCA4.Mut.223s_landscape_top30.pdf',height=15,width=30)
oncoplot(maf = smarca4,top=30, altered =TRUE,colors = col,clinicalFeatures = c('Agegroup','Gender','msi_status'),annoBorderCol='white',bgCol = "#F8F8F7",legend_height=1.5,anno_height=1.0,annotationColor = annocolors,barcode_mar=10,logColBar=TRUE,fontSize = 1.5,SampleNamefontSize = 1.8,titleFontSize = 1.8,legendFontSize = 1.8,annotationFontSize = 1.8,removeNonMutated=FALSE,fill=TRUE,keepGeneOrder = TRUE,sortByMutation =FALSE,GeneOrderSort=TRUE,sortByAnnotation = TRUE,gene_mar=20,drawRowBar=TRUE,drawColBar=TRUE,showTumorSampleBarcodes = FALSE)
dev.off()

#############
pdf('simceredx_SMARCA4.Mut.223s_lollipopPlot.pdf',height=3,width=12)
lollipopPlot(
  maf = smarca4,
  gene = 'SMARCA4',
  AACol = 'HGVSp',
  showMutationRate = TRUE,
  showDomainLabel=TRUE,labelPos = c(381,672)
)
dev.off()
################################################### simceredx somatic
res <-somaticInteractions(maf = smarca4, top = 25, pvalue = c(0.05, 0.1))
write.table(res,file='simceredx_SMARCA4.Mut.223s_Somatic_interactions.xls',sep='\t',quote=FALSE)

pdf('simceredx_SMARCA4.Mut.223s_Somatic_interactions.pdf',height=6,width=6)
somaticInteractions(maf = smarca4, top = 25, pvalue = c(0.05, 0.01),leftMar = 6,topMar = 6)
dev.off()


res <-somaticInteractions(maf = msk_smarca4, top = 25, pvalue = c(0.05, 0.1))
write.table(res,file='MSK_SMARCA4.Mut.324s_Somatic_interactions.xls',sep='\t',quote=FALSE)

pdf('MSK_SMARCA4.Mut.324s_Somatic_interactions.pdf',height=6,width=6)
somaticInteractions(maf = msk_smarca4, top = 25, pvalue = c(0.05, 0.01),leftMar = 6,topMar = 6)
dev.off()
###################################################pathway
###################################################################################
##################################################################################### panel gene list overlap
p551 <-read.table("panel551_gene_list.xls",header = F)

p341 <-read.table("Impact_341_gene_list.xls",header = F)
p410 <-read.table("Impact_410_gene_list.xls",header = F)
p468 <-read.table("Impact_468_gene_list.xls",header = F)

intersect2 <- function(...) Reduce(intersect, list(...))
overlap_gene<-intersect2(p551, p341, p410, p468)

write.table(overlap_gene,file='overlap_551_341_410_468.xls',sep='\t',quote=FALSE,row.names = FALSE)

gene_list <-read.table('overlap_551_341_410_468.xls',sep='\t',header=T)
gene_list <-as.list(gene_list$Gene)
######################################################################################

smarca4_overlap <-subsetMaf(maf=smarca4,genes=gene_list)
msk_smarca4_overlap <-subsetMaf(maf=msk_smarca4,genes=gene_list)

comp <-mafCompare(m1 =smarca4_overlap , m2 =msk_smarca4_overlap , m1Name = 'Our Cohort', m2Name = 'MSK cohort', minMut = 1)
write.table(comp$results,file='simceredx.223s_vs_MSK.324s_somatic_WES_GeneFreq.xls',sep='\t',quote=FALSE,row.names = FALSE)

pdf('simceredx.223s_vs_MSK.324s_somatic_WES_Barplot_top30_signif.gene.pdf',height=10,width=15)
coBarplot(m1 =smarca4 , m2 =msk_smarca4  , m1Name = 'Our Cohort', m2Name = 'MSK cohort',genes=c("NPM1","NRAS","IFNGR1","RAD51B","MAPK1","IL7R","REL","EZH2","SDHA","INPP4B","GNAQ","FUBP1","JAK2","PMS1","CDC73","STAG2","XPO1","RPTOR","MET","KDM6A","RB1","CTNNB1","PIK3C3","RNF43","BLM","ATR","RAD50","GNAS","CIC","EP300"))
dev.off()
