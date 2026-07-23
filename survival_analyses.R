library("survival")
library("survminer")
setwd('./SMARCA4_Data/6.MSK_survival/')
msk_clin <- read.table('msk_chord_2024_clinical_info_addSMARCA4_addTMB.xls', sep='\t',header=T)
msk_clin$Group <- factor(msk_clin$Group,levels=c("Wild","Mut"))

###################### MSI
msk_clin <- msk_clin %>% 
  mutate_if(is.character, ~ ifelse(is.na(.), "unKnown", .))

msk_clin$MSI <- ifelse(msk_clin$MSI.Type=="Instable","MSI-H",ifelse(msk_clin$MSI.Type=="Stable","MSS",(ifelse(msk_clin$MSI.Type=='Do not report',"unKnown",ifelse(msk_clin$MSI.Type=="Indeterminate","MSI-L","unKnown")))))

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msk_clin)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msk_clin)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=msk_clin,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 100, y = 0.9,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 100, y = 0.8,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 100, y = 0.7,label = paste("P","<0.001"))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)


fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Current.Age+Sex+Stage+Smoking.History+Group,data = msk_clin)
p<-ggforest(fit1,data=msk_clin,main="Hazard ratio",fontsize=1.1)
ggsave(p,file="MSK_SMARCA4_Mut_vs_Wild_OS_ForestPlot_MutiCox.pdf",width=12,height = 9)

####################################### MSKCC cohort
fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Current.Age+Sex+Stage+Smoking.History+MSI+Group,data = msk_clin)
p<-ggforest(fit1,data=msk_clin,main="Hazard ratio",fontsize=1.1)
ggsave(p,file="MSK_SMARCA4_Mut_vs_Wild_OS_ForestPlot_MutiCox_addMSI.pdf",width=12,height = 9)

###############################################################################################################
setwd('/home/data/t020412/Qinzhang_data/SMARCA4_Data/3.chinese_vs_western/')
data <-read.table("chinese_vs_western_SMARCA4.Mut_Freq.xls",sep='\t',header=T)

cols=c("#F08080","#87CEEB")
pdf("chinese_cohort_VS_MSKCC_smarca4.freq.pdf",height = 4,width = 5)
p <- ggplot(data, aes(x=Group, y=Frequency, fill=Group))+geom_bar(position=position_dodge(), stat="identity",color="white")+scale_fill_manual(values=cols)+ylim(0,10)+theme_bw()+labs(y="SMARCA4 mutation frequency(%)")
p
dev.off()

########################################## 
chemo_data <-read.table("MSK_CHORD_CRC_Chemo_total4510pts.txt",sep='\t',header=T)

chemo <-subset(msk_clin,msk_clin$Tumor_Sample_Barcode %in% chemo_data$Sample.ID)
chemo$Group <- factor(chemo$Group,levels=c("Wild","Mut"))

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = chemo)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = chemo)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=chemo,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),##
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 100, y = 0.9,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 100, y = 0.8,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 100, y = 0.7,label = paste("P","<0.001"))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="Chemo_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)


########################################## 
immune_data <-read.table("MSK_CHORD_CRC_immune_total226.txt",sep='\t',header=T)

immune <-subset(msk_clin,msk_clin$Tumor_Sample_Barcode %in% immune_data$Sample.ID)
immune$Group <- factor(immune$Group,levels=c("Wild","Mut"))

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = immune)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = immune)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=immune,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),##
                        surv.median.line = "hv",#
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 80, y = 0.9,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 80, y = 0.8,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 80, y = 0.7,label = paste("P=",pvalue))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="immune_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)


###############################################################################################################################

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Current.Age+Sex+Stage+Smoking.History+Group,data = msk_clin)
p<-ggforest(fit1,data=msk_clin,main="Hazard ratio",fontsize=1.1)
ggsave(p,file="MSK_SMARCA4_Mut_vs_Wild_OS_ForestPlot_MutiCox.pdf",width=12,height = 9)

################################################################################################################# 

mss_chemo <- subset(chemo,chemo$MSI != "MSI-H")
msih_chemo <- subset(chemo,chemo$MSI == "MSI-H")

######################################################################################################## MSS

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = mss_chemo)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = mss_chemo)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=mss_chemo,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 80, y = 0.9,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 80, y = 0.8,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 80, y = 0.7,label = paste("P=",pvalue))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="MSS_Chemo_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)

#########################################################################################  MSI-H 

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msih_chemo)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msih_chemo)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=msih_chemo,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 85, y = 1.0,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 85, y = 0.9,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 85, y = 0.8,label = paste("P=",pvalue))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="MSI-H_Chemo_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)

################################################################################################################# 

mss_all <- subset(msk_clin,msk_clin$MSI != "MSI-H")
msih_all <- subset(msk_clin,msk_clin$MSI == "MSI-H")

######################################################################################################## MSS 

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = mss_all)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = mss_all)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=mss_all,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),##
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 80, y = 0.9,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 80, y = 0.8,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 80, y = 0.7,label = paste("P=",pvalue))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="MSS_All_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)

#########################################################################################  MSI-H

fit <- surv_fit(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msih_all)

fit1 <- coxph(Surv(Overall.Survival, Overall.Survival.Status)~ Group,data = msih_all)
sum <- summary(fit1)
CI <- paste0(round(sum$conf.int[,3:4],3),collapse='-')
HR <- round(sum$coefficients[,2],3)
pvalue <- round(sum$sctest[3][[1]],3)

splots <- list()
splots[[1]]<-ggsurvplot(fit,data=msih_all,
                        xlab = "Time(Months)",
                        ylab="Overall Survival",
                        pval = F,
                        conf.int= F,
                        risk.table = T,
                        legend.title = "OS",
                        #legend.labs = levels(Data$v),
                        surv.median.line = "hv",
                        palette="lancet")

splots[[1]]$plot = splots[[1]]$plot + ggplot2::annotate("text",x = 85, y = 1.0,label = paste("HR :",round(HR,3))) + ggplot2::annotate("text",x = 85, y = 0.9,label = paste("(","95%CI:",CI,")",sep = ""))+ggplot2::annotate("text",x = 85, y = 0.8,label = paste("P=",pvalue))
res<-arrange_ggsurvplots(splots, print = F,ncol = 1, nrow = 1, risk.table.height = 0.25)

ggsave(file="MSI-H_All_MSK_SMARCA4_Mut_vs_Wild_OS.pdf",res,width=8,height = 6)
