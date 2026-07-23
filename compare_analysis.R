################################################################# SMARCA4 突变人群 拆分 MSI-H 和MSS
smarca4_mut_clin <-subset(all, all$Group=="Mut")

smarca4_mut_clin_msih <-subset(smarca4_mut_clin, smarca4_mut_clin$msi_status=="MSI-H")
smarca4_mut_clin_mss <-subset(smarca4_mut_clin, smarca4_mut_clin$msi_status=="MSS")

smarca4_maf_msih <-subsetMaf(smarca4,tsb=smarca4_mut_clin_msih$Tumor_Sample_Barcode)
smarca4_maf_mss <-subsetMaf(smarca4,tsb=smarca4_mut_clin_mss$Tumor_Sample_Barcode)

################################################################## MSK 拆分 MSI-H 和MSS
msk_smarca4_mut_clin_msih <- subset(msk_smarca4_clin, msk_smarca4_clin$MSI=="MSI-H")
msk_smarca4_mut_clin_mss <- subset(msk_smarca4_clin, msk_smarca4_clin$MSI=="MSS")

msk_smarca4_mut_maf_msih <-subsetMaf(msk_smarca4,tsb=msk_smarca4_mut_clin_msih$Tumor_Sample_Barcode)
msk_smarca4_mut_maf_mss <-subsetMaf(msk_smarca4,tsb=msk_smarca4_mut_clin_mss$Tumor_Sample_Barcode)

#####################################################################  比较 figure1B SMARCA4 Mut人群，MSI-h 的突变频率比较
smarca4_maf_mss_overlap <-subsetMaf(maf=smarca4_maf_msih,genes=gene_list)
msk_smarca4_mut_maf_mss_overlap <-subsetMaf(maf=msk_smarca4_mut_maf_msih,genes=gene_list)

comp <-mafCompare(m1 =smarca4_maf_msih_overlap , m2 =msk_smarca4_mut_maf_msih_overlap , m1Name = 'Our Cohort MSI-H SMARCA4-Mut', m2Name = 'MSK cohort MSI-H SMARCA4-Mut', minMut = 1)

res <-comp$results
write.xlsx(res, "Simceredx.MSI-H.SMARCA4-Mut.106s_vs_MSK.MSI-H.SMARCA4-Mut.199s_somatic_panel.overlap_all.compare.xlsx", rowNames = FALSE)

# 筛选显著基因，排序
res_sig <- res[res$pval < 0.05 & res$adjPval < 0.05, ]
#res_sig <- res_sig[order(res_sig$or), ]
res_sig <- res_sig[order(res_sig$`Our Cohort MSI-H SMARCA4-Mut`, decreasing = TRUE), ]
res_sig$gene <- factor(res_sig$Hugo_Symbol, levels = (res_sig$Hugo_Symbol))


top_30 <-res_sig[c(1:30),]
top_30 <- top_30 [order(top_30$`Our Cohort MSI-H SMARCA4-Mut`, decreasing = FALSE), ]

pdf('Simceredx.MSI-H.SMARCA4-Mut.106s_vs_MSK.MSI-H.SMARCA4-Mut.199s_Barplot_top30_signif.gene.pdf',height=10,width=15)
coBarplot(m1 =smarca4_maf_msih_overlap, m2 =msk_smarca4_mut_maf_msih_overlap, m1Name = 'Our Cohort MSI-H SMARCA4-Mut', m2Name = 'MSK cohort MSI-H SMARCA4-Mut',genes=top_30$Hugo_Symbol,geneMar=8)
dev.off()
##################################################################### 比较 figure1B SMARCA4 Mut人群，MSS 的突变频率比较
smarca4_maf_mss_overlap <-subsetMaf(maf=smarca4_maf_mss,genes=gene_list)
msk_smarca4_mut_maf_mss_overlap <-subsetMaf(maf=msk_smarca4_mut_maf_mss,genes=gene_list)

comp <-mafCompare(m1 =smarca4_maf_mss_overlap , m2 =msk_smarca4_mut_maf_mss_overlap , m1Name = 'Our Cohort MSS SMARCA4-Mut', m2Name = 'MSK cohort MSS SMARCA4-Mut', minMut = 1)

res <-comp$results
write.xlsx(res, "Simceredx.MSS.SMARCA4-Mut.117s_vs_MSK.MSS.SMARCA4-Mut.101s_somatic_panel.overlap_all.compare.xlsx", rowNames = FALSE)

# 筛选显著基因，排序
res_sig <- res[res$pval < 0.05 , ]
#res_sig <- res_sig[order(res_sig$or), ]
res_sig <- res_sig[order(res_sig$`Our Cohort MSS SMARCA4-Mut`, decreasing = TRUE), ]
res_sig$gene <- factor(res_sig$Hugo_Symbol, levels = (res_sig$Hugo_Symbol))


top_30 <-res_sig[c(1:30),]
top_30 <- top_30 [order(top_30$`Our Cohort MSS SMARCA4-Mut`, decreasing = FALSE), ]

pdf('Simceredx.MSS.SMARCA4-Mut.117s_vs_MSK.MSS.SMARCA4-Mut.101s_Barplot_top30_signif.gene.pdf',height=10,width=15)
coBarplot(m1 =smarca4_maf_mss_overlap , m2 =msk_smarca4_mut_maf_mss_overlap , m1Name = 'Our Cohort MSS SMARCA4-Mut', m2Name = 'MSK cohort MSS SMARCA4-Mut',genes=top_30$Hugo_Symbol,geneMar=8)
dev.off()
#################################################################### 比较 figure3c MSI-h 的SWI/SNF突变频率比较

################################################################# all人群 拆分 MSI-H 和MSS
sim_msih <-subset(all,all$msi_status=="MSI-H")
sim_mss <-subset(all,all$msi_status=="MSS")

msk_msih <- subset(all_msk, all_msk$MSI=="MSI-H")
msk_mss <- subset(all_msk, all_msk$MSI=="MSS")

sim_msih_maf_swi <-subsetMaf(maf=colon,genes=swi_gene_list,tsb =sim_msih$Tumor_Sample_Barcode )
sim_mss_maf_swi <-subsetMaf(maf=colon,genes=swi_gene_list,tsb =sim_mss$Tumor_Sample_Barcode)

write.mafSummary(maf = sim_msih_maf_swi, basename = 'sim_msih')
write.mafSummary(maf = sim_mss_maf_swi, basename = 'sim_mss')

msk_msih_maf_swi <-subsetMaf(maf=msk,genes=swi_gene_list,tsb =msk_msih$Tumor_Sample_Barcode)
msk_mss_maf_swi <-subsetMaf(maf=msk,genes=swi_gene_list,tsb =msk_mss$Tumor_Sample_Barcode)

write.mafSummary(maf = msk_msih_maf_swi, basename = 'msk_msih')
write.mafSummary(maf = msk_mss_maf_swi, basename = 'msk_mss')

sim_msih_maf_swi <-read.maf(maf='sim_msih_maftools.maf',clinicalData = sim_msih)
sim_mss_maf_swi <-read.maf(maf='sim_mss_maftools.maf',clinicalData = sim_mss)

msk_msih_maf_swi <-read.maf(maf='msk_msih_maftools.maf',clinicalData = msk_msih)
msk_mss_maf_swi <-read.maf(maf='msk_mss_maftools.maf',clinicalData = msk_mss)
###########################################################################################  MSI-h 的SWI/SNF突变频率比较

comp <-mafCompare(m1 =sim_msih_maf_swi , m2 =msk_msih_maf_swi, m1Name = 'Our Cohort MSI-H SWI/SNF', m2Name = 'MSK cohort MSI-H SWI/SNF', minMut = 1)

res <-comp$results

total_sim <- 219
total_msk <- 586

# 3. 计算每个基因在各自队列中未发生突变的样本数

res$m1 <- res$`Our Cohort MSI-H SWI/SNF`
res$m2 <-res$`MSK cohort MSI-H SWI/SNF`

res$m1_no_mut <- total_sim - res$`Our Cohort MSI-H SWI/SNF`
res$m2_no_mut <- total_msk - res$`MSK cohort MSI-H SWI/SNF`

# 4. 对每个基因重新进行 Fisher 精确检验
#    使用 apply 函数逐行处理
res$pval_corrected <- apply(res, 1, function(row) {
  contingency_table <- matrix(
    c(as.numeric(row['Our Cohort MSI-H SWI/SNF']),  as.numeric(row['m1_no_mut']),
      as.numeric(row['MSK cohort MSI-H SWI/SNF']),  as.numeric(row['m2_no_mut'])),
    nrow = 2
  )
  fisher.test(contingency_table)$p.value
})

# 5. 对校正后的 p 值进行多重检验校正 (例如 Benjamini-Hochberg)
res$adjPval_corrected <- p.adjust(res$pval_corrected, method = "BH")

write.xlsx(res, "Simceredx.MSI-H.SWI.SNF.219s_vs_MSK.MSI-H.SWI.SNF.586s_somatic_panel.overlap_all.compare.xlsx", rowNames = FALSE)

# 1. 计算突变频率（百分比）
res <- res %>%
  mutate(
    freq_sim = (m1 / total_sim) * 100,
    freq_msk = (m2 / total_msk) * 100,
    # 添加显著性标记（示例：adjP < 0.05 标记为 "*"）
    sig_label = ifelse(adjPval_corrected < 0.001, "***", ifelse(adjPval_corrected < 0.01,"**",ifelse(adjPval_corrected < 0.05,"*","ns")))
  )

res_sig <- res
# 3. 将数据转换为长格式，方便ggplot分组
# 2. 准备绘图数据（双向条形）
plot_data_bi <- res_sig %>%
  select(Hugo_Symbol, freq_sim, freq_msk, sig_label) %>%
  pivot_longer(cols = c(freq_sim, freq_msk),
               names_to = "Cohort",
               values_to = "Frequency") %>%
  mutate(
    Frequency_plot = ifelse(Cohort == "freq_sim", Frequency, -Frequency),  # Our 为正，MSK 为负
    Cohort_label = ifelse(Cohort == "freq_sim", "Our Cohort", "MSK Cohort"),
    Hugo_Symbol = factor(Hugo_Symbol, levels = rev(res_sig$Hugo_Symbol)),
    # 用于显示的数字标签（绝对值）
    label_display = paste0(round(Frequency, 1), "%")
  )

# 3. 计算 y 轴范围及标签位置
y_max <- max(plot_data_bi$Frequency, na.rm = TRUE)
y_min <- max(plot_data_bi$Frequency, na.rm = TRUE)
y_label_pos <- y_max + 25          # Our 标签在正方向最高处上方
y_msk_label_pos <- -y_min - 10      # MSK 标签在负方向最低处下方

# 显著性标记（星号放在 Our 柱子顶端上方）
sig_annot <- res_sig %>%
  mutate(y_pos = y_label_pos - 18)

# 4. 绘制双向条形图
p <- ggplot(plot_data_bi, aes(x = Hugo_Symbol, y = Frequency_plot, fill = Cohort_label)) +
  geom_col(position = "identity", width = 0.7) +   # ① 对齐，不错开
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_y_continuous(labels = abs, breaks = seq(-80, 80, by = 10)) +
  labs(title = "Mutation Frequency Comparison (including negatives)",
       x = "Gene", y = "Mutation Frequency (%)", fill = "Cohort") +
  scale_fill_manual(values = c("Our Cohort" = "#E41A1C", "MSK Cohort" = "#377EB8")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        legend.position = "bottom") +
  # 组别标签（② MSK 标签置于 y 轴最下端）
  annotate("text", x = -Inf, y = y_label_pos, label = "Our Cohort", 
           hjust = -0.2, vjust = 1.5, size = 4) +
  annotate("text", x = -Inf, y = y_msk_label_pos, label = "MSK Cohort", 
           hjust = -0.2, vjust = -0.5, size = 4)

# ③ 添加频率数字（正数在上方，负数在下方）
p <- p + 
  geom_text(data = plot_data_bi,
            aes(x = Hugo_Symbol, y = Frequency_plot, label = label_display,
                vjust = ifelse(Frequency_plot > 0, -0.5, 1.5)),  # 正数上移，负数下移
            position = "identity", size = 3, color = "black") +
  # 添加显著性星号（可选）
  geom_text(data = sig_annot, aes(x = Hugo_Symbol, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 4, vjust = 0, color = "red")

# 6. 保存
ggsave("Simceredx.MSI-H.SWI.SNF.219s_vs_MSK.MSI-H.SWI.SNF.586s_Frequency.pdf", plot = p, width = 6, height = 7)

##############################################################  MSS 的SWI/SNF突变频率比较

comp <-mafCompare(m1 =sim_mss_maf_swi , m2 =msk_mss_maf_swi, m1Name = 'Our Cohort MSS SWI/SNF', m2Name = 'MSK cohort MSS SWI/SNF', minMut = 1)

res <-comp$results

total_sim <- 3049
total_msk <- 4682

# 3. 计算每个基因在各自队列中未发生突变的样本数

res$m1 <- res$`Our Cohort MSS SWI/SNF`
res$m2 <-res$`MSK cohort MSS SWI/SNF`

res$m1_no_mut <- total_sim - res$`Our Cohort MSS SWI/SNF`
res$m2_no_mut <- total_msk - res$`MSK cohort MSS SWI/SNF`

# 4. 对每个基因重新进行 Fisher 精确检验
#    使用 apply 函数逐行处理
res$pval_corrected <- apply(res, 1, function(row) {
  contingency_table <- matrix(
    c(as.numeric(row['Our Cohort MSS SWI/SNF']),  as.numeric(row['m1_no_mut']),
      as.numeric(row['MSK cohort MSS SWI/SNF']),  as.numeric(row['m2_no_mut'])),
    nrow = 2
  )
  fisher.test(contingency_table)$p.value
})

# 5. 对校正后的 p 值进行多重检验校正 (例如 Benjamini-Hochberg)
res$adjPval_corrected <- p.adjust(res$pval_corrected, method = "BH")

write.xlsx(res, "Simceredx.MSS.SWI.SNF.3049s_vs_MSK.MSS.SWI.SNF.4682s_somatic_panel.overlap_all.compare.xlsx", rowNames = FALSE)

# 1. 计算突变频率（百分比）
res <- res %>%
  mutate(
    freq_sim = (m1 / total_sim) * 100,
    freq_msk = (m2 / total_msk) * 100,
    # 添加显著性标记（示例：adjP < 0.05 标记为 "*"）
    sig_label = ifelse(adjPval_corrected < 0.001, "***", ifelse(adjPval_corrected < 0.01,"**",ifelse(adjPval_corrected < 0.05,"*","ns")))
  )

res_sig <- res
# 3. 将数据转换为长格式，方便ggplot分组
# 2. 准备绘图数据（双向条形）
plot_data_bi <- res_sig %>%
  select(Hugo_Symbol, freq_sim, freq_msk, sig_label) %>%
  pivot_longer(cols = c(freq_sim, freq_msk),
               names_to = "Cohort",
               values_to = "Frequency") %>%
  mutate(
    Frequency_plot = ifelse(Cohort == "freq_sim", Frequency, -Frequency),  # Our 为正，MSK 为负
    Cohort_label = ifelse(Cohort == "freq_sim", "Our Cohort", "MSK Cohort"),
    Hugo_Symbol = factor(Hugo_Symbol, levels = rev(res_sig$Hugo_Symbol)),
    # 用于显示的数字标签（绝对值）
    label_display = paste0(round(Frequency, 1), "%")
  )

# 3. 计算 y 轴范围及标签位置
y_max <- max(plot_data_bi$Frequency, na.rm = TRUE)
y_min <- max(plot_data_bi$Frequency, na.rm = TRUE)
y_label_pos <- y_max + 10          # Our 标签在正方向最高处上方
y_msk_label_pos <- -y_min - 5      # MSK 标签在负方向最低处下方

# 显著性标记（星号放在 Our 柱子顶端上方）
sig_annot <- res_sig %>%
  mutate(y_pos = y_label_pos - 5)

# 4. 绘制双向条形图
p <- ggplot(plot_data_bi, aes(x = Hugo_Symbol, y = Frequency_plot, fill = Cohort_label)) +
  geom_col(position = "identity", width = 0.7) +   # ① 对齐，不错开
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_y_continuous(labels = abs, breaks = seq(-80, 80, by = 10)) +
  labs(title = "Mutation Frequency Comparison (including negatives)",
       x = "Gene", y = "Mutation Frequency (%)", fill = "Cohort") +
  scale_fill_manual(values = c("Our Cohort" = "#E41A1C", "MSK Cohort" = "#377EB8")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        legend.position = "bottom") +
  # 组别标签（② MSK 标签置于 y 轴最下端）
  annotate("text", x = -Inf, y = y_label_pos, label = "Our Cohort MSS", 
           hjust = -0.2, vjust = 1.5, size = 4) +
  annotate("text", x = -Inf, y = y_msk_label_pos, label = "MSK Cohort MSS", 
           hjust = -0.2, vjust = -0.5, size = 4)

# ③ 添加频率数字（正数在上方，负数在下方）
p <- p + 
  geom_text(data = plot_data_bi,
            aes(x = Hugo_Symbol, y = Frequency_plot, label = label_display,
                vjust = ifelse(Frequency_plot > 0, -0.5, 1.5)),  # 正数上移，负数下移
            position = "identity", size = 3, color = "black") +
  # 添加显著性星号（可选）
  geom_text(data = sig_annot, aes(x = Hugo_Symbol, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 4, vjust = 0, color = "red")

# 6. 保存
ggsave("Simceredx.MSS.SWI.SNF.3049s_vs_MSK.MSS.SWI.SNF.4682s_Frequency.pdf", plot = p, width = 6, height = 7)
