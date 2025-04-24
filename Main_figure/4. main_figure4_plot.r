load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v5.RData")
ARGO_sdiused_label2[ARGO_sdiused_label2$RNAseqID == "14TS0Z0425TR1", "rfs.event"] <- 0
ARGO_sdiused_label2[ARGO_sdiused_label2$rfs.delay > 80, "rfs.delay"] <- 80
MidStageSamples <- ARGO_sdiused_label2
#CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4')
##
CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4')
MidStageSamples <- ARGO_sdiused_label2
#MidStageSamples <- ARGO_sdiused_label2[grep("II", ARGO_sdiused_label2$TNM分期), ] ## 274
MidStageSamples$chemotherapy <- factor(MidStageSamples$chemotherapy, levels=c("Surgery alone", "Adjuvant chemotherapy"))

### CMSEMT
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamples, "rfs.delay", "rfs.event", "Immune.Subtype", Color=CMSTMESubtype, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

########################################################## IF
MidStageSamplesCMS1 <- MidStageSamples[MidStageSamples$Immune.Subtype == "CMS1" & !is.na(MidStageSamples$chemotherapy), ]  ## 81 32
MidStageSamplesCMS2 <- MidStageSamples[MidStageSamples$Immune.Subtype == "CMS2" & !is.na(MidStageSamples$chemotherapy), ]  ## 152
MidStageSamplesCMS3 <- MidStageSamples[MidStageSamples$Immune.Subtype == "CMS3" & !is.na(MidStageSamples$chemotherapy), ]  ## 54
MidStageSamplesCMS4_HighIF <- MidStageSamples[MidStageSamples$Immune.Subtype == "CMS4-IF+" & !is.na(MidStageSamples$chemotherapy), ]  ## 74
MidStageSamplesCMS4_LowIF <- MidStageSamples[MidStageSamples$Immune.Subtype == "CMS4-IF-" & !is.na(MidStageSamples$chemotherapy), ]  ## 61
MidStageSamplesCMS4 <- MidStageSamples[intersect(grep("CMS4", MidStageSamples$CMSTME), which(!is.na(MidStageSamples$chemotherapy))), ]  ## 122 114


####### dfs
### CMS1
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS1, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#E69F24", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS2, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#0273B3", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS3, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#CC79A7", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4_HighIF, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#3C5488", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4_LowIF, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#8491B4", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MidStageSamplesCMS4, "rfs.delay", "rfs.event", "chemotherapy", Color=c("#009F73", "#6D3E33"), type="mothes")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))


###########immunotherapy

load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315_2/1.cluster/Normalized_Count_TCGA622.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230310/kmeans/ARGO2_counts_normalized.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v5.RData")
checkpointgene <- c("PDCD1", "CD274","PAF1", "PDCD1LG2")
cpexp <- t(scale(t(ARGO2_counts_normalized[checkpointgene, ])))
rownames(cpexp) <- c("PD1", "PDL1", "PD2", "PDL2")

###### plot
CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4')
cpexp.pst <- PST(as.matrix(cpexp))
cpexp.pst <- data.frame(Check_Point_Genes=cpexp.pst[, 1], Sample=cpexp.pst[, 2], Expression=as.numeric(cpexp.pst[, 3]))
cpexp.pst$CMSTME <- ARGO_sdiused_label2$Immune.Subtype[match(cpexp.pst$Sample, ARGO_sdiused_label2$RNAseqID)]
#
anno_cpexp <- compare_means(Expression ~ CMSTME, method="t.test", group.by = "Check_Point_Genes", data = cpexp.pst, p.adjust.method = "BH") %>%
 mutate(y_pos = 5, p.adj = format.pval(p.adj, digits = 2))
anno_cpexp <- anno_cpexp[anno_cpexp$group1 %in% c("CMS4-IF+") | anno_cpexp$group2 %in% c("CMS4-IF+"), ]
anno_cpexp$y_pos <- rep(c(3, 4, 5, 6), 4)
#
cpexp.pst$Check_Point_Genes <- factor(cpexp.pst$Check_Point_Genes, levels=c("PD1", "PDL1", "PD2", "PDL2"))
ggboxplot(cpexp.pst, x = "CMSTME", y = "Expression", fill="CMSTME", ggtheme = theme_bw(), palette=CMSTMESubtype) +
  facet_wrap(~Check_Point_Genes, ncol = 4) +
  geom_signif(
    data=anno_cpexp,
    aes(xmin = group1, xmax = group2, annotations = formatC(p, format = "e", digits = 2), y_position = y_pos),
    manual= TRUE
  )+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

PD1GeneSig <- c("CD2", "CTSW", "GZMB", "HLA-E", "HSPA8", "ID2", "ITM2B", "LSP1", "MS4A4A", "MS4A4E", "S100A10", "SERINC3", "SERPINE2", "SHISA")
cpexp2 <- t(scale(t(ARGO2_counts_normalized[rownames(ARGO2_counts_normalized) %in% PD1GeneSig, ]))) ## 12

CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4')
cpexp2.pst <- PST(as.matrix(cpexp2))
cpexp2.pst <- data.frame(Check_Point_Genes=cpexp2.pst[, 1], Sample=cpexp2.pst[, 2], Expression=as.numeric(cpexp2.pst[, 3]))
cpexp2.pst$CMSTME <- ARGO_sdiused_label2$Immune.Subtype[match(cpexp2.pst$Sample, ARGO_sdiused_label2$RNAseqID)]


my_comparisons = list( c("CMS4-IF+", "CMS1"), c("CMS4-IF+", "CMS2"), c("CMS4-IF+", "CMS3"), c("CMS4-IF+", "CMS4-IF-"))
bp <- ggplot(cpexp2.pst, aes(x=CMSTME, y=Expression, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(x="Subtype", y = "PD1 response signature score") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


library(IOBR)
GCD_8_T_effector <- lapply(intersect(signature_tme$CD_8_T_effector, rownames(ARGO2_counts_normalized)), function(x){
	TFExpDistr <- data.frame(Expression=as.numeric(as.matrix(ARGO2_counts_normalized[x, ])), Subtype=ARGO_sdiused_label2$Immune.Subtype)
	my_comparisons = list( c("CMS4-IF+", "CMS1"), c("CMS4-IF+", "CMS2"), c("CMS4-IF+", "CMS3"), c("CMS4-IF+", "CMS4-IF-"))
	bp <- ggplot(TFExpDistr, aes(x=Subtype, y=Expression, fill=Subtype)) + geom_boxplot(width=0.8) + 
      labs(title=x,x="Subtype", y = "Expression") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)
	bp
})
PCD_8_T_effector <- ggarrange(GCD_8_T_effector[[1]], GCD_8_T_effector[[2]], GCD_8_T_effector[[3]], GCD_8_T_effector[[4]], GCD_8_T_effector[[5]], GCD_8_T_effector[[6]], GCD_8_T_effector[[7]], GCD_8_T_effector[[8]],
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                    ncol = 4, nrow = 2)
GImmune_Checkpoint <- lapply(intersect(signature_tme$Immune_Checkpoint, rownames(ARGO2_counts_normalized)), function(x){
	TFExpDistr <- data.frame(Expression=as.numeric(as.matrix(ARGO2_counts_normalized[x, ])), Subtype=ARGO_sdiused_label2$Immune.Subtype)
	my_comparisons = list( c("CMS4-IF+", "CMS1"), c("CMS4-IF+", "CMS2"), c("CMS4-IF+", "CMS3"), c("CMS4-IF+", "CMS4-IF-"))
	bp <- ggplot(TFExpDistr, aes(x=Subtype, y=Expression, fill=Subtype)) + geom_boxplot(width=0.8) + 
      labs(title=x,x="Subtype", y = "Expression") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)
	bp
})
PImmune_Checkpoint <- ggarrange(GImmune_Checkpoint[[1]], GImmune_Checkpoint[[2]], GImmune_Checkpoint[[3]], GImmune_Checkpoint[[4]], GImmune_Checkpoint[[5]], GImmune_Checkpoint[[6]],
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 4, nrow = 2)

### heatmap
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230310/kmeans/ARGO2_counts_normalized.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v5.RData")
library(IOBR)
PD1Gene <- c("PDCD1", "CD274","PAF1", "PDCD1LG2")
#PD1GeneSig <- c("CD2", "CTSW", "GZMB", "HLA-E", "HSPA8", "ID2", "ITM2B", "LSP1", "MS4A4A", "MS4A4E", "S100A10", "SERINC3", "SERPINE2", "SHISA")
PD1GeneSig <- c("CD2", "CTSW", "GZMB", "HLA-E", "ITM2B", "LSP1", "MS4A4A", "SERPINE2")
cpexp3 <- ARGO2_counts_normalized[rownames(ARGO2_counts_normalized) %in% PD1GeneSig, ] ## 12
cpexp4 <- ARGO2_counts_normalized[rownames(ARGO2_counts_normalized) %in% signature_tme$CD_8_T_effector, ]
cpexp34 <- t(scale(t(rbind(cpexp3, cpexp4))))
cpexp35 <- t(scale(t(rbind(cpexp4, rbind(cpexp3, ARGO2_counts_normalized[rownames(ARGO2_counts_normalized) %in% PD1Gene, ])))))

col_fun3 = colorRamp2(unique(c(seq(-2, 0, length.out=18), seq(0, 2, length.out=18))), colorRampPalette(c("#3660A5", "#000000", "#CB1A1E"))(35))
dend34 = cluster_within_group(cpexp34, ARGO_sdiused_label2[, "Immune.Subtype"])
colannoA34 = HeatmapAnnotation(
		CMSPlus=ARGO_sdiused_label2$Immune.Subtype,
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4'))
)
row_split34 = data.frame(rep(c("CD8 T effector", "PD1 Gene", "Check point"), c(dim(cpexp4)[1], dim(cpexp3)[1], 4)))
PGARGOIT <- Heatmap(cpexp35, col = col_fun3, cluster_rows = FALSE, cluster_columns = dend34, top_annotation = colannoA34, column_labels=rep("", dim(ARGO_sdiused_label2)[1]), row_split = row_split34)

columnOrder34 <- c(intersect(column_order(PGARGOIT), grep("CMS1", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(column_order(PGARGOIT), grep("CMS2", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(column_order(PGARGOIT), grep("CMS3", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(column_order(PGARGOIT), grep("CMS4-IF\\+", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(column_order(PGARGOIT), grep("CMS4-IF-", ARGO_sdiused_label2[, "Immune.Subtype"])))  ## 
set.seed(1234)
columnOrder34 <- c(intersect(column_order(PGARGOIT), grep("CMS1", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(sample(column_order(PGARGOIT)), grep("CMS2", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(sample(column_order(PGARGOIT)), grep("CMS3", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(sample(column_order(PGARGOIT)), grep("CMS4-IF\\+", ARGO_sdiused_label2[, "Immune.Subtype"])), intersect(sample(column_order(PGARGOIT)), grep("CMS4-IF-", ARGO_sdiused_label2[, "Immune.Subtype"])))

colannoA342 = HeatmapAnnotation(
		CMSPlus=ARGO_sdiused_label2$Immune.Subtype[columnOrder34],
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4'))
)
col.subsections34 <- as.numeric(table(ARGO_sdiused_label2[, "Immune.Subtype"]))
col_split34 = data.frame(rep(names(table(ARGO_sdiused_label2[, "Immune.Subtype"])), col.subsections34))
PGARGOIT2 <- Heatmap(cpexp35[, columnOrder34], col = col_fun3, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = colannoA342, column_labels=rep("", dim(ARGO_sdiused_label2)[1]), row_split = row_split34, column_split = col_split34)
PGARGOIT2
