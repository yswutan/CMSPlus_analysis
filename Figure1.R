### Figure1A
library(ComplexHeatmap)
library(circlize)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/model/intoeder.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")
Datasetcolor=RColorBrewer::brewer.pal(10, "Set3")[-9]
names(Datasetcolor) <- names(table(CMSclin2$dataset))
col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))
colannoA4 = HeatmapAnnotation(
		CMSPlus=CMSclin2$NewCluster,
		CMS= CMSclin2$cms_label, 
		Dataset=CMSclin2$dataset,
		MSI=CMSclin2$msi, 
		CIMP=CMSclin2$cimp, 
		KRAS=CMSclin2$kras_mut, 
		BRAF=CMSclin2$braf_mut, 
		Age=CMSclin2$age,
		Gender=CMSclin2$gender,
		Grade=CMSclin2$grade,
		Pathologic_Tstage=CMSclin2$pt, 
		Pathologic_Nstage=CMSclin2$pn, 
		Pathologic_Mstage=CMSclin2$pm, 
		pathologic_TNMstage=CMSclin2$tnm, 
    col = list(CMSPlus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', `CMS4-TME+`='#3C5488', `CMS4-TME-`='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20'), Dataset=Datasetcolor, Gender=c(female="#E69F00", male="#56B4E9"), Grade=c('1'=RColorBrewer::brewer.pal(10, "Blues")[2], '2'=RColorBrewer::brewer.pal(10, "Blues")[4], '3'=RColorBrewer::brewer.pal(10, "Blues")[6], '4'=RColorBrewer::brewer.pal(10, "Blues")[8]), Pathologic_Tstage=c('1'=RColorBrewer::brewer.pal(10, "YlGn")[2], '2'=RColorBrewer::brewer.pal(10, "YlGn")[4], '3'=RColorBrewer::brewer.pal(10, "YlGn")[6], '4'=RColorBrewer::brewer.pal(10, "YlGn")[8]), Pathologic_Nstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(10, "Paired")[2]), Pathologic_Mstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(12, "Paired")[12]), pathologic_TNMstage=c("I"=RColorBrewer::brewer.pal(10, "BuPu")[2], "II"=RColorBrewer::brewer.pal(10, "BuPu")[4], "III"=RColorBrewer::brewer.pal(10, "BuPu")[6], "IV"=RColorBrewer::brewer.pal(10, "BuPu")[8]), Age=colorRamp2(c(20, 100), c("white", "red")))
)
row.subsections1 <- c(3, 4, 5, 3)
row_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4"), row.subsections1))
col.subsections1 <- as.numeric(table(CMSclin2$NewCluster))
col_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4-TME+", "CMS4-TME-"), col.subsections1))
PPCRCSC2 <- Heatmap(gsvaall2[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = colannoA4, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split1, column_split = col_split1)
row.subsections2 <- c(8, 19)
row_split2 = data.frame(factor(rep(c("Fibroblastes", "Immune"), row.subsections2), levels=c("Immune", "Fibroblastes")))
PPCRCSC3 <- Heatmap(gsvaall2[16:42, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split2, column_split = col_split1)
ht_list = PPCRCSC2%v%PPCRCSC3
draw(ht_list)

# Figure1C
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/immune/A3EstimateScores.RData")
CMSclinallE <- merge(CMSclin2, A3EstimateScores, by.x="sample", by.y="row.names")
CMSPlusSubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
### Stromal Score
bpc1p<-formatC(t.test(CMSclinallE$StromalScore[which(CMSclinallE$CMSTME %in% "CMS4-TME+")], CMSclinallE$StromalScore[which(CMSclinallE$CMSTME %in% "CMS4-TME-")])$p.value, format = "e", digits = 2)
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
bpc1 <- ggplot(CMSclinallE, aes(x=CMSTME, y=StromalScore, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Stromal Score",x="Subtype", y = "Stromal Score") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc1p), y_position = c(1800), tip_length = 0, vjust=0)
### Immune Score
bpc2p<-formatC(t.test(CMSclinallE$ImmuneScore[which(CMSclinallE$CMSTME %in% "CMS4-TME+")], CMSclinallE$ImmuneScore[which(CMSclinallE$CMSTME %in% "CMS4-TME-")])$p.value, format = "e", digits = 2)
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
bpc2 <- ggplot(CMSclinallE, aes(x=CMSTME, y=ImmuneScore, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Immune Score",x="Subtype", y = "Immune Score") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc2p), y_position = c(2400), tip_length = 0, vjust=0)
### Tumor purity
bpc3p<-formatC(t.test(CMSclinallE$TumourPurity[which(CMSclinallE$CMSTME %in% "CMS4-TME+")], CMSclinallE$TumourPurity[which(CMSclinallE$CMSTME %in% "CMS4-TME-")])$p.value, format = "e", digits = 2)
my_comparisons = list( c("CMS4-TME+", "CMS4-TME-") )
bpc3 <- ggplot(CMSclinallE, aes(x=CMSTME, y=TumourPurity, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Tumor Purity",x="Subtype", y = "Tumor Purity") + 
      scale_fill_manual(values=CMSPlusSubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc3p), y_position = c(0.9), tip_length = 0, vjust=0)
figure <- ggarrange(bpc1, bpc2, bpc3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure

# Figure1D
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/immune")
load("/data0/tan/Task/9.CMS.SDI/data/5.PublicExpSurData/CRCSC/CRCSC.expression.clindata.RData")  ## expall.corr, clindata
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")
expcrctcga <- expall.corr[, CMSclin2$sample]
library(MCPcounter)
library(ggplot2)
CRCMCP <- MCPcounter.estimate(expcrctcga, featuresType='HUGO_symbols')
CRCMCP.pst <- PST(as.matrix(CRCMCP))
CRCMCP.pst <- data.frame(ImmuneInfiltration=CRCMCP.pst[, 1], Sample=CRCMCP.pst[, 2], Score=as.numeric(CRCMCP.pst[, 3]))
CRCMCP.pst$`CMSPlus` <- CMSclinall$CMSTME[match(CRCMCP.pst$Sample, CMSclinall$sample)]
#
anno_CRCMCP <- compare_means(Score ~ CMSPlus, method="t.test", group.by = "ImmuneInfiltration", data = CRCMCP.pst, p.adjust.method = "BH") %>%
 mutate(y_pos = 12.5, p.adj = format.pval(p.adj, digits = 2))
anno_CRCMCP <- anno_CRCMCP[anno_CRCMCP$group1 %in% c("CMS4-TME+", "CMS4-TME-") & anno_CRCMCP$group2 %in% c("CMS4-TME+", "CMS4-TME-"), ]
CMSPlusSubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-TME+'='#3C5488', 'CMS4-TME-'='#8491B4')
#
ggboxplot(CRCMCP.pst, x = "CMSPlus", y = "Score", fill="CMSPlus", ggtheme = theme_bw(), palette=CMSPlusSubtype) +
  facet_wrap(~ImmuneInfiltration, ncol = 5) + 
  geom_signif(
    data=anno_CRCMCP, 
    aes(xmin = group1, xmax = group2, annotations = formatC(p, format = "e", digits = 2), y_position = y_pos), 
    manual= TRUE
  )+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

# Figure1E
setwd("/mnt/wutan/data/9.CMSPlus/Figure3/model")
load("CMSclin2.gsvaall2.RData")
load("CRC.clin.exp.he_v2.RData")
library(survival)
library(gfplot)
SurvivalFGPlotSimplify <- function(DataFrame, Time, Status, Lable, Color, type="days"){
	if(type=="days"){
		DataFrame[, Time] <- DataFrame[, Time]/30
	}
	clin <- Surv(DataFrame[, Time], DataFrame[, Status])
	labs <- factor(DataFrame[, Lable])
	ylable <- unlist(sapply(strsplit(Time, "\\_"), function(x) x[1]))
	parameters <- list(clin=clin, labs=labs, ylable=ylable, Color=Color)
	return(parameters)
}
CRC.clin.exp.he_v2$CMSPlus <- CMSclin2$NewCluster[match(CRC.clin.exp.he_v2$sample, CMSclin2$sample)]
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "old.rfs.delay", "old.rfs.event", "CMSPlus", Color=c('#E69F24','#0273B3','#CC79A7', '#3C5488', '#8491B4'), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
