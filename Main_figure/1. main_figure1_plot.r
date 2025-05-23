library(ComplexHeatmap)
library(circlize)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/model/intoeder.RData")
###
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/CMSclin.v1.RData")
CMSclin$NewCluster <- as.character(CMSclin$NewCluster)
CMSclin$NewCluster[is.na(CMSclin$NewCluster)] <- "NOLBL"
CMSclin$gender <- tolower(CMSclin$gender)
CMSclin$tnm <- gsub("Stage ", "", CMSclin$tnm)
CMSclin$tnm <- gsub("B/IIC", "", CMSclin$tnm)
CMSclin$tnm <- gsub("A", "", CMSclin$tnm)
CMSclin$tnm <- gsub("B", "", CMSclin$tnm)
CMSclin$tnm <- gsub("C", "", CMSclin$tnm)
CMSclin$grade <- gsub("G", "", CMSclin$grade)
CMSclin$grade[CMSclin$grade %in% "Not Available"] <- NA
CMSclin$msi[CMSclin$msi %in% "Not Available"] <- NA
CMSclin$dataset <- toupper(CMSclin$dataset)
CMSclin$pn <- ifelse(CMSclin$pn == 0, 0, 1)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.all.42.v2.RData")
CMSclin2 <- CMSclin[!CMSclin$dataset %in% c("GSE17536", "GSE17537", "GSE14333"), ]  ## 2006   20
CMSclin2 <- CMSclin2[match(intersect(intoeder, CMSclin2$sample), CMSclin2$sample), ]
TMIGOrder <- read.table("/data0/tan/Task/9.CMS.SDI/data/4.CancerIntrinsicSignature/TMI_gene_signatures_order.tsv", header=F, sep="\t", quote="")
gsva.all.42 <- gsva.all.45[1:42, ]
gsvaall2 <- gsva.all.42[c(1:4, 6:7, 5, 8:15, match(TMIGOrder[, 1], rownames(gsva.all.42))[-(28:29)]), CMSclin2$sample]
gsvaall3 <- rbind(gsvaall2, gsva.all.45[43:45, CMSclin2$sample])

Datasetcolor=RColorBrewer::brewer.pal(10, "Set3")[-9]
names(Datasetcolor) <- names(table(CMSclin2$dataset))

col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))

colannoA4 = HeatmapAnnotation(
		NewCluster=CMSclin2$NewCluster,
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
    col = list(NewCluster=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20'), Dataset=Datasetcolor, Gender=c(female="#E69F00", male="#56B4E9"), Grade=c('1'=RColorBrewer::brewer.pal(10, "Blues")[2], '2'=RColorBrewer::brewer.pal(10, "Blues")[4], '3'=RColorBrewer::brewer.pal(10, "Blues")[6], '4'=RColorBrewer::brewer.pal(10, "Blues")[8]), Pathologic_Tstage=c('1'=RColorBrewer::brewer.pal(10, "YlGn")[2], '2'=RColorBrewer::brewer.pal(10, "YlGn")[4], '3'=RColorBrewer::brewer.pal(10, "YlGn")[6], '4'=RColorBrewer::brewer.pal(10, "YlGn")[8]), Pathologic_Nstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(10, "Paired")[2]), Pathologic_Mstage=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(12, "Paired")[12]), pathologic_TNMstage=c("I"=RColorBrewer::brewer.pal(10, "BuPu")[2], "II"=RColorBrewer::brewer.pal(10, "BuPu")[4], "III"=RColorBrewer::brewer.pal(10, "BuPu")[6], "IV"=RColorBrewer::brewer.pal(10, "BuPu")[8]), Age=colorRamp2(c(20, 100), c("white", "red")))
)
row.subsections1 <- c(3, 4, 5, 3)
row_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4"), row.subsections1))
col.subsections1 <- as.numeric(table(CMSclin2$NewCluster))
col_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4_HighIF", "CMS4_LowIF"), col.subsections1))
PPCRCSC2 <- Heatmap(gsvaall2[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = colannoA4, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split1, column_split = col_split1)
row.subsections2 <- c(8, 8, 11, 3)
row_split2 = data.frame(factor(rep(c("Fibroblastes", "ProTumor", "AntiTumor", "TMEscore"), row.subsections2), levels=c("ProTumor", "AntiTumor", "Fibroblastes", "TMEscore")))
PPCRCSC3 <- Heatmap(gsvaall3[16:45, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, column_labels=rep("", dim(CMSclin2)[1]), row_split = row_split2, column_split = col_split1)
ht_list = PPCRCSC2%v%PPCRCSC3
draw(ht_list)

load("/data0/tan/Task/9.CMS.SDI/data/5.PublicExpSurData/CRCSC/CRCSC.expression.clindata.RData")  ## expall.corr, clindata
library(GSVA)
gsva.all.42 <- gsva(expall.corr, CMSsignatures42.v2, method=c("gsva"), kcdf = "Gaussian", min.sz = 1, parallel.sz=40)
save(gsva.all.42, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.all.42.v2.RData")

cmssdirfs <- SurvivalFGPlotSimplify(CMSclin2[CMSclin2$dataset %in% "TCGA", ], "rfsMo", "rfsStat", "CMSTME", Color=c('#E69F24','#0273B3', '#CC79A7', '#3C5488', '#8491B4'), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

###########################################
################# ESTIMATE 
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/immune")
load("/data0/tan/Task/9.CMS.SDI/data/5.PublicExpSurData/CRCSC/CRCSC.expression.clindata.RData")  ## expall.corr, clindata
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")

expcrctcga <- expall.corr[, CMSclin2$sample]

library(estimate)
A3Estimate <- estimate(dat=expcrctcga, pro="CRCSCTCGA")
#[1] "Merged dataset includes 9913 genes (499 mismatched)."
#[1] "1 gene set: StromalSignature  overlap= 136"
#[1] "2 gene set: ImmuneSignature  overlap= 140"
A3EstimateScores <- lapply(list(A3=A3Estimate), function(x) {df <- as.data.frame(x); df$TumourPurity <- apply(x, 1, function(y) cos(0.6049872018+0.0001467884 * y[3])); df})[[1]]  ## ref: https://www.nature.com/articles/ncomms3612
rownames(A3EstimateScores) <- gsub("\\.", "-", rownames(A3EstimateScores))
save(A3EstimateScores, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/immune/A3EstimateScores.RData")

###########################################
########## CRCSCTCGA
CMSclinallE <- merge(CMSclin2, A3EstimateScores, by.x="sample", by.y="row.names")
#CMSSDIcolor <- c('#E69F24','#0273B3', '#CC79A7', '#2B2F81', '#060606')
CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4')

### Stromal Score
bpc1p<-formatC(t.test(CMSclinallE$StromalScore[which(CMSclinallE$CMSTME %in% "CMS4-IF+")], CMSclinallE$StromalScore[which(CMSclinallE$CMSTME %in% "CMS4-IF-")])$p.value, format = "e", digits = 2)
# p-value = 2.2e-16
my_comparisons = list( c("CMS4-IF+", "CMS4-IF-") )
bpc1 <- ggplot(CMSclinallE, aes(x=CMSTME, y=StromalScore, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Stromal Score",x="Subtype", y = "Stromal Score") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc1p), y_position = c(1800), tip_length = 0, vjust=0)
### Immune Score
bpc2p<-formatC(t.test(CMSclinallE$ImmuneScore[which(CMSclinallE$CMSTME %in% "CMS4-IF+")], CMSclinallE$ImmuneScore[which(CMSclinallE$CMSTME %in% "CMS4-IF-")])$p.value, format = "e", digits = 2)
# p-value = 2.2e-16
my_comparisons = list( c("CMS4-IF+", "CMS4-IF-") )
bpc2 <- ggplot(CMSclinallE, aes(x=CMSTME, y=ImmuneScore, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Immune Score",x="Subtype", y = "Immune Score") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc2p), y_position = c(2400), tip_length = 0, vjust=0)
### Tumor purity
bpc3p<-formatC(t.test(CMSclinallE$TumourPurity[which(CMSclinallE$CMSTME %in% "CMS4-IF+")], CMSclinallE$TumourPurity[which(CMSclinallE$CMSTME %in% "CMS4-IF-")])$p.value, format = "e", digits = 2)
# p-value = 2.2e-16
my_comparisons = list( c("CMS4-IF+", "CMS4-IF-") )
bpc3 <- ggplot(CMSclinallE, aes(x=CMSTME, y=TumourPurity, fill=CMSTME)) + geom_boxplot(width=0.8) + 
      labs(title="CRC Tumor Purity",x="Subtype", y = "Tumor Purity") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, annotations=c(bpc3p), y_position = c(0.9), tip_length = 0, vjust=0)


figure <- ggarrange(bpc1, bpc2, bpc3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure

