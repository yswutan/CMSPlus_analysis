
################ data
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.all.42.v2.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/CMSclin.v1.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.mvrm.v2.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.tcga.42.v2.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/model/CRC.clin.exp.he_v2.RData")
load("/data/home2/wutan/Task/MultiomicsSubtyping/data/ICGC/clin/20220211/Clin1000Focus.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.argo.42.v2.RData")#
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v4.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")
##
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model")
#### train
TrainData <- gsvaall2  ##  
TrainTrait <- CMSclin2
Y <- CMSclin2$CMSTME
#Compute weights to balance the RF
w <- 1/table(Y)
w <- w/sum(w)
weights <- rep(0, length(Y))
weights[Y == 'CMS1'] <- w['CMS1']
weights[Y == 'CMS2'] <- w['CMS2']
weights[Y == 'CMS3'] <- w['CMS3']
weights[Y == 'CMS4-IF+'] <- w['CMS4-IF+']
weights[Y == 'CMS4-IF-'] <- w['CMS4-IF-']
table(weights, Y)
#
train.vals <- as.data.frame(t(TrainData))
colnames(train.vals) <- gsub("-", "_", colnames(train.vals))
colnames(train.vals) <- gsub(" ", "_", colnames(train.vals))
train.vals$Subtype <- factor(Y)
set.seed(666)
epitme.rf.model <- ranger(Subtype~., train.vals, case.weights=weights)
save(epitme.rf.model, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/epitme.rf.model.Rdata")
confusionMatrix(factor(epitme.rf.model$predictions), factor(Y)) ## 0.8346
##
cmargomodel <- as.data.frame(table(epitme.rf.model$predictions, factor(Y)))
colnames(cmargomodel) <- c("x", "y", "n")
pcm2 <- plot_conf_mtx(
  cmargomodel, 
  xlab = "Prediction", 
  ylab = "Reference", 
  y_order = rev(c("CMS1", "CMS2", "CMS3", "CMS4-IF+", "CMS4-IF-"))
)
pcm2

# heatmap
PlotHeatMap <- function(ClinData, GSVAData, CMSTME, CMS, SDI, random=F, seed=6){
	columnAnno = HeatmapAnnotation(
		CMSTME = ClinData[, CMSTME], 
		CMS=ClinData[, CMS], 
		col = list(CMSTME =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-IF+"='#3C5488', "CMS4-IF-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"))
	)
	col_fun = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#368ABF", "#68C3A7", "#AED8A3", "#E3EA9B", "white", "#FDE18B", "#F9AF62", "#F36D46", "#D54150"))(35))
	col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))
	dendtcga = cluster_within_group(GSVAData, ClinData[, CMSTME])
	P1 <- Heatmap(GSVAData, col = col_fun, cluster_rows = FALSE, cluster_columns = dendtcga, top_annotation = columnAnno, column_labels=rep("", dim(GSVAData)[2]))
	if(!random){
		columnOrder <- c(intersect(column_order(P1), grep("CMS1", ClinData[, CMSTME])), intersect(column_order(P1), grep("CMS2", ClinData[, CMSTME])), intersect(column_order(P1), grep("CMS3", ClinData[, CMSTME])), intersect(column_order(P1), grep("CMS4-IF\\+", ClinData[, CMSTME])), intersect(column_order(P1), grep("CMS4-IF-", ClinData[, CMSTME])))  ## 
	}
	if(random){
		set.seed(seed)
		columnOrder <- c(sample(intersect(column_order(P1), grep("CMS1", ClinData[, CMSTME])), length(grep("CMS1", ClinData[, CMSTME]))), sample(intersect(column_order(P1), grep("CMS2", ClinData[, CMSTME])), length(grep("CMS2", ClinData[, CMSTME]))), sample(intersect(column_order(P1), grep("CMS3", ClinData[, CMSTME])), length(grep("CMS3", ClinData[, CMSTME]))), sample(intersect(column_order(P1), grep("CMS4-IF\\+", ClinData[, CMSTME])), length(grep("CMS4-IF\\+", ClinData[, CMSTME]))), sample(intersect(column_order(P1), grep("CMS4-IF-", ClinData[, CMSTME])), length(grep("CMS4-IF\\-", ClinData[, CMSTME]))))  ## 
	}
	if(length(SDI)>0){
		ClinData[, SDI] <- gsub("-", "_", ClinData[, SDI])
		columnAnno2 = HeatmapAnnotation(
			CMSTME = ClinData[, CMSTME][columnOrder], 
			CMS=ClinData[, CMS][columnOrder], 
			CMSSDI=ClinData[, SDI][columnOrder],
			col = list(CMSTME =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-IF+"='#3C5488', "CMS4-IF-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"), CMSSDI=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighSDI='#C2CB1E', CMS4_LowSDI='#537A34')
		))
	} else {
		columnAnno2 = HeatmapAnnotation(
			CMSTME = ClinData[, CMSTME][columnOrder], 
			CMS=ClinData[, CMS][columnOrder], 
			#CMSSDI=ClinData[, SDI][columnOrder],
			col = list(CMSTME =c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-IF+"='#3C5488', "CMS4-IF-"='#8491B4'), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73"), CMSSDI=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighSDI='#C2CB1E', CMS4_LowSDI='#537A34')
		))
	}
	TMIGOrder <- read.table("/data0/tan/Task/9.CMS.SDI/data/4.CancerIntrinsicSignature/TMI_gene_signatures_order.tsv", header=F, sep="\t", quote="")
	GSVAData2 <- GSVAData[c(1:4, 6:7, 5, 8:15, match(TMIGOrder[, 1], rownames(GSVAData))[-(28:29)]), columnOrder]
	row.subsections1 <- c(3, 4, 5, 3)
	row_split1 = data.frame(rep(c("CMS1", "CMS2", "CMS3", "CMS4"), row.subsections1))
	col.subsections1 <- as.numeric(table(ClinData[, CMSTME]))
	col_split1 = data.frame(rep(names(table(ClinData[, CMSTME])), col.subsections1))
	P2 <- Heatmap(GSVAData2[1:15, ], col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = columnAnno2, column_labels=rep("", dim(GSVAData2)[2]), row_split = row_split1, column_split = col_split1)
	row.subsections2 <- c(8, 8, 11)
	row_split2 = data.frame(factor(rep(c("Fibroblastes", "ProTumor", "AntiTumor"), row.subsections2), levels=c("ProTumor", "AntiTumor", "Fibroblastes")))
	P3 <- Heatmap(GSVAData2[16:42, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, column_labels=rep("", dim(GSVAData2)[2]), row_split = row_split2, column_split = col_split1)
	ht_list = P2%v%P3
	draw(ht_list)
}

##### predicted res

######### TCGA
load("/data0/tan/Task/9.CMS.SDI/data/1.TCGA/CRC.clin.exp.he_v2.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.tcga.42.v2.RData")

##
tcga.vals <- as.data.frame(t(gsva.tcga.42))
colnames(tcga.vals) <- gsub("-", "_", colnames(tcga.vals))
colnames(tcga.vals) <- gsub(" ", "_", colnames(tcga.vals))
probabilities1 = predict(epitme.rf.model, data = tcga.vals)$predictions
table(probabilities1)
#probabilities1 <- factor(probabilities1, levels=c("CMS1", "CMS2", "CMS3", "CMS4-IF+", "CMS4-IF-"))
#probabilities1
#    CMS1     CMS2     CMS3 CMS4-IF+ CMS4-IF- 
#      97      222       89      130       84

CRC.clin.exp.he_v2$Immune.Subtype <- probabilities1
table(CRC.clin.exp.he_v2$Immune.Subtype, CRC.clin.exp.he_v2$RF.nearestCMS.mod) 
table(CRC.clin.exp.he_v2$CMSSDI, CRC.clin.exp.he_v2$Immune.Subtype)
GetHGPlot <- function(atable, colname1, colname2, samplename, nrow=4){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow)
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3))
}
GetHGPlot(atable=CRC.clin.exp.he_v2, colname1="RF.nearestCMS.mod", colname2="Immune.Subtype", samplename="sample")
GetHGPlot(atable=CRC.clin.exp.he_v2[!is.na(CRC.clin.exp.he_v2$CMSSDI), ], colname1="CMSSDI", colname2="Immune.Subtype", samplename="sample", nrow=5)
save(CRC.clin.exp.he_v2, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")

CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', "CMS4-IF+"='#3C5488', "CMS4-IF-"='#8491B4')## new follow up (should)
### ImmuneSubtype
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "rfs.delay", "rfs.event", "Immune.Subtype", Color=CMSTMESubtype, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "old.rfs.delay", "old.rfs.event", "Immune.Subtype", Color=CMSTMESubtype, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73")
cmssdirfs <- SurvivalFGPlotSimplify(CRC.clin.exp.he_v2, "rfs.delay", "rfs.event", "RF.nearestCMS.mod", Color=CMS, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

# heatmap
PlotHeatMap(ClinData=CRC.clin.exp.he_v2, GSVAData=gsva.tcga.42, CMSTME="Immune.Subtype", CMS="RF.nearestCMS.mod", SDI="CMSSDI")




### MVRM
load("/data/home2/wutan/Task/9.CMS.SDI/data/5.PublicExpSurData/public.datasets.allgenes.RData")  ## MVRM, CIT, Jorissen, Laibe
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.mvrm.v2.RData")
#
mvrm.vals <- as.data.frame(t(gsva.mvrm))
colnames(mvrm.vals) <- gsub("-", "_", colnames(mvrm.vals))
colnames(mvrm.vals) <- gsub(" ", "_", colnames(mvrm.vals))
probabilities3 = predict(epitme.rf.model, data = mvrm.vals)$predictions
probabilities3
#    CMS1     CMS2     CMS3 CMS4-IF+ CMS4-IF- 
#      87      195       72       75       93
MVRM$clindata$Immune.Subtype <- probabilities3
table(MVRM$clindata$Immune.Subtype, MVRM$clindata$cms) 
table(MVRM$clindata$Immune.Subtype, MVRM$clindata$cms2) 

cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "rfs.delay", "rfs.event", "Immune.Subtype", Color=CMSTMESubtype, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(MVRM$clindata, "rfs.delay", "rfs.event", "cms", Color=CMS, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

MVRM.clindata <- MVRM$clindata
MVRM.clindata$Immune.Subtype <- factor(MVRM.clindata$Immune.Subtype)
GetHGPlot(atable=MVRM.clindata, colname1="cms2", colname2="Immune.Subtype", samplename="sample")
save(MVRM.clindata, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/MVRM.clindata.RData")

# heatmap
PlotHeatMap(ClinData=MVRM$clindata, GSVAData=gsva.mvrm, CMSTME="Immune.Subtype", CMS="cms", SDI=NULL)
PlotHeatMap(ClinData=MVRM.clindata[!is.na(MVRM.clindata$cms), ], GSVAData=gsva.mvrm[, !is.na(MVRM.clindata$cms)], CMSTME="Immune.Subtype", CMS="cms", SDI=NULL)
PlotHeatMap(ClinData=MVRM.clindata, GSVAData=gsva.mvrm, CMSTME="Immune.Subtype", CMS="cms", SDI=NULL, random=T, seed=6)

##
cmmvrmtes <- as.data.frame(table(substr(MVRM.clindata$Immune.Subtype, 1, 4), MVRM.clindata$cms))
colnames(cmmvrmtes) <- c("x", "y", "n")
pcm3 <- plot_conf_mtx(
  cmmvrmtes, 
  xlab = "Prediction", 
  ylab = "Reference", 
  y_order = rev(c("CMS1", "CMS2", "CMS3", "CMS4"))
)
pcm3
confusionMatrix(factor(substr(MVRM.clindata$Immune.Subtype, 1, 4)), factor(MVRM.clindata$cms2))

load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.argo.42.v2.RData")## gsva.argo.42, file=
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v4.RData")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/epitme.rf.model.Rdata")
##
argo.vals <- as.data.frame(t(gsva.argo.42))
colnames(argo.vals) <- gsub("-", "_", colnames(argo.vals))
colnames(argo.vals) <- gsub(" ", "_", colnames(argo.vals))
probabilitiesargo = predict(epitme.rf.model, data = argo.vals)$predictions
table(probabilitiesargo)
#probabilitiesargo
#    CMS1     CMS2     CMS3 CMS4-IF+ CMS4-IF- 
#      56      142       72       86       66

ARGO_sdiused_label2$Immune.Subtype <- probabilitiesargo
table(ARGO_sdiused_label2$Immune.Subtype, ARGO_sdiused_label2$CMSTME)
#           CMS1 CMS2 CMS3 CMS4_HighIF CMS4_LowIF
#  CMS1       33    5   11           7          0
#  CMS2        2  123   15           0          2
#  CMS3        4    6   61           1          0
#  CMS4-IF+    5    4    0          61         16
#  CMS4-IF-    7    6    2           2         49
table(ARGO_sdiused_label2$Immune.Subtype, ARGO_sdiused_label2$self_all_inAverSDI_5)
#           CMS1 CMS2 CMS3 CMS4_HighSDI CMS4_LowSDI
#  CMS1       33    5   11            5           2
#  CMS2        2  123   15            1           1
#  CMS3        4    6   61            1           0
#  CMS4-IF+    5    4    0           46          31
#  CMS4-IF-    7    6    2           21          30

######
GetHGPlot <- function(atable, colname1, colname2, samplename, nrow=5){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow)
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3, 4))
}
GetHGPlot(atable=ARGO_sdiused_label2, colname1="self_all_inAverSDI_5", colname2="Immune.Subtype", samplename="Sample")
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label2, "rfs.delay", "rfs.event", "Immune.Subtype", Color=CMSTMESubtype, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))

#######
####### heatmap
PlotHeatMap(ClinData=ARGO_sdiused_label2, GSVAData=gsva.argo.42, CMSTME="Immune.Subtype", CMS="cms", SDI="self_all_inAverSDI_5")
ARGO_sdiused_label2$CMSTME2 <- ARGO_sdiused_label2$CMSTME
ARGO_sdiused_label2$CMSTME2 <- gsub("CMS4_HighIF", "CMS4-IF+", ARGO_sdiused_label2$CMSTME2)
ARGO_sdiused_label2$CMSTME2 <- gsub("CMS4_LowIF", "CMS4-IF-", ARGO_sdiused_label2$CMSTME2)
ARGO_sdiused_label2$CMSTME2 <- factor(ARGO_sdiused_label2$CMSTME2, levels=c("CMS1", "CMS2", "CMS3", "CMS4-IF+", "CMS4-IF-"))
PlotHeatMap(ClinData=ARGO_sdiused_label2, GSVAData=gsva.argo.42, CMSTME="CMSTME2", CMS="cms", SDI="self_all_inAverSDI_5")
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label2, "rfs.delay", "rfs.event", "CMSTME2", Color=CMSTMESubtype, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label2, "rfs.delay", "rfs.event", "cms", Color=CMS, type="month")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
save(ARGO_sdiused_label2, file="/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v5.RData")

PlotHeatMap(ClinData=ARGO_sdiused_label2, GSVAData=gsva.argo.42, CMSTME="Immune.Subtype", CMS="cms", SDI="self_all_inAverSDI_5", random=T, seed=6)

###
cmargotes <- as.data.frame(table(substr(ARGO_sdiused_label2$Immune.Subtype, 1, 4), ARGO_sdiused_label2$cms))
colnames(cmargotes) <- c("x", "y", "n")
pcm4 <- plot_conf_mtx(
  cmargotes, 
  xlab = "Prediction", 
  ylab = "Reference", 
  y_order = rev(c("CMS1", "CMS2", "CMS3", "CMS4"))
)
pcm4
confusionMatrix(factor(substr(ARGO_sdiused_label2$Immune.Subtype, 1, 4)), factor(ARGO_sdiused_label2$cms))
