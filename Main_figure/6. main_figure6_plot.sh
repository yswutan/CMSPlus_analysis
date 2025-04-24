load("/data0/tan/Task/9.CMS.SDI/data/2.ARGO/600samples/C_ARGO_clinical_SDI.rdata")  ## "AMC_grant_lable_use"  "ARGO_sdiused_label"  "TCGA_grant_label_use"
# survival
#CMSSDIcolor <- c('#E69F24','#0273B3', '#CC79A7', '#2B2F81', '#060606')
CMScolor <- c('#E69F24','#0273B3', '#CC79A7', '#009F73')
CMSTMESubtype=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4')
CMSSDIcolor =c('#E69F24','#0273B3', '#CC79A7', '#C2CB1E', '#537A34')
ARGO_sdiused_label[ARGO_sdiused_label$sample == "14TS0Z0425TR1", "rfs.event"] <- 0
ARGO_sdiused_label[ARGO_sdiused_label$rfs.delay > 80, "rfs.delay"] <- 80
#####################################################################
################### ARGO
ARGO_sdiused_label$rfs.event[grep("复发", ARGO_sdiused_label$生存现况)] <- 1
ARGOcutoff <- median(ARGO_sdiused_label$Gau_innerAverSDI)
#[1] 0.737966
#ARGOcutoff <- 0.7476921
ARGO_sdiused_label$SDI <- "HighSDI"
ARGO_sdiused_label$SDI[ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "LowSDI"
## CMS1
ARGO_sdiused_label$CMS1.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS1.SDI[ARGO_sdiused_label$cms == "CMS1" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS1_HighSDI"
ARGO_sdiused_label$CMS1.SDI[ARGO_sdiused_label$cms == "CMS1" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS1_LowSDI"
## CMS2
ARGO_sdiused_label$CMS2.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS2.SDI[ARGO_sdiused_label$cms == "CMS2" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS2_HighSDI"
ARGO_sdiused_label$CMS2.SDI[ARGO_sdiused_label$cms == "CMS2" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS2_LowSDI"
## CMS3
ARGO_sdiused_label$CMS3.SDI <- as.character(ARGO_sdiused_label$cms)
ARGO_sdiused_label$CMS3.SDI[ARGO_sdiused_label$cms == "CMS3" & ARGO_sdiused_label$Gau_innerAverSDI > ARGOcutoff] <- "CMS3_HighSDI"
ARGO_sdiused_label$CMS3.SDI[ARGO_sdiused_label$cms == "CMS3" & ARGO_sdiused_label$Gau_innerAverSDI <= ARGOcutoff] <- "CMS3_LowSDI"
## CMS4

####### two groups
### SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "SDI", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS1.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS1", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS1.SDI", Color=c("#FCD7A2", "#DE9B2A"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS2", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS2.SDI", Color=c("#A6D2E8", "#106DA8"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS3", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "CMS3.SDI", Color=c("#E8B5D3", "#C2759F"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4.SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label[grep("CMS4", ARGO_sdiused_label$cms), ], "rfs.delay", "rfs.event", "self_all_inAverSDI_5", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
######## more than two groups
### CMS-SDI
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "self_all_inAverSDI_5", Color=CMSSDIcolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS
cmssdirfs <- SurvivalFGPlotSimplify(ARGO_sdiused_label, "rfs.delay", "rfs.event", "cms", Color=CMScolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))


#####################################################################
################### TCGA
TCGAcutoff <- median(TCGA_grant_label_use$entropy)
#[1] 0.737966
TCGA_grant_label_use$SDI <- "HighSDI"
TCGA_grant_label_use$SDI[TCGA_grant_label_use$entropy <= TCGAcutoff] <- "LowSDI"
## CMS1
TCGA_grant_label_use$CMS1.SDI <- as.character(TCGA_grant_label_use$cms)
TCGA_grant_label_use$CMS1.SDI[TCGA_grant_label_use$cms == "CMS1" & TCGA_grant_label_use$entropy > TCGAcutoff] <- "CMS1_HighSDI"
TCGA_grant_label_use$CMS1.SDI[TCGA_grant_label_use$cms == "CMS1" & TCGA_grant_label_use$entropy <= TCGAcutoff] <- "CMS1_LowSDI"
## CMS2
TCGA_grant_label_use$CMS2.SDI <- as.character(TCGA_grant_label_use$cms)
TCGA_grant_label_use$CMS2.SDI[TCGA_grant_label_use$cms == "CMS2" & TCGA_grant_label_use$entropy > TCGAcutoff] <- "CMS2_HighSDI"
TCGA_grant_label_use$CMS2.SDI[TCGA_grant_label_use$cms == "CMS2" & TCGA_grant_label_use$entropy <= TCGAcutoff] <- "CMS2_LowSDI"
## CMS3
TCGA_grant_label_use$CMS3.SDI <- as.character(TCGA_grant_label_use$cms)
TCGA_grant_label_use$CMS3.SDI[TCGA_grant_label_use$cms == "CMS3" & TCGA_grant_label_use$entropy > TCGAcutoff] <- "CMS3_HighSDI"
TCGA_grant_label_use$CMS3.SDI[TCGA_grant_label_use$cms == "CMS3" & TCGA_grant_label_use$entropy <= TCGAcutoff] <- "CMS3_LowSDI"
## CMS4

####### two groups
### SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use, "rfs.delay", "rfs.event", "SDI", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS1.SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use[grep("CMS1", TCGA_grant_label_use$cms), ], "rfs.delay", "rfs.event", "CMS1.SDI", Color=c("#FCD7A2", "#DE9B2A"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2.SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use[grep("CMS2", TCGA_grant_label_use$cms), ], "rfs.delay", "rfs.event", "CMS2.SDI", Color=c("#A6D2E8", "#106DA8"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3.SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use[grep("CMS3", TCGA_grant_label_use$cms), ], "rfs.delay", "rfs.event", "CMS3.SDI", Color=c("#E8B5D3", "#C2759F"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4.SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use[grep("CMS4", TCGA_grant_label_use$cms), ], "rfs.delay", "rfs.event", "all_inSDI_5", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
######## more than two groups
### CMS-SDI
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use, "rfs.delay", "rfs.event", "all_inSDI_5", Color=CMSSDIcolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS
cmssdirfs <- SurvivalFGPlotSimplify(TCGA_grant_label_use, "rfs.delay", "rfs.event", "cms", Color=CMScolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))


#####################################################################
################### AMC
AMCcutoff <- median(AMC_grant_lable_use$entropy)
#[1] 0.737966
AMC_grant_lable_use$SDI <- "HighSDI"
AMC_grant_lable_use$SDI[AMC_grant_lable_use$entropy <= AMCcutoff] <- "LowSDI"
## CMS1
AMC_grant_lable_use$CMS1.SDI <- as.character(AMC_grant_lable_use$cms)
AMC_grant_lable_use$CMS1.SDI[AMC_grant_lable_use$cms == "CMS1" & AMC_grant_lable_use$entropy > AMCcutoff] <- "CMS1_HighSDI"
AMC_grant_lable_use$CMS1.SDI[AMC_grant_lable_use$cms == "CMS1" & AMC_grant_lable_use$entropy <= AMCcutoff] <- "CMS1_LowSDI"
## CMS2
AMC_grant_lable_use$CMS2.SDI <- as.character(AMC_grant_lable_use$cms)
AMC_grant_lable_use$CMS2.SDI[AMC_grant_lable_use$cms == "CMS2" & AMC_grant_lable_use$entropy > AMCcutoff] <- "CMS2_HighSDI"
AMC_grant_lable_use$CMS2.SDI[AMC_grant_lable_use$cms == "CMS2" & AMC_grant_lable_use$entropy <= AMCcutoff] <- "CMS2_LowSDI"
## CMS3
AMC_grant_lable_use$CMS3.SDI <- as.character(AMC_grant_lable_use$cms)
AMC_grant_lable_use$CMS3.SDI[AMC_grant_lable_use$cms == "CMS3" & AMC_grant_lable_use$entropy > AMCcutoff] <- "CMS3_HighSDI"
AMC_grant_lable_use$CMS3.SDI[AMC_grant_lable_use$cms == "CMS3" & AMC_grant_lable_use$entropy <= AMCcutoff] <- "CMS3_LowSDI"
## CMS4

####### two groups
### SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use, "rfs.delay", "rfs.event", "SDI", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS1.SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use[grep("CMS1", AMC_grant_lable_use$cms), ], "rfs.delay", "rfs.event", "CMS1.SDI", Color=c("#FCD7A2", "#DE9B2A"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS2.SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use[grep("CMS2", AMC_grant_lable_use$cms), ], "rfs.delay", "rfs.event", "CMS2.SDI", Color=c("#A6D2E8", "#106DA8"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS3.SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use[grep("CMS3", AMC_grant_lable_use$cms), ], "rfs.delay", "rfs.event", "CMS3.SDI", Color=c("#E8B5D3", "#C2759F"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS4.SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use[grep("CMS4", AMC_grant_lable_use$cms), ], "rfs.delay", "rfs.event", "self_all_inSDI_5", Color=c("#C2CB1E", "#537A34"), type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
######## more than two groups
### CMS-SDI
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use, "rfs.delay", "rfs.event", "self_all_inSDI_5", Color=CMSSDIcolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))
### CMS
cmssdirfs <- SurvivalFGPlotSimplify(AMC_grant_lable_use, "rfs.delay", "rfs.event", "cms", Color=CMScolor, type="months")
plot_KMCurve(cmssdirfs$clin, cmssdirfs$labs, color=cmssdirfs$Color, font="Helvetica", xlab = "Follow-up (Months)", ylab = paste(cmssdirfs$ylable, "(prob.)", sep=" "))


