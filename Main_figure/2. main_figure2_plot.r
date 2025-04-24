## pathway analysis in 5 subtype
options(stringsAsFactors=FALSE)
source("~/R/functions/GSEAHeatmap.R")
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/3.differential/exp")
library(HTSanalyzeR2)

############## GSEA
load("DEList.RData")
load("/data0/tan/Task/9.CMS.SDI/data/4.CancerIntrinsicSignature/ListGSC_20211028.RData")
GCsignatures <- openxlsx::read.xlsx("/data0/tan/Task/14.amc/4.SIA/20221215/GCsig/41591_2015_BFnm3850_MOESM38_ESM.xlsx")
Ensembl2Name2NCBI <- read.table(paste0("~/Task/Tools/GeneConversion/", "Ensembl2Name2NCBI_20210903", ".txt"), header=T, sep="\t", quote="")
GCSignatures <- lapply(colnames(GCsignatures), function(x){
	genelists <- GCsignatures[, x]
	geneids <- Ensembl2Name2NCBI$NCBI.gene..formerly.Entrezgene..ID[Ensembl2Name2NCBI$Gene.name %in% genelists]
	geneids <- geneids[!is.na(geneids)]
	return(geneids)
})
names(GCSignatures) <- gsub("\\(", "", gsub("\\)", "", colnames(GCsignatures)))
ListGSC[[13]] <- GCSignatures
names(ListGSC)[13] <- "GCSignatures"
load("~/Task/Tools/Ref/CMS.gs.rdata")
load("~/Task/9.CMS.SDI/data/4.CancerIntrinsicSignature/CRC.Cell.gs.Entrez.RData")
ListGSC[[14]] <- CMS.gs
names(ListGSC)[14] <- "CMS.gs"
ListGSC2 <- ListGSC[c(2, 4, 7, 8, 9, 10, 11, 12, 13, 14)]
###
DEIS.de.rslt <- lapply(DEList, function(x){
	x$table
})
DEIS.GSEA.rslt <- CCSEnrichmentHeatmapGSEA(de.rslt=DEIS.de.rslt, ListGSC=ListGSC2, log2FoldChange="logFC", cores=50, filename="CMS.DEIF")
load("CMS.DEIF.GSEA.result.RData")
DEIS.GSEA.rslt <- GSEA.rslt
DEIS.heat.dat.p <- CCSEnrichmentHeatmapData(DEIS.GSEA.rslt, P="Pvalue", terms=TRUE)
DEIS.heat.dat.fdr <- CCSEnrichmentHeatmapData(DEIS.GSEA.rslt, P="Adjusted.Pvalue", terms=TRUE)
###

############## Selection
DEIS.heat.dat.fdr$Genesetname <- sapply(strsplit(rownames(DEIS.heat.dat.fdr), "\\."), function(x) x[1])
table(DEIS.heat.dat.fdr$Genesetname)
DEIS.heat.dat.p$Genesetname <- sapply(strsplit(rownames(DEIS.heat.dat.p), "\\."), function(x) x[1])
table(DEIS.heat.dat.p$Genesetname)

labelcolors <- list(Subtype = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4', 'CMS4-IF- vs CMS4-IF+'='#1E8BA9'))

######  CRCSignature
CRCSignaturedata <- PrepareForHeatmap(heat.dat=DEIS.heat.dat.p, class="CRCSignature", labelcolors=labelcolors)
breaks <- c(seq(-4, 4, length.out = 100))
GeneSets.p <- pheatmap(CRCSignaturedata$heat.dat,
	       #filename = "CRCSignature.pdf",
         fontsize =10, breaks = breaks,
         show_rownames=T,annotation_col = CRCSignaturedata$annotation_col,
				 annotation_colors = CRCSignaturedata$ann_colors,
         color = colorRampPalette(c("#4A7F4B", "#FFFFFF", "#B41600"))(100),
         # border_color = F,
         annotation_names_col=F,show_colnames=T,
         cluster_cols =F,cluster_rows=F,fontsize_col=12,
         cellwidth = 15,cellheight=15,
         main = "CRC Signatures")

load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")  ## CMSclin2, gsvaall2
############################################################ age
CMSclin2$Age <- "50-"
CMSclin2$Age[CMSclin2$age > 50 & CMSclin2$age<= 65] <- "51-65"
CMSclin2$Age[CMSclin2$age > 65 & CMSclin2$age <= 75] <- "66-75"
CMSclin2$Age[CMSclin2$age > 75] <- "76+"
table(CMSclin2$Age)
#  50- 51-65 66-75   76+ 
#  724   350   380   348 

## distribution
### 1.1
fisher.test(table(CMSclin2$CMSTME, CMSclin2$Age), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.0004998
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$Age), workspace = 2e7, simulate.p.value=TRUE)
#0.2839
SDI1.1_Age <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="Age", Sample="sample", adjust=FALSE)
SDI1.1_Age
#            50-   51-65 66-75      76+
#CMS1     0.1610 1.00000 0.947 4.73e-06
#CMS2     0.9840 0.01030 0.105 8.85e-01
#CMS3     0.3880 0.94100 0.171 5.46e-01
#CMS4-IF+ 0.6740 0.00511 0.819 9.30e-01
#CMS4-IF- 0.0252 0.51600 0.688 9.90e-01
SDI1.1_Age <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="Age", Sample="sample", adjust=TRUE)
SDI1.1_Age
#           50-  51-65 66-75      76+
#CMS1     0.489 1.0000 1.000 9.46e-05
#CMS2     1.000 0.0684 0.420 1.00e+00
#CMS3     0.969 1.0000 0.489 1.00e+00
#CMS4-IF+ 1.000 0.0511 1.000 1.00e+00
#CMS4-IF- 0.126 1.0000 1.000 1.00e+00

age1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$Age))
pp1 <- ggplot(age1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("#1E8BA9", "#FECB24", "#DA1F27", "#704122")) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Age") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) # + scale_x_discrete(guide = guide_axis(n.dodge = 2))


############################################################ Gender 

### 1.1
fisher.test(table(CMSclin2$CMSTME, CMSclin2$gender), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.004498
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$gender), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.146
SDI1.1_Gender <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="gender", Sample="sample", adjust=FALSE)
SDI1.1_Gender 
#         female   male
#CMS1     0.0197 1.0000
#CMS2     0.9760 0.0141
#CMS3     0.9570 0.0419
#CMS4-IF+ 0.0607 0.7000
#CMS4-IF- 0.3600 0.5740
SDI1.1_Gender <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="gender", Sample="sample", adjust=TRUE)
SDI1.1_Gender 
#         female   male
#CMS1     0.0987 1.0000
#CMS2     1.0000 0.0987
#CMS3     1.0000 0.1400
#CMS4-IF+ 0.1520 1.0000
#CMS4-IF- 0.7190 0.9560
Gender1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$gender))
pp2 <- ggplot(Gender1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Gender") + ggtitle("p-value = 0.005") + theme(axis.text.x=element_text(angle=60,hjust=1)) # + scale_x_discrete(guide = guide_axis(n.dodge = 2))


############################################################ tnm
CMSclin2$TNM <- CMSclin2$tnm
CMSclin2$TNM <- gsub("A", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("B", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("C", "", CMSclin2$TNM)
CMSclin2$TNM <- gsub("/II", "", CMSclin2$TNM)
table(CMSclin2$TNM, useNA="always")

### 1.1
fisher.test(table(CMSclin2$CMSTME, CMSclin2$TNM), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.0004998
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$TNM), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.2984
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="TNM", Sample="sample", adjust=FALSE)
SDI1.1_TNM 
#              I      II   III      IV
#CMS1     0.5110 0.00323 0.912 1.00000
#CMS2     0.0105 0.59800 0.235 0.21100
#CMS3     0.0153 0.62900 0.942 0.97800
#CMS4-IF+ 0.9970 0.85100 0.273 0.00279
#CMS4-IF- 1.0000 0.96000 0.158 0.19600
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="TNM", Sample="sample", adjust=TRUE)
SDI1.1_TNM 
#              I     II   III     IV
#CMS1     1.0000 0.0323 1.000 1.0000
#CMS2     0.0703 1.0000 0.587 0.5870
#CMS3     0.0764 1.0000 1.000 1.0000
#CMS4-IF+ 1.0000 1.0000 0.606 0.0323
#CMS4-IF- 1.0000 1.0000 0.587 0.5870
TNM1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$TNM))
pp3 <- ggplot(TNM1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("I"=RColorBrewer::brewer.pal(9, "BuPu")[2], "II"=RColorBrewer::brewer.pal(9, "BuPu")[4], "III"=RColorBrewer::brewer.pal(9, "BuPu")[6], "IV"=RColorBrewer::brewer.pal(9, "BuPu")[8])) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "TNM") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 


############################################################ stage

### 1.1
fisher.test(table(CMSclin2$CMSTME, CMSclin2$stage), workspace = 2e7, simulate.p.value=TRUE)
#p-value = p-value = 0.0004998
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$stage), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.5937
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="stage", Sample="sample", adjust=FALSE)
SDI1.1_TumorLocation1 
#              1     2       3        4
#CMS1     0.3010 0.025 0.96900 1.000000
#CMS2     0.5240 0.175 0.84900 0.260000
#CMS3     0.0173 0.727 0.87800 0.981000
#CMS4-IF+ 0.9580 0.962 0.03340 0.000849
#CMS4-IF- 0.9350 0.932 0.00691 0.305000
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="stage", Sample="sample", adjust=TRUE)
SDI1.1_TumorLocation1 
#             1     2      3     4
#CMS1     0.677 0.125 1.0000 1.000
#CMS2     1.000 0.582 1.0000 0.677
#CMS3     0.115 1.000 1.0000 1.000
#CMS4-IF+ 1.000 1.000 0.1340 0.017
#CMS4-IF- 1.000 1.000 0.0691 0.677
TumorLocation11.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$stage))
pp4 <- ggplot(TumorLocation11.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c("1"=RColorBrewer::brewer.pal(9, "YlGn")[2], "2"=RColorBrewer::brewer.pal(9, "YlGn")[4], "3"=RColorBrewer::brewer.pal(9, "YlGn")[6], "4"=RColorBrewer::brewer.pal(9, "YlGn")[8])) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Stage") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1))  #+
      #theme(plot.title = element_text(size = 9, face = "bold"))

########
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
ccc <- CRC.clin.exp.he_v2
colnames(ccc)[colnames(ccc) %in% "Immune.Subtype"] <- "Subtype"
colnames(ccc)[colnames(ccc) %in% "tumor.location"] <- "TumorLocation"

### 1.1
fisher.test(table(ccc$Subtype, ccc$TumorLocation), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.0004998
ccc_CMS4 <- ccc[grep("CMS4", ccc$Subtype), ]
fisher.test(table(ccc_CMS4$Subtype, ccc_CMS4$TumorLocation), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.02399
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=ccc, GroupName1="Subtype", GroupName2="TumorLocation", Sample="sample", adjust=FALSE)
SDI1.1_TumorLocation1 
#                L        R
#CMS1     1.00e+00 9.82e-20
#CMS2     7.63e-07 1.00e+00
#CMS3     9.78e-01 3.25e-02
#CMS4-IF+ 2.67e-01 7.70e-01
#CMS4-IF- 1.29e-04 1.00e+00
SDI1.1_TumorLocation1 <- HyperGeoForTwoGroups(T90Clin=ccc, GroupName1="Subtype", GroupName2="TumorLocation", Sample="sample", adjust=TRUE)
SDI1.1_TumorLocation1 
#                L        R
#CMS1     1.00e+00 9.82e-19
#CMS2     3.81e-06 1.00e+00
#CMS3     1.00e+00 8.12e-02
#CMS4-IF+ 5.33e-01 1.00e+00
#CMS4-IF- 4.29e-04 1.00e+00
TumorLocation11.1box <- as.data.frame(table(ccc$Subtype, ccc$TumorLocation))
pp5 <- ggplot(TumorLocation11.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Dark2")[3:4]) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "Tumor Location") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1))  #+
      #theme(plot.title = element_text(size = 9, face = "bold"))

### pn
fisher.test(table(CMSclin2$CMSTME, CMSclin2$pn), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.0004998
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$pn), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.8841
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pn", Sample="sample", adjust=FALSE)
SDI1.1_TNM 
#              0       1
#CMS1     0.0131 0.99900
#CMS2     0.1180 0.19000
#CMS3     0.1160 0.99000
#CMS4-IF+ 0.9910 0.00606
#CMS4-IF- 1.0000 0.10300
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pn", Sample="sample", adjust=TRUE)
SDI1.1_TNM 
#              0      1
#CMS1     0.0655 1.0000
#CMS2     0.2350 0.3170
#CMS3     0.2350 1.0000
#CMS4-IF+ 1.0000 0.0606
#CMS4-IF- 1.0000 0.2350
TNM1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$pn))
ppn1 <- ggplot(TNM1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(10, "Paired")[2])) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "pn") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1)) 

### pm
fisher.test(table(CMSclin2$CMSTME, CMSclin2$pm), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.0004998
CMSclin2_CMS4 <- CMSclin2[grep("CMS4", CMSclin2$CMSTME), ]
fisher.test(table(CMSclin2_CMS4$CMSTME, CMSclin2_CMS4$pm), workspace = 2e7, simulate.p.value=TRUE)
#p-value = 0.4289
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pm", Sample="sample", adjust=FALSE)
SDI1.1_TNM 
#              0       1
#CMS1     0.0595 1.00000
#CMS2     0.0907 0.21100
#CMS3     0.4840 0.97800
#CMS4-IF+ 0.9840 0.00279
#CMS4-IF- 0.9710 0.19600
SDI1.1_TNM <- HyperGeoForTwoGroups(T90Clin=CMSclin2, GroupName1="CMSTME", GroupName2="pm", Sample="sample", adjust=TRUE)
SDI1.1_TNM 
#             0      1
#CMS1     0.298 1.0000
#CMS2     0.302 0.4220
#CMS3     0.806 1.0000
#CMS4-IF+ 1.000 0.0279
#CMS4-IF- 1.000 0.4220
TNM1.1box <- as.data.frame(table(CMSclin2$CMSTME, CMSclin2$pm))
ppm1 <- ggplot(TNM1.1box, aes(fill=Var2, y=Freq, x=Var1)) + 
    geom_bar(position="fill", stat="identity", width=0.5) + scale_fill_manual(values=c('0'='#E2E2E2','1'=RColorBrewer::brewer.pal(12, "Paired")[12])) + # scale_fill_viridis(discrete = T) +
		xlab("Subtype") + ylab("Percentage") + labs(fill = "pm") + ggtitle("p-value = 0.0005") + theme(axis.text.x=element_text(angle=60,hjust=1))

figure <- ggarrange(pp1, pp2, pp3, pp4, pp5, ppn1, ppm1,
                    labels = c("A", "B", "C", "D", "E", "F", "G"),
                    ncol = 4, nrow = 2)
figure

setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/5.otherclassifier")
CCCRCids <- read.table("CCCRC_offical_label.txt", header=T, sep="\t", quote="")
CRISids <- read.table("CRIS.label.txt", header=T, sep="\t", quote="", fill=T)
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CMSclin2.gsvaall2.RData")  ## CMSclin2, gsvaall2
CCCRCids$id <- gsub("\\.", "-", CCCRCids$id)
length(intersect(CMSclin2$sample, CCCRCids$id))  ## 1244
length(intersect(CMSclin2$sample, CRISids$Sample.ID)) ## 1351
SampleInt <- Reduce(intersect, list(CMSclin2$sample, CCCRCids$id, CRISids$Sample.ID))  ## 1239
CMSclin3 <- CMSclin2[CMSclin2$sample %in% SampleInt, ]
CMSclin3$CCCRC <- CCCRCids$CCCRC[match(CMSclin3$sample, CCCRCids$id)]
CMSclin3$CRIS.TSP <- CRISids$CRIS.TSP.Class.assignment[match(CMSclin3$sample, CRISids$Sample.ID)]
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230514/gsva.all.45.v2.RData")
gsva.all.45.clin3 <- gsva.all.45[, CMSclin3$sample]




##
labelcolors <- list(Subtype = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', 'CMS4-IF+'='#3C5488', 'CMS4-IF-'='#8491B4'))
col_fun2 = colorRamp2(unique(c(seq(-1, 0, length.out=18), seq(0, 1, length.out=18))), colorRampPalette(c("#1283E6", "#FFFFFF", "#FF8B00"))(35))

######
CMSplusha = HeatmapAnnotation(
    CMSplus = CMSclin3$NewCluster, 
    MSI = CMSclin3$msi,
    CIMP = CMSclin3$cimp,
    KRAS = CMSclin3$kras_mut,
    BRAF = CMSclin3$braf_mut,
    col = list(CMSplus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), CCCRC =c(C1="#72598D", C2="#FB9B03", C3="#FB0206", C4="#0765F7"), CRIS =c(Nolabel="#DEDEDE", CRISA="#EF4C2A", CRISB="#D42127", CRISC="#1C275C", CRISD="#018647", CRISE="#00AD9B"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20')),
    na_col = "white"
)
CMSplushaP <- Heatmap(gsva.all.45.clin3[43:45, ], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = CMSplusha, column_labels=rep("", 1239))
#draw(CMSplusha)

CMSclin32 <- CMSclin3[order(CMSclin3$CCCRC), ]
CCCRCha = HeatmapAnnotation(
    CCCRC = CMSclin32$CCCRC, 
    MSI = CMSclin32$msi,
    CIMP = CMSclin32$cimp,
    KRAS = CMSclin32$kras_mut,
    BRAF = CMSclin32$braf_mut,
    col = list(CMSplus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), CCCRC =c(C1="#72598D", C2="#FB9B03", C3="#FB0206", C4="#0765F7"), CRIS =c(Nolabel="#DEDEDE", CRISA="#EF4C2A", CRISB="#D42127", CRISC="#1C275C", CRISD="#018647", CRISE="#00AD9B"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20')),
    na_col = "white"
)
CCCRChaP <- Heatmap(gsva.all.45.clin3[43:45, CMSclin32$sample], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = CCCRCha, column_labels=rep("", 1239))
#draw(CCCRCha)

CMSclin33 <- CMSclin3[order(CMSclin3$CRIS.TSP), ]
CMSclin33 <- CMSclin33[CMSclin33$CRIS.TSP != "" & !is.na(CMSclin33$CRIS.TSP), ]
CRISha = HeatmapAnnotation(
    CRIS = CMSclin33$CRIS.TSP, 
    MSI = CMSclin33$msi,
    CIMP = CMSclin33$cimp,
    KRAS = CMSclin33$kras_mut,
    BRAF = CMSclin33$braf_mut,
    col = list(CMSplus=c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4_HighIF='#3C5488', CMS4_LowIF='#8491B4', NOLBL="#E2E2E2"), CMS = c(CMS1='#E69F24',CMS2='#0273B3', CMS3='#CC79A7', CMS4="#009F73", NOLBL="#E2E2E2"), CCCRC =c(C1="#72598D", C2="#FB9B03", C3="#FB0206", C4="#0765F7"), CRIS =c(Nolabel="#DEDEDE", CRISA="#EF4C2A", CRISB="#D42127", CRISC="#1C275C", CRISD="#018647", CRISE="#00AD9B"), MSI = c(msi='#782069',"mss"='#606EAE'), CIMP = c("CIMP.High"='#37AEB5', "CIMP.Low"='#A0D3B0', 'CIMP.Neg'='#E2E2E2'), KRAS=c('0'='#E2E2E2','1'='#851B20'), BRAF=c('0'='#E2E2E2','1'='#851B20')),
    na_col = "white"
)
CRIShaP <- Heatmap(gsva.all.45.clin3[43:45, CMSclin33$sample], col = col_fun2, cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = CRISha, column_labels=rep("", 1115))
#draw(CRISha)


GetHGPlot <- function(atable, colname1, colname2, samplename){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow=length((table(atable[, colname1]))))
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3, 4))
}
GetHGPlot(atable=CMSclin33, colname1="cms_label", colname2="CRIS.TSP", samplename="sample")
GetHGPlot(atable=CMSclin33, colname1="CMSTME", colname2="CRIS.TSP", samplename="sample")
##

GetHGPlot <- function(atable, colname1, colname2, samplename){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 0.05, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow=length((table(atable[, colname1]))))
	bk = c(seq(0, 0.05, length.out=20), seq(0.051, 1, length.out=20))
	my_palette <- c(rev(colorpanel(n=20, low="#BFD8EA",high="#F0F3FA")), rev(colorpanel(n=20, low="#02519D",high="#BFD8EA")))
	pheatmap(SIAGCadj, scale="none", cluster_row = FALSE, cluster_col = FALSE, color = rev(my_palette), breaks = bk, display_numbers = SIAGCshow, fontsize =13, cellwidth = 30, cellheight=30, gaps_col=c(1, 2, 3, 4), gaps_row=c(1, 2, 3, 4))
}
GetHGPlot(atable=CMSclin32, colname1="cms_label", colname2="CCCRC", samplename="sample")
GetHGPlot(atable=CMSclin32, colname1="CMSTME", colname2="CCCRC", samplename="sample")
#########
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/5.otherclassifier")
GetHGMatrix <- function(atable, colname1, colname2, samplename, filename, nrow=4, adjust=TRUE){
	SIAGCadj <- HyperGeoForTwoGroups(T90Clin=atable, GroupName1=colname1, GroupName2=colname2, Sample=samplename, adjust=adjust)
	SIAGCtable <- as.data.frame.matrix(table(atable[, colname1], atable[, colname2]))
	SIAGCshow <- matrix(ifelse(c(SIAGCadj) < 2, paste0(c(as.matrix(SIAGCtable)), " (", formatC(c(SIAGCadj), format = "e", digits = 2), ")"), c(as.matrix(SIAGCtable))), nrow=nrow)
	SIAGCshow <- as.data.frame(SIAGCshow)
	colnames(SIAGCshow) <- names(table(atable[, colname2]))
	rownames(SIAGCshow) <- names(table(atable[, colname1]))
	print(SIAGCshow)
	write.table(SIAGCshow, paste0(filename, ".txt"), row.names=T, col.names=T, sep="\t", quote=F)
}
GetHGMatrix(atable=CMSclin32, colname1="msi", colname2="CCCRC", samplename="sample", filename="msi2CCCRC", nrow=2, adjust=TRUE)
GetHGMatrix(atable=CMSclin32, colname1="cimp", colname2="CCCRC", samplename="sample", filename="cimp2CCCRC", nrow=3, adjust=TRUE)
GetHGMatrix(atable=CMSclin32, colname1="kras_mut", colname2="CCCRC", samplename="sample", filename="kras2CCCRC", nrow=2, adjust=TRUE)
GetHGMatrix(atable=CMSclin32, colname1="braf_mut", colname2="CCCRC", samplename="sample", filename="braf2CCCRC", nrow=2, adjust=TRUE)

GetHGMatrix(atable=CMSclin3, colname1="msi", colname2="NewCluster", samplename="sample", filename="msi2CMSPlus", nrow=2)
GetHGMatrix(atable=CMSclin3, colname1="cimp", colname2="NewCluster", samplename="sample", filename="cimp2CMSPlus", nrow=3)
GetHGMatrix(atable=CMSclin3, colname1="kras_mut", colname2="NewCluster", samplename="sample", filename="kras2CMSPlus", nrow=2)
GetHGMatrix(atable=CMSclin3, colname1="braf_mut", colname2="NewCluster", samplename="sample", filename="braf2CMSPlus", nrow=2)

GetHGMatrix(atable=CMSclin33, colname1="msi", colname2="CRIS.TSP", samplename="sample", filename="msi2CRIS", nrow=2)
GetHGMatrix(atable=CMSclin33, colname1="cimp", colname2="CRIS.TSP", samplename="sample", filename="cimp2CRIS", nrow=3)
GetHGMatrix(atable=CMSclin33, colname1="kras_mut", colname2="CRIS.TSP", samplename="sample", filename="kras2CRIS", nrow=2)
GetHGMatrix(atable=CMSclin33, colname1="braf_mut", colname2="CRIS.TSP", samplename="sample", filename="braf2CRIS", nrow=2)

CMSclin32$immune <- gsva.all.45[43, CMSclin32$sample]
CMSclin32$stroma <- gsva.all.45[44, CMSclin32$sample]
one.way <- aov(immune~ CCCRC, data = CMSclin32)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 1.829335e-212
one.way <- aov(stroma~ CCCRC, data = CMSclin32)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 2.730782e-156


CMSclin3$immune <- gsva.all.45[43, CMSclin3$sample]
CMSclin3$stroma <- gsva.all.45[44, CMSclin3$sample]
one.way <- aov(immune~ NewCluster, data = CMSclin3)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 4.018496e-91
one.way <- aov(stroma~ NewCluster, data = CMSclin3)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 7.367882e-216

CMSclin33$immune <- gsva.all.45[43, CMSclin33$sample]
CMSclin33$stroma <- gsva.all.45[44, CMSclin33$sample]
one.way <- aov(immune~ CRIS.TSP, data = CMSclin33)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 9.158263e-27
one.way <- aov(stroma~ CRIS.TSP, data = CMSclin33)
summary(one.way)[[1]][["Pr(>F)"]][1]
[1] 2.487168e-30



