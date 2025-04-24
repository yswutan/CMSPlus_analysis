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


