## pathway analysis in CMS4 subtype
source("~/R/functions/GSEAHeatmap.R")
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20211014/1.1/CMS4.ThreeCorhot.GSEA.result.RData")
CMS4.GSEA.rslt <- GSEA.rslt
CMS4.heat.dat.p <- CCSEnrichmentHeatmapData(CMS4.GSEA.rslt, P="Pvalue", terms=TRUE)
CMS4.heat.dat <- CMS4.heat.dat.p[, c(2, 1, 3)]

####################################### heatmap
## CMS4.ThreeCorhot
setwd("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/1.picture/CMS4SDIpathway")
pretty_colours <- c('#1E8BA9', '#DA1F27', '#704122')
labelcolors <- list(Corhot = c(ARGO=pretty_colours[1], TCGA=pretty_colours[2], AMC=pretty_colours[3]))
CCSEnrichmentHeatmapRename(CMS4.heat.dat, labelcolors=labelcolors, colorheatmap=colorRampPalette(c("#532788", "#FFFFFF", "#B35807"))(100), cancerintrinsic=TRUE, breaks = c(seq(-2, 2, length.out = 100)))
load("Figure.RData")
pdf("pathway.pdf", height = 9.1, width = 10)
	r1 <- plot_grid(Signature.p$gtable, Pathways.p$gtable, labels = NULL,
	                nrow = 1, align = "h")
	r2 <- plot_grid(Immune.p$gtable, Estimate.p$gtable, labels = NULL,
	                nrow = 2, align = "v")
	r3 <- plot_grid(Metabolism.p$gtable, r2, labels = NULL,
	                nrow = 1, align = "h")
	plot_grid(r1, r3, labels = NULL,
	          nrow = 2, align = "none")
	dev.off()

 ####### ARGO
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230315/1.cluster/ARGO_sdiused_label2.v5.RData")
my_comparisons = list( c("CMS4-IF+", "CMS1"), c("CMS4-IF+", "CMS2"), c("CMS4-IF+", "CMS3"), c("CMS4-IF+", "CMS4-IF-"))
bpi11 <- ggplot(ARGO_sdiused_label2, aes(x=Immune.Subtype, y=Gau_innerAverSDI, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="SDI",x="CMS-EMT", y = "SDI") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


my_comparisons = list( c("CMS4-IF+", "CMS1"), c("CMS4-IF+", "CMS2"), c("CMS4-IF+", "CMS3"), c("CMS4-IF+", "CMS4-IF-"))
bpi22 <- ggplot(ARGO_sdiused_label2, aes(x=Immune.Subtype, y=Infiltrating_STR, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="STR",x="Subtype", y = "STR") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


bpi33 <- ggplot(ARGO_sdiused_label2, aes(x=Immune.Subtype, y=Infiltrating_LYM, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="LYM",x="Subtype", y = "LYM") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


bpi44 <- ggplot(ARGO_sdiused_label2, aes(x=Immune.Subtype, y=TUM, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="TUM",x="Subtype", y = "TUM") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


figure <- ggarrange(bpi11, bpi22, bpi33, bpi44,
                    labels = c("A", "B", "C", "D"),
                    ncol = 4, nrow = 1)
figure

####### TCGA
load("/data0/tan/Task/9.CMS.SDI/Procedures/StartFrom20230606/2.model/CRC.clin.exp.he_v2.RData")
bpi1 <- ggplot(CRC.clin.exp.he_v2, aes(x=Immune.Subtype, y=entropy, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="SDI",x="Subtype", y = "SDI") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


bpi2 <- ggplot(CRC.clin.exp.he_v2, aes(x=Immune.Subtype, y=STR, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="STR",x="Subtype", y = "STR") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


bpi3 <- ggplot(CRC.clin.exp.he_v2, aes(x=Immune.Subtype, y=LYM, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="LYM",x="Subtype", y = "LYM") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


bpi4 <- ggplot(CRC.clin.exp.he_v2, aes(x=Immune.Subtype, y=TUM, fill=Immune.Subtype)) + geom_boxplot(width=0.8) + labs(title="TUM",x="Subtype", y = "TUM") + 
      scale_fill_manual(values=CMSTMESubtype) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      geom_signif(comparisons=my_comparisons, test = "t.test", step_increase = .1, tip_length = 0, vjust=0)


figure2 <- ggarrange(bpi1, bpi2, bpi3, bpi4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 4, nrow = 1)
figure2

