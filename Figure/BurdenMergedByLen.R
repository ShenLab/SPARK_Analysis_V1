library(ggpubr)
library(ggh4x)
library(ggtext)

# Add 95 CI for contra table
add_poci<-function(X, level=0.95) {
	X$RateRatioHi<-array()
	X$RateRatioLo<-array()
	X$Pvalue<-array()
	for(ii in 1:nrow(X)) {
		res<-poisson.test(c(X$Count1[ii],X$Count2[ii]), c(X$N_Samps1[ii], X$N_Samps2[ii]), conf.level=level)
		X$RateRatioLo[ii]<-res$conf.int[1]
		X$RateRatioHi[ii]<-res$conf.int[2]
		X$Pvalue[ii]<-res$p.value
	}
	return(X)
}

# Read contra table for given sample group and gene set
# Note that we calculated SE for rate diff
read_tab <- function(file, gs=c("ConsLong", "NonConsLong")) {
	X<-read_tsv(file, na=".") %>%
		filter( SampGroup1=="Affected" & SampGroup2=="Unaffected" & GeneSet %in% gs & 
				(VarClass=="LoF" | VarClass=="Dmis_REVEL0.5")) %>%
		mutate(GeneSet = factor(GeneSet, levels=gs),
			   VarClass = factor(VarClass, levels=c("LoF", "Dmis_REVEL0.5")),
			   RateDiff = Rate1-Rate2,
		   	   RateDiffSE = sqrt(Count1)/N_Samps1+sqrt(Count2)/N_Samps2 ) %>%
		select(GeneSet, VarClass, SampGroup1, N_Samps1, Count1, SampGroup2, N_Samps2, Count2, RateRatio, RateDiff, RateDiffSE)
	return(add_poci(X))
}

prep_dat <- function(...) {
	X<-read_tab("DNVBurdenPlotStratGeneLen.contra.txt", ...)
	Y<-read_tab("DNVBurdenPlotStratGeneLenExclKnown.contra.txt", ...)
	levels(Y$VarClass) <- c("LoF_2", "Dmis_REVEL0.5_2")
	Z<-rbind(X, Y)
	Z$VarClass<-with(Z, factor(VarClass, levels=c("LoF", "LoF_2", "Dmis_REVEL0.5", "Dmis_REVEL0.5_2")))
	return(Z)
}

varcols<-c('#ca0020', '#f4a582', '#0571b0', '#92c5de')
varlabs<-c("*De novo* LoFs", "... Excl known genes", "*De novo* D-mis (REVEL>=0.5)", "... Excl known genes")


W<-prep_dat(gs=c("ConsLong", "NonConsLong"))
FELong<-ggplot(W, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin= RateRatioLo, ymax= RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) +
		scale_color_manual(values=varcols, labels=varlabs) +
		coord_cartesian(clip="off", ylim=c(0.5, 4)) + labs(y="Rate ratio vs unaffected") + theme_classic()  + 
		theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank())

RTLong<-ggplot(W, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		geom_hline(yintercept=0, color="darkgray") + labs(y="Excess per case") + 
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) + 
		scale_fill_manual(values=varcols, labels=varlabs) +
		guides(fill = guide_legend(override.aes = list(size=1))) + coord_cartesian(clip="off", ylim=c(-0.005, 0.04)) +
		theme_classic() + theme(legend.title=element_blank(), legend.position="none")

Z<-prep_dat(gs=c("ConsShort", "NonConsShort"))
FEShort<-ggplot(Z, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin= RateRatioLo, ymax= RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) +
		scale_color_manual(values=varcols, labels=varlabs) +
		coord_cartesian(clip="off", ylim=c(0.5, 4)) + labs(y="Rate ratio vs unaffected") + theme_classic()  + 
		theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank())

RTShort<-ggplot(Z, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		geom_hline(yintercept=0, color="darkgray") + labs(y="Excess per case") + 
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) + 
		scale_fill_manual(values=varcols, labels=varlabs) +
		guides(fill = guide_legend(override.aes = list(size=1))) + coord_cartesian(clip="off", ylim=c(-0.005, 0.04)) +
		theme_classic() + theme(legend.title=element_blank(), legend.position="none")


X<-prep_dat(gs=c("LOEUF_D1+D2_Long", "LOEUF_D3+D4_Long", "LOEUF_D5+D6_Long", "LOEUF_D7-D10_Long"))
FELong2<-ggplot(X, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin=RateRatioLo, ymax=RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 20)","[20, 40)", "[40, 60)", "[60, 100)")) + 
		scale_color_manual(values=varcols, labels=varlabs) +
		theme_classic() + labs(y="Rate ratio") + guides(color=guide_legend(nrow=2)) + coord_cartesian(clip="off", ylim=c(0.5, 4)) + 
		theme(legend.title=element_blank(), legend.position=c(0.6, 0.9), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.text=element_markdown())

RTLong2<-ggplot(X, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_hline(yintercept=0, color="dimgray") + 
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 20)","[20, 40)", "[40, 60)", "[60, 100)")) +
		scale_fill_manual(values=varcols, labels=varlabs) +
		labs(y="Rate difference") + guides(fill=guide_legend(nrow=2, override.aes = list(size = 3))) + 
		coord_cartesian(ylim=c(-0.005, 0.04)) + theme_classic() + theme(legend.title=element_blank(), axis.title.y=element_blank(), legend.position=c(0.6, 0.9), legend.text=element_markdown())

Y<-prep_dat(gs=c("LOEUF_D1+D2_Short", "LOEUF_D3+D4_Short", "LOEUF_D5+D6_Short", "LOEUF_D7-D10_Short"))
FEShort2<-ggplot(Y, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin=RateRatioLo, ymax=RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 20)","[20, 40)", "[40, 60)", "[60, 100)")) + 
		scale_color_manual(values=varcols, labels=varlabs) +
		theme_classic() + labs(y="Rate ratio") + guides(color=guide_legend(nrow=2)) + coord_cartesian(clip="off", ylim=c(0.5, 4.5)) + 
		theme(legend.title=element_blank(), legend.position=c(0.6, 0.9), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.text=element_markdown())

RTShort2<-ggplot(Y, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_hline(yintercept=0, color="dimgray") + 
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 20)","[20, 40)", "[40, 60)", "[60, 100)")) +
		scale_fill_manual(values=varcols, labels=varlabs) +
		labs(y="Rate difference") + guides(fill=guide_legend(nrow=2, override.aes = list(size = 3))) + 
		coord_cartesian(ylim=c(-0.005, 0.04)) + theme_classic() + theme(legend.title=element_blank(), axis.title.y=element_blank(), legend.position=c(0.6, 0.9), legend.text=element_markdown())



CombLong<-ggarrange(FELong, FELong2, RTLong, RTLong2, ncol=2, nrow=2, widths=c(1,2), heights=c(3,2), align="hv")
ggsave("FoldEnrichRate_pLI+LOEUF_Long.png", CombLong, height=5, width=8, unit="in")


CombShort<-ggarrange(FEShort, FEShort2, RTShort, RTShort2, ncol=2, nrow=2, widths=c(1,2), heights=c(3,2), align="hv")
ggsave("FoldEnrichRate_pLI+LOEUF_Short.png", CombShort, height=5, width=8, unit="in")


CombTab<-rbind(X, Y, W, Z)
write_tsv(CombTab, "BurdenMergedByLen.txt")

