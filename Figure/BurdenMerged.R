library(ggpubr)
library(ggh4x)
library(ggtext)

# Add 95 CI for contra table
add_poci<-function(X, level=0.95) {
	X$RateRatioHi<-array()
	X$RateRatioLo<-array()
	for(ii in 1:nrow(X)) {
		res<-poisson.test(c(X$Count1[ii],X$Count2[ii]), c(X$N_Samps1[ii], X$N_Samps2[ii]), conf.level=level)
		X$RateRatioLo[ii]<-res$conf.int[1]
		X$RateRatioHi[ii]<-res$conf.int[2]
	}
	return(X)
}

# Read contra table for given sample group and gene set
# Note that we calculated SE for rate diff
read_tab <- function(file, sp=c("Affected"), gs=paste0("LOEUF_D", seq(1,10))) {
	X<-read_tsv(file, na=".") %>%
		filter( SampGroup1 %in% sp & SampGroup2=="Unaffected" & GeneSet %in% gs & 
				(VarClass=="LoF" | VarClass=="Dmis_REVEL0.5")) %>%
		mutate(GeneSet = factor(GeneSet, levels=gs),
			   VarClass = factor(VarClass, levels=c("LoF", "Dmis_REVEL0.5")),
			   RateDiff = Rate1-Rate2,
		   	   RateDiffSE = sqrt(Count1)/N_Samps1+sqrt(Count2)/N_Samps2 ) %>%
		select(GeneSet, VarClass, SampGroup1, N_Samps1, Count1, SampGroup2, N_Samps2, Count2, RateRatio, RateDiff, RateDiffSE)
	return(add_poci(X))
}

read_tab2 <- function(file, sp=c("Affected"), gs=c("Constrained", "NonCons")) {
	X<-read_tsv(file, na=".") %>% 
		filter( SampGroup1 %in% sp & SampGroup2=="Unaffected" & GeneSet %in% gs & 
				(VarClass=="LoF" | VarClass=="Dmis_REVEL0.5")) %>%
		mutate(GeneSet = factor(GeneSet, levels=gs),
			   VarClass = factor(VarClass, levels=c("LoF", "Dmis_REVEL0.5")),
			   RateDiff = Rate1-Rate2,
		   	   RateDiffSE = sqrt(Count1)/N_Samps1+sqrt(Count2)/N_Samps2 ) %>%
		select(GeneSet, VarClass, SampGroup1, N_Samps1, Count1, SampGroup2, N_Samps2, Count2, RateRatio, RateDiff, RateDiffSE)
	return(add_poci(X))
}


prep_dat <- function(reader, ...) {
	X<-reader("DNVBurdenPlot.contra.txt", ...)
	Y<-reader("DNVBurdenPlotExclKnown.contra.txt", ...)
	levels(Y$VarClass) <- c("LoF_2", "Dmis_REVEL0.5_2")
	Z<-rbind(X, Y)
	Z$VarClass<-with(Z, factor(VarClass, levels=c("LoF", "LoF_2", "Dmis_REVEL0.5", "Dmis_REVEL0.5_2")))
	return(Z)
}

#
# Fold enrichment and rate-diff across LOEUF bins	
# 

# Colors and labels for diff types of variants
varcols<-c('#ca0020', '#f4a582', '#0571b0', '#92c5de')
varlabs<-c("*De novo* LoFs", "... Excl known genes", "*De novo* D-mis (REVEL>=0.5)", "... Excl known genes")

Z<-prep_dat(read_tab)
W<-prep_dat(read_tab2)

FE<-ggplot(Z, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin=RateRatioLo, ymax=RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 10)","[10, 20)", "[20, 30)", "[30, 40)", "[40, 50)", "[50, 60)", "[60, 70)", "[70, 80)", "[80, 90)", "[90, 100)")) + 
		scale_color_manual(values=varcols, labels=varlabs) +
		theme_classic() + labs(y="Rate ratio") + guides(color=guide_legend(nrow=2)) + coord_cartesian(clip="off", ylim=c(0.5, 4.5)) + 
		theme(legend.title=element_blank(), legend.position=c(0.6, 0.85), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.text=element_markdown())

RT<-ggplot(Z, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_hline(yintercept=0, color="dimgray") + 
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		scale_x_discrete(name="gnomAD LOEUF decile (%)", labels=c("[0, 10)","[10, 20)", "[20, 30)", "[30, 40)", "[40, 50)", "[50, 60)", "[60, 70)", "[70, 80)", "[80, 90)", "[90, 100)")) +
		scale_fill_manual(values=varcols, labels=varlabs) +
		labs(y="Rate difference") + guides(fill=guide_legend(nrow=2, override.aes = list(size = 3))) + 
		coord_cartesian(ylim=c(-0.005, 0.05)) + theme_classic() + theme(legend.title=element_blank(), axis.title.y=element_blank(), legend.position=c(0.6, 0.8), legend.text=element_markdown())

FE3<-ggplot(W, aes(x=GeneSet, y=RateRatio, color=VarClass)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") +
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin= RateRatioLo, ymax= RateRatioHi), position=position_dodge(0.75), width=0.3) +
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) +
		scale_color_manual(values=varcols, labels=varlabs) +
		coord_cartesian(clip="off", ylim=c(0.5, 4.5)) + labs(y="Rate ratio vs unaffected") + theme_classic()  + 
		theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_blank())

RT3<-ggplot(W, aes(x=GeneSet, y=RateDiff, fill=VarClass)) +
		geom_bar(position=position_dodge(), width=0.8, stat="identity") +
		geom_errorbar(aes(ymin=RateDiff-RateDiffSE, ymax=RateDiff+RateDiffSE), position=position_dodge(0.8), width=0.25, color="darkgray") +
		geom_hline(yintercept=0, color="darkgray") + labs(y="Excess per case") + 
		scale_x_discrete(name="ExAC pLI", labels=c(">=0.5", "<0.5")) + 
		scale_fill_manual(values=varcols, labels=varlabs) +
		guides(fill = guide_legend(override.aes = list(size=1))) + coord_cartesian(clip="off", ylim=c(-0.005, 0.05)) +
		theme_classic() + theme(legend.title=element_blank(), legend.position="none")

Comb<-ggarrange(FE3, FE, RT3, RT, ncol=2, nrow=2, widths=c(1,3), heights=c(3,2), align="hv")
ggsave("FoldEnrichRate_pLI+LOEUF.png", Comb, height=5, width=10, unit="in")


Comb<-ggarrange(FE3, FE, RT3, RT, ncol=2, nrow=2, widths=c(1,3), heights=c(2,2), align="hv")
ggsave("FoldEnrichRate_pLI+LOEUF2.png", Comb, height=4.4, width=10, unit="in")

#
# Further compare low and high functional groups
#

Z2 <-prep_dat(read_tab, sp=c("AffLowFunc", "AffHighFunc"), gs=c("LOEUF_D1", "OtherCons") ) %>%
	 mutate(SampGroup1=factor(SampGroup1, levels=c("AffLowFunc", "AffHighFunc")),
	 		Rate1=Count1/N_Samps1,
	 		Rate1SE=sqrt(Count1)/N_Samps1,
	 		VarGroup=ifelse(VarClass %in% c("LoF_2", "Dmis_REVEL0.5_2"), "After", "Before"),
	 		VarClass2=ifelse(VarClass %in% c("LoF", "LoF_2"), "LoF", "Dmis_REVEL0.5")) %>%
	 mutate(VarGroup=factor(VarGroup, levels=c("Before", "After")))


LOEUFbin<-c("Top 10% LOEUF", "Other")
names(LOEUFbin)<-c("LOEUF_D1", "OtherCons")
Genebin<-c("ExAC pLI>=0.5 or top 20% LOEUF", "... After excluding known genes")
names(Genebin)<-c("Before", "After")

LOEUFbin2<-c("Top 10% LOEUF", "Other")
names(LOEUFbin2)<-c("LOEUF_D1", "OtherCons")
Genebin2<-c("", "")
names(Genebin2)<-c("Before", "After")

#Z2Top <-Z2 %>% filter(GeneSet=="LOEUF_D1" & !str_detect(VarClass, "_2")) %>% 
#		select(SampGroup1, VarClass, RateRatio) %>% spread(SampGroup1, RateRatio)
Z2$VarClass2<-with(Z2, factor(VarClass2, levels=c("LoF", "Dmis_REVEL0.5")))

FE2<-ggplot(Z2, aes(x=SampGroup1, y=RateRatio, color=VarClass2)) + 
		geom_hline(yintercept=1, color="dimgray", linetype="dashed") + 
		geom_point(position=position_dodge(0.75), size=3) + 
		geom_errorbar(aes(ymin=RateRatioLo, ymax=RateRatioHi), position=position_dodge(0.75), width=0.3) +
		labs(y="Rate ratio") + guides(color=guide_legend(nrow=2)) + coord_cartesian(ylim=c(0, 8), clip="off") + 
		scale_color_manual(values=varcols[c(1,3)], labels=varlabs[c(1,3)]) +
		scale_x_discrete(labels=c("LoFunc", "HiFunc"))  +
		facet_nested(.~VarGroup*GeneSet, scale="free_x", space="free", labeller=labeller(GeneSet=LOEUFbin, VarGroup=Genebin2), switch="both") + 
		theme_classic() + 
		theme(legend.title=element_blank(), legend.position=c(0.75, 0.9), axis.title.x=element_blank(), 
			  strip.background=element_rect(color="white"), strip.placement="outside", legend.text=element_markdown())

RT2<-ggplot(Z2, aes(x=SampGroup1, y=Rate1, fill=VarClass2)) + 
	geom_hline(yintercept=0, color="dimgray") + 
	geom_bar(position=position_dodge(), width=0.8, stat="identity") + 
	geom_errorbar(aes(ymin=Rate1-Rate1SE, ymax=Rate1+Rate1SE), position=position_dodge(0.8), width=0.25, color="darkgray") +
	facet_nested(.~VarGroup*GeneSet, scale="free_x", space="free", labeller=labeller(GeneSet=LOEUFbin, VarGroup=Genebin), switch="both") +
	scale_fill_manual(values=varcols[c(1,3)], labels=varlabs[c(1,3)]) +
	scale_x_discrete(labels=c("LoFunc", "HiFunc")) +
	labs(y="DNV per case") + guides(fill=guide_legend(nrow=2, override.aes = list(size = 3))) + 
	theme_classic() + 
	theme(strip.placement="outside", strip.background=element_rect(color="white"), 
		  axis.title.x=element_blank(), legend.title=element_blank(), legend.position=c(0.75, 0.95), legend.text=element_markdown())
		
Comb<-ggarrange(FE2, RT2, nrow=2, align="hv", heights=c(2.5, 3))
ggsave("HiLowFunc_AllCons.png", Comb, height=5.5, width=5, unit="in")

	



