#
# Plot gene set enrichment results
#
library(ggtext)

# Read func set for plot
blacklist<-c("CELF4Targets_OC0.3", "HighForecASD_0.25", "ASDdnEnrich_TADAPP0.5",
			 "KnownNDDFull", "KnownNDDStrict", "NDDdnEnrich_TADAPP0.5", "HighsHet_0.15")

deltaLoF<-0.0890-0.0349
deltaDmis<-0.1029-0.0666
X<-read_tsv("GeneSetsDnEnrichSummary.txt") %>% 
	filter(VarClass=="Damaging" & NGenesInSet>140 & NGenesInSet<2000 & !GeneSet %in% blacklist ) %>%
	mutate(RateDiffLoF=RateDiffLoF/deltaLoF, RateDiffLoFSE=RateDiffLoFSE/deltaLoF, 
		   RateDiffDmis=RateDiffDmis/deltaDmis, RateDiffDmisSE=RateDiffDmisSE/deltaDmis) %>%
	mutate(RateDiffLoF=ifelse(RateDiffLoF>1, 1, RateDiffLoF),
		   RateDiffDmis=ifelse(RateDiffDmis>1, 1, RateDiffDmis))

Y<-read_tsv("GeneSets_NumGenesInSet.txt") %>%
	filter(!GeneSet %in% blacklist) %>%
	mutate(Class=ifelse(Class=="Transcriptome" | Class=="Proteome", "Transcriptome_Proteome", Class))

#column_to_rownames(Y, var="GeneSet")

Z<-inner_join(X, Y, by="GeneSet") %>% 
	mutate(Background=factor(Background, levels=c("ExclKnown", "AllCons")),
		   Class=factor(Class, levels=c("Transcriptome_Proteome", "Regulome", "Prediction", "Genetic", "Constraint")),
		   Signif=ifelse(Pval<5e-5, "***", ifelse(Pval<5e-4, "**", ifelse(Pval<5e-3, "*", ""))))


# Plot fold enrichment
# Should indicate gene sets that are significantly enriched
classlab<-c("Genetic\nevidence", "Gene constraint", "ASD gene\nprediction", "Neuronal regulome", "Transcriptome &\nproteome")
names(classlab)<-c("Genetic", "Constraint", "Prediction", "Regulome", "Transcriptome_Proteome")

FoldEnrich<-ggplot(Z, aes(x=FoldEnrich, y=Description, color=Background)) + 
	geom_vline(xintercept=1.0, color="salmon", linetype="dashed") +
	geom_errorbar(aes(xmin=CI95Lower, xmax=CI95Upper), position=position_dodge(0.75), width=0) + 
	geom_point(size=2, position=position_dodge(0.75)) + 
	geom_text(aes(x=FoldEnrich, label=Signif), show_guide=FALSE, size=3, position=position_dodge(0.75), vjust="left", hjust="middle") +
	scale_color_manual(values=c("gray", "black"), labels=c("... Excl known genes", "All constrained genes")) +
	facet_grid(Class~., scales="free", space="free", labeller=labeller(Class=classlab), switch="both") + 
	guides(color=guide_legend(reverse=TRUE, nrow=2)) + labs(x="Fold enrichment") +
	coord_cartesian(xlim=c(0.95, 1.15))  + theme_classic() +
	theme(strip.placement="left", strip.text=element_text(size=11), axis.title.y=element_blank(), 
		  legend.title=element_blank(), legend.position="bottom",
		  axis.text.y=element_markdown(), strip.background=element_rect(color="white"))

ggsave("FoldEnrich.png", FoldEnrich, width=4, height=12, unit="in")

varcols<-c('#ca0020', '#f4a582', '#0571b0', '#92c5de')
varlabs<-c("*De novo* LoFs", "... Excl known genes", "*De novo* D-mis (REVEL>=0.5)", "... Excl known genes")

Z1<-Z%>%select(Class,Background,GeneSet,Description,RateDiffLoF,RateDiffLoFSE) %>% 
	rename(RateDiff=RateDiffLoF, RateDiffSE=RateDiffLoFSE) %>% 
	mutate(VarClass=ifelse(Background=="AllCons", "LoF", "LoF_2"), VarClass2="LoF")

Z2<-Z%>%select(Class,Background,GeneSet,Description,RateDiffDmis,RateDiffDmisSE) %>% 
	rename(RateDiff=RateDiffDmis, RateDiffSE=RateDiffDmisSE) %>% 
	mutate(VarClass=ifelse(Background=="AllCons", "Dmis", "Dmis_2"), VarClass2="Dmis")

Z3<-rbind(Z1,Z2) %>% 
	mutate(VarClass=factor(VarClass, levels=c("LoF", "LoF_2", "Dmis", "Dmis_2")),
		   VarClass2=factor(VarClass2, levels=c("LoF", "Dmis")),
		   Class=factor(Class, levels=c("Transcriptome_Proteome", "Regulome", "Prediction", "Genetic", "Constraint")))

FracDiff<-ggplot(Z3, aes(x=RateDiff, y=Description, fill=VarClass)) +
	geom_point(size=2, shape=23, position=position_dodge(0.8)) +
	scale_fill_manual(values=varcols, labels=varlabs) +
	facet_grid(Class~., scales="free", space="free", labeller=labeller(Class=classlab), switch="both") + 
	guides(fill=guide_legend(nrow=2)) + labs(x="% Excess DNVs in cases") + 
	coord_cartesian(xlim=c(0, 1))  + theme_classic() +
	theme(strip.placement="left", axis.title.y=element_blank(), legend.title=element_blank(), legend.position="bottom",
		  axis.text.y=element_markdown(), strip.background=element_rect(color="white"), legend.text=element_markdown(),
		  panel.grid.major = element_line("lightgray",0.25), panel.grid.minor = element_line("lightgray",0.25))

ggsave("FracDiff.png", FracDiff, width=6, height=12, unit="in")


FracDiff2<-ggplot(Z3, aes(x=RateDiff, y=Description, fill=VarClass2)) +
	geom_point(data=Z3 %>% filter(Background=="AllCons"), size=2, shape=23, position=position_dodge(0.75)) +	
	geom_point(data=Z3 %>% filter(Background=="ExclKnown"), aes(fill=VarClass), size=2, shape=23, position=position_dodge(0.75)) +
	scale_fill_manual(values=varcols[c(3,4,1,2)], labels=varlabs[c(3,4,1,2)]) +
	facet_grid(Class~., scales="free", space="free", labeller=labeller(Class=classlab), switch="both") + 
	guides(fill=guide_legend(nrow=2)) +	labs(x="% Excess DNVs in cases") + 
	coord_cartesian(xlim=c(0, 1))  + theme_classic() +
	theme(strip.placement="left", strip.background=element_rect(color="white"),
		  axis.title.y=element_blank(), axis.text.y=element_markdown(), axis.ticks.y=element_blank(),
		  legend.title=element_blank(), legend.position="bottom", legend.text=element_markdown(),
		  panel.grid.major = element_line("lightgray",0.25), panel.grid.minor = element_line("lightgray",0.25))

ggsave("FracDiff2.png", FracDiff2, width=4, height=12, unit="in")



