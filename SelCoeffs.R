# Plot estimated selection coefficient for known and other genes 
# based on mutation rates and curated CAFs from population data
library(ggpubr)
library(ggrepel)
# We plot mutation rate agains CAFs to get a sense of s

CAF<-read_tsv("Data/MergeCuratedRaw.txt", na=".") %>% 
	filter(LOEUFbin<=2 & str_detect(Filter, "All")) %>% 
	mutate( Known=ifelse(is.na(Known), "Other", "Known"),
			FracDenovo=dnLoFMegaCase*15857/(UnrelCaseInherit*(14208+4512/2)+dnLoFMegaCase*15857),
			nLoFwithInh=UnrelCaseInherit+dnLoFMegaCase,
			nLoFMegaCase=UnrelCase,
			CaseCAF=UnrelCase/(2*32024),
			gnomADexomeCAF=gnomADexome/(2*125748),
			gnomADexomeNonNeuroCAF=gnomADexomeNonNeuro/(2*104068),
			gnomADgenomeCAF=gnomADgenome/(2*76156),
			gnomADgenomeNonNeuroCAF=gnomADgenomeNonNeuro/(2*67442),
			TopMedCAF=TopMed/(2*132345) ) %>%
	select( Known, GeneID, HGNC, LOEUFbin, pAllEnrichMeta, gnomADexomeNonNeuropMegaCC, 
			FracDenovo, nLoFwithInh, nLoFMegaCase, CaseCAF, gnomADexomeCAF, gnomADexomeNonNeuroCAF,
			gnomADgenomeCAF, gnomADgenomeNonNeuroCAF, TopMedCAF )

Cov<-read_tsv("Data/CovPatchFinal.txt") %>% unique()

Mut<-read_tsv("Data/hg19_mutrate_7mer_SPARK30K_NonAdj.txt") %>%
		select(GeneID, HGNC, Mu_LoF)

dir.create("SelCoeffs")

# Merge CAF, Cov and Mut
X<-CAF %>% inner_join(Cov) %>% inner_join(Mut) %>%
	mutate(	gnomADexomeLoFRate=Mu_LoF*gnomADgenomeOver15, gnomADgenomeLoFRate=Mu_LoF*gnomADgenomeOver15,
		    TopMedLoFRate=Mu_LoF*TopMedOver15 )

# First plot CAF against mutation rates, adding lines to demarcating regions with different selection coefficients

# Plotting 
gExome<-ggplot(X, aes(x=gnomADexomeCAF, y=gnomADexomeLoFRate, color=Known)) +
	geom_point(alpha=0.6) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
 	scale_color_manual(values=c("hotpink", "grey40"), labels=c("Known genes", "Other genes")) +
	scale_x_continuous(name="CAF of HC LoFs in gnomAD exomes (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	theme_classic() + 
	theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

gGenome<-ggplot(X, aes(x=gnomADgenomeCAF, y=gnomADgenomeLoFRate, color=Known)) +
	geom_point(alpha=0.6) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
 	scale_color_manual(values=c("hotpink", "grey40"), labels=c("Known genes", "Other genes")) +
	scale_x_continuous(name="CAF of HC LoFs in gnomAD genomes (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	theme_classic() + 
	theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

TopMed<-ggplot(X, aes(x=TopMedCAF, y=TopMedLoFRate, color=Known)) +
	geom_point(alpha=0.6) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
 	scale_color_manual(values=c("hotpink", "grey40"), labels=c("Known genes", "Other genes")) +
	scale_x_continuous(name="CAF of HC LoFs in TopMed (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	theme_classic() + 
	theme(legend.title=element_blank(), legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

Comb<-ggarrange(gExome, gGenome, TopMed, ncol=3, common.legend=TRUE, legend="bottom")
ggsave("SelCoeffs/CAFvsMutRate.png", Comb, height=4.5, width=12, unit="in")



# Then plot subset of genes with more than 10 LoFs in familial samples
Y <- X %>% mutate(FdnBin=case_when( FracDenovo<=0.1 ~ "<=0.1", 
							   		FracDenovo>0.1 & FracDenovo<=0.5 ~ "0.1~0.5",
									FracDenovo>0.5 ~ ">0.5" )) %>%
	mutate(FdnBin=factor(FdnBin, levels=c("<=0.1", "0.1~0.5", ">0.5"))) %>%
	filter(nLoFwithInh > 10)

Highlight<-Y %>% filter(HGNC %in% c("GIGYF1", "KDM5B", "PTEN", "NF1", "ANK2"))

gExome2<-ggplot(Y, aes(x=gnomADexomeCAF, y=gnomADexomeLoFRate, color=FdnBin, shape=Known)) +
	geom_point(size=2) +
	geom_text_repel(data=Highlight, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.6) +
	geom_point(data=Highlight, size=3) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD exomes (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	scale_color_discrete(name="Fraction of de novo LoFs in ASD") +
	scale_shape_discrete(name="", labels=c("Known genes", "Other genes")) +
	guides(color=guide_legend(order=1), shape=guide_legend(order=2)) +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))


gGenome2<-ggplot(Y, aes(x=gnomADgenomeCAF, y=gnomADgenomeLoFRate, color=FdnBin, shape=Known)) +
	geom_point(size=2) +
	geom_text_repel(data=Highlight, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.6) +
	geom_point(data=Highlight, size=3) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD genomes (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	scale_color_discrete(name="Fraction of de novo LoFs in ASD") +
	scale_shape_discrete(name="", labels=c("Known genes", "Other genes")) +
	guides(color=guide_legend(order=1), shape=guide_legend(order=2)) +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))


TopMed2<-ggplot(Y, aes(x=TopMedCAF, y=TopMedLoFRate, color=FdnBin, shape=Known)) +
	geom_point(size=2) +
	geom_text_repel(data=Highlight, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_point(data=Highlight, size=3) +
	geom_abline(slope=sqrt(c(0.01, 0.1, 0.5)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in TopMed (‰)", trans="sqrt",
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4, 6e-4, 8e-4), 
					   labels=c(0.002, 0.02, 0.1, 0.2, 0.4, 0.6, 0.8)) + 
	scale_y_continuous(name="LoF mutation rate", trans="sqrt",  
					   breaks=c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 4e-5)) +
	scale_color_discrete(name="Fraction of de novo LoFs in ASD") +
	scale_shape_discrete(name="", labels=c("Known genes", "Other genes")) +
	guides(color=guide_legend(order=1), shape=guide_legend(order=2)) +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

Comb2<-ggarrange(gExome2, gGenome2, TopMed2, ncol=3, common.legend=TRUE, legend="bottom")
ggsave("SelCoeffs/SelCoefvsFracDenovo.png", Comb2, height=4.5, width=12, unit="in")

# Also try to group genes by sel coef bin, and test their overall frac de novo
Z <- X %>% mutate(SelCoeffgnomADexome =
						case_when( gnomADexomeCAF>=gnomADexomeLoFRate/0.1 ~ "<=0.1", 
								   gnomADexomeCAF<gnomADexomeLoFRate/0.5 ~ ">0.5",
								   TRUE ~ "0.1~0.5" ),
				  SelCoeffgnomADgenome = 
				  		case_when( gnomADgenomeCAF>=gnomADgenomeLoFRate/0.1 ~ "<=0.1", 
								   gnomADgenomeCAF<gnomADgenomeLoFRate/0.5 ~ ">0.5",
							   	   TRUE ~ "0.1~0.5" ),
				  SelCoeffTopMed = 
				  		case_when( TopMedCAF>=TopMedLoFRate/0.1 ~ "<=0.1", 
								   TopMedCAF<TopMedLoFRate/0.5 ~ ">0.5",
							   	   TRUE ~ "0.1~0.5" )) %>%
	mutate( SelCoeffgnomADexome=factor(SelCoeffgnomADexome, levels=c("<=0.1", "0.1~0.5", ">0.5")),
			SelCoeffgnomADgenome=factor(SelCoeffgnomADgenome, levels=c("<=0.1", "0.1~0.5", ">0.5")),	
			SelCoeffTopMed=factor(SelCoeffTopMed, levels=c("<=0.1", "0.1~0.5", ">0.5")) ) 
 
plot_fdnbars <- function(Z, groupname, title) {
	label<-ifelse(str_detect(groupname, "TopMed"), "TopMed", 
				 ifelse(str_detect(groupname, "exome"), "gnomAD exomes", "gnomAD genomes"))
	Z %>% group_by_at(groupname) %>% 
	 	summarize(nLoFs=sum(nLoFwithInh), 
	 			  Fdenovo=sum(FracDenovo*nLoFwithInh, na.rm=TRUE)/sum(nLoFwithInh, na.rm=TRUE)) %>%
	 	ggplot(aes_string(x=groupname, y="Fdenovo", fill=groupname)) + 
	 	geom_bar(stat="identity") + coord_flip(ylim=c(0, 0.8)) + 
	 	labs(y="Fraction of de novo LoFs in ASD cases", x="Selection coeff", title=title) +
	 	theme_classic() + theme(legend.position="none", axis.title.x=element_text(size=10))
} 

gExome3<-plot_fdnbars(Z, "SelCoeffgnomADexome", "gnomAD exomes")
gGenome3<-plot_fdnbars(Z, "SelCoeffgnomADgenome", "gnomAD genomes")
TopMed3<-plot_fdnbars(Z, "SelCoeffTopMed", "TopMed") 

Comb3<-ggarrange(gExome3, gGenome3, TopMed3, ncol=3)
ggsave("SelCoeffs/SelCoefBinvsFDenovo.png", Comb3, height=2, width=12, unit="in")


# For known and novel ASD genes, plot case rate vs control rate
# then adding seletion coeff bins
rr<-function(e1, n1, e2, n2) {
	irr = (e1/n1) / (e2/n2)
	irr
}

rr.lower<-function(e1, n1, e2, n2, conf=0.9) {
	lb = n2/n1 * (e1/(e2+1)) * 1/qf(2*(e2+1), 2*e1, p = (1-conf)/2, lower.tail = FALSE)
	lb
}

rr.upper<-function(e1, n1, e2, n2, conf=0.9) {
	ub = n2/n1 * ((e1+1)/e2) * qf(2*(e1+1), 2*e2, p = (1-conf)/2, lower.tail = FALSE)
	ub
}

# Select known genes with Pdenovo<1e-4, and Pmega<0.05, or novel genes
SelectedGenes <- Z %>% filter(pAllEnrichMeta<1e-4 & gnomADexomeNonNeuropMegaCC<0.05 & LOEUFbin<=2 & !is.na(Known) | 
							  HGNC %in% c("NAV3", "ITSN1", "RALGAPB", "SPTBN1", "SCAF1", "HNRNPUL2", "MARK2", "DSCAM", "MED13") ) %>%
			mutate(CaseCAFSE=sqrt(CaseCAF*(1-CaseCAF)/(2*32024)),
				   gnomADexomeCAFSE=sqrt(gnomADexomeCAF*(1-gnomADexomeCAF)/(2*125748)),
				   gnomADexomeNonNeuroCAFSE=sqrt(gnomADexomeNonNeuroCAF*(1-gnomADexomeNonNeuroCAF)/(2*104068)),
				   gnomADgenomeCAFSE=sqrt(gnomADgenomeCAF*(1-gnomADgenomeCAF)/(2*76156)),
				   gnomADgenomeNonNeuroCAFSE=sqrt(gnomADgenomeNonNeuroCAF*(1-gnomADgenomeNonNeuroCAF)/(2*67442)),
				   TopMedCAFSE=sqrt(TopMedCAF*(1-TopMedCAF)/(2*132345))) 

Highlight2<-SelectedGenes %>% filter(HGNC %in% c("GIGYF1", "KDM5B", "PTEN", "NF1"))

gExome4<-SelectedGenes %>% ggplot(aes(x=gnomADexomeCAF, y=CaseCAF, color=SelCoeffgnomADexome)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=gnomADexomeCAF-gnomADexomeCAFSE, xmax=gnomADexomeCAF+gnomADexomeCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight2, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD exomes (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimateds selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))


gGenome4<-SelectedGenes %>% ggplot(aes(x=gnomADgenomeCAF, y=CaseCAF, color=SelCoeffgnomADgenome)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=gnomADgenomeCAF-gnomADgenomeCAFSE, xmax=gnomADgenomeCAF+gnomADgenomeCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight2, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD genomes (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimated selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

TopMed4<-SelectedGenes %>% ggplot(aes(x=TopMedCAF, y=CaseCAF, color=SelCoeffTopMed)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=TopMedCAF-TopMedCAFSE, xmax=TopMedCAF+TopMedCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight2, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in TopMed (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimated selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

Comb4<-ggarrange(gExome4, gGenome4, TopMed4, ncol=3, common.legend=TRUE, legend="bottom")
ggsave("SelCoeffs/SelCoefbinvsASDRR.png", Comb4, height=4.5, width=12, unit="in")


Highlight3<-SelectedGenes %>% filter(HGNC %in% c("ITSN1", "NAV3", "SCAF1", "HNRNPUL2", "MARK2"))

gExome5<-SelectedGenes %>% ggplot(aes(x=gnomADexomeCAF, y=CaseCAF, color=SelCoeffgnomADexome)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=gnomADexomeCAF-gnomADexomeCAFSE, xmax=gnomADexomeCAF+gnomADexomeCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight3, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD exomes (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimateds selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))


gGenome5<-SelectedGenes %>% ggplot(aes(x=gnomADgenomeCAF, y=CaseCAF, color=SelCoeffgnomADgenome)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=gnomADgenomeCAF-gnomADgenomeCAFSE, xmax=gnomADgenomeCAF+gnomADgenomeCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight3, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in gnomAD genomes (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimated selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

TopMed5<-SelectedGenes %>% ggplot(aes(x=TopMedCAF, y=CaseCAF, color=SelCoeffTopMed)) + 
	geom_point() +  
	geom_errorbar(aes(xmin=TopMedCAF-TopMedCAFSE, xmax=TopMedCAF+TopMedCAFSE), alpha=0.4) +
	geom_errorbar(aes(ymin=CaseCAF-CaseCAFSE, max=CaseCAF+CaseCAFSE), alpha=0.4) +
	geom_text_repel(data=Highlight3, aes(label=HGNC), color="black", size=3, fontface="italic", box.padding=0.8) +
	geom_abline(slope=sqrt(c(1, 5, 25)), color="steelblue", linetype="dashed", intercept=0) + 
	scale_x_continuous(name="CAF of HC LoFs in TopMed (‰)", trans="sqrt", 
					   breaks=c(2e-6, 2e-5, 1e-4, 2e-4, 4e-4), labels=c(0.002, 0.02, 0.1, 0.2, 0.4)) + 
	scale_y_continuous(name="CAF of HC LoFs in ASD cases (‰)", trans="sqrt",
					   breaks=c(5e-5, 2e-4, 5e-4, 1e-3), labels=c(0.05, 0.2, 0.5, 1.0)) +
	scale_color_discrete(name="Estimated selection coeff s ") +
	theme_classic() + 
	theme(legend.position="bottom", legend.text=element_text(size=11), axis.title.x=element_text(size=10))

Comb5<-ggarrange(gExome5, gGenome5, TopMed5, ncol=3, common.legend=TRUE, legend="bottom")
ggsave("SelCoeffs/SelCoefbinvsASDRR2.png", Comb5, height=4.5, width=12, unit="in")

