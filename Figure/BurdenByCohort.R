library(ggh4x)
library(ggpubr)

X<-read.table("BurdenByCohort.txt", header=TRUE, as.is=TRUE)

# Calculate fold enrichment and confidence interval
X$RR<-X$ObsCount/X$ExpCount
X$RRlower<-vector(mode="numeric", dim(X)[1])
X$RRupper<-vector(mode="numeric", dim(X)[1])
for(ii in seq(1, dim(X)[1])) {
	pt<-poisson.test(X[ii, "ObsCount"], X[ii, "ExpCount"])
	X$RRlower[ii]<-pt$conf.int[1]
	X$RRupper[ii]<-pt$conf.int[2]
}
# Calculate SE for rates
X$RateSE<-with(X, sqrt(ObsCount)/N_Samps)

X$Cohort<-factor(X$Cohort, levels=c("ASC", "MSSNG", "SSC", "SPARK", "Control"))
X$VarClass<-factor(X$VarClass, levels=c("LoF", "Dmis_REVEL0.5", "Silent"),
				labels=c("LoF", "Dmis(REVEL>=0.5)", "Silent"))
X$FamHist<-factor(X$FamHist, levels=c("SPX", "MPX", "Control"), 
				labels=c("No or unknown\nfamily history", "Known\nhistory", ""))

Y<-subset(X, GeneSet=="All")

g<-ggplot(Y, aes(x=Cohort, y=RR, color=Cohort)) +
	geom_point(position=position_dodge(0.6), size=5) + 
	geom_errorbar(aes(ymin=RRlower, ymax=RRupper), position=position_dodge(0.6), width=0.55, size=1.2) +
	geom_hline(yintercept=1, color="darkgray", linetype="dashed") + labs(y="Fold enrichment") +
	scale_color_discrete(labels=c("ASC (N=4076)", "MSSNG (n=3363,\n35% with family history)", 
								  "SSC (n=2654)", "SPARK (n=7015,\n25% with family history)", "Control (n=5764)")) +
	facet_nested(.~VarClass*FamHist, scale="free_x", space="free") + theme_classic() +
	theme(strip.text.x = element_text(size = 11), strip.background=element_rect(color="white"),
		  axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank())

ggsave("FoldEnrich.png", g, height=4, width=10.5, unit="in")


f<-ggplot(Y, aes(x=Cohort, y=ObsRate, fill=Cohort)) + 
	geom_bar(stat="identity",position=position_dodge(0.8)) + 
	geom_errorbar(aes(ymin=ObsRate-RateSE, ymax=ObsRate+RateSE), color="gray", position=position_dodge(0.8), width=0.5) +
	scale_fill_discrete(labels=c("ASC (N=4076)", "MSSNG (n=3363\n35% with family history)", 
								 "SSC (n=2654)", "SPARK (n=7015\n25% with family history)", "Control (n=5764)")) +
	facet_nested(.~VarClass*FamHist, scale="free_x", space="free") + 
	labs(y="Variant rate") + theme_classic() +
	theme(strip.text.x = element_text(size = 11), strip.background=element_rect(color="white"),
		  axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank())  

ggsave("DNVRate.png", f, height=3, width=10.5, unit="in")


