# Draw cumulative distribution of LoF rate per gene
# and highlight known genes
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(readr)

# Making CDF plots for gene mut rates
MutRate<-read_tsv("Data/hg19_mutrate_7mer_SPARK30K.txt", na=".")
MutRate<-subset(MutRate, !is.na(CytoBand) & !str_detect(CytoBand, "X"))

Known<-read_tsv("Data/FullKnownGenes.txt")
KnownRate<-MutRate %>% filter(GeneID %in% Known$GeneID)

Selected<-read_tsv("Data/CandGenes.txt")
SelectedRate<-MutRate %>% filter(GeneID %in% Selected$GeneID)

dir.create("GeneLen")

#Marks<-c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
Marks<-c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

Quant<-quantile(MutRate$Mu_LoF, Marks)
Seg<-data.frame(x=c(rep(10^-8, length(Marks)), Quant), 
				xend=rep(Quant, 2), 
				y=rep(Marks, 2), 
				yend=c(Marks, rep(-0.05, length(Marks))))

# CDF of LoF rate, showing that known genes are biased toward large genes
LoF<-ggplot(MutRate, aes(Mu_LoF)) + stat_ecdf(geom="step") + 
	stat_ecdf(data=KnownRate, color="salmon", geom="step") + 
	geom_segment(data=Seg, aes(x=x, xend=xend, y=y, yend=yend), linetype="dotted") +
#	scale_x_continuous(trans="log10", breaks=c(2.0e-7, 4.9e-7, 6.8e-7, 1.2e-6, 2.0e-6, 3.4e-6, 5.4e-6, 7.2e-6, 1.3e-5)) +
#	scale_y_continuous(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) +
	scale_x_continuous(trans="log10", breaks=c(2.0e-7, 4.9e-7, 1.2e-6, 2.0e-6, 3.4e-6, 7.2e-6, 1.3e-5),
					   labels=c("2e-7", "4.9e-7", "1.2e-6", "2e-6", "3.4e-6", "7.2e-6", "1.3e-5")) +
	scale_y_continuous(breaks=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) +
	coord_cartesian(xlim=c(10^-7, 10^-4.6), ylim=c(0,1)) +
	labs(x="LoF mutation rate", y="Quantile") +
	theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

MutRate2 <- MutRate %>% filter(Mu_Dmis_REVEL0.5>1e-9)
Marks2<-c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
Quant2<-quantile(MutRate2$Mu_Dmis_REVEL0.5, Marks2)
Seg2<-data.frame(x=c(rep(10^-9, length(Marks2)), Quant2), 
				 xend=rep(Quant2, 2), 
				 y=rep(Marks2, 2), 
				 yend=c(Marks2, rep(-0.05, length(Marks2))) )

Dmis<-ggplot(MutRate, aes(Mu_Dmis_REVEL0.5)) + stat_ecdf(geom="step") + 
	stat_ecdf(data=KnownRate, color="salmon", geom="step") + 
	scale_x_continuous(trans="log10", breaks=c(6.9e-8, 9.9e-7, 2.6e-6, 5.6e-6, 1.4e-5, 2.9e-5)) +
	scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) +
	geom_segment(data=Seg2, aes(x=x,xend=xend,y=y,yend=yend), linetype="dotted") +
	labs(x="Dmis(REVEL>=0.5) rate", y="Quantile") + 
	coord_cartesian(xlim=c(10^-8.1, 10^-4.2), ylim=c(0,1)) + 
	labs(x="Dmis (REVEL>=0.5) mutation rate", y="Quantile") +
	theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Comb<-ggarrange(LoF, Dmis, ncol=2)

ggsave("GeneLen/MutRateCDFs.png", Comb, width=9, height=4.5, unit="in")

LoF2<-ggplot(MutRate, aes(Mu_LoF)) + stat_ecdf(geom="step") + 
	stat_ecdf(data=SelectedRate, color="orange", geom="step") + 
	geom_segment(data=Seg, aes(x=x, xend=xend, y=y, yend=yend), linetype="dotted") +
#	scale_x_continuous(trans="log10", breaks=c(2.0e-7, 4.9e-7, 6.8e-7, 1.2e-6, 2.0e-6, 3.4e-6, 5.4e-6, 7.2e-6, 1.3e-5)) +
#	scale_y_continuous(breaks=c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) +
	scale_x_continuous(trans="log10", breaks=c(2.0e-7, 4.9e-7, 1.2e-6, 2.0e-6, 3.4e-6, 7.2e-6, 1.3e-5),
					   labels=c("2e-7", "4.9e-7", "1.2e-6", "2e-6", "3.4e-6", "7.2e-6", "1.3e-5")) +
	scale_y_continuous(breaks=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) +
	coord_cartesian(xlim=c(10^-7, 10^-4.6), ylim=c(0,1)) +
	labs(x="LoF mutation rate", y="Quantile") +
	theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Dmis2<-ggplot(MutRate, aes(Mu_Dmis_REVEL0.5)) + stat_ecdf(geom="step") + 
	stat_ecdf(data=SelectedRate, color="orange", geom="step") + 
	scale_x_continuous(trans="log10", breaks=c(6.9e-8, 9.9e-7, 2.6e-6, 5.6e-6, 1.4e-5, 2.9e-5)) +
	scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) +
	geom_segment(data=Seg2, aes(x=x,xend=xend,y=y,yend=yend), linetype="dotted") +
	labs(x="Dmis(REVEL>=0.5) rate", y="Quantile") + 
	coord_cartesian(xlim=c(10^-8.1, 10^-4.2), ylim=c(0,1)) + 
	labs(x="Dmis (REVEL>=0.5) mutation rate", y="Quantile") +
	theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Comb2<-ggarrange(LoF2, Dmis2, ncol=2)

ggsave("GeneLen/MutRateCDFs2.png", Comb2, width=9, height=4.5, unit="in")


