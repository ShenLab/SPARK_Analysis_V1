#
# Power caluclation by Zuk et al. 2014
#
library(ggtext)
library(RColorBrewer)

dir.create("Heatmap")

# Calculate NPC
# n: sample size, lambda: effect size RR-1, f: CAF
NCP<- function(n, lambda, f) {
  return(4*n*((1+lambda)*f*log(1+lambda) + (1 - (1 + lambda)*f)*log((1-(1+lambda)*f)/(1-f))))
}

# Calculate samp size given NPC and lambda/f
# by reversing the above formula
SampSizeByNCP<-function(ncp, lambda, f) {
	return( ncp/(4*((1+lambda)*f*log(1+lambda) + (1 - (1 + lambda)*f)*log((1-(1+lambda)*f)/(1-f)))) )
}

# Power given NCP and signif
PowerByNCP <-function(signif, ncp){
  return(1 - pchisq(qchisq(1 - signif, 1), 1, ncp))
}


# Power given on n, lambda, and f
PowerByParam <- function(n, signif, lambda, f) {
	ncp<-NCP(n, lambda, f)
	return(PowerByNCP(signif, ncp))
}

# Power given genetic parameters
# Based on n (samp size), signif (significance level), RR (relative risk), scoeff (selection coeff), and mu (mutation rate mu)
# under the assumption that a gene does not have pleiotropic effect, scoeff = (RR-1)*prev*sD under-estimated for most genes
PowerByArch<-function(n, signif, RR, scoeff, mu) {
	# First convert mu to f
	f<-mu/scoeff
	# Calculate NPC from RR and f
	ncp<-NCP(n, RR-1, f)
	# Then calculate power
	return(PowerByNCP(signif, ncp))
}

# Solve NPC based on desired power and signif level
SolveNCP <- function(signif, power) {
  return( uniroot(function(ncp) PowerByNCP(signif,ncp)-power, lower = 0, upper = 500, tol=.001)[[1]])
}

# Sample size given genetic parameters
SampSizeByArch <- function(power, signif, RR, scoeff, mu) {
	ncp<-SolveNCP(signif, power)
	return(SampSizeByNCP(ncp, RR-1, mu/scoeff))
}


# Given sample size 32,024, signifance level of 9e-6 (5400 autosomal constrained genes)
# Calculate power as a function of RR (1+lambda), and selection coefficient (scoeff) for genes of different size
# Note that RR should be coupled with scoeff as shown from known and novel ASD genes
# In particular, for RR>10, we cannot expect s to be smaller than 0.1 due to reduece reproductive fitness
# For selection coeff of 0.1, max relative risk is 7 by assuming reduced reprod fitness of 0.71
plot_power<-function(mu=1e-5, RR=seq(1, 20, 0.01), scoeff=seq(0.01, 0.5, 0.001)) {
	sampsize<-32024
	signif<-9e-6
	param<-expand_grid(RR=RR, scoeff=scoeff)
	param$power<-PowerByArch(sampsize, signif, param$RR, param$scoeff, mu)
	g<-param %>% filter(scoeff > RR*1/54*0.71) %>% 
		ggplot(aes(x=RR, y=scoeff, fill=power)) + geom_tile() + 
		scale_fill_distiller(name="Power", palette="Spectral", limits=c(0,1), 
							 breaks=c(1e-3, 0.2, 0.4, 0.6, 0.8, 0.999),
							 labels=c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
		scale_x_continuous(name="Relative risk", breaks=c(1, 5, 10, 15, 20, 25, 30)) + 
		scale_y_continuous(name="Selection coefficient", 
						   breaks=c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
		labs(title=paste0("mu = ", mu)) +
		coord_cartesian(expand=0, ylim=c(0, max(scoeff))) + theme_classic() 
	ggsave(paste0("Heatmap/PowerMu", mu, ".png"), g, height=4, width=4.5)
}

plot_sampsize<-function(mu=1e-5, power=0.8, RR=seq(1, 20, 0.01), scoeff=seq(0.01, 0.5, 0.001)) {
	ncurrent<-32024
	signif<-9e-6
	myPalette <- colorRampPalette(rev(brewer.pal(7, "Spectral")))
	param<-expand_grid(RR=RR, scoeff=scoeff)
	param$sampsize<-SampSizeByArch(power, signif, param$RR, param$scoeff, mu)
	param$sampsize<-ifelse(param$sampsize<0.1*ncurrent, 0.1*ncurrent, param$sampsize)
	f<-param %>% filter(scoeff > RR*1/54*0.71) %>% 
		ggplot(aes(x=RR, y=scoeff, fill=sampsize/ncurrent)) + geom_tile() +
		scale_fill_gradientn(name="Sample size\n(N/32,024)", trans="log2", 
							 limit=c(0.1, 10), breaks=c(1/10, 1/5, 1/2, 1, 2, 5, 10),
							 labels=c("0.1", "0.2", "0.5", "1", "2", "5", "10"), colours = myPalette(7)) + 
		stat_contour(aes(z=log2(sampsize/32024)), color="gray50", breaks=log2(c(1,2,5)), linetype="dotted") +
		scale_x_continuous(name="Relative risk", breaks=c(1, 5, 10, 15, 20, 25, 30)) + 
		scale_y_continuous(name="Selection coefficient", 
						   breaks=c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
		labs(title=paste0("mu = ", mu)) +
		coord_cartesian(expand=0, ylim=c(0, max(scoeff))) + theme_classic()
	ggsave(paste0("Heatmap/SampSizeMu", mu, ".png"), f, height=4, width=4.75)	
}

# Making power plots and sample size plots for genes of different sizes
for(ii in c(seq(1, 12), 15, 20)) {
	plot_power(mu=1e-6*ii)
	plot_sampsize(mu=1e-6*ii, power=0.9)
}

