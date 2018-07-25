#!/usr/bin/r
pdf("figures/R_plots/plot_degrees.pdf")
tab = read.table("figures/R_plots/files/degrees.txt")
plot(tab$V1,tab$V4,type="l",col="black",lwd=3,xlab="Degree",ylab="Proportion of genes")
points(tab$V1,tab$V3,type="l",col="red",lwd=6)
points(tab$V1,tab$V2,type="l",col="blue",lwd=6)
legend(2.3,0.6,c("Theoretical ideal genomes","PNJ","RAW"),lwd=6,col=c("black","red","blue"),bg="white")
dev.off()

pdf("figures/R_plots/boxplot_content.pdf")
r=read.table("figures/R_plots/files/content.txt")
boxplot(r$V1,r$V2,r$V3,xlab="Gene sets",ylab="Number of genes",names=c("Extant","anc_RAW","anc_PNJ"))
dev.off()


# boxplot_scaffolds:
d<-read.table(file="figures/R_plots/files/boxplot_scaffolds.txt",header=FALSE, sep="\t", dec=",", comment.char="!",na.string="-99")
x1= factor(d[,1], levels=c("XNS ext", "XNS anc", "X ext", "X anc", "WG ext", "WG anc"))
pdf("figures/R_plots/boxplot_scaffolds.pdf")
boxplot(d[,2]~x1, plot=TRUE, main="Number of scaffolds in extant and ancestral genomes")
dev.off()


# boxplot_rearrangements:
d<-read.table(file="figures/R_plots/files/boxplot_rearrangements.txt",header=FALSE, sep="\t", dec=",", comment.char="!",na.string="-99")
x1= factor(d[,1], levels=c("X phylogeny, non scaffold", "X phylogeny, scaffold"))
pdf("figures/R_plots/boxplot_rearrangements.pdf")
boxplot(d[,2]~x1, plot=TRUE, main="Number of rearrangements per branch")
dev.off()
