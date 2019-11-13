#!/anaconda2/envs/reseq/bin/Rscript
args <- commandArgs(T)

if (length(args) != 5) {
	print("Usage: Rscript Fst_estimate.r plink_prefix (bed, bim.fam files needed) group_file group_col (specify which columns containing group information,  Default <Group>) output_prefix generateplot <Default 0 means no plot, acceplt 1 to do manhatan plot>")
	q()
}

plink <- args[1]
group <- args[2]
group_name <- args[3]
outf <- args[4]
generateplot <- args[5]

library("snpStats")
snp.test<-read.plink(paste(plink,".bed",sep=""),paste(plink,".bim",sep=""),paste(plink,".fam",sep=""),na.strings = c("0", "-9"), sep = "." )

clade<- read.table(group, header = T, row.names=1, sep="\t")

clade <- clade[match(rownames(snp.test$genotypes), rownames(clade)),, drop=F]

fac<-factor(clade[[group_name]])
fst.test<-Fst(snp.test$genotypes,fac)
weighted.mean(fst.test$Fst, fst.test$weight)
X<-cbind(snp.test$map[,0],snp.test$map[,1],snp.test$map[,4],fst.test$Fst,fst.test$weight)
colnames(X)<-c("Chr","Pos","Fst","Weight")
write.table(X,paste(outf,"snpStats.Fst.out",sep="."),sep = "\t",col.names = T,quote=F)

if(generateplot == 1) {
	library(qqman)
	pdf(paste(outf, "snpStats.FST_chr.pdf", sep="."))
	nn<-which(X$Fst>0)
	data<-X[nn,]
	manhattan(data, chr = "Chr", bp = "Pos", p = "Fst",
			  col = c("blue", "green4"), chrlabs = c(1:16),
			  suggestiveline = F, genomewideline = F,
			  cex = 0.8,cex.axis = 0.9,
			  highlight = NULL, logp = F, annotatePval = NULL,
			  annotateTop = F,ylab="Fst")
	dev.off()
	q()
}
