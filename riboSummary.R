
# posterior median and HPD interval summaries for each strain X condition

load("PA33988_alkanes_totalRNA_param_infs_100516.RData")
params <- c("lambda","alpha","prob")

for (param in params) {
	cur.sample <- fun.jags.bv$BUGSoutput$sims.matrix[,param]
	print(param)
	if (param == "lambda") {
		print(median(1/cur.sample))
		print(HPDinterval(as.mcmc(1/cur.sample)))
	}
	else if (param == "alpha") {
		beta.sample <- fun.jags.bv$BUGSoutput$sims.matrix[,"beta"]
		mean.sample <- cur.sample / beta.sample
		print(median(mean.sample))
		print(HPDinterval(as.mcmc(mean.sample)))
		print(median(sqrt(mean.sample / beta.sample)))
		print(HPDinterval(as.mcmc(sqrt(mean.sample / beta.sample))))
	}
	else {
		print(median(cur.sample))
		print(HPDinterval(as.mcmc(cur.sample)))
	}

}

# determine labels of each gene ("low" or "high")
library(coda)
library(VennDiagram)

strains <- c("PA33988","PAO1")
conditions <- c("alkanes","glycerol")
dataTypes <- c("totalRNA","riboseq")

for (i in 1:length(strains)) {
	for (j in 1:length(conditions)) {

		lowTrans <- c()
		lowOcc <- c()
		highTrans <- c()
		highOcc <- c()

		# open full ID list (with gene duplicates)
		fullIDs <- read.delim(paste("/Users/mayhew5/Desktop/riboProfiling/",strains[i],"_",conditions[j],"_fullGeneList_withDups.txt",sep=""),sep="\t",header=TRUE)

		genomeCategories <- matrix(NA,nrow=length(fullIDs[,1]),ncol=2)
		rownames(genomeCategories) <- fullIDs[,1]
		colnames(genomeCategories) <- c("High Tx","High Occ.")		

		for (k in 1:length(dataTypes)) {

			print("STRAIN")
			print(strains[i])
			print("CONDITION")
			print(conditions[j])
			print("DATATYPE")
			print(dataTypes[k])

			load(paste("/Users/mayhew5/Desktop/riboProfiling/",strains[i],"_",conditions[j],"_",dataTypes[k],"_param_infs_100516.RData",sep=""))
			params <- c("lambda","alpha","prob")

			# for (param in params) {
			# 	cur.sample <- fun.jags.bv$BUGSoutput$sims.matrix[,param]
			# 	print(param)
			# 	if (param == "lambda") {
			# 		print(median(1/cur.sample))
			# 		print(HPDinterval(as.mcmc(1/cur.sample)))
			# 	}
			# 	else if (param == "alpha") {
			# 		beta.sample <- fun.jags.bv$BUGSoutput$sims.matrix[,"beta"]
			# 		mean.sample <- cur.sample / beta.sample
			# 		print(median(mean.sample))
			# 		print(HPDinterval(as.mcmc(mean.sample)))
			# 		print(median(sqrt(mean.sample / beta.sample)))
			# 		print(HPDinterval(as.mcmc(sqrt(mean.sample / beta.sample))))
			# 	}
			# 	else {
			# 		print(median(cur.sample))
			# 		print(HPDinterval(as.mcmc(cur.sample)))
			# 	}

			# }

			readData <- read.table(paste("/Users/mayhew5/Desktop/riboProfiling/",strains[i],"_",conditions[j],"_",dataTypes[k],".txt",sep=""),sep="\t",header=TRUE)

			# check for and remove rows with duplicated ids
			readData <- readData[!duplicated(readData[,1]),]
			rownames(readData) <- readData[,1]
			readData <- readData[,2:4]

			paramInfs <- fun.jags.bv$BUGSoutput

			# if PAO1 in glycerol, generate "model fit" figure
			if ((strains[i] == "PAO1") & (conditions[j] == "glycerol") & (dataTypes[k] == "riboseq")) {
				pdf(paste(strains[i],"_",conditions[j],"_modelDemo.pdf",sep=""))
				hist(c(readData[,1],readData[,2],readData[,3]),col="blue",freq=FALSE,main=paste(strains[i],conditions[j],sep=" "),xlab="RPKM",xlim=c(0,1000),breaks=10000)

				curLambda <- median(paramInfs$sims.matrix[,"lambda"])
				curAlpha <- median(paramInfs$sims.matrix[,"alpha"])
				curBeta <- median(paramInfs$sims.matrix[,"beta"])

				curve(dexp(x,curLambda),lwd=3,col="green",add=TRUE)
				curve(dgamma(x,curAlpha,curBeta),lwd=3,col="orange",add=TRUE)
				dev.off()
			}

			matR <- as.matrix(paramInfs)

			lowlist <- c()
			highlist <- c()

			for (r in 1:nrow(readData)) {
				clusterPred <- mean(paramInfs$sims.matrix[,paste0("labels[",r,"]")])
				if (clusterPred > 0.5) {
					lowlist[length(lowlist)+1] <- rownames(readData)[r]
				}
				else {
					highlist[length(highlist)+1] <- rownames(readData)[r]
				}
			}

			if (dataTypes[k] == "totalRNA") {
				lowTrans <- lowlist
				highTrans <- highlist

				for (m in 1:length(rownames(genomeCategories))) {
					if (rownames(genomeCategories)[m] %in% lowTrans) {
						genomeCategories[m,1] <- 0
					}
					else {
						genomeCategories[m,1] <- 1
					}
				}
			}
			else {
				lowOcc <- lowlist
				highOcc <- highlist

				for (m in 1:length(rownames(genomeCategories))) {
					if (rownames(genomeCategories)[m] %in% lowOcc) {
						genomeCategories[m,2] <- 0
					}
					else {
						genomeCategories[m,2] <- 1
					}
				}
			}

			write(lowlist,file=paste(strains[i],"_",conditions[j],"_low_",dataTypes[k],"_Genes_101216.txt",sep=""))
			write(highlist,file=paste(strains[i],"_",conditions[j],"_high_",dataTypes[k],"_Genes_101216.txt",sep=""))


		}

		# write genome category matrix out
		write.table(genomeCategories,paste("/Users/mayhew5/Desktop/riboProfiling/",strains[i],"_",conditions[j],"_genomeCategories_101216.txt",sep=""),quote=FALSE,sep="\t")

		# generate pie chart for breakdown of different types of genes
		highTransHighOcc <- length(intersect(highTrans,highOcc))
		highTransLowOcc <- length(intersect(highTrans,lowOcc))
		lowTransHighOcc <- length(intersect(lowTrans,highOcc))
		lowTransLowOcc <- length(intersect(lowTrans,lowOcc))

		slices <- c(highTransHighOcc, highTransLowOcc, lowTransHighOcc, lowTransLowOcc)
		print(slices)
		lbls <- c("High Tx, High Occ.", "High Tx, Low Occ.", "Low Tx, High Occ.", "Low Tx, Low Occ.")
		pct <- round(slices/sum(slices)*100)
		lbls <- paste(lbls, pct)
		lbls <- paste(lbls,"%",sep="")
		pdf(paste(strains[i],"_",conditions[j],"_funcBreakdown_101216.pdf",sep=""))
		pie(slices, labels = lbls, col=rainbow(length(slices)), main=paste(strains[i],conditions[j],sep=" "),radius=0.6)
		dev.off()

	}
}

# generate Venn diagram within the strain-condition for breakdown (e.g. what proportion of genes are high transcribed, high occ.?)



# Pr(pi_{PAO1,gly} > pi_{PAO1,alk}) and Pr(pi_{33988,gly} > pi_{33988,alk})


# Pr(pi_{PAO1,gly} > pi_{33988,gly}) and Pr(pi_{PAO1,alk} > pi_{33988,alk})


# Single-gene inferences of LT, HT membership if Pr(labels_i == 1 > 0.5), set to LT, generate gene lists
library(VennDiagram)

group1Data <- read.delim("/Users/mayhew5/Desktop/riboProfiling/PAO1_glycerol.csv",sep=",",header=TRUE,row.names=1)
load("/Users/mayhew5/Desktop/riboProfiling/PAO1_glycerol_param_infs_070516.RData")

group1Inf <- fun.jags.bv$BUGSoutput

group2Data <- read.delim("/Users/mayhew5/Desktop/riboProfiling/PAO1_alkanes.csv",sep=",",header=TRUE,row.names=1)
load("/Users/mayhew5/Desktop/riboProfiling/PAO1_alkanes_param_infs_070516.RData")

group2Inf <- fun.jags.bv$BUGSoutput

R1 <- as.matrix(group1Data)
R2 <- as.matrix(group2Data)

R1means <- apply(R1,1,mean)
R2means <- apply(R2,2,mean)

group1LTlist <- c()
group1HTlist <- c()

for (i in 1:nrow(R1)) {
	clusterPred <- mean(group1Inf$sims.matrix[,paste0("labels[",i,"]")])
	if (clusterPred > 0.5) {
		group1LTlist[length(group1LTlist)+1] <- rownames(group1Data)[i]
	}
	else {
		group1HTlist[length(group1HTlist)+1] <- rownames(group1Data)[i]
	}
}

write(group1LTlist,file="PAO1_glycerol_lowTransGenes_072616.txt")
write(group1HTlist,file="PAO1_glycerol_highTransGenes_072616.txt")

group2LTlist <- c()
group2HTlist <- c()

for (i in 1:nrow(R2)) {
	clusterPred <- mean(group2Inf$sims.matrix[,paste0("labels[",i,"]")])
	if (clusterPred > 0.5) {
		group2LTlist[length(group2LTlist)+1] <- rownames(group2Data)[i]
	}
	else {
		group2HTlist[length(group2HTlist)+1] <- rownames(group2Data)[i]
	}
}

write(group2LTlist,file="PAO1_alkanes_lowTransGenes_072616.txt")
write(group2HTlist,file="PAO1_alkanes_highTransGenes_072616.txt")



# 1) Generate Venn diagram of within strain and across strain comparisons

group1HTonly <- intersect(group2LTlist,group1HTlist)
group2HTonly <- intersect(group1LTlist,group2HTlist)

print(length(group1LTlist))
print(length(group2LTlist))

par(cex=1.5)
grid.newpage()
draw.pairwise.venn(length(group2LTlist), length(group1LTlist), length(intersect(group1LTlist,group2LTlist)), category = c("33988-gly (LT)","33988-alk (LT)"), lty = rep("blank", 
    2), fill = c("blue", "light blue"), alpha = rep(0.5, 2), cat.pos = c(45, 
    -45), cat.dist = rep(0.05, 2),cex=2.5,cat.cex=c(1.7,1.7),ext.text=FALSE)

group1LTgroup2HT <- setdiff(group1LTlist,group2LTlist)
group2LTgroup1HT <- setdiff(group2LTlist,group1LTlist)

print(length(group1LTgroup2HT))
print(length(group2LTgroup1HT))

#write.table(group1LTgroup2HT,"PAO1_lowalk_PA33988_highalk",quote=FALSE,col.names=FALSE,row.names=FALSE)
#write.table(group2LTgroup1HT,"PAO1_highalk_PA33988_lowalk",quote=FALSE,col.names=FALSE,row.names=FALSE)


# 2) Pick out candidate genes that "flip" LT-HT status across conditions, strains (rank by log(avg RPKM)) - ASK Sarah and Stephanie which comparison most meaningful