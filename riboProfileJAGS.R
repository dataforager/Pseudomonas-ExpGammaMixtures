library(coda)
library(rjags)
library(R2WinBUGS)
library(R2jags)

write.readProfile.model <- function(fileName,alphaEst,betaEst) {
    fileObj <- file(fileName)



    fileLines=c("model{","for(i in 1:N){","# assuming indepdence between observations for each gene given cluster label",
        "ld_exp[i] <- logdensity.exp(R[i,1],lambda) + logdensity.exp(R[i,2],lambda) + logdensity.exp(R[i,3],lambda)",
        "ld_gamma[i] <- logdensity.gamma(R[i,1],alpha,beta) + logdensity.gamma(R[i,2],alpha,beta) + logdensity.gamma(R[i,3],alpha,beta)",
        #"density[i] <- exp(ld[i,labels[i]] - log(exp(ld[i,1])+exp(ld[i,2])))",
        "density[i] <- exp(ld_exp[i]*labels[i] + ld_gamma[i]*(1-labels[i]))",
        "# 'ones trick' from Matthew Denwood, Jens Preussner",
        "# http://jenzopr.github.io/stats/2016/04/15/jags-finite-component-mixture.html",
        "ones[i] ~ dbern(density[i])",
        "labels[i] ~ dbern(prob)",
        "}",

    #"muStar <- 2*log(mu) - (1/2)*log( (1/tau) + (mu^2) )",
    #"tauStar <- 1/( log( ( (1/tau) + (mu^2) )/(mu^2) ) )",

    "# elicit priors here",
    #"pis ~ ddirch(c(1,1))",
    "prob ~ dbeta(1,1)",
    "# choose prior params such that prior mean = 0.4605237, corresponds to Pr(R > 5.0) ~ 0.1",
    "# also subject to Pr(lambda < 0.2302344) ~ 0.1, corresponds to Pr(R > 10.0) ~ 0.1",
    "lambda ~ dgamma(5.317999,11.546939)",
    #"lambda ~ dunif(10e-6,10000)",
    "# centered at mean of all data - precision chosen s.t. Pr(mu < 5) ~ 0.1",
    #paste0("mu ~ dnorm(",muMean,",",1/(muSD^2),")"),
    #"tau ~ dunif(1/(5000^2),1)",
    paste0("alpha ~ dunif(1/(10000^2),",10000*(alphaEst-(1/(10000^2))),")"),
    paste0("beta ~ dunif(1/(10000^2),",10000*(betaEst-(1/(10000^2))),")"),
    #"alpha ~ dunif(10e-6,10000)",
    #"beta ~ dunif(10e-6,10000)",
    "}")

    writeLines(fileLines, fileObj)
    close(fileObj)
}

readData <- read.table("/Users/mayhew5/Desktop/riboProfiling/PA33988_glycerol_riboseq.txt",sep="\t",header=TRUE)

# check for and remove rows with duplicated ids
readData <- readData[!duplicated(readData[,1]),]
rownames(readData) <- readData[,1]
readData <- readData[,2:4]
print(readData)

# Check on this with Sarah & Stephanie
readData[which(is.na(readData))] <- 0.0

allReadData <- c(readData[,1],readData[,2],readData[,3])
allReadData <- allReadData[which(allReadData > 0.0)]

# No transcript can be completely unoccupied - arbitrary to help fitting (0.1)
readData[which(readData[,1] == 0.0),1] <- min(allReadData)/2.0
readData[which(readData[,2] == 0.0),2] <- min(allReadData)/2.0
readData[which(readData[,3] == 0.0),3] <- min(allReadData)/2.0

numComps <- 2

R <- as.matrix(readData)
N <- nrow(R)
ones <- rep(1,N)
#labels <- sample(c(1,2),N,replace=TRUE)
#pis <- rep(1/numComps,numComps)

dat = list(R=R,N=N,ones=ones)

initVals<-vector(mode="list",length=1)

# compute init estimates
#start.est <- c()
#start.est[1] <- mean(c(R))
#start.est[2] <- 1/(sd(c(R))^2)
beta.init <- mean(c(R))/(sd(c(R))^2)
alpha.init <- mean(c(R))*beta.init

initVals[[1]]<-list(lambda=0.4605237, alpha=alpha.init, beta=beta.init)

paramsToSave <- c("lambda","alpha","beta","prob","labels","density")

modelFileName <- "readProfile.model"
#muSD <- (5.0 - start.est[1])/(-1.281552)
write.readProfile.model(modelFileName,alpha.init,beta.init)

fun.jags.bv<-jags(data=dat,inits=initVals,parameters.to.save=paramsToSave,
                  model.file=modelFileName,n.chains=1, n.iter=10000,
                  n.burnin=1000,n.thin=1,DIC=F,progress.bar="text")

save(fun.jags.bv,file="PA33988_glycerol_riboseq_param_infs_100516.RData")