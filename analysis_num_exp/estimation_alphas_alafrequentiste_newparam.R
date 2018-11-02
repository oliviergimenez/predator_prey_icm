## FB 23/08/2018 - Analysis of numerical experiment
# OG 24/08/2018 - added the plots a la frequentiste
# FB 02/09/2018 - New set of parameters
rm(list=ls())
# get parameters used to simulate data
alpha <- rep(NA,8)
alpha[1] <- 0.5
alpha[2] <- -0.01 # slopes
alpha[3] <- 0.5
alpha[4] <- -0.025 # slopes
alpha[5] <- 0.5
alpha[6] <- 0.01 # slopes
alpha[7] <- 1.5 #2.5 too much
alpha[8] <- -0.005 # slopes

# get all estimates for T = 10
file_path = "../numerical_experiment/new_param/T=10/"
scenar = c("1_all_data_types/","2_only_counts/","3_all_predator_cmr_prey/","4_all_predator_repro_prey/")

estim_T10 <- vector("list",4)
nsims <- 100
for (ks in 1:4){ # for each scenario  
	estim_mean <- matrix(NA,nsims,8)
	for (i in 1:nsims){ # for each simulation  	
  	# load each simulation
  	load(paste(file_path,scenar[ks],'fitted/EstimIPM_',i,'.RData',sep=""))
  	# get alpha estimates for each simulation
  	estim_mean[i,] <- ipm.fit$BUGSoutput$mean$alpha
  }
  # get quantiles
  alpha_quantile <- apply(estim_mean,2,quantile, probs=c(.025,.25,.5,.75,.975))
  estim_T10[[ks]] <- alpha_quantile
  estim_T10[[ks]][3,] <- apply(estim_mean,2,mean) ## replacing the median by the mean
}
  
# # get all estimates for T = 30
# file_path = "../numerical_experiment/T=30/"
# scenar = c("1_all_data_types/","2_only_counts/","3_all_predator_cmr_prey/","4_all_predator_repro_prey/")
# 
# estim_T30 <- vector("list",4)
# nsims <- 75
# for (ks in 1:4){ # for each scenario  
# 	estim_mean <- matrix(NA,nsims,8)
# 	for (i in 1:nsims){ # for each simulation  	
#   	# load each simulation
#   	load(paste(file_path,scenar[ks],'fitted/EstimIPM_',i,'.RData',sep=""))
#   	# get alpha estimates for each simulation
#   	estim_mean[i,] <- ipm.fit$BUGSoutput$mean$alpha
#   }
#   # get quantiles
#   alpha_quantile <- apply(estim_mean,2,quantile, probs=c(.025,.25,.5,.75,.975))
#   estim_T30[[ks]] <- alpha_quantile
# }

# unlist and rbind matrices, needed to set ylim
temp <- do.call(rbind, estim_T10)

#png("T10_newparam.png",res=300,width=14,height=8,unit='in') 

pdf("alphas_T10_newparam.pdf",width = 10,height = 8)
par(mfrow=c(2,4))
for (j in 1:8){ #for each alpha, it plots in a new panel
  for (k in 1:4){
    if (k==1) plot(1, estim_T10[[k]][3,j],type="n",xlab="",ylab="",main = substitute(paste(alpha[index]), list(index = j)), 
                   xlim=c(1,4), ylim=c(apply(temp,2,min)[j]-apply(temp,2,sd)[j],apply(temp,2,max)[j]+apply(temp,2,sd)[j]),axes=FALSE)
    
    lines(rep(k,2), estim_T10[[k]][c(1,5),j],lwd=3)
    lines(rep(k,2), estim_T10[[k]][c(2,4),j],lwd=6)
    points(k, estim_T10[[k]][3,j],pch=19,cex=1.5,col="lightblue")
    ticks <- 1:4     
    axis(side=1,at=ticks, labels=c('S 1','S 2','S 3','S 4'))
    #axis(side=1,at=ticks, labels=c('scenario 1','scenario 2','scenario 3','scenario 4'))
    axis(side=2)
  }
  abline(h=alpha[j],col='red')
  box()
}

dev.off()
