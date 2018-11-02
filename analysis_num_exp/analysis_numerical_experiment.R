## FB 23/08/2018 - Analysis of numerical experiment

file_path = "../numerical_experiment/T=10/"
scenar = c("1_all_data_types/","2_only_counts/","3_all_predator_cmr_prey/","4_all_predator_repro_prey/")

### Alpha parameters
n.occasions = 10
alpha <- rep(NA,8)
N <- seq(10,500,length=n.occasions)

# intra-species DD - juvenile predator survival fn of adult predator abundance
alpha[1] <- 0.5
alpha[2] <- -0.01 # slopes
1/(1+exp(-(alpha[1] + alpha[2] * N)))

# inter-species DD - juvenile prey survival fn of adult predator abundance
alpha[3] <- 0.5
alpha[4] <- -0.025 # slopes
1/(1+exp(-(alpha[3] + alpha[4] * N)))

# inter-species DD - predator fecundity fn of juvenile prey abundance
alpha[5] <- 0
alpha[6] <- 0.004 # slopes
exp(alpha[5] + alpha[6] * N)
# if juv prey gets above 500 we have a problem, otherwise OK

# intra-species DD - prey fecundity fn of adult prey abundance
alpha[7] <- 2 #1.5 #2.5 too much
alpha[8] <- -0.005 # slopes
exp(alpha[7] + alpha[8] * N)


nsims = 100

for (ks in 1:4)
  { # for each scenario
  
  estim_mean <- NULL

        #plotting
        par(mfrow=c(2,4))
        height = c(4,2,3,2,1,1,1,1)

  for (i in 1:nsims)
      {
        load(paste(file_path,scenar[ks],'fitted/EstimIPM_',i,'.RData',sep=""))
        print(i)
        # Estimate bias
        estim_mean <- rbind(estim_mean,ipm.fit$BUGSoutput$mean$alpha)
        
        # Plot figures
        for (j in 1:8){ #for each alpha, it plots in a new panel. 
        #Because there are eight panels, at the end of the loop we're back to panel 1      
        if (i==1){
          y=density(ipm.fit$BUGSoutput$sims.list$alpha[,j])
          plot(y,ylim=c(0,height[j]),lwd=0.5,xlim=c(-3,3),ylab="Pr(alpha|data)", xlab=paste("alpha",j,sep=""),main="",col="grey")
          abline(v=alpha[j],col="red",lwd=3)
          } else {
          lines(density(ipm.fit$BUGSoutput$sims.list$alpha[,j]),lwd=0.5,col="grey")
        }
      }
}
bias_alphas <- (apply(estim_mean,2,mean) - alpha) 

}

