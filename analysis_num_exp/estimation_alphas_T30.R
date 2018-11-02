## FB 23/08/2018 - Analysis of numerical experiment
# T=30

rm(list=ls())
graphics.off()

file_path = "../numerical_experiment/T=30/"
scenar = c("1_all_data_types/","2_only_counts/","3_all_predator_cmr_prey/","4_all_predator_repro_prey/")

### Alpha parameters
n.occasions = 30
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


nsims = 75

for (ks in 1:4)
  { # for each scenario
  
  estim_mean <-matrix(NA,nsims,8)

        #plotting
        dev.new()
        par(mfrow=c(2,4))
        height = c(1,100,2.5,100,5,1000,5,1000)
        #height = c(1,10,2.5,10,5,1000,5,1000)
        width = c(3,0.1,3,0.1,3,0.01,1,0.01)
        #width = c(3,1,3,1,3,0.01,1,0.01)
        
     
# Plot figures
for (j in 1:8){ #for each alpha, it plots in a new panel. 
        
  for (i in 1:nsims)
      {
        load(paste(file_path,scenar[ks],'fitted/EstimIPM_',i,'.RData',sep=""))
        #print(i)
        # Estimate bias
        estim_mean[i,j] <- ipm.fit$BUGSoutput$mean$alpha[j]
        
        if (i==1){
          y=density(ipm.fit$BUGSoutput$sims.list$alpha[,j])
          plot(y,ylim=c(0,height[j]),lwd=0.5,xlim=c(alpha[j]-width[j],alpha[j]+width[j]),ylab="Pr(alpha|data)", xlab=paste("alpha",j,sep=""),main="",col="grey")
          abline(v=alpha[j],col="red",lwd=3)
          } else {
          lines(density(ipm.fit$BUGSoutput$sims.list$alpha[,j]),lwd=0.5,col="grey")
        }
      }
}

bias_alphas <- (apply(estim_mean,2,mean) - alpha) 
print(bias_alphas)
}

