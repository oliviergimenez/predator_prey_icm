## FB 03/09/2018 - redo plots for Fig 1 and 2
library(mcmcplots)

source("figlabel.R")

# get parameters used to simulate data
alpha <- rep(NA,8)
alpha[1] <- 0.5
alpha[2] <- -0.01 # slopes
alpha[3] <- 0.5
alpha[4] <- -0.025 # slopes
alpha[5] <- 0
alpha[6] <- 0.004 # slopes
alpha[7] <- 2 #1.5 #2.5 too much
alpha[8] <- -0.005 # slopes


# get all estimates for T = 10
file_path = "../numerical_experiment/T=10/"
scenar = c("1_all_data_types/","2_only_counts/","3_all_predator_cmr_prey/","4_all_predator_repro_prey/","0_simulation/")

ks=1 #scenario
i = 1 # first simulation

# retrieve info
load(paste(file_path,scenar[ks],'fitted/EstimIPM_',i,'.RData',sep=""))

# summarize posteriors
print(ipm.fit, digits = 2)

# traceplots
traplot(ipm.fit,c("mean.pP", "alpha", "mean.pV"))

# posterior densities
denplot(ipm.fit,c("mean.pP", "alpha", "mean.pV"))

# compare actual vs estimated abundance
n.occasions = 10

# load simulated data
load(paste(file_path,scenar[5],'sim_data/SimulationIPM_',i,'.RData',sep=""))

pdf(file="Fig1.pdf",width=8,height=8)
par(mfrow=c(2,1))
# real vs estimated predator abundance
plot(1:n.occasions,jags.data.out $countsP,type='l',lwd=3,col='blue',ylim=c(0,max(jags.data.out $countsP)),ylab='counts',xlab='years',main='predator')
lines(1:n.occasions, ipm.fit$BUGSoutput$mean$NtotP,type='l',lwd=3,col='green')
legend('bottomright',col=c('blue','green'),legend=c('truth','estimated'),lty=1,lwd=3)
fig_label("A",cex=2)
# real vs estimated prey abundance
plot(1:n.occasions,jags.data.out $countsV,type='l',lwd=3,col='blue',ylim=c(0,max(jags.data.out $countsV)),ylab='counts',xlab='years',main='prey')
lines(1:n.occasions, ipm.fit$BUGSoutput$mean$NtotV,type='l',lwd=3,col='green')
legend('bottomright',col=c('blue','green'),legend=c('truth','estimated'),lty=1,lwd=3)
fig_label("B",cex=2)
dev.off()

### Figure 2

# real vs estimated parameters involved in DD relationships 
alpha_est <- ipm.fit$BUGSoutput$mean$alpha
data.frame(real=alpha,estimated=alpha_est )

# plot DD relationships (actual and estimated)
N <- seq(10,500,length=n.occasions) #density index

# intra-species DD - juvenile predator survival fn of adult predator abundance
surv_juvP_intrasp = 1/(1+exp(-(alpha[1] + alpha[2] * N)))
surv_juvP_intrasp_est = 1/(1+exp(-(alpha_est[1] + alpha_est[2] * N)))

# inter-species DD - juvenile prey survival fn of adult predator abundance
surv_juvV_intersp = 1/(1+exp(-(alpha[3] + alpha[4] * N)))
surv_juvV_intersp_est = 1/(1+exp(-(alpha_est[3] + alpha_est[4] * N)))

# inter-species DD - predator fecundity fn of juvenile prey abundance
fecP_intersp = exp(alpha[5] + alpha[6] * N)
fecP_intersp_est = exp(alpha_est[5] + alpha_est[6] * N)

# intra-species DD - prey fecundity fn of adult prey abundance
fecV_intrasp =exp(alpha[7] + alpha[8] * N)
fecV_intrasp_est =exp(alpha_est[7] + alpha_est[8] * N)


pdf(file="Fig2_v1.pdf",width=8,height=8)

par(mfrow=c(2,2))

plot(N,surv_juvP_intrasp,type='l',lwd=3,col='blue',ylab='Juvenile P survival',ylim=c(0,1),xlab='P abundance',main='INTRA-DD')
lines(N,surv_juvP_intrasp_est,type='l',lwd=3,col='green')
legend('topright',col=c('blue','green'),legend=c('actual','estimated'),lty=1,lwd=3)
fig_label("A",cex=2)

plot(N,surv_juvV_intersp,type='l',lwd=3,col='blue',ylab='Juvenile V survival',ylim=c(0,1),xlab='P abundance',main='INTER-DD')
lines(N,surv_juvV_intersp_est,type='l',lwd=3,col='green')
legend('topright',col=c('blue','green'),legend=c('actual','estimated'),lty=1,lwd=3)
fig_label("B",cex=2)

plot(N,fecP_intersp,type='l',lwd=3,col='blue',ylab='P fecundity',ylim=c(0,10),xlab=' Juv V abundance',main='INTER-DD')
lines(N,fecP_intersp_est,type='l',lwd=3,col='green')
legend('topright',col=c('blue','green'),legend=c('actual','estimated'),lty=1,lwd=3)
fig_label("C",cex=2)

plot(N,fecV_intrasp,type='l',lwd=3,col='blue',ylab='V fecundity',ylim=c(0,10),xlab='Adult V abundance',main='INTRA-DD')
lines(N,fecV_intrasp_est,type='l',lwd=3,col='green')
legend('topright',col=c('blue','green'),legend=c('actual','estimated'),lty=1,lwd=3)
fig_label("D",cex=2)

dev.off()