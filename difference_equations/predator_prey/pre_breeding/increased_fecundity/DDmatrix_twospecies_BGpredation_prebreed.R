#####################################################################################################################################
### FBarraquand & OGimenez 14/11/2017 -- Two species density-dependent matrix model // predator -prey version
### 'New' version with density-dependence intra and inter affecting four vital rates, each vital rate affected once
### 18/03/2018 New set of parameters, slightly more realistic
### 07/09/2018 New version of pre-breeding model that we use statistically
### 09/09/2018 Shut off intraspecific interactions to demonstrate the blow-up. 
#####################################################################################################################################

library('Matrix')

####################### Definition of functions #####################################################

logistic=function(x){
return(1/(1+exp(-x)))
}


#logistic=function(x,x_thresh=1,steepness=1){
#return(1/(1+exp(-steepness*(x-x_thresh))))
#}

# Computes vital rates from the state vector N 
vital_rates=function(params,N){

############ Parameters used as inputs #################################################################
### similar to statistical model
#########################################################################################################

VR=matrix(NA,nrow=4,ncol=4)
phiA = params[1:2] # 
alphas = params[3:10] #

f=rep(NA,4)
### same order as in statistical model (not much logic)
f[1] = logistic(alpha[1] + alpha[2] *N[4]) # intra-species DD - juvenile predator survival fn of adult predator abundance
f[2] = logistic(alpha[3] + alpha[4] * N[4]) # inter-species DD - juvenile prey survival fn of adult predator abundance
f[3] = exp(alpha[5] + alpha[6] * N[1]) # inter-species DD - predator fecundity fn of juvenile prey abundance
f[4] = exp(alpha[7] + alpha[8] * N[2]) # intra-species DD - prey fecundity fn of adult prey abundance


VR_prey = matrix(c(0,0.5 * f[4]* f[2], 
            phiA[1],  phiA[1]),nrow = 2, ncol = 2,byrow = TRUE)

VR_pred =  matrix(c(0,0.5* f[3] * f[1], 
             phiA[2],phiA[2]),nrow = 2, ncol = 2,byrow = TRUE)

VR=bdiag(VR_prey,VR_pred)

return(VR)

}

############################################# end of functions def. ############################################################


################### Main program ############################################################################################### 

### Initialization
set.seed(42) 

# Define parameters
tmax = 500 #max time for sim

### Max rate parameters
# Prey (max) adult survival (v for victim)
phiA_V = 0.6
# Predator adult survival 
phiA_P = 0.7

# Other parameters 

alpha <- rep(NA,8)
N <- seq(10,500,length=10)
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
alpha[7] <- 5 #Here 2 shows easily the blowup, 1.5 takes up too much time
alpha[8] <- -0.005 # slopes
exp(alpha[7] + alpha[8] * N)


# Bundle the whole shabang into a single vector, which we then feed to the function
params = c(phiA_V,phiA_P,alpha)

# Initialize the state vector
N=matrix(NA,4,tmax)
N[,1] = N[,1] = c(100,100,20,20) #runif(4,50,100) ### we state 4 values for 2 species and 2 stages

############################# Loop over time #################################################################

for (t in 1:(tmax-1)){
	## Compute vital rates as a function of the initial state
	vital_rate_matrix = as.matrix(vital_rates(params,N[,t])) 
	N[,t+1] = vital_rate_matrix %*% N[,t] 
}
N[,tmax]
matplot(t(N[,1:tmax]))
matlines(t(N[,1:tmax]))
############################################ end of simulation ##################################################################
nplotstart = 100
nplotend = 200 #tmax
pdf(file="TimePlot_BGpredation_prebreed_highfec.pdf",width=12,height=8)
par(cex=1.2)
matplot(t(N[,nplotstart:nplotend]),xlab="Time",ylab="Densities")
matlines(t(N[,nplotstart:nplotend]),lwd=2)
dev.off()

pdf(file="Log_TimePlot_BGpredation_prebreed_highfec.pdf",width=12,height=8)
par(cex=1.2)
matplot(t(log(N[,nplotstart:nplotend])),xlab="Time",ylab="Densities")
matlines(t(log(N[,nplotstart:nplotend])),lwd=2)
dev.off()
