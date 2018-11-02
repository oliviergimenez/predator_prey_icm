#####################################################################################################################################
### FBarraquand & OGimenez 14/11/2017 -- Two species density-dependent matrix model // predator -prey version
### 'New' version with density-dependence intra and inter affecting four vital rates, each vital rate affected once
### 18/03/2018 New set of parameters, slightly more realistic
#####################################################################################################################################

library('Matrix')

####################### Definition of functions #####################################################

logistic=function(x,x_thresh=1,steepness=1){
return(1/(1+exp(-steepness*(x-x_thresh))))
}

# Computes vital rates from the state vector N 
vital_rates=function(params,N){

############ Parameters used as inputs #################################################################
### we have params = c(phiJ_V,phiA_V,phiJ_P,phiA_P,gamma_V,gamma_P,b_P,b_V,beta,eta,V_thresh,P_thresh)
### The vector of state variables is N = c(N^J_V,N^A_V,N^J_P,N^A_P)
#########################################################################################################

###################################### Projection matrix ####################################################
# we have a projection matrix of the form 
# for prey: [(1-gamma_V)*phiJ_V, b_V, 
#             gamma_V*phiJ_V,   phiA_V]
### We could eventually add a + bJ_V juvenile fertility to the 1st line first column for e.g. passerine birds
# for predator: [(1-gamma_P)*phiJ_P, b_P, 
#             gamma_P*phiJ_P,   phiA_P]
##################################################################################################################

########### Density-dependence ##############################################################################
# Here we use 4 different vital rates upon which the densities impose some density-dependence
# Each of these 4 VRs is affected by different element of the state vector N   
# Now beta[i] models threshold value for the density-dependent term i
# and eta[i] controls the steepness of the logistic function employed
# Thus we lose 2 parameters compared with the previous version
#############################################################################################################

VR=matrix(NA,nrow=4,ncol=4)
beta = params[9:12] # beta = c(15,10,15,5) # predation kicks in at higher densities than prey DD
eta = params[13:16] # eta = c(-1,-1,1,-1)
### Sidenote: it looks like the logistic allows two implicitly connect better than other functions
### the prey mortality and predator reproduction process because we can specify one as the inverse
### of the other, and with similar threshold of juv prey density for both processes

f=rep(NA,4)
###for (i in 1:4){f[i] = logistic(N[i],beta[i],eta[i])} ### was completely wrong
f[1] = logistic(N[2],beta[1],eta[1])
f[2] = logistic(N[4],beta[2],eta[2])
f[3] = logistic(N[4],beta[3],eta[3])
f[4] = logistic(N[1],beta[4],eta[4])


VR_prey = matrix(c((1-params[5])*params[1]*f[2], params[7]*f[1], 
            params[5]*params[1]*f[2],  params[2]),nrow = 2, ncol = 2, byrow = TRUE)

VR_pred =  matrix(c((1-params[6])*params[3]*f[3], params[8]*f[4], 
             params[6]*params[3]*f[3], params[4]),nrow = 2, ncol = 2, byrow = TRUE)

VR=bdiag(VR_prey,VR_pred)

return(VR)

### Eventually, the O blocks of this matrix could be used to model interactions later on 
### (i.e., matter flow leaving the prey for the predator) with appropriate density-dependencies
}

############################################# end of functions def. ############################################################


################### Main program ############################################################################################### 

### Initialization
set.seed(42) 

# Define parameters
tmax = 100 #max time for sim

### Max rate parameters
# Prey (max) survival (v for victim)
phiJ_V = 0.4
phiA_V = 0.6
# Predator survival 
phiJ_P = 0.7
phiA_P = 0.9
# Transition rate from juvenile to adulthood
gamma_V = 1
gamma_P = 0.5 #from 1/4 to 1
# Reproduction (b or pi or f) max rate // only adults reproduce
b_P = 3
b_V = 7 ### previously 5

###################### Density-dependence parameters #####################################
### Now beta is the vector of threshold values and eta the vector of steepness
#beta = c(15,10,15,5) # predation kicks in at higher densities than prey DD
beta = c(15,30,5,10)
eta = c(-1,-1,-1,1)
### Sidenote: it looks like the logistic allows two implicitly connect better than other functions
### the prey mortality and predator reproduction process because we can specify one as the inverse
### of the other, and with similar threshold of juv prey density for both processes

### --- would be easier to compute as fraction of equilibrium densities... 

# Bundle the whole shabang into a single vector, which we then feed to the function
params = c(phiJ_V,phiA_V,phiJ_P,phiA_P,gamma_V,gamma_P,b_P,b_V,beta,eta)

# Initialize the state vector
N=matrix(NA,4,tmax)
N[,1] = runif(4,50,100) ### we state 4 values for 2 species and 2 stages

############################# Loop over time #################################################################

for (t in 1:(tmax-1)){
	## Compute vital rates as a function of the initial state
	vital_rate_matrix = as.matrix(vital_rates(params,N[,t])) 
	N[,t+1] = vital_rate_matrix %*% N[,t] 
}
N[,tmax]
matplot(t(N[,20:tmax]))
matlines(t(N[,20:tmax]))
############################################ end of simulation ##################################################################

pdf(file="TimePlot_BGpredation_logistic.pdf",width=12,height=8)
par(cex=1.2)
matplot(t(N[,20:tmax]),xlab="Time",ylab="Densities")
matlines(t(N[,20:tmax]),lwd=2)
dev.off()

pdf(file="TimePlot_BGpredation_logistic_log.pdf",width=12,height=8)
par(cex=1.2)
matplot(t(log(N[,20:tmax])),xlab="Time",ylab="Densities")
matlines(t(log(N[,20:tmax])),lwd=2)
dev.off()

### Note: I can make a program to see where we are in parameter space with the density-independent model
# We want to have the largest eigenvalue close to 1 in terms of modulus
### or can't we look at it this way? 

# Computes vital rates from the state vector N 
max_modulusEigen_DImodel=function(params, factorCorrect = 1 ){
  
  VR=matrix(NA,nrow=4,ncol=4)
  # factorCorrect // May be used to lower the vital rates based on the fact that their are lower in the DD model
  VR_prey = matrix(c((1-params[5])*params[1]*factorCorrect, params[7], 
                     params[5]*params[1]*factorCorrect,  params[2]),nrow = 2, ncol = 2, byrow = TRUE)
  
  VR_pred =  matrix(c((1-params[6])*params[3], params[8]*factorCorrect, 
                      params[6]*params[3], params[4]),nrow = 2, ncol = 2, byrow = TRUE)
  
  VR=bdiag(VR_prey,VR_pred)
  max_eig = max(abs(eigen(VR)$values))
  return(max_eig)
}

max_modulusEigen_DImodel(params)
max_modulusEigen_DImodel(params,factorCorrect = 0.5)
max_modulusEigen_DImodel(params,factorCorrect = 0.1)

## This parameter set seems to highlight interesting 4-point cycles for the prey 
## coexisting with longer cycles for the predator, which is super interesting. 

