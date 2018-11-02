############################################################################################
# Two-species integrated population model
# Combining counts, repro and capture-recapture data
# Data simulation and model fitting in Jags
# O. Gimenez and F. Barraquand 
# Notes:
#		We relied on code from Kery and Schaub book (chapter 11) initially written
# 		for analyses that were carried out in Schaub, Gimenez, Sierro, Arlettaz (2007)
#		Use of integrated modeling to enhance estimates of population dynamics
# 		obtained from limited data. Conservation Biology 21: 945–955. 
# 		We were also inspired by some code from Abadi, Gimenez, Jakobere, Stauberf, 
#		Arlettaz, Schaub (2012). Estimating the strength of density dependence in the 
# 		presence of observation errors using integrated population models. Ecological 
#		Modelling 242: 1-9.
############################################################################################

# COMMENTS:
# O. Gimenez 20/11/2017: draft code for simulations and model fitting
# O. Gimenez 15/03/2018: change the code to have individuals marked in each year 
# F. Barraquand 25/07/2018: change to formulation with predation on prey juveniles & non-zero top-down effect
# O. Gimenez 26/07/2018: add productivity likelihood component
# F. Barraquand 31/07/2018: add informative priors and investigated scaling (/100) on all DD relationships
# O. Gimenez 02/08/2018: removed scaling and tweaked the prior on alpha a little
# O. Gimenez 03/08/2018: got rid of the measurement error in the count model by using very small 
# 						 values for the variance of the obs error, got rid of informative priors, 
# 						 and cleaned up comments
# F. Barraquand 13/08/2018: added nonlinear (1-exp(-x)) predator productivity

# MODEL:
# Consider 2 species with age-structured demography 
# 2 age-classes, Besbeas et al 2003 formulation

# with matrix pop model for predator: [                0, phi^J_t(predator) * f^A(prey),
#                                      phi^J_t(predator),             phi^A_t(predator)]

# and matrix pop model for prey: [                        0, phi^J_t(predator) * f^A(prey),
#                                             phi^A_t(prey),             	 phi^A_t(prey)]

# DENSITY-DEPENDENCE (DD) RELATIONSHIPS (time is implicit in notation)

# Inter-DD
# logit(phi^J_V) = beta_0 + beta_1 * N^A_P 	# effect of adult predator on juvenile prey survival
# logit(f_P) = eta_0 + eta_1 * N^J_V   		# effect of juvenile prey on the fecundity of the predator

# Intra-DD
# log(f_V) = eta_0 + eta_1 * N^A_V 			# effect of prey on prey fecundity
# logit(phi^J_P) = beta_0 + beta_1 * N^A_P 	# effect of adult predator on juvenile predator survival

############################################################################################
######################################## SIMULATION ########################################
############################################################################################

# simulations are performed in Jags 
# (yes, that's possible, check out https://oliviergimenez.github.io/post/sim_with_jags/)

# load packages
library(R2jags)
library(mcmcplots)

# code to simulate data
sink("icm.jags")
cat("
data{
# predator = P, prey = V
#-- Priors and constraints
# Survival and recapture probabilities, as well as detection
for (i in 1:nindP){
   for (t in fP[i]:(n.occasions-1)){
      phiP[i,t] <- betaP[xP[i,t],t]
      pP[i,t] <- mean.pP # constant predator detection
      } #t
   } #i
for (t in 1:(n.occasions-1)){
	logit(betaP[1,t]) <- alpha[1] + alpha[2] * NadP[t] # intra-species DD - juvenile predator survival fn of adult predator abundance
	betaP[2,t] <- mean.sadP # constant adult predator survival 
} #t

for (i in 1:nindV){
   for (t in fV[i]:(n.occasions-1)){
      phiV[i,t] <- betaV[xV[i,t],t]
      pV[i,t] <- mean.pV  # constant prey detection
      } #t
   } #i
for (t in 1:(n.occasions-1)){
	logit(betaV[1,t]) <- alpha[3] + alpha[4] * NadP[t] # inter-species DD - juvenile prey survival fn of adult predator abundance 
	betaV[2,t] <- mean.sadV  # constant adult prey survival -- previously mean.sjuvV was constant
} #t
for (t in 1:(n.occasions-1)){
   sjuvP[t] <- betaP[1,t] # juvenile predator survival
   sadP[t] <- betaP[2,t] # adult predator survival
   fecP[t] <-exp(alpha[5])*N1V[t]/(exp(alpha[6])+N1V[t]) # inter-species DD - predator fecundity fn of juvenile prey abundance
   sjuvV[t] <- betaV[1,t] # juvenile prey survival
   sadV[t] <- betaV[2,t] # adult prey survival
   log(fecV[t]) <- alpha[7] + alpha[8] * NadV[t] # intra-species DD - prey fecundity fn of adult prey abundance
   }

#-- Likelihoods
# 1. Likelihood for population population count data (state-space model)
# System process
N1P[1] <- N1initP # init juv predator abundance
NadP[1] <- NadinitP # init adult predator abundance
N1V[1] <- N1initV # init juv prey abundance
NadV[1] <- NadinitV # init adult prey abundance
for(t in 2: n.occasions){
	mean1P[t] <- fecP[t-1] / 2 * sjuvP[t-1] * NadP[t-1]
	N1P[t] ~ dpois(mean1P[t])
	NadP[t] ~ dbin(sadP[t-1], NtotP[t-1])
	mean1V[t] <- fecV[t-1] / 2 * sjuvV[t-1] * NadV[t-1]
	N1V[t] ~ dpois(mean1V[t])
	NadV[t] ~ dbin(sadV[t-1], NtotV[t-1])
}
for(t in 1:n.occasions){
	NtotP[t] <- NadP[t] + N1P[t] # total predator abundance
	NtotV[t] <- NadV[t] + N1V[t] # total prey abundance
}
# Observation process
for(t in 1:n.occasions){
	countsP[t] ~ dnorm(NtotP[t], tauyP) 
	countsV[t] ~ dnorm(NtotV[t], tauyV)
}
# 2. Likelihood for capture-recapture data (state-space model)
# predator
for (i in 1:nindP){
   # Define latent state at first capture
   zP[i,fP[i]] <- 1
   mu2P[i,fP[i]] <- 1 * zP[i,fP[i]] # detection is 1 at first capture (conditional on first capture)
   chP[i,fP[i]] ~ dbern(mu2P[i,fP[i]])
   for (t in (fP[i]+1):n.occasions){
      # State process
      zP[i,t] ~ dbern(mu1P[i,t])
      mu1P[i,t] <- phiP[i,t-1] * zP[i,t-1]
      # Observation process
      chP[i,t] ~ dbern(mu2P[i,t])
      mu2P[i,t] <- pP[i,t-1] * zP[i,t]
      } #t
   } #i
# prey
for (i in 1:nindV){
   # Define latent state at first capture
   zV[i,fV[i]] <- 1
   mu2V[i,fV[i]] <- 1 * zV[i,fV[i]] # detection is 1 at first capture (conditional on first capture)
   chV[i,fV[i]] ~ dbern(mu2V[i,fV[i]])
   for (t in (fV[i]+1):n.occasions){
      # State process
      zV[i,t] ~ dbern(mu1V[i,t])
      mu1V[i,t] <- phiV[i,t-1] * zV[i,t-1]
      # Observation process
      chV[i,t] ~ dbern(mu2V[i,t])
      mu2V[i,t] <- pV[i,t-1] * zV[i,t]
      } #t
   } #i
# 3. Likelihood for productivity data (Poisson regression)
# total number of nestlings counted in year t (J) is a Poisson
# with rate = product of nb of surveyed broods (R) and productivity (f)
for(t in 1:(n.occasions-1)){
	JP[t] ~ dpois(rhoP[t]) # predator prod
	JV[t] ~ dpois(rhoV[t]) # prey prod
	rhoP[t] <- RP[t] * fecP[t]
	rhoV[t] <- RV[t] * fecV[t]
}
}
model{
fake <- 0
}
",fill = TRUE)
sink()

# specify parameters for simulations 

# number of survey occasions (same for predator/prey, could be changed though)
n.occasions <- 10

# parameters for predator; for simplicity, all individuals are marked as young 
# and released at each sampling occasion (but the last one)
nindj <- 100 # nb of individuals marked as young per year
fj <- as.numeric(gl(n.occasions-1,nindj))
fP <- fj
nindP <- length(fP)

# create matrices X indicating age classes
xj <- matrix(NA, ncol = n.occasions-1, nrow = nindP)
for (i in 1: nindP){
   for (t in fj[i]:(n.occasions-1)){
   	  #if (fj[i]>(n.occasions-1)) next
      xj[i,t] <- 2
      xj[i,fj[i]] <- 1   
      } #t
   } #i
xP <- xj

# demog parameters
mean.pP <- 0.7 # detection is the same for youngs and adults
mean.sadP <- 0.7 # constant adult survival
sigma.y <- 0.001 # basically no observation error
tauyP <- 1/(sigma.y*sigma.y)

# init abundance
N1initP <- 20
NadinitP <- 20

# annual number of surveyed broods 
RP <- rep(20,n.occasions-1)

# parameters for prey; for simplicity, all individuals are marked as young 
# and released at each sampling occasion (but the last one)
nindj <- 100 # nb of individuals marked as young per year
fj <- as.numeric(gl(n.occasions-1,nindj))
fV <- fj
nindV <- length(fV)

# create matrices X indicating age classes
xj <- matrix(NA, ncol = n.occasions-1, nrow = nindV)
for (i in 1: nindP){
   for (t in fj[i]:(n.occasions-1)){
   	  #if (fj[i]>(n.occasions-1)) next
      xj[i,t] <- 2
      xj[i,fj[i]] <- 1   
      } #t
   } #i
xV <- xj

# demo parameters
mean.pV <- 0.7 # detection is the same for youngs and adults
mean.sadV <- 0.6 # constant adult survival ##previously juvenile
sigma.y <- 0.001 # observation error
tauyV <- 1/(sigma.y*sigma.y)

# init abundance
N1initV <- 100
NadinitV <- 100

# annual number of surveyed broods 
RV <- rep(50,n.occasions-1)

# Specify parameters for DD relationships
# Plot relationships to check they're biologically relevant

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
alpha[5] <- 1.5 #log(max fecundity)
alpha[6] <- 4 # slopes
exp(alpha[5])*N/(exp(alpha[6])+N) 

# intra-species DD - prey fecundity fn of adult prey abundance
alpha[7] <- 2 #1.5 #2.5 too much
alpha[8] <- -0.005 # slopes
exp(alpha[7] + alpha[8] * N)

# build list of data to be passed in jags
jags.data <- list(n.occasions = n.occasions, alpha = alpha, tauyP = tauyP, N1initP = N1initP, NadinitP = NadinitP, xP = xP, mean.pP = mean.pP, fP = fP, nindP = nindP, tauyV = tauyV, N1initV = N1initV, NadinitV = NadinitV, xV = xV, mean.pV = mean.pV, fV = fV, nindV = nindV, mean.sadV = mean.sadV, mean.sadP = mean.sadP, RP = RP, RV = RV)

# simulate data
jmodel <- jags.model("icm.jags", jags.data, inits=NULL,n.chains = 1,n.adapt = 1)
jsample <- coda.samples(jmodel, c("NadP","N1P","NadV","N1V","chP","chV","countsP","countsV","JP","JV"), n.iter=1, thin = 1)

# format data to be used in the model fitting step
Simulated <- coda::as.mcmc(jsample)
Simulated
dim(Simulated)
colnames(Simulated)

# get number of produced offspring
mask <- grep("JP",colnames(Simulated))
JP <- as.vector(Simulated[mask])
JP
mask <- grep("JV",colnames(Simulated))
JV <- as.vector(Simulated[mask])
JV

# get number of counted animals
mask <- grep("countsP",colnames(Simulated))
countsP <- as.vector(Simulated[mask])
countsP
mask <- grep("countsV",colnames(Simulated))
countsV <- as.vector(Simulated[mask])
countsV

# plot latent abundance and counts; should be indistinguishable 
# because we specified almost negligeable obs error on counts
plot(1:n.occasions, countsP, type='l', lwd=3, ylim=c(0,max(countsV,countsP)), col='red', ylab='counts', xlab='years')
lines(1:n.occasions, countsV, type='l', lwd=3, col='blue')
legend('topright', col=c('red','blue'), legend=c('predator','prey'), lty=1, lwd=3)

# get encounter histories for predators
mask <- grep("chP",colnames(Simulated)) # get CR histories for predators
names_colP <- colnames(Simulated)[mask] # get names of columns
names_colP_col <- sub('.*,\\s*','', names_colP) # extract column index (1)
names_colP_col <- gsub("\\]", "", names_colP_col) # extract column index (2)
names_colP_col <- as.numeric(names_colP_col) # extract column index (3)
names_colP_col
names_colP_row <- sub('\\s*,.*','', names_colP) # extract row index (1)
names_colP_row <- gsub("chP", "", names_colP_row) # extract row index (2)
names_colP_row <- gsub("\\[", "", names_colP_row) # extract row index (3)
names_colP_row <- as.numeric(names_colP_row) # extract row index (4)
names_colP_row
chP <- matrix(0,nrow= nindP,ncol=n.occasions)
index <- 1
for (i in 1:length(names_colP_col)){
	chP[names_colP_row[i], names_colP_col[i]] <- Simulated[mask][i]
	index <- index + 1
}
head(chP)

# get encounter histories for preys
mask <- grep("chV",colnames(Simulated))
names_colV <- colnames(Simulated)[mask] # get names of columns
names_colV_col <- sub('.*,\\s*','', names_colV) # extract column index (1)
names_colV_col <- gsub("\\]", "", names_colV_col) # extract column index (2)
names_colV_col <- as.numeric(names_colV_col) # extract column index (3)
names_colV_col
names_colV_row <- sub('\\s*,.*','', names_colV) # extract row index (1)
names_colV_row <- gsub("chV", "", names_colV_row) # extract row index (2)
names_colV_row <- gsub("\\[", "", names_colV_row) # extract row index (3)
names_colV_row <- as.numeric(names_colV_row) # extract row index (4)
names_colV_row
chV <- matrix(0,nrow= nindV,ncol=n.occasions)
index <- 1
for (i in 1:length(names_colV_col)){
	chV[names_colV_row[i], names_colV_col[i]] <- Simulated[mask][i]
	index <- index + 1
}
head(chV)

############################################################################################
######################################## ESTIMATION ########################################
############################################################################################

# specify model
sink("ipm-prod_2species.jags")
cat("
model {
#-------------------------------------------------
#  2-species integrated population model (Besbeas et al. formulation)
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant, except when intra/inter DD
#  - N is for prey, P is for predator
# INTER DENSITY DEPENDENCE
# logit(phi^A_V) = beta_0 + beta_1 * N^A_P
# log(f_P) = eta_0 + eta_1 * N^J_V
# INTRA DENSITY DEPENDENCE
# log(f_V) = eta_0 + eta_1 * N^A_V
# logit(phi^J_P) = beta_0 + beta_1 * N^A_P
# Note: We use huge precision = tiny variance 
# 		and neglect obs error
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------

# Initial population sizes
n1P ~ dnorm(25, 100000)T(0,)     # 1-year predators
nadP ~ dnorm(25, 100000)T(0,)    # adult predators
N1P[1] <- round(n1P)
NadP[1] <- round(nadP)
n1V ~ dnorm(25, 100000)T(0,)     # 1-year preys
nadV ~ dnorm(25, 100000)T(0,)    # adult preys
N1V[1] <- round(n1V)
NadV[1] <- round(nadV)

# Survival and recapture probabilities, as well as productivity
for (i in 1:nindP){
   for (t in fP[i]:(nyears-1)){
      phiP[i,t] <- betaP[xP[i,t],t]
      pP[i,t] <- mean.pP # constant predator detection
      } #t
   } #i
for (t in 1:(nyears-1)){
	logit(betaP[1,t]) <- alpha[1] + alpha[2] * NadP[t] # intra-species DD - juvenile predator survival fn of adult predator abundance	
	betaP[2,t] <- mean.sadP # constant adult predator survival
} #t

for (i in 1:nindV){
   for (t in fV[i]:(nyears-1)){
      phiV[i,t] <- betaV[xV[i,t],t]
      pV[i,t] <- mean.pV # constant prey detection
      } #t
   } #i
for (t in 1:(nyears-1)){
	logit(betaV[1,t]) <- alpha[3] + alpha[4] * NadP[t] # fn of adult predator abundance, previously mean.sjuvV, constant juvenile prey survival
  betaV[2,t] <- mean.sadV #constant adult prey survival 
} #t

mean.sadV ~ dunif(0, 1) # prior for juvenile prey survival -- info prior
mean.sadP ~ dunif(0, 1) # prior for adult predator survival
mean.pP ~ dunif(0,1) # prior for mean recapture 
mean.pV ~ dunif(0,1) # prior for mean recapture
for (t in 1:(nyears-1)){
   sjuvP[t] <- betaP[1,t] # juvenile predator survival
   sadP[t] <- betaP[2,t] # adult predator survival
   fecP[t] <- exp(alpha[5])*N1V[t]/(exp(alpha[6])+N1V[t]) # inter-species DD - predator fecundity fn of juvenile prey abundance
   sjuvV[t] <- betaV[1,t] # juvenile prey survival
   sadV[t] <- betaV[2,t] # adult prey survival
   log(fecV[t]) <- alpha[7] + alpha[8] * NadV[t] # intra-species DD - prey fecundity fn of adult prey abundance
   }

for (j in 1:4){
		alpha[2*j-1] ~ dnorm(0,1) # priors on density-dependence parameters
  # also tried sd changed from 1 to 2 because of alpha[5] that is about 5
  alpha[2*j] ~ dnorm(0,1) # same scale here
}
#alpha[5] ~ dlnorm(0,1) # strictly positive -- does not improve things

#-------------------------------------------------
# 2. The likelihoods of the isolated data sets
#-------------------------------------------------
# 2.1. Likelihood for population population count data (state-space model)
# 2.1.1 System process
for (t in 2:nyears){
	mean1P[t] <- fecP[t-1] / 2 * sjuvP[t-1] * NadP[t-1]
	N1P[t] ~ dpois(mean1P[t])
	NadP[t] ~ dbin(sadP[t-1], NtotP[t-1])
	mean1V[t] <- fecV[t-1] / 2 * sjuvV[t-1] * NadV[t-1]
	N1V[t] ~ dpois(mean1V[t])
	NadV[t] ~ dbin(sadV[t-1], NtotV[t-1])
}
for (t in 1:nyears){
	NtotP[t] <- NadP[t] + N1P[t] # total predator abundance
	NtotV[t] <- NadV[t] + N1V[t] # total prey abundance
}
# 2.1.2 Observation process
for (t in 1:nyears){
countsP[t] ~ dnorm(NtotP[t], 100000)
countsV[t] ~ dnorm(NtotV[t], 100000)
}

# 2.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# State-space modelling formulation

# predator
for (i in 1:nindP){
	# Define latent state at first capture
	zP[i,fP[i]] <- 1
	for (t in (fP[i]+1): nyears){
		# State process
		zP[i,t] ~ dbern(mu1P[i,t])
		mu1P[i,t] <- phiP[i,t-1] * zP[i,t-1]
		# Observation process
		chP[i,t] ~ dbern(mu2P[i,t])
		mu2P[i,t] <- pP[i,t-1] * zP[i,t]
	} #t
} #i

# prey
for (i in 1:nindV){
	# Define latent state at first capture
	zV[i,fV[i]] <- 1
	for (t in (fV[i]+1): nyears){
		# State process
		zV[i,t] ~ dbern(mu1V[i,t])
		mu1V[i,t] <- phiV[i,t-1] * zV[i,t-1]
		# Observation process
		chV[i,t] ~ dbern(mu2V[i,t])
		mu2V[i,t] <- pV[i,t-1] * zV[i,t]
	} #t
} #i

# 2.3 Likelihood for productivity data: Poisson regression
# total number of nestlings counted in year t (J) is a Poisson
# with rate = product of nb of surveyed broods (R) and productivity (f)
for(t in 1:(nyears-1)){
	JP[t] ~ dpois(rhoP[t]) # predator prod
	JV[t] ~ dpois(rhoV[t]) # prey prod
	rhoP[t] <- RP[t] * fecP[t]
	rhoV[t] <- RV[t] * fecV[t]
}

}
",fill = TRUE)
sink()

# build list of data to be passed in jags (use info for capture-recapture models, see Kery and Schaub book)
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }
jags.data <- list(nyears = ncol(chP), countsP = countsP, chP = chP, nindP = nrow(chP), fP = fP, xP = xP, zP = known.state.cjs(chP), countsV = countsV, chV = chV, nindV = nrow(chV), fV = fV, xV = xV, zV = known.state.cjs(chV), JP = JP,  JV = JV, RP = RP, RV = RV)

# specify initial values: we assigned the values we used to simulated data as initial values, 
# other parameters (latent states in particular) are assigned initial values by Jags 
inits <- function(){list(
mean.pP = mean.pP, n1P = N1initP, nadP = NadinitP, mean.sadP = mean.sadP, 
mean.sadV = mean.sadV, mean.pV = mean.pV, n1V = N1initV, nadV = NadinitV,
alpha = alpha)}  

# specify which parameters we'd like to monitor
parameters <- c("mean.pP", "N1P", "NadP", "NtotP","alpha","mean.pV", "N1V", "NadV", "NtotV")

# MCMC settings
ni <- 20000
nb <- 10000
nc <- 2

# run jags
ipm.prod <- jags(jags.data, inits, parameters, "ipm-prod_2species.jags", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())

# summarize posteriors
print(ipm.prod, digits = 2)

# traceplots
traplot(ipm.prod,c("mean.pP","alpha", "mean.pV"))

# posterior densities
denplot(ipm.prod,c("mean.pP","alpha", "mean.pV"))

# correlation between posterior densities 
parcorplot(ipm.prod,parms = "alpha")

# compare actual vs estimated abundance

par(mfrow=c(1,2))
# real vs estimated predator abundance
plot(1:n.occasions,countsP,type='l',lwd=3,col='blue',ylim=c(0,max(countsP)),ylab='counts',xlab='years',main='predator')
lines(1:n.occasions, ipm.prod$BUGSoutput$mean$NtotP,type='l',lwd=3,col='green')
legend('bottomright',col=c('blue','green'),legend=c('truth','estimated'),lty=1,lwd=3)

# real vs estimated prey abundance
plot(1:n.occasions,countsV,type='l',lwd=3,col='blue',ylim=c(0,max(countsV)),ylab='counts',xlab='years',main='prey')
lines(1:n.occasions, ipm.prod$BUGSoutput$mean$NtotV,type='l',lwd=3,col='green')
legend('bottomright',col=c('blue','green'),legend=c('truth','estimated'),lty=1,lwd=3)

# real vs estimated parameters involved in DD relationships 
alpha_est <- ipm.prod$BUGSoutput$mean$alpha
data.frame(real=alpha,estimated=alpha_est )

# plot DD relationships (actual and estimated)
N <- seq(10,1000,length=n.occasions) #density index

# intra-species DD - juvenile predator survival fn of adult predator abundance
surv_juvP_intrasp = 1/(1+exp(-(alpha[1] + alpha[2] * N)))
surv_juvP_intrasp_est = 1/(1+exp(-(alpha_est[1] + alpha_est[2] * N)))

# inter-species DD - juvenile prey survival fn of adult predator abundance
surv_juvV_intersp = 1/(1+exp(-(alpha[3] + alpha[4] * N)))
surv_juvV_intersp_est = 1/(1+exp(-(alpha_est[3] + alpha_est[4] * N)))

# inter-species DD - predator fecundity fn of juvenile prey abundance
fecP_intersp = exp(alpha[5])*N/(exp(alpha[6])+N)
fecP_intersp_est = exp(alpha_est[5])*N/(exp(alpha_est[6])+N) 

# intra-species DD - prey fecundity fn of adult prey abundance
fecV_intrasp =exp(alpha[7] + alpha[8] * N)
fecV_intrasp_est =exp(alpha_est[7] + alpha_est[8] * N)

pdf(file="DD_saturating_michaelis.pdf",width=8,height=8)

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
