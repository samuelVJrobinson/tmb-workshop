#Time-series / state-space models in TMB
#Cahill 27 Nov 2019 
rm(list=ls(all=TRUE))
#Simulate some fake data, analyze it, make lots of plots
require(TMB)
require(tidyverse)
theme_set(theme_classic())
require(TMBhelper)
require(mgcv)
Nyears = 50 #number of years of data
B_init = 20 #initial population size
mu_lambda = 1.05 #Mean annual population growth rate
SigP = 0.15 #Process (temporal) variation in growth rate lambda
SigO = 7 #Variance of the observation term 

y_obs = B = Nyears
B[1] = B_init #initialize  the N vector

#Simulate lambdas with process error
set.seed(42)
lambda_t = rnorm(Nyears-1, mu_lambda, SigP)

#Simulate the exponential growth model:
for(t in 1:(Nyears-1)){
  B[t+1] = B[t]*lambda_t[t]
}

#Generate the observed data | True Biomass
for(t in 1:Nyears){
  y_obs[t] = rnorm(1, B[t], SigO) 
}

#Fit a smoother across time (not the same model, but important alternative)
#could also do a loess, but the Galpern posse seems to have a thing for Simon's wiggle-lines
Gam = gam( y_obs ~ s(I(1:Nyears)) )
ypred_t = predict(Gam, se=TRUE)
gam.check(Gam)

# Compile model
setwd( "~/Documents/tmb-workshop/week 4" )
Version = "state_space_exponential"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs ) #Years 
Parameters = list( "logB0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears-1) ) #Parameter starting values
Random = c("lambda_t") #Random effect: actual population growth rate

Use_REML = TRUE #If REML, marginalize across lambda_t as well as logB0, mu_lambda
if( Use_REML==TRUE ) Random = union( Random, c("logB0","mu_lambda") )

# Build function
dyn.load( dynlib("state_space_exponential") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

#Run function through optimizer
Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj ) #Get standard errors 
#SEHat  = as.list( Opt$SD, "Std. Error" )

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

# Get reporting and SEs
Report = Obj$report() # Estimates of parameters (not mu_lambda?   I think this is where REPORT or ADREPORT are returning values inside .cpp function)
ParHat = as.list( Opt$SD, "Estimate" ) #Get ML estimates from Opt - I think this is getting the working parameter estimates (e.g. log_sigma0) rather than derived parameters like sigma0 (as above) 
ParHat[["biomass_t"]] = unname(SD$value[names(SD$value)=="biomass_t"]) #Extract biomass from SD$value

ParHat

#Plot stuff
plot.data = data.frame(true_biomass=B, true_lambda=c(lambda_t, NA), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = c(ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, NA),  
                       tmb_lambda_mu = c(ParHat$lambda_t, NA), 
                       tmb_lambda_upper = c(ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96, NA)
                       )
plot.data
png( file="results.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

plot(0, 0, ylim = c(min(plot.data$tmb_lower), max(plot.data$tmb_upper)), 
     xlim = c(0.5, Nyears), ylab = "Population
     Biomass (tons)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lower, plot.data$tmb_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
points(plot.data$y_obs, type="p", pch=16, col="black")
lines(plot.data$true_biomass, type="l", col="darkorange")# truth
lines(plot.data$tmb_mu, type="l", col="blue") #estimated
lines(plot.data$gam_mu, type="l", col="red") #estimated gam

legend(x = 30, y = 45, legend = c("True", "Observed", "Estimated (TMB)", 
                                                       "Estimated (Gam)"),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black", "blue", "red"),
       bty = "n", cex = 1)
dev.off()

#Here's a thing gam cannot give you (to my knowledge):

#Plot stuff
plot.data = data.frame(true_lambda=lambda_t,
                       #Get lambda_t from ParHat, SD for lambda_t from SD object 
                       tmb_lambda_lower = ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, 
                       tmb_lambda_mu = ParHat$lambda_t, 
                       tmb_lambda_upper = ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96
)

png( file="lambda.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

plot(0, 0, ylim = c(0.5, 1.5), 
     xlim = c(0.5, Nyears-1), ylab = "Lambda", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:nrow(plot.data), (nrow(plot.data):1)), 
        y = c(plot.data$tmb_lambda_lower, plot.data$tmb_lambda_upper[(Nyears-1):1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
lines(ParHat$lambda_t)
lines(plot.data$true_lambda, pch=16, type="b", col="darkorange")
abline(h=1.0, col="blue", lty=2)
legend('bottomright', 
       legend = c(expression(paste('True ',lambda)), expression(paste('Estimated ',lambda))),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black"),
       bty = "n", cex = 1)

dev.off()

#---------------------------------------------------------
#Do it again but predict left out years this time
dyn.unload( dynlib("state_space_exponential") )

Data = list( "Nyears"=Nyears, "Y_obs_t"=y_obs )
Data$Y_obs_t[47:50]=NA
Parameters = list( "logB0"=0, "log_sigmaP"=1, "log_sigmaO"=1, "mu_lambda"=1, "lambda_t"=rep(0,Nyears-1) )

compile( paste0(Version,".cpp") )
dyn.load( dynlib("state_space_exponential") )

Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

Opt = fit_tmb( obj=Obj, newtonsteps=1 )
Opt

SD = sdreport( Obj )
SD

#Did it converge
final_gradient = Obj$gr( Opt$par )
if( any(abs(final_gradient)>0.0001) | SD$pdHess==FALSE ) stop("Not converged")

ParHat = as.list( Opt$SD, "Estimate" )
ParHat$lambda_t
ParHat[["biomass_t"]] = SD$value[names(SD$value)=="biomass_t"]
Report = Obj$report()

plot.data = data.frame(true_biomass=B, true_lambda=c(lambda_t, NA), 
                       y_obs = y_obs, gam_lower = ypred_t$fit - ypred_t$se.fit*1.96,gam_mu = ypred_t$fit,  
                       gam_upper = ypred_t$se.fit*1.96+ypred_t$fit, 
                       tmb_lower = ParHat$biomass_t - SD$sd[names(SD$value)=="biomass_t"]*1.96 , tmb_mu = ParHat$biomass_t,  
                       tmb_upper = SD$sd[names(SD$value)=="biomass_t"]*1.96 + ParHat$biomass_t, 
                       tmb_lambda_lower = c(ParHat$lambda_t - SD$sd[names(SD$value)=="lambda_t"]*1.96, NA),  
                       tmb_lambda_mu = c(ParHat$lambda_t, NA), 
                       tmb_lambda_upper = c(ParHat$lambda_t + SD$sd[names(SD$value)=="lambda_t"]*1.96, NA), 
                       year=1:Nyears
)

png( file="predict.png", width=8, height=5, res=800, units="in" )
par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )

plot(0, 0, ylim = c(min(plot.data$tmb_lower), max(plot.data$tmb_upper)), 
     xlim = c(0.5, Nyears), ylab = "Population
     Biomass (tons)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2,
     axes = TRUE)
polygon(x = c(1:Nyears, Nyears:1), y = c(plot.data$tmb_lower, plot.data$tmb_upper[Nyears:1]),
        col = "gray90", border = "gray90") #confidence intervals tmb
points(plot.data$y_obs, type="p", pch=16, col="black")
lines(plot.data$true_biomass, type="l", col="darkorange")# truth
lines(plot.data$tmb_mu, type="l", col="blue") #estimated
points(plot.data$y_obs[47:50]~plot.data$year[47:50], col="red", pch=16)

legend(x = 30, y = 40, legend = c("True", "Observed", "Estimated (TMB)", "Predicted"),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("darkorange","black", "blue", "red"),
       bty = "n", cex = 1)

dev.off()
#---------------------------------------------------------
#Redo it with 
#SigP = SigO = 1; B0 = 20;
#for a 'fun' learning experience
#i.e., see Auger-methe et al. 2016

# Same thing, but in Stan---------------------------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Set up data and initial values in chain
dat <- list(N=Nyears,y_obs=y_obs)
startVals <- function() {
  list(logB0 = 3, sigmaProc = 0.15, sigmaObs = 7, mu_lambda = 1, 
       lambda_t = rep(1,dat$N-1))
}

#true values: logB0 = exp(20) = 2.99, mu_lambda = 1.05, sigmaProc = 0.15, sigmaObs = 7

stanMod1 <- stan(file='state_space_exponential.stan',data=dat,chains=4,init=startVals,
                 iter=5000)

#Parameters
p <- c('logB0','sigmaProc','sigmaObs','mu_lambda')

stan_trace(stanMod1,pars=p) #Looks OK, but lower than estimate

mod1Res <- extract(stanMod1)

#Doing a bad job at estimating sigmaObs/sigmaProc. Bad priors?

data.frame(draw=unlist(mod1Res))

as.data.frame(t(sapply(mod1Res[1:4],function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rownames_to_column(var='var') %>% 
  mutate(actual=c(log(B_init),SigP,SigO,mu_lambda)) %>% 
  rename(med=2,lwr=3,upr=4) %>% 
  ggplot(aes(x=1))+facet_wrap(~var,scales='free_y')+
  geom_pointrange(aes(y=med,ymax=upr,ymin=lwr),col='red')+
  geom_point(aes(y=actual))
  






