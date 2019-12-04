// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y_obs;
}

parameters {
  real logB0; //Initial biomass estimate
  real<lower=0> sigmaProc; //Process SD
  real<lower=0> sigmaObs; //Observation SD
  real mu_lambda; //Generating parameter for growth rate
  vector[N-1] lambda_t; //Actual population growth rate at time t
}

transformed parameters{
  vector[N] biomass_t; //Latent biomass variable
  biomass_t[1] = exp(logB0);
  
  for(i in 1:(N-1)){ //Generating process
    biomass_t[i+1] = lambda_t[i]*biomass_t[i]; //biomass(t+1) = lambda(t)*biomass(t)
  }
}

model {
  lambda_t ~ normal(mu_lambda,sigmaProc); //generating process for lambda_t
  y_obs ~ normal(biomass_t,sigmaObs); //Observation process
  
  //Priors
  logB0 ~ normal(0,20);
  sigmaProc ~ gamma(1,1);
  sigmaObs ~ gamma(1,1);
  mu_lambda ~ normal(0,20);
}

