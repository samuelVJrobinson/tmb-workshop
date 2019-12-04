#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt ); //Number of observations
  DATA_VECTOR( y_t ); //Observations
  
  // Parameters
  PARAMETER( x0 ); //Initial value?
  PARAMETER( log_sigmaP ); //Sigma for process model
  PARAMETER( log_sigmaO ); //Sigma for observation model
  PARAMETER( alpha ); //Amout of drift
  PARAMETER_VECTOR( x_t ); //Process
  
  // Objective 
  Type jnll = 0;
  
  //random effect
  jnll -= dnorm( x_t(0), x0, exp(log_sigmaP), true ); 
  
  //Process model -> x_t(i) ~ Normal(x_t(i-1) + alpha, sigmaP)
  for( int t=1; t<nt; t++){ 
    jnll -= dnorm( x_t(t), x_t(t-1) + alpha, exp(log_sigmaP), true );
  }
  
  //Observation model -> y_t ~ Normal(x_t(t), sigmaO)
  for( int t=0; t<nt; t++){ // Probability of observations conditional on fixed and random effect values
    jnll -= dnorm( y_t(t), x_t(t), exp(log_sigmaO), true );
  }
  
  // Reporting:
  
  //Back-transform log-sigmas to sigmas
  Type sigmaP = exp(log_sigmaP); 
  Type sigmaO = exp(log_sigmaO);
  
  //Things we want reported
  REPORT( sigmaP ); 
  REPORT( sigmaO );
  REPORT( x_t );
  
  //Things we want SEs on
  ADREPORT( sigmaP );
  ADREPORT( sigmaO );
  ADREPORT( x_t );
  
  return jnll;
}
