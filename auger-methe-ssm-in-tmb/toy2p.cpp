/*----------------------- SECTION A --------------------------*/
// Link to the TMB package
#include <TMB.hpp>

/*----------------------- SECTION B --------------------------*/
// Define main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  /*----------------------- SECTION C --------------------------*/
  // Specify the input data
  DATA_VECTOR(y);
  
  // For one-step-ahead residuals
  DATA_VECTOR_INDICATOR(keep, y);
  
  // Specify the parameters
  PARAMETER(logSdP); // Log of st. dev. for the process variation
  PARAMETER(logSdO); // Log of st. dev. for the observation error
  
  // Specify the random effect/states
  PARAMETER_VECTOR(z);
  
  /*----------------------- SECTION D --------------------------*/
  // Transform standard deviations
  // exp(par) is a trick to make sure that the estimated sd > 0
  Type sdp = exp(logSdP);
  Type sdo = exp(logSdO);
  
  /*----------------------- SECTION E --------------------------*/
  // Define the variable that will keep track of
  // the negative log-likelihood (nll)
  Type nll = 0.0;
  
  /*----------------------- SECTION F --------------------------*/
  // Calculate the contribution to the negative log-likelihood
  // of the process equation for t=1,...,T
  // Remember that we fixed z_0 = 0
  for(int i = 1; i < z.size(); ++i){
    nll -= dnorm(z(i), z(i-1), sdp, true);
    
    //*----------------------- SECTION G --------------------------*/
    // Simulation block for process equation
    SIMULATE {
      z(i) = rnorm(z(i-1), sdp);
    }
  }
  
  /*----------------------- SECTION H --------------------------*/
  // Calculate the contribution to the negative log-likelihood
  // of the observation equation for t=1,...,T
  // Remember, the first element of z is at t=0,
  // while the first element of y is at t=1
  for(int i = 0; i < y.size(); ++i){
    nll -= keep(i)*dnorm(y(i), z(i+1), sdo, true);
    
    //*----------------------- SECTION I --------------------------*/
    // Simulation block for observation equation
    SIMULATE {
      y(i) = rnorm(z(i+1), sdo);
    }
  }
  
  
  /*----------------------- SECTION J --------------------------*/
  // State the transformed parameters to report
  // Using ADREPORT will return the point values and the standard errors
  // Note that we only need to specify this for parameters
  // we transformed, see section D above
  // The other parameters, including the random effects (states),
  // will be returned automatically
  ADREPORT(sdp);
  ADREPORT(sdo);
  
  /*----------------------- SECTION K --------------------------*/
  // Report simulated values
  SIMULATE{
    REPORT(z);
    REPORT(y);
  }
  
  
  /*----------------------- SECTION L --------------------------*/
  // State that we want the negative log-likelihood to be returned
  // This is what the optimizer in R will minimize
  return nll;
  
}
