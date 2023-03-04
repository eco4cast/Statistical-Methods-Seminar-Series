/*----------------------- SECTION A --------------------------*/
// State that we need the TMB package
#include <TMB.hpp>
// Needed for the multivariate normal function: MVNORM_t
using namespace density;

/*----------------------- SECTION B --------------------------*/
// Define main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  /*----------------------- SECTION C --------------------------*/
  // Specify the input data
  DATA_ARRAY(y); // Longitude and latitude of locations
  DATA_VECTOR(idx); // Index linking irregular obs. time i to state time t
  DATA_VECTOR(ji); // When is obs. i relative to state t-1 and t
  DATA_MATRIX(ac); // Argos category
  
  // For one-step-ahead residuals
  DATA_ARRAY_INDICATOR(keep, y);
  
  // Input parameters - i.e. parameters to estimate.
  PARAMETER(logitGamma); // Autocorrelation
  PARAMETER(logSdLon); // Process standard deviation in lon
  PARAMETER(logSdLat); // Process standard deviation in lat
  PARAMETER(logPsiLon); // Scaling parameter for error values
  PARAMETER(logPsiLat); // Scaling parameter for error values
  
  // The true/unobserved locations of the animal, i.e states
  PARAMETER_MATRIX(z);
  
  /*----------------------- SECTION D --------------------------*/
  // Transformation of the input parameters to model format
  /* These transformations are made to insured that
   the parameters have sensical values.
   They do not change the model, they are only a computational trick. */
  Type gamma = 1.0/(1.0+exp(-logitGamma)); // b/c we want 0 < gamma < 1
  Type sdLon = exp(logSdLon); // b/c we want sd > 0
  Type sdLat = exp(logSdLat); // b/c we want sd > 0
  Type psiLon = exp(logPsiLon); // b/c we want psi > 0
  Type psiLat = exp(logPsiLat); // b/c we want psi > 0
  
  /*----------------------- SECTION E --------------------------*/
  // Setting the bivariate normal
  // Variance-covariance matrix for the process equation.
  // We are using a bivariate normal even if there is no correlation.
  matrix<Type> covs(2,2);
  covs << sdLon*sdLon, 0, 0, sdLat*sdLat;
  
  // Notes:
  // - the mean of MVNORM_t is fixed at 0,
  //   thus we only specify the variance-covariance matrix.
  // - MVNORM_t evaluates the negative log density
  MVNORM_t<Type> nll_dens(covs);
  
  /*----------------------- SECTION F --------------------------*/
  // Creating a variable that keeps track of the negative log likelihood
  Type nll = 0.0;
  
  /*----------------------- SECTION G --------------------------*/
  // Creating a variable for the value used as data in bivariate normal
  vector<Type> tmp(2);
  
  /*----------------------- SECTION H --------------------------*/
  // Initializing the model
  // For the 2nd time step, we assume that we have a simple random walk:
  // z_1 = z_0 + eta
  // Here tmp is simply the difference between the first two states:
  // d_1 = z_1 - z_0
  tmp = z.row(1) - z.row(0);
  nll += nll_dens(tmp);
  
  /*----------------------- SECTION I --------------------------*/
  SIMULATE{
    z.row(1) = vector<Type>(z.row(0)) + nll_dens.simulate();
  }
  
  /*----------------------- SECTION J --------------------------*/
  // nll contribution of the process equation after the 2nd time step
  // Notes:
  // - .row(i) gives the ith row of parameter or data matrix,
  //   so locations at time i.
  // - loop here is the size of the states vector z,
  //   thus, if we have missing data this loop will be longer
  //   than the observation equation
  //   and will estimate the state value even for
  //   time steps where we have no observation.
  for(int i = 2; i < z.rows(); ++i){
    tmp = (z.row(i) - z.row(i-1)) - (z.row(i-1) - z.row(i-2)) * gamma;
    nll += nll_dens(tmp);
    
    /*----------------------- SECTION K --------------------------*/
    SIMULATE{
      z.row(i) = vector<Type>((1+gamma)*z.row(i-1) - gamma*z.row(i-2)) + 
        nll_dens.simulate();
    }
  }
  
  /*----------------------- SECTION L --------------------------*/
  // nll contribution of the observation equation
  // The loop here is just the size of the observation vector.
  // We use the time index found in idx to relate
  // the appropriate state to the observation.
  for(int i = 0; i < y.matrix().rows(); ++i){
    // Interpolate the value at time of observation
    // CppAD::Integer because the index used to get value
    // in a vector in c++ needs to be a integer to work.
    tmp = y.matrix().row(i) -
      ((1.0-ji(i))*z.row(CppAD::Integer(idx(i)-1)) +
      ji(i)*z.row(CppAD::Integer(idx(i))));
    
    // We are using the non-standardised t distribution,
    // that is why we are correcting the pdf.
    // Longitude
    nll -= keep(i,0) * (log(1/(psiLon*ac(i,0))) +
      dt(tmp(0)/(psiLon*ac(i,0)),ac(i,1),true));
    // Latitude
    nll -= keep(i,1) * (log(1/(psiLat*ac(i,2))) +
      dt(tmp(1)/(psiLat*ac(i,2)),ac(i,3),true));
    
    
    /*----------------------- SECTION M --------------------------*/
    // Simulation of observations
    y(i,0) = psiLon * ac(i,0) * rt(ac(i,1)) + 
    ((1.0-ji(i))*z(CppAD::Integer(idx(i)-1),0) +
    ji(i)*z(CppAD::Integer(idx(i)),0));
    y(i,1) = psiLat * ac(i,2) * rt(ac(i,3)) + 
      ((1.0-ji(i))*z(CppAD::Integer(idx(i)-1),1) +
      ji(i)*z(CppAD::Integer(idx(i)),1));
  }
  
  /*----------------------- SECTION N --------------------------*/
  // Report the parameters and their standard errors in their model format
  ADREPORT(gamma);
  ADREPORT(sdLon);
  ADREPORT(sdLat);
  ADREPORT(psiLon);
  ADREPORT(psiLat);
  
  /*----------------------- SECTION O --------------------------*/
  // Report simulated values
  SIMULATE{
    REPORT(z)
    REPORT(y)
  }
  
  return nll;
}
