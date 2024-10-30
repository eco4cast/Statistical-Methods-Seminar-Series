##Read in the data for our model

mod_dat = readRDS("./mod_dat.rds")

##Load the RTMB library
library(RTMB)

##Number of years in the model
Y = 10

##Adjust the model data to line up with 1 to 10 which what will be used in the model
POPs = mod_dat$POPs
POPs$byear_J = POPs$byear_J - 20


##The model negative log likelihood function
model <- function(parameters){

  ##Our parameters
  ##----------------

  ##Initial adult population
  ## We pass it in on the log scale since we want it to be positive and
  ##from the optimizers perspective between -Inf and Inf
  init_N_a = exp(parameters$log_init_N_a)

  ##Rate of growth of the (adult) population
  r_of_growth = exp(parameters$log_r_of_growth)

  ##The population dynamics
  ##----------------------

  ##The vector  containing our adult abundance
  N_a = numeric(Y)
  ##The first year is our init_N_a parameter
  N_a[1] = init_N_a
  ##Now we calculate the rest using the rate of growth
  for(i in 2:Y){
    N_a[i] = N_a[i-1]*r_of_growth
  }

  ##The expected POP probability
  ##------------------------------
  ##Since fecundity is constant among all adults it's just 1/N_a in a given year

  e_POP_prob = numeric(Y)
  for(i in 1:Y){
    e_POP_prob[i] = 1/N_a[i]
  }



  ##The negative log-likelihood
  ##-----------------------------
  nll = 0

  ##We go through the data.frame of POPs subtracting each element to the nll
  for(i in 1:nrow(POPs)){
    nll = nll - dbinom(POPs$n_POP[i],POPs$n_comp[i],e_POP_prob[i],TRUE) 
  }

  ##output some quantities of interest
  REPORT(N_a)
  REPORT(e_POP_prob)

  ##Get the standard errors of some quantities
  ADREPORT(N_a)

  ##Return the negative log-likelihood
  nll
}

##Create the list of initial parameters
params = list(
  log_init_N_a = log(300),
  log_r_of_growth = log(1)
)

##Create the RTMB model object from our model
obj = MakeADFun(model,params)

##optimize the model to find the best set of parameters that explain this data

##obj$par is the intial parameters we gave the model
##obj$fn is our negative log likelihood function
##obj$gr is the gradient or first derivative of the negative log likelihood function created for us by RTMB
opt = nlminb(obj$par,obj$fn,obj$gr)

##Look at the objects we wanted output
report = obj$report()
report

##Have RTMB find the standard errors of the parameters and the adult population size
sdrep = sdreport(obj)
sdrep
ssdr = summary(sdrep)

##Create 95% Confidence interval for pop size
N_a_full = ssdr[rownames(ssdr) == "N_a",]
N_a_lower = N_a_full[,1]-1.96*N_a_full[,2]
N_a_upper = N_a_full[,1]+1.96*N_a_full[,2]

##Compare the size of the adult population 
plot(1:10,report$N_a,type="l",col="red",ylim=c(min(N_a_lower)*.90,max(N_a_upper)*1.10))
lines(1:10,sim$adult_pop[21:30])

##Add the CIs
lines(1:10,N_a_lower,col="red",lty=2)
lines(1:10,N_a_upper,col="red",lty=2)
