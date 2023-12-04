# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
# INPUT. Fitted models.

#	OUTPUT. Model fits computed by the cross-validation, with fitting
# (which is part of cross-validation) done for multiple RUNs:
# first short MCMC chains (to provide some results fast), and then with increasingly long MCMC chains
# (up to the longest run performed in S2). The results are stored in the files
# "models/MF_thin_1_samples_5_chains_4.Rdata" (RUN 0),
# "models/MF_thin_1_samples_250_chains_4.Rdata" (RUN 1),
# "models/MF_thin_10_samples_250_chains_4.Rdata" (RUN 2), 
# "models/MF_thin_100_samples_250_chains_4.Rdata" (RUN 3), and so on.
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
nfolds = NULL #Default: two-fold cross-validation
cv.level = "site"
#cv.level = "id"
nParallel = NULL #Default: nParallel = nChains
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
#nfolds = 2 #change the number of CV-folds
#nParallel = 1 #Set to 1 to disable parallel computing
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(nfolds)) nfolds = 2

library(Hmsc)
samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nChains = 4
if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename.in = file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
  filename.out = file.path(modelDir,paste("MF_thin_", as.character(thin),
                                          "_samples_", as.character(samples),
                                          "_chains_",as.character(nChains),
                                          "_nfolds_", as.character(nfolds),
                                          "_cvlevel_",cv.level,
                                          ".Rdata",sep = ""))
  if(file.exists(filename.out)){
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    print("model fit had been computed already")
  } else {
    if(file.exists(filename.in)){
      print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
      print(date())
      load(file = filename.in) #models
      nm = length(models)
      
      MF = list()
      MFCV = list()
      WAIC = list()
      
      for(mi in 1:nm){
        print(paste0("model = ",names(models)[mi]))
        m = models[[mi]]
        preds = computePredictedValues(m)
        MF[[mi]] = evaluateModelFit(hM=m, predY=preds)
        partition = createPartition(m, nfolds = nfolds, column = cv.level) #USE column = ...
        preds = computePredictedValues(m,partition=partition, nParallel = nParallel)
        MFCV[[mi]] = evaluateModelFit(hM=m, predY=preds)
        WAIC[[mi]] = computeWAIC(m)
      }
      names(MF)=names(models)
      names(MFCV)=names(models)
      names(WAIC)=names(models)
      
      save(MF,MFCV,WAIC,file = filename.out)
    }
  }
  Lst = Lst + 1
}
