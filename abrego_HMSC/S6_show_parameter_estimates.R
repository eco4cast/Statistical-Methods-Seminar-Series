# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. the Fitted models.

# OUTPUT: Parameter estimates illustrated (for highest RUN of S2) in the file
# "results/parameter_estimates.pdf", the text file "results/parameter_estimates.txt", 
# as well as given numerically in multiple csv files (one per parameter type) named 
# "results/parameter_estimates_[parameter_name].csv".
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
support.level.beta = NULL #Default: 0.95
support.level.gamma = NULL #Default: 0.95
support.level.omega = NULL #Default: 0.9
var.part.order.explained = NULL #Default: in variance partitioning of explained variance, species are shown in the order they are in the model
var.part.order.raw = NULL #Default: in variance partitioning of raw variance, species are shown in the order they are in the model
show.sp.names.beta = NULL #Default: species names shown in beta plot if there are at most 30 species and no phylogeny
plotTree = NULL #Default: tree is plotted in Beta plot if the model includes it
omega.order = NULL #Default: species shown in the order they are in the model
show.sp.names.omega = NULL #Default: species names shown in beta plot if there are at most 30 species
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
#use support levels to selected the level of statistical support shown
#support.level.beta = 0.8
#support.level.gamma = 0.8
#support.level.omega = 0.8

#use var.part.order.explained to select which order species are shown in the raw variance partitioning
#var.part.order.raw should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "decreasing" if you wish to order according to explanatory power
#var.part.order.explained = list()
#var.part.order.explained[[1]] = 0
#var.part.order.explained[[2]] = c(2,1)

#use var.part.order.raw to select which order species are shown in the explained variance partitioning
#same options apply as for var.part.order.explained
#var.part.order.raw = list()
#var.part.order.raw[[1]] = "decreasing"
#var.part.order.raw[[2]] = c(1,2)

#use show.sp.names.beta to choose to show / not show species names in betaPlot
#if given, show.sp.names.beta should be a vector with length equalling number of models
#show.sp.names.beta = TRUE

#use plotTree to choose to plot / not plot the tree in betaPlot
#if given, plotTree should be a vector with length equalling number of models
#plotTree = c(FALSE,FALSE)

#use omega.order to select which order species are shown in omega plots
#omega.order should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "AOE" if you wish to use the angular order of the eigenvectors.
omega.order = list()
omega.order[[1]] = "AOE"
omega.order[[2]] = "AOE"
#Default: species shown in the order they are in the model
#show.sp.names.omega = c(TRUE,FALSE) #Default: species names shown in beta plot if there are at most 30 species
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(support.level.beta)) support.level.beta = 0.95
if(is.null(support.level.gamma)) support.level.gamma =  0.95
if(is.null(support.level.omega)) support.level.omega =  0.9

library(Hmsc)
library(colorspace)
library(corrplot)
library(writexl)

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nst = length(thin_list)
nChains = 4

text.file = file.path(resultDir,"/parameter_estimates.txt")
cat(c("This file contains additional information regarding parameter estimates.","\n","\n",sep=""),file=text.file)

for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){break}
}
if(file.exists(filename)){
  load(filename)
  cat(c("\n",filename,"\n","\n"),file=text.file,sep="",append=TRUE)
  nm = length(models)
  if(is.null(var.part.order.explained)){
    var.part.order.explained = list()
    for(j in 1:nm) var.part.order.explained[[j]] = 0
  }
  if(is.null(var.part.order.raw)){
    var.part.order.raw = list()
    for(j in 1:nm) var.part.order.raw[[j]] = 0
  }
  if(is.null(omega.order)){
    omega.order = list()
    for(j in 1:nm) omega.order[[j]] = 0
  }
  
  modelnames = names(models)
  
  pdf(file= file.path(resultDir,"parameter_estimates.pdf"))
  for(j in 1:nm){
    cat(c("\n",names(models)[j],"\n","\n"),file=text.file,sep="",append=TRUE)
    m = models[[j]]
    if(m$XFormula=="~."){
      covariates = colnames(m$X)[-1]
    } else {
      covariates = attr(terms(m$XFormula),"term.labels")
    }
    if(m$nr+length(covariates)>1 & m$ns>1){
      preds = computePredictedValues(m)
      VP = computeVariancePartitioning(m)
      vals = VP$vals
      mycols = rainbow(nrow(VP$vals))
      MF = evaluateModelFit(hM=m, predY=preds)
      R2 = NULL
      if(!is.null(MF$TjurR2)){
        TjurR2 = MF$TjurR2
        vals = rbind(vals,TjurR2)
        R2=TjurR2
      }
      if(!is.null(MF$R2)){
        R2=MF$R2
        vals = rbind(vals,R2)
      }
      if(!is.null(MF$SR2)){
        R2=MF$SR2
        vals = rbind(vals,R2)
      }
      filename = file.path(resultDir, paste("parameter_estimates_VP_",modelnames[j],".csv"))
      write.csv(vals,file=filename)
      if(!is.null(VP$R2T$Beta)){
        filename = file.path(resultDir,paste("parameter_estimates_VP_R2T_Beta",modelnames[j],".csv"))
        write.csv(VP$R2T$Beta,file=filename)
      }
      if(!is.null(VP$R2T$Y)){
        filename = file.path(resultDir, paste("parameter_estimates_VP_R2T_Y",modelnames[j],".csv"))
        write.csv(VP$R2T$Y,file=filename)
      }
      if(all(var.part.order.explained[[j]]==0)){
        c.var.part.order.explained = 1:m$ns
      } else {
        if(all(var.part.order.explained[[j]]=="decreasing")){
          c.var.part.order.explained = order(R2, decreasing = TRUE)
        }
        else {
          c.var.part.order.explained  = var.part.order.explained[[j]]
        }
      }
      VP$vals = VP$vals[,c.var.part.order.explained]
      cat(c("\n","var.part.order.explained","\n","\n"),file=text.file,sep="",append=TRUE)
      cat(m$spNames[c.var.part.order.explained],file=text.file,sep="\n",append=TRUE)
      plotVariancePartitioning(hM=m, VP=VP, main = paste0("Proportion of explained variance, ",modelnames[j]), cex.main=0.8, cols = mycols, args.leg=list(bg="white",cex=0.7))
      if(all(var.part.order.raw[[j]]==0)){
        c.var.part.order.raw = 1:m$ns
      } else {
        if(all(var.part.order.raw[[j]]=="decreasing")){
          c.var.part.order.raw = order(R2, decreasing = TRUE)
        }
        else {
          c.var.part.order.raw  = var.part.order.raw[[j]]
        }
      }
      if(!is.null(R2)){
        VPr = VP
        for(k in 1:m$ns){
          VPr$vals[,k] = R2[k]*VPr$vals[,k]
        }
        VPr$vals = VPr$vals[,c.var.part.order.raw]
        cat(c("\n","var.part.order.raw","\n","\n"),file=text.file,sep="",append=TRUE)
        cat(m$spNames[c.var.part.order.raw],file=text.file,sep="\n",append=TRUE)
        plotVariancePartitioning(hM=m, VP=VPr,main=paste0("Proportion of raw variance, ",modelnames[j]),cex.main=0.8, cols = mycols, args.leg=list(bg="white",cex=0.7),ylim=c(0,1))
      }
    }
  }
  for(j in 1:nm){
    m = models[[j]]
    if(m$nc>1){
      postBeta = getPostEstimate(m, parName="Beta")
      filename = file.path(resultDir, paste("parameter_estimates_Beta_",modelnames[j],".xlsx"))
      me = as.data.frame(t(postBeta$mean))
      me = cbind(m$spNames,me)
      colnames(me) = c("Species",m$covNames)
      po = as.data.frame(t(postBeta$support))
      po = cbind(m$spNames,po)
      colnames(po) = c("Species",m$covNames)
      ne = as.data.frame(t(postBeta$supportNeg))
      ne = cbind(m$spNames,ne)
      colnames(ne) = c("Species",m$covNames)
      vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
      writexl::write_xlsx(vals,path = filename)
      if(is.null(show.sp.names.beta)){
        c.show.sp.names = (is.null(m$phyloTree) && m$ns<=30) 
      } else {
        c.show.sp.names = show.sp.names.beta[j]
      }
      c.plotTree = !is.null(m$phyloTree)
      if(!is.null(plotTree)){
        c.plotTree = c.plotTree & plotTree[j]
      }
      plotBeta(m, post=postBeta, supportLevel = support.level.beta, param="Sign",
               plotTree = c.plotTree,
               covNamesNumbers = c(TRUE,FALSE),
               spNamesNumbers=c(c.show.sp.names,FALSE),
               cex=c(0.6,0.6,0.8))
      mymain = paste0("BetaPlot, ",modelnames[j])
      if(!is.null(m$phyloTree)){
        mpost = convertToCodaObject(m)
        rhovals = unlist(poolMcmcChains(mpost$Rho))
        mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",round(mean(rhovals>0),2))
      }
      title(main=mymain, line=2.5, cex.main=0.8)
    }
  }
  for(j in 1:nm){
    m = models[[j]]      
    if(m$nt>1 & m$nc>1){
      postGamma = getPostEstimate(m, parName="Gamma")
      plotGamma(m, post=postGamma, supportLevel = support.level.gamma, param="Sign",
                covNamesNumbers = c(TRUE,FALSE),
                trNamesNumbers=c(m$nt<21,FALSE),
                cex=c(0.6,0.6,0.8))
      title(main=paste0("GammaPlot ",modelnames[j]), line=2.5,cex.main=0.8)
    }
  }
  for(j in 1:nm){
    m = models[[j]]
    if(m$nr>0 & m$ns>1){
      OmegaCor = computeAssociations(m)
      for (r in 1:m$nr){
        toPlot = ((OmegaCor[[r]]$support>support.level.omega) + (OmegaCor[[r]]$support<(1-support.level.omega))>0)*sign(OmegaCor[[r]]$mean)
        if(is.null(show.sp.names.omega)){
          c.show.sp.names = (m$ns<=30) 
        } else {
          c.show.sp.names = show.sp.names.omega[j]
        }
        if(!c.show.sp.names){
          colnames(toPlot)=rep("",m$ns)
          rownames(toPlot)=rep("",m$ns)
        }
        if(all(omega.order[[j]]==0)){
          plotOrder = 1:m$ns
        } else {
          if(all(omega.order[[j]]=="AOE")){
            plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
          } else {
            plotOrder = omega.order[[j]]
          }
        }
        cat(c("\n","omega.order","\n","\n"),file=text.file,sep="",append=TRUE)
        cat(m$spNames[plotOrder],file=text.file,sep="\n",append=TRUE)
        mymain = paste0("Associations, ",modelnames[j], ": ",names(m$ranLevels)[[r]])
        if(m$ranLevels[[r]]$sDim>0){
          mpost = convertToCodaObject(m)
          alphavals = unlist(poolMcmcChains(mpost$Alpha[[1]][,1]))
          mymain = paste0(mymain,", E[alpha1] = ",round(mean(alphavals),2),", Pr[alpha1>0] = ",round(mean(alphavals>0),2))
        }
        corrplot(toPlot[plotOrder,plotOrder], method = "color",
                 col=colorRampPalette(c("blue","white","red"))(3),
                 mar=c(0,0,1,0),
                 main=mymain,cex.main=0.8)
        me = as.data.frame(OmegaCor[[r]]$mean)
        me = cbind(m$spNames,me)
        colnames(me)[1] = ""
        po = as.data.frame(OmegaCor[[r]]$support)
        po = cbind(m$spNames,po)
        colnames(po)[1] = ""
        ne = as.data.frame(1-OmegaCor[[r]]$support)
        ne = cbind(m$spNames,ne)
        colnames(ne)[1] = ""
        vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
        filename = file.path(resultDir, paste("parameter_estimates_Omega_",modelnames[j],"_",names(m$ranLevels)[[r]],".xlsx"))
        writexl::write_xlsx(vals,path = filename)
      }
    }
  }
  dev.off()
}

