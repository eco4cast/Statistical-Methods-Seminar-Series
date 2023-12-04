# Set the base directory
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Original datafiles of the case study, placed in the data folder.

#	OUTPUT. Unfitted models, i.e., the list of Hmsc model(s) that have been defined but not fitted yet,
# stored in the file "models/unfitted_models.RData".
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
# LOAD PACKAGES (BEGINNING)
##################################################################################################
library(Hmsc)
library(ape) # we need this to construct a taxonomic tree
##################################################################################################
# LOAD PACKAGES (END)
##################################################################################################


##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


##################################################################################################
# READ AND SELECT SPECIES DATA (BEGINNING)
##################################################################################################
data = read.csv(file.path(dataDir, "species.csv"))
Y = as.matrix(data[,2:658])
hist(Y)
##################################################################################################
# READ AND SELECT SPECIES DATA (END)
##################################################################################################


##################################################################################################
# READ AND MODIFY ENVIRONMENTAL DATA (BEGINNING)
##################################################################################################
XData = read.csv(file.path(dataDir, "environment.csv"), as.is = FALSE)
head(XData)
XData$id = as.factor(XData$id)
plot(XData)
hist(XData$volume)
XData$volume = log(XData$volume)
hist(XData$volume)
hist(XData$decay)
XData$decay[XData$decay==4]=3
hist(XData$decay)
##################################################################################################
# READ AND MODIFY ENVIRONMENTAL DATA (BEGINNING)
##################################################################################################


##################################################################################################
# READ AND MODIFY TRAIT DATA (BEGINNING)
##################################################################################################
TrData = read.csv(file.path(dataDir, "traits.csv"), as.is = FALSE)
head(TrData)
colnames(Y)[1:10]
TrData$species[1:10]
all(colnames(Y) == TrData$species)
rownames(TrData) = TrData$species
plot(TrData)
hist(TrData$volume)
TrData$volume = log(TrData$volume)
hist(TrData$volume)
hist(TrData$shape)
TrData$shape = log(TrData$shape)
hist(TrData$shape)
##################################################################################################
# READ AND MODIFY TRAIT DATA (END)
##################################################################################################


##################################################################################################
# SELECT COMMON SPECIES (BEGINNING)
##################################################################################################
prev = colSums(Y)
hist(prev)
sum(prev>=10)
sum(prev>=20)
sel.sp = (prev>=20)
Y = Y[,sel.sp] #presence-absence data for selected species
TrData = droplevels(TrData[sel.sp,])
dim(TrData)
##################################################################################################
# SELECT COMMON SPECIES (END)
##################################################################################################


##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# STUDY DESIGN
studyDesign = data.frame(site=XData$site, id=XData$id)
# RANDOM EFFECT STRUCTURE, HERE Site (hierarchical study design)
rL.site = HmscRandomLevel(units = levels(studyDesign$site))
# and optionally id, if we are interested in species associations at that level
rL.id = HmscRandomLevel(units = levels(studyDesign$id))
# REGRESSION MODEL FOR ENVIRONMENTAL COVARIATES.
XFormula = ~ tree+volume+decay+index
# REGRESSION MODEL FOR TRAITS
TrFormula = ~ fb+orn+shape+volume
# CONSTRUCT TAXONOMICAL TREE TO BE USED AS PROXY FOR PHYLOGENETIC TREE
taxonomicTree = as.phylo(~phylum/class/order/family/genus/species,data = TrData, collapse = FALSE)
taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))
plot(taxonomicTree,cex=0.5)

# CONSTRUCT THE MODELS.

# PRESENCE-ABSENCE MODEL FOR INDIVIDUAL SPECIES (COMMON ONLY)
m = Hmsc(Y=Y, XData = XData,  XFormula = XFormula,
         TrData = TrData, TrFormula = TrFormula,
         phyloTree = taxonomicTree,
         distr="probit",
         studyDesign = studyDesign, ranLevels=list(site=rL.site, id=rL.id))
##################################################################################################
# SET UP THE MODEL (END)
##################################################################################################


##################################################################################################
# COMBINING AND SAVING MODELS (START)
##################################################################################################
models = list(m)
names(models) = c("presence-absence model")
save(models, file = file.path(modelDir, "unfitted_models.RData"))
##################################################################################################
# COMBINING AND SAVING MODELS (END)
##################################################################################################


##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
##################################################################################################
