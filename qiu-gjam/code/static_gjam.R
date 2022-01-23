library(gjam) ##versoin 2.5
library(magrittr)

#############################################
###########code for static gjam #############
#############################################

### load in data
dataList <- readRDS('data/beetles_static.rds')
xdata    <- dataList$xdata
ydata    <- dataList$ydata
edata    <- dataList$edata

### optional, trim Y to have a lower numbers of species
# species  <- gjamTrimY(ydata, 100, OTHER = FALSE)$y %>% colnames()
# ydata    <- ydata[, species] %>% as.data.frame()

### model details
effort   <- list(columns=1:ncol(ydata), values = edata$values)
rl       <- list(r = 8, N = 25) ##dimension reduction
ml       <- list(ng = 200, burnin = 100, 
                 typeNames = 'DA',  ###specify type
                 effort = effort, 
                 reductList = rl)

### model formula
fstring  <- '~def.JJA + tmean.JJA + cec30 + CwdVolume.sqrt + u1 + u2 + u3 + land'

### fit gjam
out      <- gjam(as.formula(fstring), 
                 xdata = xdata, 
                 ydata = ydata, 
                 modelList = ml)

### visualization parameters
pl              <- list(SMALLPLOTS = T, GRIDPLOTS=T, 
                        SAVEPLOTS = T, PLOTALLY = T, 
                        outFolder = 'static_output/')

### visualization
vis_out         <- gjamPlot(output = out, 
                            plotPars = pl)

#############################################
###########code for future projection #######
#############################################

### load future environment
future_data     <- readRDS('data/beetles_static_future.rds')
future_xdata    <- future_data$future_xdata
gjamOut         <- future_data$gjamOut
new_effort      <- future_data$new_effort

### prediction under RCP scenarios
new             <- list(xdata = future_xdata, 
                        nsim = 200,
                        effort = new_effort)

new_p           <- gjamPredict(output = gjamOut, newdata = new)


#############################################
###########code for dynamic gjam ############
#############################################
source('code/gjamTimeFunctions.R') ### function of gjamTimePrior
dataList          <- readRDS('data/beetles_dynamic.rds')
xdata             <- dataList$xdata
ydata             <- dataList$ydata
edata             <- dataList$edata
timeList          <- dataList$timeList ###required for gjamTime
gtmp              <- dataList$gtmp     ###required for gjamTime
zeroAlpha         <- gtmp$zeroAlpha    ###required for gjamTime

###set up beta
betaPrior         <- list(lo = list(intercept = -Inf), 
                          hi = list(intercept = Inf))
formulaBeta       <- as.formula('~anom.def.JJA + anom.tmmn.DJF')

###set up rho
formulaRho        <- as.formula( '~anom.def.JJA + anom.tmmn.DJF')
lo                <- list(intercept = -.01)
hi                <- list(intercept = .3)
rhoPrior          <- list(lo = lo, hi = hi)

###set up alpha
alphaSign              <- matrix(-1, ncol(ydata), ncol(ydata))
colnames(alphaSign)    <- rownames(alphaSign) <- colnames(ydata)
alphaSign[ zeroAlpha ] <- 0

###effort
edata                  <- edata[,colnames(ydata)] %>% as.data.frame()
effort                 <- list(columns = 1:ncol(ydata), values = edata)

###priors
priorList              <- list( formulaBeta = formulaBeta, 
                                formulaRho = formulaRho,
                                betaPrior = betaPrior, 
                                rhoPrior = rhoPrior, 
                                alphaSign = alphaSign)

tmp                    <- gjamTimePrior( xdata, ydata, edata, priorList)
timeList               <- mergeList(timeList, tmp)

###model list
modelList              <- list(typeNames = 'DA', ng = 100, burnin=50,
                               timeList=timeList, effort = effort)

###fit gjamTime
out2                   <- gjam(formulaBeta, xdata=xdata, 
                               ydata=ydata, modelList=modelList)

###visualization
pl                     <- list(SMALLPLOTS = T, GRIDPLOTS=T, 
                               SAVEPLOTS = T, PLOTALLY = T, 
                               outFolder = 'dynamic_output/')
vis_out                <- gjamPlot(output = out2, 
                                   plotPars = pl)
