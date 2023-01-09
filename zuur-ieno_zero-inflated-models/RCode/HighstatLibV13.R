#Library files for courses provided by: Highland Statistics Ltd.
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

#Copyright Highland Statistics LTD.

#####################################################################
#VIF FUNCTION.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS





##################################################################
##################################################################
#Here are some functions that we took from the pairs help file and
#modified, or wrote ourselves. To cite these, use the r citation: citation()

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
     cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}

##################################################################
panel.smooth2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.smooth = "black", span = 2/3, iter = 3, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = 1, ...)
}

##################################################################
panel.lines2=function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                       cex = 1, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}
}

##################################################################
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}
##################################################################
##################################################################



##################################################################
##################################################################
#Functions for variograms
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
#Make a variogram for one variable
#To use, type:  MyVariogram(XUTM, YUTM, E , MyDistance=10)
# XUTM is x coordinates
# XUTM is y coordinates
# E is variable used in sample variogram
# MyDistance is the cutoff value for the distances

MyVariogram <- function(x,y,z, MyDistance) {
  library(gstat)
  mydata      <- data.frame(z, x, y)
  coordinates(mydata)    <- c("x", "y")  
  Var <- variogram(z ~ 1, mydata, cutoff = MyDistance)
  data.frame(Var$np, Var$dist, Var$gamma)
}

#Function for making multiple variograms in an xyplot
#To use, type:  MultiVariogram(Z, MyVar,XUTM, YUTM, MyDistance=10)
# Z is a data frame with all the data 
# Character string with variable names that will be used in the xyplot
# XUTM is x coordinates
# XUTM is y coordinates
# MyDistance is the cutoff value for the distances

MultiVariogram <- function(Z, MyVar, x, y, MyDistance) {
  #Z is the data frame with data
  #MyVar is a list of variables for for which variograms are calculated
  #x, y: spatial coordinates
  #MyDistance: limit for distances in the variogram
  
  library(lattice)
  VarAll<- c(NA,NA,NA,NA)
  for (i in MyVar){
    vi <- MyVariogram(x,y,Z[,i], MyDistance)
    vii <- cbind(vi, i)
    VarAll <- rbind(VarAll,vii)
  }
  VarAll <- VarAll[-1,]
  
  P <- xyplot(Var.gamma ~ Var.dist | factor(i), col = 1, type = "p", pch = 16,
              data = VarAll,
              xlab = "Distance",
              ylab = "Semi-variogram",
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = TRUE,
                            x = list(relation = "same"),
                            y = list(relation = "same"))
              )
  
  print(P)
}


MultiVariogramVector <- function(Z, MyVar, x, y, MyDistance, ID) {
  #Z is the data frame with data
  #MyVar is the variable for for which variograms are calculated
  #x, y: spatial coordinates
  #MyDistance: limit for distances in the variogram
  #ID is the grouping variabl
  library(ggplot2)
  VarAll<- c(NA,NA,NA,NA, NA)
  
  for (i in levels(ID)){
    Zi <- subset(Z, ID == i)
    Zi <- droplevels(Zi)
    vi <- MyVariogram(Zi[,x],Zi[,y],Zi[,MyVar], MyDistance)
    vi$ID <- i
    VarAll <- rbind(VarAll,vi)
  }
  VarAll <- VarAll[-1,]
  
  p <- ggplot(data = VarAll,
              aes(x = Var.dist, 
                  y = Var.gamma))
  p <- p + geom_point(size = 0.5, alpha = 0.5)
  p <- p + geom_smooth(se = FALSE,
                       col = "black",
                       span = 0.9)
  p <- p + xlab("Distance") + ylab("Semi-variogram")
  p <- p + theme(text = element_text(size = 7))
  
  p <- p + theme(legend.position="none") 
  p <- p + facet_wrap(~ID, scale = "free_y")
  print(p)
}

#End variogram code
##########################################################

#Function for multi-panel Cleveland dotplot.
#The input file must contain no categorical variables
Mydotplot <- function(DataSelected){

P <- dotplot(as.matrix(as.matrix(DataSelected)),
          groups=FALSE,
          strip = strip.custom(bg = 'white',
                               par.strip.text = list(cex = 1.2)),
          scales = list(x = list(relation = "free", draw = TRUE),
                        y = list(relation = "free", draw = FALSE)),
          col=1, cex  = 0.5, pch = 16,
          xlab = list(label = "Value of the variable", cex = 1.5),
          ylab = list(label = "Order of the data from text file", cex = 1.5))
  
print(P)  
  }


#Add more code here:


Mybwplot <- function(Z, MyVar, TargetVar){
#Multipanel boxplots
#Z: data set
#MyVar: character string
#TargetVar: variable for the x-axis..must be a factor
  
  AllY <- as.vector(as.matrix(Z[,MyVar]))
  AllX <- rep(Z[,TargetVar], length(MyVar))
  ID <- rep(MyVar, each = nrow(Z))
  
P <- bwplot(AllY ~ factor(AllX) | ID, horizontal = FALSE,
         ylab = "", xlab = "",
         scales = list(alternating = TRUE,cex.lab = 1.5,
                       x = list(relation = "same",rot =90, abbreviate = TRUE, cex = 1.5),
                       y = list(relation = "free", draw = FALSE)),
         strip = strip.custom(bg = 'white',
                              par.strip.text = list(cex = 1.2)),
         cex = .5,
         par.settings = list(
           box.rectangle = list(col = 1),
           box.umbrella  = list(col = 1),
           plot.symbol   = list(cex = .5, col = 1)))
print(P)
  }



#######################################################
MyxyplotBin <- function(Z, MyV, NameY1) {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  
  
  library(mgcv)
  library(lattice)
  
  P <- xyplot(AllY ~ AllX | factor(AllID), col = 1,
              strip = function(bg='white', ...) strip.default(bg='white', ...),
              scales = list(alternating = TRUE, 
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              xlab = "Covariate",
              ylab = "Probability of presence",
              panel=function(x,y){
                panel.grid(h=-1, v= 2)
                panel.points(x,y,col=1)
                tmp<-gam(y~s(x, k = 4), family = binomial)
                MyData <- data.frame(x = seq(min(x), max(x), length = 25))
                p1 <- predict(tmp, newdata = MyData, type ="response")
                panel.lines(MyData$x,p1, col = 1, lwd = 3)
              })
  
  print(P)
}
#######################################################

#######################################################
Myxyplot <- function(Z, MyV, NameY1, MyXlab = "", MyYlab="") {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  
  
  library(mgcv)
  library(lattice)
  
  P <- xyplot(AllY ~ AllX|factor(AllID), col = 1,
              xlab = list(MyXlab, cex = 1.5),
              #ylab = list("Response variable", cex = 1.5),
              #ylab = list("Pearson residuals", cex = 1.5),
              ylab = list(MyYlab, cex = 1.5),
              #layout = c(2,2),   #Modify
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = TRUE,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = 1)
                panel.loess(x, y, span = 0.8,col = 1, lwd = 2)
                })
  
  print(P)
}
#######################################################

MyxyplotPolygon <- function(Z, MyV, NameY1) {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  
  
  library(mgcv)
  library(lattice)
  Z <- xyplot(AllY ~ AllX|factor(AllID), col = 1,
              xlab = list(label = "Explanatory variables", cex = 1.5),
              ylab = "",
              strip = function(bg='white',cex.lab = 1.5,...)
                strip.default(bg='white', ...),
              scales = list(alternating = TRUE,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                t1 <- gam(y~s(x))
                MD1 <- data.frame(x=seq(from = min(x, na.rm = TRUE),
                                        to = max(x, na.rm = TRUE),
                                        length = 100))
                P1 <- predict(t1,   se.fit = TRUE)
                I1 <- order(x)
                xs <- sort(x)
                panel.lines(xs, P1$fit[I1], col = 1)
                panel.polygon(c(xs, rev(xs)),
                              c(P1$fit[I1]-2*P1$se.fit[I1],
                                rev(P1$fit[I1]+2*P1$se.fit[I1])),
                              col = gray(0.7),
                              density = 10 )
                panel.grid(h=-1, v= 2)
                panel.abline(0,0)
                panel.points(x, y, col = 1)
                
              })
  #Because the xyplot is inside a function you need to print 
  #construction below
  print(Z)
}

################################################
#Mypairs
#Make fancy pair plots
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
      cex.labels =  2,
      lower.panel = function(x, y, digits=2, prefix="", cex.cor = 6) {
        panel.cor(x, y, digits, prefix, cex.cor)}, 
      upper.panel =  function(x, y) points(x, y, 
                                           pch = 16, cex = 0.8, 
                                           col = gray(0.1)))
 #print(P)
}


###########################################################
#Show double zeros in a row-by-column matrix
#To run: ShowZeros(Species)
#Better use short names for the row labels

ShowZeros <- function(XX) {
  #To run this function type: 
  #  ShowZeros(Species)
  	
  #Load packages	
  library(ggplot2)
  library(reshape2)
  library(zoo)
	
    N <- ncol(XX)
    X <- matrix(nrow = N, ncol = N)
    for (i in 1:N) {
	  for (j in i:N){
       X[i,j] <- sum(XX[,i] == 0 & XX[,j] == 0, na.rm = TRUE)
	   X[j,i] <- X[i,j]
	  }
    }
    Zeros <- X / nrow(XX)
    colnames(Zeros) <- colnames(XX)
    rownames(Zeros) <- colnames(XX)




    #Taken from: http://www.r-bloggers.com/controlling-heatmap-colors-with-ggplot2/
    #Set color representation for specific values of the data distribution
    BreakPoints <- c(0, 0.25, 0.5, 0.75, 1)
    ColourQuantiles <- quantile(Zeros, probs = BreakPoints)

    ## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
    color_palette <- colorRampPalette(c("#edf8b1", "#7fcdbb", "#2c7fb8"))(length(ColourQuantiles) - 1)
  
    #Or change to:
    #color_palette[4] <- "red"
    #color_palette[3] <- "orange"
    #color_palette[2] <- "yellow"
    #color_palette[1] <- "white"

 
    # prepare label text (use two adjacent values for range text)
    label_text <- rollapply(round(ColourQuantiles, 2), 
                            width = 2, 
                            by = 1, 
                            FUN = function(i) paste(i, collapse = " - "))
 
    # discretize matrix; this is the most important step, 
    # where for each value we find category of predefined 
    # ranges (modify probs argument of quantile to detail the colors)
    mod_mat <- matrix(findInterval(Zeros, ColourQuantiles, all.inside = TRUE), 
                      nrow = nrow(Zeros))
 
 
    #Prepare the graph
    Zeros.m          <- melt(Zeros)
    Zeros.m$NewValue <- as.vector(mod_mat)
    p <- ggplot(Zeros.m, aes(x = Var1, y = Var2, fill = factor(NewValue))) 
    p <- p + geom_tile(color = "black")
    p <- p + scale_fill_manual(values = color_palette, name = "", labels = label_text) 
    p <- p + scale_x_discrete(expand = c(0, 0))
    p <- p + scale_y_discrete(expand = c(0, 0))
    p <- p + theme(axis.text.x = element_text(angle=45, hjust = 1, size = 5))
    p <- p + xlab("") + ylab("") 
    p <- p + labs(title = "Percentage zeros in common")
    print(p)
}
################################################################################


MyDotplot.ggp2 <- function(Z, varx, Ncol = 5, TextSize = 15, PointSize = 1,
                           DropLabels = FALSE, MyTitle = "") {
  library(ggplot2)	
  K     <- length(varx)
  MyData <- data.frame(Y = rep(1:nrow(Z), K),
                       X = as.vector(as.matrix(Z[, varx])),
                       Var = rep(varx, each = nrow(Z))) 

  p <- ggplot(MyData, aes(y = Y, x = X))
  p <- p + geom_point(size= PointSize) + ylab("Order of the data") + xlab("Range of the data")
  if (DropLabels) { 
  	   p <- p + theme(axis.text = element_blank(),
  	                  strip.text = element_text(size=TextSize),
  	                  axis.ticks = element_blank() )} else {
       p <- p + theme(strip.text = element_text(size=TextSize))	    	
  	    }
  p <- p + facet_wrap(~ Var, scales = "free_x", ncol=Ncol)
  p <- p + labs(title = MyTitle)
  p	
}

Myfapply <- function(Zfactor, INDEX) {
  INDEX <- factor(INDEX)	
  AllLevels <- levels(INDEX)
  MyNewFactor <- NULL
  for (i in AllLevels) {
	ThisLevel <- as.character(Zfactor[INDEX == i] [1])
	MyNewFactor <- c(MyNewFactor, ThisLevel)
  }
  MyNewFactor <- factor(MyNewFactor)
  MyNewFactor	
}


######################################################
MyMultipanel.ggp2 <- function(Z, varx, vary, 
                              ylab = "Response variable",
                              addSmoother = FALSE,
                              addRegressionLine = FALSE,
                              addHorizontalLine = FALSE) {
  K <- length(varx)
  MyData <- data.frame(Y = rep(as.vector(as.matrix(Z[,vary])), K),
                       X = as.vector(as.matrix(Z[, varx])),
                       Var = rep(varx, each = nrow(Z))) 
  library(ggplot2)
  p <- ggplot(MyData, aes(y = Y, x = X))
  p <- p + geom_point(size = 0.5) + ylab(ylab) + xlab("Covariates")
  p <- p + theme(text = element_text(size=15))
  if (addSmoother == TRUE) {
  	 p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1)
  }
  if (addRegressionLine == TRUE) {
  	 p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1, method = "lm")
  }
  if (addRegressionLine == TRUE) {
  	 p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1, method = "lm")
  }
  if (addHorizontalLine == TRUE) {
  	 p <- p + geom_hline(yintercept = 0)
  }
  p <- p + facet_wrap(~ Var, scales = "free_x")
  suppressMessages(print(p)) 	
}
######################################################


#Standardize the continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}



# Support functions for Volume II of the INLA
# book
tentf <- function(xCov, xKnots, j){
  dj <- xKnots * 0
  dj[j] <- 1
  approx(xKnots, dj, xCov)$y
}
tentf.X <- function(xCov, xKnots){
  nk <- length(xKnots)
  n  <- length(xCov)
  X  <- matrix(NA, n, nk)
  for (j in 1:nk){ X[,j]<- tentf(xCov, xKnots,j)}
  X
}



LongLatToUTM <- function(x,y,zone, Hemisphere = "north"){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")  
  if (Hemisphere == "north"){
     res <- spTransform(xy, 
                        CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  }
  
  if (Hemisphere == "south"){
    res <- spTransform(xy, 
                       CRS(paste("+proj=utm +zone=",zone,
                     " +ellps=WGS84 +datum=WGS84 +units=m +no_defs +south",sep='')))
  }
  
  return(as.data.frame(res))
}



# Function to plot the betas and CIs of different models next to each other
MyCompareBetasofModels <- function(AllModels, ModelNames ){
  # Determine the number of models:
  NumModels <- length(AllModels)
  
  # Get the model results out of the list and rbind the values
  Combined <- NULL
  for (i in 1:NumModels){
    Combined <- rbind(Combined, AllModels[[i]])
  }
  
  #What are the names of the betas?
  BetaNames <- rownames(AllModels[[i]])
  
  #How many betas in each model?
  K <- length(BetaNames)
  
  #Add the names of the models:
  Combined$ModelNames <- factor(rep(ModelNames, each = K),
                                levels = ModelNames)
  
  #Add panel identifier  
  Combined$PanelID <- rep(BetaNames, length(ModelNames))
  colnames(Combined) <- c("Mean", "Lo", "Up", "ModelNames", "PanelID")
  
  # Plot the graph
  p <- ggplot()
  p <- p + geom_point(data = Combined,
                      aes(x = ModelNames, 
                          y = Mean))
  p <- p + geom_errorbar(data = Combined,
                         aes(x = ModelNames, 
                             ymax = Up, 
                             ymin = Lo), 
                         width = 0.2)
  p <- p + xlab("Parameters") + ylab("Values")
  p <- p + theme(text = element_text(size=15)) 
  p <- p + facet_wrap( ~ PanelID, scales = "free_y")
  p <- p + theme(legend.position="none", axis.text.x = element_text(angle = 60, hjust = 1)) 
  print(p)
}



# Get the spatial parameters of an INLA model
MySpatialParams <- function(Model, ThisSpde) {
  SpFi <- inla.spde2.result(inla = Model, name = "w", spde = ThisSpde, do.transfer = TRUE) 
  Kappa  <- inla.emarginal(function(x) x, SpFi$marginals.kappa[[1]] )
  sigmau <- inla.emarginal(function(x) sqrt(x),SpFi$marginals.variance.nominal[[1]] )
  Range <- inla.emarginal(function(x) x, SpFi$marginals.range.nominal[[1]] )
  Out <- c(Kappa, sigmau, Range)
  names(Out) <- c("Kappa", "Sigma_u", "Range")
  Out
}


# Copied from Haakon Bakka
# Support function to plot a SRF
PlotField2 <- function(field, mesh, ContourMap, xlim, ylim, Add=FALSE, MyTitle = "",...){
  #stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = xlim, 
                              ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # And plot the whole thing
  image.plot(list(x = proj$x, 
                  y = proj$y,
                  z = field.proj), 
             xlim = xlim, 
             ylim = ylim,
             asp = 1,
             add = Add,
             main = MyTitle,
             ...)  
}


# This function is taken from:
# https://haakonbakka.bitbucket.io/btopic105.html
local.find.correlation = function(Q, location, mesh) {
  sd = sqrt(diag(inla.qinv(Q)))
  # - the marginal standard deviations
  A.tmp = inla.spde.make.A(mesh=mesh, loc = matrix(c(location[1], location[2]),1,2))
  # - create a fake A matrix, to extract the closest mesh node index
  id.node = which.max(A.tmp[1, ])
  # - index of the closest node
  print(paste('The location used was c(', 
              round(mesh$loc[id.node, 1], 4), ', ', 
              round(mesh$loc[id.node, 2], 4), ')' ))
  # - location of the closest node
  # - should be close to the location input
  # - sometimes used to plot a black dot
  ## Solve a matrix system to find the column of the covariance matrix
  Inode = rep(0, dim(Q)[1]); Inode[id.node] = 1
  covar.column = solve(Q, Inode)
  corr = drop(matrix(covar.column)) / (sd*sd[id.node])
  return(corr)
}



bri.hyperpar.summary <- function (r) 
{
  irp = r$internal.marginals.hyperpar
  hrp = r$marginals.hyperpar
  hypnames = names(irp)
  iip = grep("precision", hypnames)
  for (i in 1:length(irp)) {
    if (i %in% iip) {
      irp[[i]] = bri.hyper.sd(irp[[i]], internal = TRUE)
    }
    else {
      irp[[i]] = hrp[[i]]
      hypnames[i] = names(hrp)[i]
    }
  }
  ts = t(sapply(irp, bri.density.summary))
  hypnames = sub("Log precision", "SD", hypnames)
  row.names(ts) = hypnames
  ts
}


bri.hyper.sd <- function (prec, internal = FALSE) 
{
  if (internal) {
    inla.tmarginal(function(x) 1/sqrt(exp(x)), prec)
  }
  else {
    inla.tmarginal(function(x) 1/sqrt(x), prec)
  }
}


bri.density.summary <- function (dens) 
{
  m = inla.emarginal(function(xx) c(xx, xx^2), dens)
  q = inla.qmarginal(c(0.025, 0.5, 0.975), dens)
  s = sqrt(max(0, m[2] - m[1]^2))
  md = inla.mmarginal(dens)
  c(mean = m[1], sd = s, q0.025 = q[1], q0.5 = q[2], q0.975 = q[3], 
    mode = md)
}





#' This is a support file for making a PCA biplot
MyBiplot <- function(x, MySTD, MyAlpha) {
  ZZ <- x
  # get the names for the loadings from the data file
  VarNames <- colnames(ZZ)       
  
  # calculation of new matrix based upon choice 
  # between correlation ("std") and covariance ("cen") relationships
  if (MySTD == "std") {ZZnc <- decostand(ZZ, "standardize")} 
  if (MySTD == "cen") {ZZnc <- scale(ZZ, center = TRUE, scale = FALSE)}
  
  # compute the singular value decomposition of the rectangular matrix
  a <- svd(ZZnc)
  
  # calculation of information needed to amke a biplot based upon choice 
  # between correlation ("alpha=0") and distance ("alpha=1") biplot
  if (MyAlpha == 0) {
    G <- a$u %*% diag(a$d) 
    G12 <- G[,1:2]
    H <-   t(a$v) 
    H12 <- t(H[1:2,])
  } 
  
  if (MyAlpha == 1) {
    G <- a$u  
    G12 <- G[,1:2]
    H <-   diag(a$d) %*% t(a$v) 
    H12 <- t(H[1:2,])
  }
  
  # PLOTTING A BIPLOT
  # start to build a biplot by framing a square plot
  par(pty = "s")                       
  
  # calculate the maximum value tod etermine the scale of your axes
  MaxLoadings <- max(abs(H12))
  
  # create an empty box
  plot(0,0, type = "n", axes= F,
       xlim = c(-1.3 * MaxLoadings, 1.3 * MaxLoadings),
       ylim = c(-1.3 * MaxLoadings, 1.3 * MaxLoadings),
       xlab = "Axis 1",
       ylab = "Axis 2")                
  
  # determine the number of variables from the number of variable names
  N <- length(VarNames)
  # draw a line for each variable
  # with a particular color (col), line width (lwd), etc.
  # adding variable names at a set distance from the line
  # adding labels to the axes 1 and 2
  for (i in 1:N){
    p1 <- c(0, H12[i, 1])
    p2 <- c(0, H12[i, 2])
    lines(p1, p2, lwd = 2)            
    text(1.3 * H12[i,1],  1.3 * H12[i,2], 
         VarNames[i], cex = 1)
  }
  axis(1)                             
  axis(2)                             
  
  #Add scores for observations (sites)
  MaxScores <- max(abs(G12))          #Ensure that both axes have same range
  par(pty = "s", new = TRUE)          #The new = TRUE does the trick (superimposed this graph on previous graph)
  plot(G12[,1], G[,2], axes = FALSE, 
       xlab = "", ylab ="", 
       xlim = c(-MaxScores, MaxScores), 
       ylim = c(-MaxScores,MaxScores),
       pch = 16, cex = 0.5)          #Play with different col, cex, pch
  axis(3)                            #Add axis
  axis(4)                            #Add axis
  box()                              #Add box to make it look normal
}





#' Support file for plotting a SRF
PlotField3 <- function(field, mesh, ContourMap, 
                       xlim, ylim, Add=FALSE, MyTitle = "",
                       ColourScheme = 1, Xlabel = "", Ylabel = "", ...){
  #stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, xlim = xlim, ylim = ylim, 
                              dims = c(300, 300))
  # The function inla.mesh.project can then be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # Put all info together
  df <- expand.grid(x = proj$x, y = proj$y)
  df$w.pm <- as.vector(field.proj)
  
  w <- ggplot()
  w <- w + geom_raster(data = df, aes(x = x, y = y, fill = w.pm))
  w <- w + coord_fixed(ratio = 1) + theme_bw() + ggtitle(MyTitle)
  w <- w + xlab(Xlabel) + ylab(Ylabel)
  
  if (ColourScheme == 1){
    w <- w + scale_fill_gradient2(name = "",
                                  #limits = c(-1, 0.8),
                                  midpoint = 0, 
                                  low = "blue", 
                                  mid = "white", 
                                  high = "red",
                                  na.value = "transparent")}
  
  if (ColourScheme == 2){
    w <- w + scale_fill_gradient(name = "",
                                 low = "blue", 
                                 high = "orange",
                                 na.value = "transparent")}
  if (ColourScheme == 3){
    w <- w + scale_fill_viridis(name = "",
                                na.value = "transparent")}
  if (ColourScheme == 4){
    w <- w + scale_fill_gradientn(name = "",
                                  colours = rainbow(5),
                                  na.value = "transparent")}
  
  w
}






# ACF
# https://rh8liuqy.github.io/ACF_PACF_by_ggplot2.html
ggplot.corr <- function(data, lag.max = 24, ci = 0.95, large.sample.size = TRUE, horizontal = TRUE,...) {
  
  require(ggplot2)
  require(dplyr)
  require(cowplot)
  
  if(horizontal == TRUE) {numofrow <- 1} else {numofrow <- 2}
  
  list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
  N <- as.numeric(list.acf$n.used)
  df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
  df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
  df1$lag.acf[2] <- 0
  df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
  df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
  df1$acfstd[1] <- 0
  df1 <- select(df1, lag, acf, acfstd)
  
  list.pacf <- acf(data, lag.max = lag.max, type = "partial", plot = FALSE)
  df2 <- data.frame(lag = list.pacf$lag,pacf = list.pacf$acf)
  df2$pacfstd <- sqrt(1/N)
  
  if(large.sample.size == TRUE) {
    plot.acf <- ggplot(data = df1, aes( x = lag, y = acf)) +
      geom_area(aes(x = lag, y = qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*acfstd), fill = "#B9CFE7") +
      geom_col(fill = "#4373B6", width = 0.1) +
      scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
      scale_y_continuous(name = element_blank(), 
                         limits = c(min(df1$acf,df2$pacf),1)) +
      ggtitle("ACF") +
      theme_bw()
    
    plot.pacf <- ggplot(data = df2, aes(x = lag, y = pacf)) +
      geom_area(aes(x = lag, y = qnorm((1+ci)/2)*pacfstd), fill = "#B9CFE7") +
      geom_area(aes(x = lag, y = -qnorm((1+ci)/2)*pacfstd), fill = "#B9CFE7") +
      geom_col(fill = "#4373B6", width = 0.7) +
      scale_x_continuous(breaks = seq(0,max(df2$lag, na.rm = TRUE),6)) +
      scale_y_continuous(name = element_blank(),
                         limits = c(min(df1$acf,df2$pacf),1)) +
      ggtitle("PACF") +
      theme_bw()
  }
  else {
    plot.acf <- ggplot(data = df1, aes( x = lag, y = acf)) +
      geom_col(fill = "#4373B6", width = 0.7) +
      geom_hline(yintercept = qnorm((1+ci)/2)/sqrt(N), 
                 colour = "sandybrown",
                 linetype = "dashed") + 
      geom_hline(yintercept = - qnorm((1+ci)/2)/sqrt(N), 
                 colour = "sandybrown",
                 linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
      scale_y_continuous(name = element_blank(), 
                         limits = c(min(df1$acf,df2$pacf),1)) +
      ggtitle("ACF") +
      theme_bw()
    
    plot.pacf <- ggplot(data = df2, aes(x = lag, y = pacf)) +
      geom_col(fill = "#4373B6", width = 0.7) +
      geom_hline(yintercept = qnorm((1+ci)/2)/sqrt(N), 
                 colour = "sandybrown",
                 linetype = "dashed") + 
      geom_hline(yintercept = - qnorm((1+ci)/2)/sqrt(N), 
                 colour = "sandybrown",
                 linetype = "dashed") + 
      scale_x_continuous(breaks = seq(0,max(df2$lag, na.rm = TRUE),6)) +
      scale_y_continuous(name = element_blank(),
                         limits = c(min(df1$acf,df2$pacf),1)) +
      ggtitle("PACF") +
      theme_bw()
  }
  cowplot::plot_grid(plot.acf, plot.pacf, nrow = numofrow)
}



