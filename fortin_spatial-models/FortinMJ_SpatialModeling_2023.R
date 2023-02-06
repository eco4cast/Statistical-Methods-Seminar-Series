## Install librairies
#require(ggplot2)
require(MASS)
require(nlme)
require(lme4)
require(MuMIn) 
require(spdep)
require(spgwr)
require(vegan)

## Read Birdbird - Forest data
bird = read.csv("bird_forest.csv", header=TRUE)
str(bird)
plot(bird$xUTM,bird$yUTM,pch=20)

#ggplot(bird, aes(xUTM, yUTM, color=Bird)) + 
#  ggtitle("Overbird") + geom_point(show.legend=T)

## OSL - Linear Regression
plot(bird$Forest,bird$Bird)
Bird.lm=lm(Bird ~ Forest, data=bird)
summary(Bird.lm)

## Residual analysis
Fitted.Bird=fitted.values(Bird.lm)
Resid.Bird=residuals(Bird.lm)
plot(bird$Bird,Fitted.Bird)

## Linear regression - Trend Suface (x-y coordinates)
Bird.lm.xy=lm(Bird ~ Forest + xUTM + yUTM, data=bird)
summary(Bird.lm.xy)

## Generalized Linear Mixed Model
# Fit model including random effect
Bird.GLMMr = lmer(Bird ~ (1 | Zones), data=bird,REML=FALSE)
summary(Bird.GLMMr)
# Nakagawa & Schielzeth's: R2m=fixed effects, R2c=fixed and random effects
#r.squaredGLMM(Bird.GLMMr)

# Fit model including fixed and random effects
Bird.GLMM = lmer(Bird ~ Forest + (1 | Zones), data=bird,REML=FALSE)
summary(Bird.GLMM)
# Nakagawa & Schielzeth's: R2m=fixed effects, R2c=fixed and random effects
#r.squaredGLMM(Bird.GLMM)

## GLS
Bird.GLSx = gls(Bird ~ Forest,
   correlation = corAR1(form = ~ 1 | xUTM), data=bird)
summary(Bird.GLSx)

Bird.GLSy = gls(Bird ~ Forest,
   correlation = corAR1(form = ~ 1 | yUTM), data=bird)
summary(Bird.GLSy)
  
Bird.corLin <- gls(Bird ~ Forest, correlation = corLin(form = ~ xUTM + yUTM), data = bird)
summary(Bird.corLin)

Bird.corSpher <- gls(Bird ~ Forest, correlation = corSpher(form = ~ xUTM + yUTM, nugget = TRUE), data = bird)
summary(Bird.corSpher)

## make a listw object
xy=cbind(bird$xUTM,bird$yUTM)
# Delaunay Tessellation
xy.delo=tri2nb(xy)
xy.W=nb2listw(xy.delo, glist=NULL, style="W", zero.policy=NULL)
fdist <- lapply(xy.delo, function(x) 1 - x/max(dist(xy)))
xy.listw <- nb2listw(xy.delo, glist = fdist)

## Spatial Lag Regression
#Bird.lag=lagsarlm(Bird ~ Forest, data = bird, listw=xy.listw) 
#summary(Bird.lag)

## Spatial Error Regression
#Bird.err = errorsarlm(Bird ~ Forest, data = bird, listw = xy.listw)
#summary(Bird.err)

## Comparing models
#AIC(Bird.lm, Bird.lm.xy, Bird.GLMMr, Bird.GLMM, Bird.GLSx, Bird.GLSy, Bird.corLin, Bird.corSpher,Bird.lag,Bird.err)
AIC(Bird.lm, Bird.lm.xy, Bird.GLMMr, Bird.GLMM, Bird.GLSx, Bird.GLSy, Bird.corLin, Bird.corSpher)

## Spatial Filtering - dbMEMs / Multiple Scales
d.xy=dist(xy)
pcnm.xy=pcnm(d.xy)

Bird.dbMEMs=lm(Bird ~ Forest + scores(pcnm.xy), data=bird)
summary(dbBird.MEMs)

# Keep only 10 highly significant MEMs
pcnm.10=cbind(pcnm.xy$vectors[,1],pcnm.xy$vectors[,19],pcnm.xy$vectors[,31],
              pcnm.xy$vectors[,45],pcnm.xy$vectors[,94],pcnm.xy$vectors[,163],
              pcnm.xy$vectors[,170],pcnm.xy$vectors[,232],pcnm.xy$vectors[,300],
              pcnm.xy$vectors[,378])
Bird.dbMEM10=lm(Bird ~ Forest + pcnm.10, data=bird)
summary(dbBird.MEM10)

## GWR
BirdGWRbandwidth <- gwr.sel(Bird ~ Forest, data=bird, coords=cbind(bird$xUTM,bird$yUTM),adapt=T)
Bird.gwr = gwr(Bird ~ Forest, data=bird, coords=cbind(bird$xUTM,bird$yUTM),
               adapt=BirdGWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 
Bird.gwr
Bird.results<-as.data.frame(Bird.gwr$SDF)
head(Bird.results)
bird$grw.R2=Bird.results$localR2
summary(bird$grw.R2)
boxplot(bird$grw.R2)
