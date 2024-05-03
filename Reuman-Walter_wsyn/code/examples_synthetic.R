#This script implements a bunch of synthetic examples, mostly taken from the
#wsyn vignette, for use in the talk.
#
#Reuman

#***setup

resloc<-"../results/synthetic_results/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

#***wavelet example

set.seed(101)

#Set up an artificial time series that changes timescale of oscillation in the 
#middle, and has some white noise
time1<-1:100; time2<-101:200; times<-c(time1,time2)

ts1<-c(sin(2*pi*time1/15),0*time2) #period 15 initially
ts2<-c(0*time1,sin(2*pi*time2/8)) #then period 8
ts3<-rnorm(200,mean=0,sd=0.5) #add some white noise

ts<-ts1+ts2+ts3 

#Now apply the wavelet transform, obtaining an object of class `wt`. Default 
#parameter values for `scale.min`, `scale.max.input`, `sigma` and `f0` are usually 
#good enough forinitial data exploration.
ts<-wsyn::cleandat(ts,times,clev=1)$cdat
wtres<-wsyn::wt(ts,times)

#now plot the magnitude and phase and verify that they make sense given how the 
#data were constructed
grDevices::pdf(file=paste0(resloc,"WaveletExample_magnitude.pdf"))
wsyn::plotmag(wtres)
grDevices::dev.off()

grDevices::pdf(file=paste0(resloc,"WaveletExample_phase.pdf"))
wsyn::plotphase(wtres)
grDevices::dev.off()

#***wavelet phasor mean field example

#Set up multi-location data that all has three components: 1) a sin wave of period 
#10 for the first half of the data and then 5 for the second half, present in all\
#locations; 2) a sin wave of period 3, phase randomized independently in all 
#locations; 3) white noise. 
times1<-0:50; times2<-51:100; times<-c(times1,times2)
ts1<-c(sin(2*pi*times1/10),sin(2*pi*times2/5))+1.1 

dat<-matrix(NA,11,length(times))
for (counter in 1:dim(dat)[1])
{
  ts2<-3*sin(2*pi*times/3+2*pi*runif(1))+3.1
  ts3<-rnorm(length(times),0,1.5)
  dat[counter,]<-ts1+ts2+ts3    
}
dat<-wsyn::cleandat(dat,times,1)$cdat

#synchrony is not visually obvious
grDevices::pdf(file=paste0(resloc,"WPMFExample_timeseries.pdf"))
plot(times,dat[1,]/10+1,type='l',xlab="Time",ylab="Time series index",ylim=c(0,12))
for (counter in 2:dim(dat)[1])
{
  lines(times,dat[counter,]/10+counter)
}
grDevices::dev.off()

#Nor can synchrony be readily detected by examining the $55$ pairwise correlation coefficients
cmat<-cor(t(dat))
diag(cmat)<-NA
cmat<-as.vector(cmat)
cmat<-cmat[!is.na(cmat)]
grDevices::pdf(file=paste0(resloc,"WPMFExample_PairwiseCorrelation.pdf"))
hist(cmat,30,xlab="Pearson correlation",ylab="Count")
grDevices::dev.off()

#so now apply the wpmf
wpmfres<-wsyn::wpmf(dat,times,sigmethod="quick")
grDevices::pdf(file=paste0(resloc,"WPMFExample_wpmf.pdf"))
wsyn::plotmag(wpmfres,sigthresh=0.95)
grDevices::dev.off()
#you can see what we know is there because we built it in, so this demonstrates the utility of the method

#***wavelet coherence

#make some environmental data, 2 superimposed sin waves in all locations and white noise
times<-(-3:100); lt<-length(times)
ts1<-sin(2*pi*times/10)
ts2<-5*sin(2*pi*times/3)
x<-matrix(ts1+ts2,11,lt,byrow=TRUE)+
  matrix(rnorm(11*lt,0,1.5),11,lt)

#make biological data - the population is a moving average of the environment plus white noise
times<-0:100; lt<-length(times)
y<-matrix(NA,11,lt) #the driven (biological) variable
for (i in 1:101) {
  y[,i]<-apply(FUN=mean,X=x[,i:(i+2)],MARGIN=1) }
y<-y+matrix(rnorm(11*lt,mean=0,sd=3),11,lt)

x<-x[,4:104]
x<-wsyn::cleandat(x,times,1)$cdat 
y<-wsyn::cleandat(y,times,1)$cdat

#The relationship between the environmental and biological variables cannot 
#readily be detected using ordinary correlation methods 
allcors<-c()
for (counter in 1:dim(x)[1])
{
  allcors[counter]<-cor(x[counter,],y[counter,])
}
grDevices::pdf(file=paste0(resloc,"CoherenceExample_LocalBiolEnvCorrs.pdf"))
hist(allcors,xlab="Biol/env correlation",ylab="Count")
grDevices::dev.off()

#However, the function `coh` can be used to compute the coherence, which reveals a 
#timescale-specific relationship
if (file.exists(paste0(resloc,"CoherenceResults.Rds")))
{
  cohres<-readRDS(file=paste0(resloc,"CoherenceResults.Rds"))
} else
{
  cohres<-wsyn::coh(dat1=x,dat2=y,times=times,norm="powall",
         sigmethod="fftsurrog1",nrand=1000,
         f0=0.5,scale.max.input=28)
  saveRDS(cohres,file=paste0(resloc,"CoherenceResults.Rds"))
}

#now display what you get
grDevices::pdf(file=paste0(resloc,"CoherenceExample_PlotCoh1.pdf"))
wsyn::plotmag(cohres)
grDevices::dev.off()

grDevices::pdf(file=paste0(resloc,"CoherenceExample_PlotRanks.pdf"))
wsyn::plotrank(cohres)
grDevices::dev.off()

#now add p-values for a couple bands
cohres<-wsyn::bandtest(cohres,c(8,12))
cohres<-wsyn::bandtest(cohres,c(2,4))

#now display what you get again
grDevices::pdf(file=paste0(resloc,"CoherenceExample_PlotCoh2.pdf"))
wsyn::plotmag(cohres)
grDevices::dev.off()

#***wavelet linear model example

set.seed(3221) 

#First create a diver variable composed of an oscillation of period $12$ years and
#an oscillation of period $3$ years, and normally
#distributed white noise of mean $0$ and standard deviation $1.5$.
lts<-12
sts<-3
mats<-3
times<-seq(from=-mats,to=100)
ts1<-sin(2*pi*times/lts)
ts2<-sin(2*pi*times/sts)
numlocs<-10
d1<-matrix(NA,numlocs,length(times)) #the first driver
for (counter in 1:numlocs)
{
  d1[counter,]<-ts1+ts2+rnorm(length(times),mean=0,sd=1.5)
}

#Next create a second driver, again composed of an oscillation of period $12$ years and
#an oscillation of period $3$ years, and normally
#distributed white noise of mean $0$ and standard deviation $1.5$.
ts1<-sin(2*pi*times/lts)
ts2<-sin(2*pi*times/sts)
d2<-matrix(NA,numlocs,length(times)) #the second driver
for (counter in 1:numlocs)
{
  d2[counter,]<-ts1+ts2+rnorm(length(times),mean=0,sd=1.5)
}

#Next create an irrelevant environmental variable. With real data, of course,
#one will not necessarily know in advance whether an environmental
#variable is irrelevant to a population system. But, for the purpose of 
#demonstrating the methods, we are playing the dual role of data creator and 
#analyst.
dirrel<-matrix(NA,numlocs,length(times)) #the irrelevant env var
for (counter in 1:numlocs)
{
  dirrel[counter,]<-rnorm(length(times),mean=0,sd=1.5)
}

#The population in each location is a combination of the two drivers,
#plus local variability. Driver 1 is averaged over 3 time steps
#in its influence on the populations, so
#only the period-12 variability in driver 1 influences the populations. 
pops<-matrix(NA,numlocs,length(times)) #the populations
for (counter in (mats+1):length(times))
{
  aff1<-apply(FUN=mean,X=d1[,(counter-mats):(counter-1)],MARGIN=1)
  aff2<-d2[,counter-1]
  pops[,counter]<-aff1+aff2+rnorm(numlocs,mean=0,sd=3)
}
pops<-pops[,times>=0]
d1<-d1[,times>=0]
d2<-d2[,times>=0]
dirrel<-dirrel[,times>=0]
times<-times[times>=0]

#If only the data were available and we were unaware of how they were generated,
#we may want to infer the causes of synchrony and its timescale-specific patterns 
#in the populations.   
#Start by fitting a model with all three predictors. 
dat<-list(pops=pops,d1=d1,d2=d2,dirrel=dirrel)
dat<-lapply(FUN=function(x){wsyn::cleandat(x,times,1)$cdat},X=dat) #prep data for call to wlm
wlm_all<-wsyn::wlm(dat,times,resp=1,pred=2:4,norm="powall",scale.max.input=28)

#We will carry out analyses for this model at long timescales 
#($11$ to $13$ years) and short timescales 
#($2$ to $4$ years) simultaneously. First test whether we can drop each variable.
wlm_all_dropi<-wsyn::wlmtest(wlm_all,drop="dirrel",sigmethod="fft",nrand=100) #few randomizations for speed, but in real life use at least 1000
wlm_all_drop1<-wsyn::wlmtest(wlm_all,drop="d1",sigmethod="fft",nrand=100)
wlm_all_drop2<-wsyn::wlmtest(wlm_all,drop="d2",sigmethod="fft",nrand=100)

#the two timescale bands we will use
blong<-c(11,13)
bshort<-c(2,4)

#tests for whether we can drop dirrel
wlm_all_dropi<-wsyn::bandtest(wlm_all_dropi,band=blong)
wlm_all_dropi<-wsyn::bandtest(wlm_all_dropi,band=bshort)

#display/save results
grDevices::pdf(file=paste0(resloc,"WLM_example_DropIrrel.pdf"))
wsyn::plotrank(wlm_all_dropi)
grDevices::dev.off()

#tests for whether we can drop d1
wlm_all_drop1<-wsyn::bandtest(wlm_all_drop1,band=blong)
wlm_all_drop1<-wsyn::bandtest(wlm_all_drop1,band=bshort)

#display/save results
grDevices::pdf(file=paste0(resloc,"WLM_example_Drop1.pdf"))
wsyn::plotrank(wlm_all_drop1)
grDevices::dev.off()

#tests for whether we can drop d2
wlm_all_drop2<-wsyn::bandtest(wlm_all_drop2,band=blong)
wlm_all_drop2<-wsyn::bandtest(wlm_all_drop2,band=bshort)

#display/save results
grDevices::pdf(file=paste0(resloc,"WLM_example_Drop2.pdf"))
wsyn::plotrank(wlm_all_drop2)
grDevices::dev.off()

#it turns out we can drop dirrel on both bands, so do it, i.e., make a new model
#where it's dropped
wlm_good<-wsyn::wlm(dat,times,resp=1,pred=2:3,norm="powall",scale.max.input=28)

#now look at how much synchrony is explained, short timescales
se<-wsyn::syncexpl(wlm_good)
se_short<-se[se$timescales>=bshort[1] & se$timescales<=bshort[2],]
round(100*colMeans(se_short[,c(3:8)])/mean(se_short$sync),4)
#almost all the synchrony that can be explained is explained by d1

#now look at how much synchrony is explained, long timescales
se_long<-se[se$timescales>=blong[1] & se$timescales<=blong[2],]
round(100*colMeans(se_long[,c(3:8)])/mean(se_long$sync),4)
#synchrony explained by d1, d2 and interactions



