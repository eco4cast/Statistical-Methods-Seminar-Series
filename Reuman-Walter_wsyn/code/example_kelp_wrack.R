## This script implements a worked empirical example based on a subset of the 
## results in Walter et al. (2024) Spatial synchrony cascades across ecosystem 
## boundaries and up food webs via resource subsidies. PNAS 121(2) e23152120.

## The data provided to accompany this example have been pre-processed to 
## (where appropriate) spatially aggregate, organize, and gap-fill.
## Original data sources are publicly accessible:

## Santa Barbara Coastal LTER, J. Dugan, SBC LTER: Beach: 
## Time-series of beach wrack cover and biomass, ongoing since 2008 ver 16. 
## Environmental Data Initiative (2021). 
## https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.40.16.

## T. Bell, SBC LTER: Reef: California kelp canopy and environmental variable dynamics ver 1. 
## Environmental Data Initiative (2023). https://doi.org/10.6073/pasta/c40db2c8629cfa3fbe80fdc9e086a9aa. 


## Examine coherence between kelp forest production and kelp wrack deposition

rm(list=ls())
graphics.off()

library(wsyn)


#This sets the working directory to the ./code subdirectory but only works if using RStudio
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 


resloc<-"../results/kelp_results/"
if (!dir.exists(resloc))
{
  dir.create(resloc,recursive=TRUE)
}

## load datasets ----------------------------------------------------------------------------------

#these files have already been organized into locations-by-time matrices and been gap-filled by
#replacing NAs with month-specific medians
wrack <- as.matrix(read.csv("../data/wrack.csv"))
beachwidth <- as.matrix(read.csv("../data/beachwidth.csv"))
wavediff <- as.matrix(read.csv("../data/wavediff.csv"))
kelp <- as.matrix(read.csv("../data/kelpproduction.csv"))

#replicate replicate the region-wide kelp time series 5x to match site dimensions
kelp <- rbind(c(kelp),c(kelp),c(kelp),c(kelp),c(kelp))

#this is the time (year and month) information corresponding to columns in the data matrices
yearmonth <- read.csv("../data/time.csv")
tt=1:132 #vector of timesteps


## apply cleandat() to normalize, detrend, and scale ----------------------------------------------

beachwidth.cln <- wsyn::cleandat(beachwidth, tt, clev=5)$cdat
wavediff.cln <- wsyn::cleandat(wavediff, tt, clev=5)$cdat
kelp.cln <- wsyn::cleandat(kelp, tt, clev=5)$cdat
wrack.cln <- wsyn::cleandat(wrack, tt, clev=5)$cdat


## Create wavelet mean fields ---------------------------------------------------------------------

wmf.wrack <- wsyn::wmf(wrack.cln, tt)
wpmf.wrack <- wsyn::wpmf(wrack.cln, tt, sigmethod="quick")

pdf(paste0(resloc,"wrackWMF.pdf"))
wsyn::plotmag(wmf.wrack)
mtext("Wavelet mean field")
dev.off()

pdf(paste0(resloc,"wrackWPMF.pdf"))
wsyn::plotmag(wpmf.wrack)
mtext("Wavelet phasor mean field")
dev.off()

#Just adding timescale bands to these figs for pedagogical reasons
pdf(paste0(resloc,"wrackWMF2.pdf"))
wsyn::plotmag(wmf.wrack)
abline(h=log2(c(8,16)), lwd=3)
mtext("Wavelet mean field")
dev.off()

pdf(paste0(resloc,"wrackWPMF2.pdf"))
wsyn::plotmag(wpmf.wrack)
mtext("Wavelet phasor mean field")
abline(h=log2(c(8,16)), lwd=3)
dev.off()



## Wavelet spatial coherence analyses -------------------------------------------------------------

#set timescale ranges
subann <- c(0,8) #months
annual <- c(8,16) #months
intann <- c(16,60) #months

#nrand is relatively small to speed computation. Often 10,000 is used for publications with the fast method
wrackXkelp <- wsyn::coh(wrack.cln, kelp.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXkelp <- wsyn::bandtest(wrackXkelp, subann)
wrackXkelp <- wsyn::bandtest(wrackXkelp, annual)
wrackXkelp <- wsyn::bandtest(wrackXkelp, intann)

pdf(paste0(resloc,"wrackXkelp_cohmag.pdf"))
wsyn::plotmag(wrackXkelp) #coherence on subannual, annual and interannual
dev.off()

pdf(paste0(resloc,"wrackXkelp_cohphase.pdf"))
wsyn::plotphase(wrackXkelp) #phase lagged on annual and interannual 
dev.off()

wrackXwidth <- wsyn::coh(wrack.cln, beachwidth.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXwidth <- wsyn::bandtest(wrackXwidth, subann)
wrackXwidth <- wsyn::bandtest(wrackXwidth, annual)
wrackXwidth <- wsyn::bandtest(wrackXwidth, intann)
wsyn::plotmag(wrackXwidth) #coherence on seasonal, annual, interannual
wsyn::plotphase(wrackXwidth) #in phase on annual, borderline in-phase or lagged on long timescales

wrackXwavediff <- wsyn::coh(wrack.cln, wavediff.cln, tt, sigmethod="fast", norm="powall", nrand=1000)
wrackXwavediff <- wsyn::bandtest(wrackXwavediff, subann)
wrackXwavediff <- wsyn::bandtest(wrackXwavediff, annual)
wrackXwavediff <- wsyn::bandtest(wrackXwavediff, intann)
wsyn::plotmag(wrackXwavediff) #coherence on seasonal, annual, interannual
wsyn::plotphase(wrackXwavediff) #in phase on annual, borderline in-phase or lagged on long timescales



## Wavelet linear model of wrack ------------------------------------------------------------------

#makes a wavelet linear model based on the set of putative predictors that were coherent with
#wrack on the annual timescale band

datlist <- list(wrack.cln, kelp.cln, wavediff.cln, beachwidth.cln)
wlm.annual <- wsyn::wlm(datlist, tt, resp=1, pred=c(2,3,4), norm="powall", f0=1)
#get whole-model p-value by 'dropping' all predictors
#note: can drop a single predictor to get a p-value for that predictor
wlmtest.annual <- wsyn::wlmtest(wlm.annual, drop=2:4, sigmethod = "fft", nrand=100)
wlmtest.annual <- wsyn::bandtest(wlmtest.annual, annual)
print(wlmtest.annual$bandp)

#output and plot predicted synchrony
predsync.annual <- wsyn::predsync(wlm.annual)

pdf(paste0(resloc,"wrackWLM_pred.pdf"))
wsyn::plotmag(predsync.annual, zlims=c(0,max(Mod(wmf.wrack$values), na.rm=TRUE)))
abline(h=log2(annual), lwd=3)
mtext("WLM predicted synchrony")
dev.off()

#look at synchrony explained across all timescales
syncexpl.annual <- wsyn::syncexpl(wlm.annual)
print(round(syncexpl.annual,3))

#subset to wavelet components in annual timescale band
syncexpl.annual <- syncexpl.annual[syncexpl.annual$timescales > min(annual) 
                                   & syncexpl.annual$timescales < max(annual),]
print(round(syncexpl.annual,3))

#average across whole annual timescale band
print(colMeans(syncexpl.annual))

#overall synchrony explained in the annual timescale band
print(mean(syncexpl.annual$syncexpl)/mean(syncexpl.annual$sync))

#overall 'crossterms' across the annual timescale band
mean(syncexpl.annual$crossterms)/mean(syncexpl.annual$syncexpl)

#crossterms are a measure of spatial independence that evaluates 
#an assumption of the synchrony explained partitioning scheme.
#We have considered values <10% of the synchrony explained acceptable.


