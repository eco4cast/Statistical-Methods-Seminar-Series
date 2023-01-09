#    Highland Statistics Ltd.
#    www.highstat.com
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.




# Section 1: Data description----

#* Subsection 1.1: Paper----
#'  The data used in this exercise are taken from:
#'   Roulin A, Bersier LF (2007). Nestling barn owls beg more 
#'   intensely in the presence of their mother than their father. 
#'   Animal Behaviour 74: 1099–1106.

#' We also used these data in:
#'   - Zuur et al. (2009). Mixed Effects Models and Extensions in Ecology 
#'     with R (2009)

#'   - Zuur and Ieno (2018). Book: Beginner's Guide to Spatial, Temporal 
#'     and Spatial-Temporal Ecological Data Analysis with R-INLA. 
#'     Volume II: GAM and Zero-Inflated Models.
   




#* Subsection 1.2: Background biology----

#' Roulin and Bersier (2007) gathered data from 27 owl nests.
#'   - Using microphones and a video recorder, they investigated vocal 
#'     behaviour of nestlings.

#' Sibling negotiation:
#'    -Is the number of calls made by all offspring, in the absence of 
#'     the parents, during 30-sec time periods recorded at 15-min intervals.


#' Biological background:
#'  - By vocalising in the absence of parents, siblings inform one another 
#'    of their relative hunger level. 
#'  - Vocal individuals are hungrier than their siblings. 
#'  - By vocalising loudly, a hungry owlet informs its siblings that it is 
#'    willing to compete to monopolise the next delivered food item. 
#'  - Being informed, its siblings temporarily retreat from the contest, 
#'    thereby reducing the level of competition. 

#' Sibling competition allows siblings to optimally adjust effort in 
#' competition for parent-provided resources.




#* Subsection 1.3: Setup of the field study----

#' 1. Food treatment: 
#'    - Satiated versus deprived.
#'    - Half of the nests were given extra prey in the morning 
#'      preceding recording.
#'    - Prey remains were removed from the other half of the nests. 
#'      This was called ‘food-satiated’ and ‘food-deprived’.
#'    - Sampling was conducted on two nights.
#'      Food treatment was reversed on the second night.


#' 2. Arrival Time:
#'    Recordings were made between 21.30 h and 05.30 h. 
#'    Arrival time reflects the time when a parent arrived at the perch 
#'    with prey.

#' 3. Sex of the parent.

#' 4. Brood size.




#* Subsection 1.4: Underlying questions----

#' Is sibling negotiation related to:
#'   - Sex of the parent.
#'   - Food treatment 
#'   - Arrival time of the parent. 

#' Broodsize can be used as an offset or as a covariate (either
#' as a linear term or as a smoother).


#' To simplify the models and graphs in this presentation, we will refrain 
#' from using interactions, and we will also drop broodsize. This is the
#' wrong thing to do..but it makes this 1-hour presentation easier.




#  Section 2: Import the data and load the packages----

setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/Data")
Owls <- read.table("Owls2.txt", 
                   header = TRUE,
                   na.strings = "NA",
                   stringsAsFactors = TRUE,
                   dec = ".")

#' Check the imported data.
str(Owls)
names(Owls)



#' We load the following packages.
library(lattice)
library(ggplot2)
library(sp)
library(ggmap)
library(MASS)
library(rgl)
library(plyr)
library(DHARMa)
library(performance)
library(glmmTMB)
library(mgcv)
library(gratia)
source("HighstatLibV13.R") #' Our support file





# Section 3: Prepare the data----


#' Rename SiblingNegotiation to a shorter name:
Owls$NCalls <- Owls$SiblingNegotiation      





#' Get the spatial coordinates. Convert from Swiss coordinates to WGS84.

#' Define coordinates as Swiss coordinates:
Owls.spdf <- Owls 
coordinates(Owls.spdf) <- ~Xcoord + Ycoord
projSWISS <- CRS("+init=epsg:21781")
proj4string(Owls.spdf) <- projSWISS

#' Convert the Swiss coordinates to wgs84.
WGS84 <- CRS("+proj=longlat +datum=WGS84")
Owls_wgs84.spdf <- spTransform(Owls.spdf, WGS84)

#' And add the wgs84 coordinates to the Owls object.
Owls$Longitude <- Owls_wgs84.spdf$Xcoord
Owls$Latitude  <- Owls_wgs84.spdf$Ycoord

#' Also get UTM coordinates.
LL <- LongLatToUTM(x = Owls$Longitude, 
                   y = Owls$Latitude, 
                   Hemisphere = "north", 
                   zone = 32)
Owls$Xkm <- LL[,"X"] / 1000
Owls$Ykm <- LL[,"Y"] / 1000




#' Put the sampling dates in chronological order.
Owls$Date2 <- factor(Owls$Date, 
                     levels = c(
                       "28/05/97",  "29/05/97",
                       "16/06/97", "17/06/97", "18/06/97", "19/06/97", "23/06/97", "24/06/97", 
                       "25/06/97", "26/06/97", "27/06/97", "28/06/97", "29/06/97", "30/06/97", 
                       "01/07/97", "02/07/97", "03/07/97", "04/07/97", "12/07/97", "13/07/97", 
                       "14/07/97", "15/07/97", "21/07/97", "22/07/97", "23/07/97", "24/07/97", 
                       "25/07/97", "26/07/97", "27/07/97", "09/08/97", "10/08/97", "11/08/97", 
                       "17/08/97",  "18/08/97"))







# Section 4: Data exploration----


#* Subsection 4.1: Spatial locations of the nests----

#' Plot the spatial locations on map      
glgmap   <- get_map(location = c(6.7, 46.7, 7.2, 47),
                    maptype= "terrain",
                    col = "color")    
p <- ggmap(glgmap)
p <- p + geom_point(aes(x = Longitude, y = Latitude),
                    data = Owls,
                    size = 4) 
p <- p + xlab("Longitude") + ylab("Latitude")  
p <- p + theme(text = element_text(size = 15)) 
p
#' How far does an owl fly?


#' If you don't have online access, then use:
xyplot(Latitude ~ Longitude,
       data = Owls,
       aspect = "iso")





#* Subsection 4.2: Outliers----

#' Make Cleveland dotplots.
MyVar <- c("ArrivalTime", "BroodSize",  
           "Longitude", "Latitude", "NCalls")
Mydotplot(Owls[,MyVar])

#' Comments:
#'  -The response variable (NCalls) has plenty of zeros.
#'  -The covariates look ok (nothing odd).
#'  -Broodsize has limited unique values. Tricky to use as a smoother?




#* Subsection 4.3: Collinearity----

#' We will skip this part during the ESA presentation.
#' In short, there is no strong collinearity.





#* Subsection 4.4: Relationships----

#' Plot the response variable NCalls vs. ArrivalTime.
p1 <- ggplot(Owls, 
            aes(x = ArrivalTime, y = NCalls)) +
      xlab("Arrival time") + 
      ylab("Number of calls") +
      geom_point() + 
      geom_smooth()
p1
#' Is this a non-linear relationship? That would make biological sense.







#* Subsection 4.5: Dependency----


#' Visualise the nest effect.
p2 <- ggplot(data = Owls, 
            aes(x = Nest, 
                y = NCalls)) +
      xlab("Nest") + ylab("Number of calls") +
      geom_boxplot() + 
      theme(text = element_text(size=12), 
            legend.position="none",
            axis.text.x = element_text(size = 5, angle =40,hjust = 0.5))
p2
#' There seems to be a nest effect.


#' How many nests and how many observations per nest do we have?
length(unique(Owls$Nest)) #Number of nests
table(Owls$Nest)          #Number of observations per nest
#' That is all ok-ish.




#' Plot the data as time series.
p <- ggplot(data = Owls,
            aes(x = ArrivalTime, 
                y = NCalls,
                group = FoodTreatment,
                col = FoodTreatment)) + 
     geom_line() + 
     xlab("Arrival time (hours)") + 
     ylab("Number of calls") + 
     facet_wrap(~ Nest, ncol = 8)
p
#' We have two short time series for each nest, and some of these
#' nests are close to each other!





#* Subsection 4.6: Zero inflation----

#' Do we have zero inflation?
100 * sum(Owls$NCalls == 0, na.rm = TRUE) / nrow(Owls)
plot(table(Owls$NCalls), type = "h", ylab = "Frequency")
# Potential trouble!





#* Subsection 4.7: Conclusions data exploration----

#'  - Plenty of zeros in NCalls.
#'  - There is dependency inside a nest (we need something with random effects?).
#'  - Is there a non-linear arrival time effect?
#'  - Is there spatial correlation between the nests?
#'  - Small number of nests.
#'  - No clear collinearity.






# Section 5: Model formulation----

#' We have count data:
#'   - A Poisson GLMM seems to be a decent starting point.
#'   - One can also argue that a Poisson GAMM (with a smoothing function of 
#'     arrival time) should be the starting point.


#' We will apply the following Poisson GLMM.
#'  NCalls_ij ~ Poisson(mu_ij)
#'  E(NCalls_ij)   = mu_ij
#'  var(NCalls_ij) = mu_ij

#'  log(mu_ij) = Intercept + SexParent + FoodTreatment + ArrivalTime + a_i

#'  a_i is the random intercept for nests.
#'  Imposes dependency between observations from the same nest.
#'  a_i ~ N(0, sigma^2_Nst)    i = 1, ..., 27




# Section 6: Poisson GLMM----



#* Subsection 6.1: Execute the Poisson GLMM----

#' To avoid numerical problems it is wise to standardize the continuous 
#' covariates.
Owls$ArrivalTimec <- MyStd(Owls$ArrivalTime)

#' Use glmmTMB to execute the Poisson GLMM.
M1 <- glmmTMB(NCalls ~ 1 + SexParent + FoodTreatment + ArrivalTimec + (1 | Nest),
              data = Owls,
              family = poisson)

summary(M1)
drop1(M1, test = "Chi")
performance(M1)
#' Looks all ok.


#' This is the fitted model.

#' For SexParent = Females and FoodTreatment = Deprived:
#' log(mu_ij) = 1.92 - 0.24 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Deprived:
#' log(mu_ij) = 1.92 + 0.04 - 0.24 * ArrivalTimec + Nest effect 

#' For SexParent = Females and FoodTreatment = Satiated:
#' log(mu_ij) = 1.92 -0.58 - 0.24 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Satiated:
#' log(mu_ij) = 1.92 + 0.04 -0.58 - 0.24 * ArrivalTimec + Nest effect 





#* Subsection 6.2: Model validation Poisson GLMM----

#' Do we have overdispersion?
testDispersion(M1)
#' Sort of.


#' Can this model cope with the 25% of zeros?
testZeroInflation(M1)
#' No.


#' Can this model cope with the non-linear arrival time effect?
#' Get scaled quantile residuals (using DHARMa) and plot these against 
#' arrival time.

#' Get scaled quantile residuals.
E1Eqr <- simulateResiduals(fittedModel = M1, plot = FALSE)

plotResiduals(E1Eqr, 
              form = Owls$ArrivalTime, 
              xlab = "Rank-transformed arrival time") 
#' Nice non-linear patterns. We can see the same with Pearson residuals.

Owls$E1.pears <- resid(M1, type = "pearson")
p3 <- ggplot(Owls, 
             aes(x = ArrivalTime, y = E1.pears)) +
        xlab("Arrival time") + 
        ylab("Pearson residuaals") +
        geom_point() + 
        geom_smooth()
p3
#' There is indeed a non-linear residual pattern!


#' Homework: 
#'   -Plot the scaled quantile residuals versus the fitted values.
#'   -Plot the scaled quantile residuals versus each covariate in the model.
#'   -Plot the scaled quantile residuals versus each covariate not in the model.
#'   -Check the scaled quantile residuals for spatial patterns.
#'   -Check the random effects for spatial patterns (spoiler alert: these
#'    are spatially correlated).




#* Subsection 6.3: Fitted values of the Poisson GLMM----

#' We will plot the fitted values using the following steps:
#'  A. Specify covariate values for predictions.
#'  B. Calculate the predicted values.
#'  C. Get a 95% CI for the predicted values.
#'  D. Plot the whole thing.


#' A. Define a grid of covariate values without extrapolation.
MyData1 <- ddply(Owls, 
                .(SexParent, FoodTreatment), 
                summarize,
                ArrivalTimec = seq(min(ArrivalTimec), 
                                   max(ArrivalTimec), 
                                   length = 50))
head(MyData1)
MyData1$Nest <- NA #' Ignore the random effects in the predictions (also 
                  #' something one can argue about).


#' B. Get the corresponding predicted values for a typical nest.
P1 <- predict(M1, 
              newdata = MyData1, 
              re.form =~ 0, 
              se = TRUE)


#' C. Calculate the predicted values and 95% CI.
MyData1$mu   <- exp(P1$fit)
MyData1$SeUp <- exp(P1$fit + 1.96 * P1$se.fit) 
MyData1$SeLo <- exp(P1$fit - 1.96 * P1$se.fit) 
head(MyData1)



#' D. And plot the whole thing

#' First back-standardize arival time.
MyData1$ArrivalTime <- MyData1$ArrivalTimec * sd(Owls$ArrivalTime) + mean(Owls$ArrivalTime)

p <- ggplot() + 
       geom_point(data = Owls, 
                  aes(x = ArrivalTime, 
                      y = NCalls)) +
       geom_line(data = MyData1, 
                 aes(x = ArrivalTime, 
                     y = mu)) + 
       geom_ribbon(data = MyData1,
                   aes(x = ArrivalTime, 
                       ymax = SeUp, 
                       ymin = SeLo), 
                   alpha = 0.5) + 
       xlab("ArrivalTime") + ylab("Number of calls") + 
       theme(text = element_text(size = 15), 
             legend.position="none") + 
       facet_grid(SexParent ~ FoodTreatment, 
                  scales = "fixed")
p







#* Subsection 6.4: What is the Poisson distribution doing?----


#' What is the Poisson distribution doing?
#' We are going to make a 3-d scatterplot. The code is not so nice on 
#' the eyes. Don't try to fully understand it. Just run it.
#' It can probably be programmed way more efficient.
#' MAC USERS: You may need to upgrade/install Quartz.



#' Focus on a specific SexParent/FoodTreatment combination.
ThisSexParent = "Female"
ThisFoodTr = "Deprived"

#' Select the data for the selected SexParent and FoodTreatment levels.
MyData.sub <- subset(MyData1, 
                     SexParent == ThisSexParent & 
                     FoodTreatment == ThisFoodTr)

Owls.sub <- subset(Owls, 
                   SexParent == ThisSexParent & 
                   FoodTreatment == ThisFoodTr)


#' Three variable for 3-d scatterplot
x      <- MyData.sub$ArrivalTime  #' Covariate
y.pois <- MyData.sub$mu           #' Fitted values
z      <- 0 * x                   #' A vector with 0s.

#' Number of Poisson density curves to plot.
K <- 6

#' We will pock the following rows from MyData.sub:
TheseRows <- round(seq(1, nrow(MyData.sub), length= K))
TheseRows


#' Plot the results.
#' First we plot the observed data for the selected SexParent/FoodTreatment
#' levels.
plot3d(x = Owls.sub$ArrivalTime,
       y = Owls.sub$NCalls,
       z = rep(0, nrow(Owls.sub)),
       type = "s",
       size = 0.8,
       lit = FALSE,
       zlim = c(0, 0.4),
       xlab = "Arrival time",
       ylab = "NCalls",
       zlab = "Density functions")

#'You can spin this graph! Don't close it!

#' Add the fitted values of the Poisson GLMM.
lines3d(x, y.pois, z, 
        col = 2,
        lwd = 3)

#' Add the Poisson density curves. 
interleave <- function(v1, v2) as.vector(rbind(v1, v2))
for (i in TheseRows){
  yseq <- seq(0, 26, by = 1)
  
  #' Calculate Poisson density values.
  dens.pois <- dpois(yseq, lambda = MyData.sub$mu[i])
  
  #' Dotted line 
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, 0)
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 1,
          lty = 2)
  
  #' Plot the Poisson density values.
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, dens.pois)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = "purple")
}

#' Now compare the 3d graph with the relevant ggplot2 panel.
#' Questions: 
#'     -What is the role of the Poisson distribution?
#'     -How can you improve the model?
#'     -What do you do after a meal?       





#* Subsection 6.5: Now what?----

#' Results presented above showed that the Poisson GLMM:
#'   - Is slightly overdispersed.
#'   - It cannot cope with the 26% of zeros.
#'   - Non-linear residual patterns with arrival time.

#' And we haven't even checked for spatial dependency yet. We should have done
#' that!


#' What shall we do next? Use arrival time as a smoother? That means a
#' GAMM. Makes biological sense!


#' I sensible plan would be:
#'   - Try a GAMM with s(ArrivalTime).
#'   - If still trouble, move on towards zero-inflated models.
#'   - If there is spatial trouble, then move on to something
#'     like a zero-inflated spatial GAM (that is fun!).





#* Subsection 6.6: The wrong approach----

#' The wrong approach would be to apply a series of models with different
#' distributions, and compare them all.

#' But that is exactly what we will do! 

#' Why? Because this is a teaching exercise.





# Section 7: Wrong approach 1: NB GLMM----


#* Subsection 7.1: Model formulation----

#' We will apply the following NB GLMM.
#'  NCalls_ij ~ NB(mu_ij, theta)
#'  E(NCalls_ij)   = mu_ij
#'  var(NCalls_ij) = mu_ij + mu_ij ^2 / theta

#'  log(mu_ij) = Intercept + SexParent + FoodTreatment + ArrivalTime + a_i
#'  a_i ~ N(0, sigma^2_Nst)    i = 1, ..., 27




#* Subsection 7.2: Execute the NB GLMM----

#' NB GLMM
M2 <- glmmTMB(NCalls ~ 1 + SexParent + FoodTreatment + ArrivalTimec + (1 | Nest),
              data = Owls,
              family = "nbinom2")


summary(M2)
performance(M2)


#' This is the fitted model.

#' For SexParent = Females and FoodTreatment = Deprived:
#' log(mu_ij) = 2.00 - 0.22 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Deprived:
#' log(mu_ij) = 2.00 + 0.07 - 0.22 * ArrivalTimec + Nest effect 

#' For SexParent = Females and FoodTreatment = Satiated:
#' log(mu_ij) = 2.00 - 0.68 - 0.22 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Satiated:
#' log(mu_ij) = 2.00 + 0.07 - 0.68 - 0.22 * ArrivalTimec + Nest effect 

#' The estimated parameters are rather similar to those of the Poisson GLMM!






#* Subsection 7.3: Model validation NB GLMM----

#' Do we have overdispersion?
testDispersion(M2)
#' Does this indicate an underdispersed NB? 
#' Try a generalised Poisson GLMM! Or stick to the Poisson distribution.


#' Can this model cope with the 26% of zeros?
testZeroInflation(M2)
#' No...it cannot. That is odd. Quite often the NB GLMM is well capable of
#' dealing with lots of zeros.



#' Can this model cope with the non-linear arrival time effect?
#' Get scaled quantile residuals and plot these against arrival time.

#' Get scaled quantile residuals.
E2Eqr <- simulateResiduals(fittedModel = M2, plot = FALSE)

plotResiduals(E2Eqr, 
              form = Owls$ArrivalTime, 
              xlab = "Rank-transformed arrival time") 
#' Nice non-linear patterns. We can see the same with Pearson residuals.

Owls$E2.pears <- resid(M2, type = "pearson")
p3 <- ggplot(Owls, 
             aes(x = ArrivalTime, y = E2.pears)) +
       xlab("Arrival time") + 
       ylab("Pearson residuaals") +
       geom_point() + 
       geom_smooth()
p3
#' There is indeed a non-linear residual pattern!

#' Homework: 
#'   -Plot the scaled quantile residuals versus the fitted values.
#'   -Plot the scaled quantile residuals versus each covariate in the model.
#'   -Plot the scaled quantile residuals versus each covariate not in the model.
#'   -Check the scaled quantile residuals for spatial patterns.
#'   -Check the random effects for spatial patterns (spoiler alert: these
#'    are spatially correlated).





#* Subsection 7.4: Fitted values of the NB GLMM----

#' We will plot the fitted values using the following steps:
#'  A. Specify covariate values for predictions
#'  B. Calculate the predicted values.
#'  C. Get a 95% CI for the predicted values.
#'  D. Plot the whole thing.


#' A. Define a grid of covariate values without extrapolation.
MyData2 <- ddply(Owls, 
                .(SexParent, FoodTreatment), 
                summarize,
                ArrivalTimec = seq(min(ArrivalTimec), 
                                   max(ArrivalTimec), 
                                   length = 50))
head(MyData2)
MyData2$Nest <- NA #' Ignore the random effects in the predictions (also 
#' something one can argue about).


#' B. Get the corresponding predicted values for a typical nest.
P2 <- predict(M2, 
              newdata = MyData2, 
              re.form =~ 0, 
              se = TRUE)


#' C. Calculate the predicted values and 95% CI.
MyData2$mu   <- exp(P2$fit)
MyData2$SeUp <- exp(P2$fit + 1.96 * P2$se.fit) 
MyData2$SeLo <- exp(P2$fit - 1.96 * P2$se.fit) 
head(MyData2)



#' D. And plot the whole thing

#' First back-standardize arival time.
MyData2$ArrivalTime <- MyData2$ArrivalTimec * sd(Owls$ArrivalTime) + mean(Owls$ArrivalTime)

p <- ggplot()
p <- p + geom_point(data = Owls, 
                    aes(x = ArrivalTime, 
                        y = NCalls))
p <- p + geom_line(data = MyData2, 
                   aes(x = ArrivalTime, 
                       y = mu))
p <- p + geom_ribbon(data = MyData2,
                     aes(x = ArrivalTime, 
                         ymax = SeUp, 
                         ymin = SeLo), 
                     alpha = 0.5)
p <- p + xlab("ArrivalTime") + ylab("Number of calls")
p <- p + theme(text = element_text(size=15)) 
p <- p + facet_grid(SexParent ~ FoodTreatment, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p







#* Subsection 7.5: What is the NB distribution doing?----


#' We are going to make a 3-d scatterplot. The code is not so nice on 
#' the eyes. Don't try to fully understand it. Just run it.
#' It can probably be programmed way more efficient.
#' MAC USERS: You may need to upgrade/install Quartz.




#' Focus on a specific SexParent/FoodTreatment combination.
ThisSexParent = "Female"
ThisFoodTr = "Deprived"

#ThisSexParent = "Male"
#ThisFoodTr = "Satiated"


#' Select the data for the selected SexParent and FoodTreatment levels.
MyData.sub <- subset(MyData2, 
                     SexParent == ThisSexParent & 
                     FoodTreatment == ThisFoodTr)

Owls.sub <- subset(Owls, 
                   SexParent == ThisSexParent & 
                   FoodTreatment == ThisFoodTr)


#' Three variable for 3-d scatterplot
x      <- MyData.sub$ArrivalTime  #' Covariate
y.nb   <- MyData.sub$mu           #' Fitted values
z      <- 0 * x                   #' A vector with 0s.

#' Number of NB density curves to plot.
K <- 6

#' We will pock the following rows from MyData.sub:
TheseRows <- round(seq(1, nrow(MyData.sub), length= K))
TheseRows


#' Plot the results.
#' First we plot the observed data for the selected SexParent/FoodTreatment
#' levels.
plot3d(x = Owls.sub$ArrivalTime,
       y = Owls.sub$NCalls,
       z = rep(0, nrow(Owls.sub)),
       type = "s",
       size = 0.8,
       lit = FALSE,
       zlim = c(0, 0.4),
       xlab = "Arrival time",
       ylab = "NCalls",
       zlab = "Density functions")

#' You can spin this graph! Don't close it!

#' Add the fitted values of the NB GLMM.
lines3d(x, y.nb, z, 
        col = 2,
        lwd = 3)

#' Add the NB density curves. 
interleave <- function(v1, v2) as.vector(rbind(v1, v2))
for (i in TheseRows){
  yseq <- seq(0, 26, by = 1)
  
  #' Calculate NB density values.
  dens.NB <-  dnbinom(yseq, size = sigma(M2), mu = MyData.sub$mu[i])
  
  #' Dotted line 
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, 0)
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 1,
          lty = 2)
  
  #' Plot the NB density values.
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, dens.NB)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = "purple")
}

#' Now compare the 3d graph with the relevant ggplot2 panel.
#' Questions: 
#'     -What is the role of the NB distribution?
#'     -How can you improve the model?
#'     -What do you do after a meal?       


#' Odd that the NB GLMM cannot cope that well with the zero
#' inflation. Double check:
simulationOutput.M1 <- simulateResiduals(fittedModel = M1, n = 1000)
simulationOutput.M2 <- simulateResiduals(fittedModel = M2, n = 1000)

par(mfrow = c(2, 1))
testZeroInflation(simulationOutput.M1)  #' Poisson GLMM
testZeroInflation(simulationOutput.M2)  #' NB GLMM
par(mfrow = c(1, 1))
#' Well..certainly the NB GLMM generates way more zeros. But not enough.









# Section 8: Wrong approach 2: ZIP GLMM----


#' Let us now do the ZIP GLMM.


#* Subsection 8.1: Model formulation----


#' We will apply the following ZIP GLMM.
# NCalls_ij      ~ ZIP(mu_ij, Pi)
# E(NCalls_ij)   = (1 - Pi) * mu_ij
# var(NCalls_ij) = (1 - Pi) * (mu + Pi * mu^2)

#' log(mu_ij) = Intercept + SexParent + FoodTreatment + ArrivalTime + a_i
#' logit(Pi) = Intercept

#'  a_i ~ N(0, sigma^2_Nest)    i = 1, ..., 27





#* Subsection 8.2: Execute the ZIP GLMM----

M3 <- glmmTMB(NCalls ~ 1 + SexParent + FoodTreatment + ArrivalTimec + (1 | Nest),
              data = Owls,
              ziformula = ~1,
              family = poisson)

drop1(M3, test = "Chi")
summary(M3)



#' We will now write down the fitted model.


#' 1. We first calculate the Pi:

#'           exp(-1.05360)
#' Pi =  -------------------   =  0.258
#'       1 + exp(-1.05360)



#' 2. Next we write down the expressions for the count part.

#' For SexParent = Females and FoodTreatment = Deprived:
#' log(mu_ij) = 2.18 - 0.15 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Deprived:
#' log(mu_ij) = 2.18 - 0.02 - 0.15 * ArrivalTimec + Nest effect 

#' For SexParent = Females and FoodTreatment = Satiated:
#' log(mu_ij) = 2.18 - 0.21 - 0.15 * ArrivalTimec + Nest effect 

#' For SexParent = Males and FoodTreatment = Satiated:
#' log(mu_ij) = 2.18 - 0.02 - 0.21 - 0.15 * ArrivalTimec + Nest effect 



#' 3. Now we combine them.
#' Recall that this was the model:

#'  E(NCalls_ij)   = (1 - Pi) * mu_ij
#'  log(mu_ij) = Intercept + SexParent + FoodTreatment + ArrivalTime + a_i
#'  logit(Pi) = Intercept



#' For the expected number of calls we therefore have:

#' For SexParent = Females and FoodTreatment = Deprived:
#' Exp(NCalls_ij) = (1 - 0.258) * exp(2.18 - 0.15 * ArrivalTimec + Nest effect) 

#' For SexParent = Males and FoodTreatment = Deprived:
#' Exp(NCalls_ij) = (1 - 0.258) * exp(2.18 - 0.02 - 0.15 * ArrivalTimec + Nest effect) 

#' For SexParent = Females and FoodTreatment = Satiated:
#' Exp(NCalls_ij) = (1 - 0.258) * exp(2.18 - 0.21 - 0.15 * ArrivalTimec + Nest effect) 

#' For SexParent = Males and FoodTreatment = Satiated:
#' Exp(NCalls_ij) = (1 - 0.258) * exp(2.18 - 0.02 - 0.21 - 0.15 * ArrivalTimec + Nest effect) 



#' Are you still with us? Better make a picture of this...in a moment.





#* Subsection 8.3: Model validation ZIP GLMM----


#' You now have to hope that a support package like DHARMa can cope with
#' more exotic distributions like a ZIP. And if it cannot, that it will then
#' tell you. We prefer to do all these DHARMa simulations outselves.

#' Do we have overdispersion?
simulationOutput.M3 <- simulateResiduals(fittedModel = M3, n = 1000)
par(mfrow = c(1,1))
testZeroInflation(simulationOutput.M3)  #' ZIP GLMM
#' That seems to be ok.

testZeroInflation(M3)  #' ZIP GLMM
#' That seems to be ok.

#' Can this model cope with the 25% of zeros?
testZeroInflation(M3)
#' Yes it can.


#' Can this model cope with the non-linear arrival time effect?
#' Get scaled quantile residuals and plot these against arrival time.

#' Get scaled quantile residuals.
E3Eqr <- simulateResiduals(fittedModel = M3, plot = FALSE)

plotResiduals(E3Eqr, 
              form = Owls$ArrivalTime, 
              xlab = "Rank-transformed arrival time") 
#' Nice non-linear patterns. We can see the same with Pearson residuals.



#' Homework: 
#'   -Plot the scaled quantile residuals versus the fitted values.
#'   -Plot the scaled quantile residuals versus each covariate in the model.
#'   -Plot the scaled quantile residuals versus each covariate not in the model.
#'   -Check the scaled quantile residuals for spatial patterns.
#'   -Check the random effects for spatial patterns (spoiler alert: these
#'    are spatially correlated).





#* Subsection 8.4: Fitted values of the ZIP GLMM----

#' We will use the principles explained in this paper:
#' Brooks et al. (2017).
#' glmmTMB Balances Speed and Flexibility Among Packages for Zero-inflated
#' Generalized Linear Mixed Modeling

#' https://journal.r-project.org/archive/2017/RJ-2017-066/RJ-2017-066.pdf
#' See page 391. Alternative prediction method.


#' Sketch the fitted values:
#'  A. Specify covariate values for predictions.
#'  B. Create X matrix with expand.grid.
#'  C. Calculate predicted values.
#'  D. Calculate standard errors (SE) for predicted values via simulation.
#'  E. Plot predicted values.
#'  F. Plot predicted values +/- 2* SE.


#' A. Define a grid of covariate values without extrapolation.
#' We will predict ZIP values for these covariate values + SEs.

#' A. Define a grid of covariate values without extrapolation.
MyData3 <- ddply(Owls, 
                .(SexParent, FoodTreatment), 
                summarize,
                ArrivalTimec = seq(min(ArrivalTimec), 
                                   max(ArrivalTimec), 
                                   length = 50))
head(MyData3)
MyData3$Nest <- NA #' Ignore the random effects in the predictions (also 
                  #' something one can argue about).



#' B. Create X matrix from MyData
Xcount <- model.matrix(~ SexParent + FoodTreatment + ArrivalTimec, 
                       data = MyData3) 



#' C. Calculate predicted values.
#' Count part.
Betas     <- fixef(M3)$cond
eta.count <- Xcount %*% Betas
mu.count  <- exp(eta.count) 

#' For the next subsection:
MyData3$mu.count <- mu.count


#' Bernoulli part.
ziformula <- M3$modelInfo$allForm$ziformula #Intercept
X.zi      <- model.matrix(lme4::nobars(ziformula), MyData3)
gammas    <- fixef(M3)$zi       
eta.zi    <- X.zi %*% gammas    
Pi        <- exp(eta.zi) / (1 + exp(eta.zi))
Pi 

#' For the next subsection:
MyData3$Pi <- Pi


#' ZIP expected values: (1-Pi) * mu
MyData3$Pred.zip <- (1 - Pi) * mu.count 
head(MyData3)




#'  D. Calculate standard errors (SE) for predicted values
#' This is now getting difficult. But it is all copy-paste from the paper.

#' Copy from the paper:
#'  For the standard errors/confidence intervals, we could use posterior 
#'  predictive simulations (i.e., draw multivariate normal samples from the 
#'  parameter for the fixed effects). This conditions on/ignores uncertainty in 
#' the random-effect parameters.

#' Another reference for this is Wood (2018). He calls this posterior sampling
#' and does it in the context of GAMs.


#' The code below is going to simulate betas from:  
#'   N(Betas, SIGMA)
#' where Betas are the estimated betas from the model and
#' SIGMA is the variance-covariance matrix of the betas.
#' And we are going to do it 1000 times.

#' Get the betas and variance-covariance matrix of model M3.
Betas <- fixef(M3)$cond
SIGMA <- vcov(M3)$cond

#' Simulate the betas. 
set.seed(101)
Betas.1000 <- mvrnorm(1000,
                      mu = Betas,
                      Sigma = SIGMA)
#' What is this?
head(Betas.1000)
dim(Betas.1000)
#' Answer: 1000 simulated sets of betas for the regression parameters 
#'         in the count part.

#' Calculate X * Betas.1000.
eta.1000 <- Xcount %*% t(Betas.1000)
mu.1000  <- exp(eta.1000)
#' What is mu.1000? 
#' Answer: That is exp(X * beta) for 1000 simulated sets of betas


#' Binary part
#' Now we simulate the gammas 1000 times.
gammas.1000 <- mvrnorm(1000,
                       mu = gammas,
                       Sigma = vcov(M3)$zi)

#' We just simulated 1000 sets of gammas from N(gammas, vcov(M3)$zi)

#' Now calculate eta.zi 1000 times.
#' And calculate Pi 1000 times.
eta.zi.1000 <- X.zi %*% t(gammas.1000)
Pi.1000 <- exp(eta.zi.1000) / (1 + exp(eta.zi.1000))
#' What is this? 
#' This is a simulated Pi 1000 times...for each row in MyData
#' Because we only use an intercept, the rows are identical.
Pi.1000[, 1:10]



#' Get the ZIP expected values.
#' We now calculate   (1-Pi) * mu   a thousand times:
Zip.1000 <- (1 - Pi.1000) * mu.1000

#' Check dimensions.
dim(Pi.1000)
dim(mu.1000)
dim(Zip.1000)
#' 50 values for arrival time, 2 SexParent combinations, 2 Treatment levels.
#' And we have 1000 times predicted values.



#' Now use quantiles to obtain the 95% confidence intervals from
#' the 1000 sets of Zip.1000
ci.zip <- t(apply(Zip.1000,
                  MARGIN = 1,     #' Work on the rows
                  FUN = quantile, #' Calculate quantiles
                  c(0.025, 0.975)))
ci.zip <- data.frame(ci.zip)
names(ci.zip) <- c("zip.low", "zip.high")
head(ci.zip)

#' And add these to Mydata3.
MyData3$SeUp <- ci.zip$zip.high
MyData3$SeLo <- ci.zip$zip.low
head(MyData3)



#' D. And plot the whole thing

#' First back-standardize arrival time.
MyData3$ArrivalTime <- MyData3$ArrivalTimec * sd(Owls$ArrivalTime) + mean(Owls$ArrivalTime)

p <- ggplot()
p <- p + geom_point(data = Owls, 
                    aes(x = ArrivalTime, 
                        y = NCalls))
p <- p + geom_line(data = MyData3, 
                   aes(x = ArrivalTime, 
                       y = Pred.zip))   #'  <--- Changed
p <- p + geom_ribbon(data = MyData3,
                     aes(x = ArrivalTime, 
                         ymax = SeUp, 
                         ymin = SeLo), 
                     alpha = 0.5)
p <- p + xlab("ArrivalTime") + ylab("Number of calls")
p <- p + theme(text = element_text(size=15)) 
p <- p + facet_grid(SexParent ~ FoodTreatment, 
                    scales = "fixed")
p <- p + theme(legend.position="none") 
p








#* Subsection 8.5: What is the ZIP distribution doing?----


#' We are going to make a 3-d scatterplot. The code is not so nice on 
#' the eyes. Don't try to fully understand it. Just run it.
#' It can probably be programmed way more efficient.
#' MAC USERS: You may need to upgrade/install Quartz.




#' Focus on a specific SexParent/FoodTreatment combination.
ThisSexParent = "Female"
ThisFoodTr = "Deprived"

#ThisSexParent = "Male"
#ThisFoodTr = "Satiated"


#' Select the data for the selected SexParent and FoodTreatment levels.
MyData.sub <- subset(MyData3, 
                     SexParent == ThisSexParent & 
                       FoodTreatment == ThisFoodTr)

Owls.sub <- subset(Owls, 
                   SexParent == ThisSexParent & 
                     FoodTreatment == ThisFoodTr)


#' Three variable for 3-d scatterplot
x      <- MyData.sub$ArrivalTime  #' Covariate
y.zip  <- MyData.sub$Pred.zip     #' Fitted values..called Pred.zip here!
z      <- 0 * x                   #' A vector with 0s.

#' Number of ZIP density curves to plot.
K <- 6

#' We will pock the following rows from MyData.sub:
TheseRows <- round(seq(1, nrow(MyData.sub), length= K))
TheseRows


#' Plot the results.
#' First we plot the observed data for the selected SexParent/FoodTreatment
#' levels.
plot3d(x = Owls.sub$ArrivalTime,
       y = Owls.sub$NCalls,
       z = rep(0, nrow(Owls.sub)),
       type = "s",
       size = 0.8,
       lit = FALSE,
       zlim = c(0, 0.4),
       xlab = "Arrival time",
       ylab = "NCalls",
       zlab = "Density functions")

#' You can spin this graph! Don't close it!

#' Add the fitted values of the ZIP GLMM.
lines3d(x, y.zip, z, 
        col = 2,
        lwd = 3)

#' Add the ZIP density curves. 
interleave <- function(v1, v2) as.vector(rbind(v1, v2))
for (i in TheseRows){
  yseq <- seq(0, 26, by = 1)
  
  #' Calculate ZIP density values.
  dens.ZIP <- VGAM::dzipois(yseq, 
                            lambda = MyData.sub$mu.count[i], 
                            pstr0 = MyData.sub$Pi[i])
  
  
  
  #' Dotted line 
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, 0)
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 1,
          lty = 2)
  
  #' Plot the ZIP density values.
  rb = cbind(MyData.sub$ArrivalTime[i], yseq, dens.ZIP)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = "purple")
}

#' Now compare the 3d graph with the relevant ggplot2 panel.
#' Questions: 
#'     -What is the role of the ZIP distribution?
#'     -How can you improve the model?
#'     -What do you do after a meal?       
#'     -How would you allow for different height in probabilities
#'      for NCalls = 0 per Treatment/SexParent level?



#' Another question: What would a ZINB distribution doing?








# Section 9: Wrong approach 3: ZINB GLMM----


#' There is actually no need to do a ZINB as the ZIP is not overdispersed.
#' But let's do it anyway!


#' ZINB GLMM
M4 <- glmmTMB(NCalls ~ 1 + SexParent + FoodTreatment + ArrivalTimec +(1 | Nest),
              data = Owls,
              ziformula = ~1,
              family = "nbinom2")

#' This model can cope with the 26% of zeros and it is not overdispersed.
#' But it has strong non-linear residual patterns wit arrival time.


#' Compare all 4 models via AICs.
AIC(M1, M2, M3, M4)


#' That is odd. Why would the ZINB GLMM be better than the ZIP GLMM?
#' But it is only slightly better.
#' We will stop at this point as it becomes a rather long presentation.
#' And we are on the wrong track anyway.


#' To add more options, we strongly suggest that you also try the
#' generalised Poisson GLMM. It also has a zero-inflation cousin.
#' See: 
#'   The World of Zero-Inflated Models (2021) 
#'   Zuur and Ieno (2021)
#'   www.highstat.com






# Section 10: Comparing the Poisson, NB and ZIP distributions----

#' We will now put all three distributions (Poisson, NB and ZIP) in
#' one 3-d graph.


#' Focus on a specific SexParent/FoodTreatment combination.
ThisSexParent = "Female"
ThisFoodTr = "Deprived"

#ThisSexParent = "Male"
#ThisFoodTr = "Satiated"


#' Select the data for the selected SexParent and FoodTreatment levels.
#' Poisson results
MyData1.sub <- subset(MyData1, 
                     SexParent == ThisSexParent & 
                     FoodTreatment == ThisFoodTr)

#' NB results
MyData2.sub <- subset(MyData2, 
                      SexParent == ThisSexParent & 
                      FoodTreatment == ThisFoodTr)

#' ZIP results
MyData3.sub <- subset(MyData3, 
                      SexParent == ThisSexParent & 
                      FoodTreatment == ThisFoodTr)

#' Observed data
Owls.sub <- subset(Owls, 
                   SexParent == ThisSexParent & 
                     FoodTreatment == ThisFoodTr)


#' Three variable for 3-d scatterplot
#' It is the same for each model.
x <- MyData1.sub$ArrivalTime  #' Covariate
z <- 0 * x                    #' A vector with 0s.


#' Fit for each model
y.pois <- MyData1.sub$mu          #' Fitted values
y.nb   <- MyData2.sub$mu           #' Fitted values
y.zip  <- MyData3.sub$Pred.zip    #' Fitted values..called Pred.zip here!


#' Number of density curves to plot.
K <- 4

#' We will pock the following rows from MyData1.sub:
TheseRows <- round(seq(1, nrow(MyData1.sub), length= K))
TheseRows


#' Plot the results.
#' First we plot the observed data for the selected SexParent/FoodTreatment
#' levels.
plot3d(x = Owls.sub$ArrivalTime,
       y = Owls.sub$NCalls,
       z = rep(0, nrow(Owls.sub)),
       type = "s",
       size = 0.8,
       lit = FALSE,
       zlim = c(0, 0.4),
       xlab = "Arrival time",
       ylab = "NCalls",
       zlab = "Density functions")

#' You can spin this graph! Don't close it!

#' Add the fitted values for each model
lines3d(x, y.pois, z, col = 3, lwd = 2)
lines3d(x, y.nb, z, col = 2, lwd = 2)
lines3d(x, y.zip, z, col = 1, lwd = 2)


#' Add the density curves. 
interleave <- function(v1, v2) as.vector(rbind(v1, v2))

Threshold <- 0.4 #' Horizontal displacement between the three density curves.
for (i in TheseRows){
  yseq <- seq(0, 26, by = 1)
  
  #' Calculate Poisson density values.
  dens.pois <- dpois(yseq, lambda = MyData1.sub$mu[i])
  
  #' Calculate NB density values.
  dens.NB <-  dnbinom(yseq, size = sigma(M2), mu = MyData2.sub$mu[i])
  
  #' Calculate ZIP density values.
  dens.ZIP <- VGAM::dzipois(yseq, 
                            lambda = MyData3.sub$mu.count[i], 
                            pstr0 = MyData3.sub$Pi[i])
  
  
  
  #' Dotted line (same for each distribution).
  rb = cbind(MyData1.sub$ArrivalTime[i], yseq, 0)
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 1,
          lty = 2)
  

  #' Plot the Poisson density values.
  rb = cbind(MyData1.sub$ArrivalTime[i] - Threshold, yseq, dens.pois)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = 3)

  
  lines3d(rb[,1], rb[,2], rb[,3], 
          col = 3,
          lty = 2)
  
  
    #' Plot the NB density values.
  rb = cbind(MyData2.sub$ArrivalTime[i], yseq, dens.NB)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = 2)
  
  
  #' Plot the ZIP density values.
  rb = cbind(MyData3.sub$ArrivalTime[i] + Threshold, yseq, dens.ZIP)
  segments3d(interleave(rb[,1], rb[,1]),
             interleave(rb[,2], rb[,2]),
             interleave(rep(0, nrow(rb)), rb[,3]),
             alpha = 0.9,
             col = 1)
}












# Section 11: Poisson GAMM----


#' This should have been the second model, or perhaps even the starting
#' point of the entire analysis.
M5 <- gam(NCalls ~ 1 + SexParent + FoodTreatment + 
                   s(ArrivalTime, bs = "cr") +
                   s(Nest, bs = "re"),
          select = TRUE,
          method = "REML",
          data = Owls,
          family = "poisson")

summary(M5)
draw(M5, select = 1)


#' Now the whole game starts over again. 
#'  -Execute a NB GAMM.
#'  -Execute a ZIP GAMM.
#'  -Execute a ZINB GAMM.

#' In our course:  
#'    Zero-inflated GLM, GAM, GLMM and GAMM for the analysis of 
#'    spatial and spatial-temporal correlated data using R-INLA
#' we apply a zero-inflated Poisson GAM with spatial correlation in 
#' these data. The reason for that is as follows.

#' We will extract the estimated random effects for the nests, obtained by
#' the GAMM.



#' We first specify a MyData data frame containing the Site identities
MyData <- data.frame(Nest = levels(Owls$Nest))

#' Next we set all other covariates to something (what is irrelevant).
MyData$ArrivalTime   <- 0
MyData$SexParent     <- factor("Female", levels = c("Female", "Male"))
MyData$FoodTreatment <- factor("Deprived", levels = c("Deprived", "Satiated"))


#' This is what we just created:
head(MyData, 50)


#' We then use the predict function with type = "terms"
P5 <- predict(M5, 
              MyData, 
              type = "terms", 
              se.fit = TRUE)
head(P5$fit)

#' If you were to add up all the values on a specific row, then
#' We extract the estimated random effects and their standard errors:
MyData$ai    <- P5[["fit"]][ , "s(Nest)"]
MyData$ai.se <- P5[["se.fit"]][ , "s(Nest)"]

#' Here is our final product. The ai column contains the estimated
#' random effects. They are the same as those in coef(M2). There is no
#' benefit in extracting them in this way.
head(MyData)

#' We put the random effects in a dotchart.
MyData$RowID <- 1:nrow(MyData)

p <- ggplot(data = MyData, 
            aes(x = ai, 
                y = RowID)) 
p <- p + geom_point(size = 1, col = grey(0.5)) 
p <- p + xlab("Estimated random effect")
p <- p + ylab("Nest identity")                       
p <- p + theme(text = element_text(size = 15)) 
p <- p + geom_vline(xintercept = 0, lty = 2, col = "red")
p

#' And add errorbars.
p2 <- p + geom_errorbarh(data = MyData,
                         aes(y = RowID, 
                             xmax = ai + 1.96 * ai.se, 
                             xmin = ai - 1.96 * ai.se),
                         col = "blue",
                         height = 0.05,
                         alpha = 0.2)
p2
#' We can't see the spatial information in this graph.
#' To do this, we first calculate average latitude and longitude 
#' values for each site.
MyData$Latitude  <- tapply(Owls$Latitude, FUN = mean, INDEX = Owls$Nest)
MyData$Longitude <- tapply(Owls$Longitude, FUN = mean, INDEX = Owls$Nest)

#' We will plot longitude versus latitude, and superimpose the 
#' random effects on this graph. We will use different colours
#' for positive and negative random effects, and a large dot
#' refers to a large (in absolute sense) random effect. There should
#' be no spatial clustering of random effects.
MyData$Size <- abs(MyData$ai)
MyData$Sign <- ifelse(MyData$ai >=0, "Positive", "Negative")

p <- ggplot(data = MyData, 
            aes(x = Longitude, 
                y = Latitude,
                size = Size,
                col = Sign)) 
p <- p + geom_point() 
p <- p + xlab("Longitude")
p <- p + ylab("Latitude")     
p <- p + coord_fixed(ratio = 1) 
p <- p + theme(text = element_text(size = 15)) 
p
#' Bugger...do we have a spatial pattern in these random effects?










# Section 12: Final comments----

# We can also ask the questions:
# 1. What is driving the chicks to make noise?
#    This is a noise/no-noise question.
# 2. And once they make noise, what is driving the 
#    amount of noise?
# These two questions lead to a ZAP GLMM.

# We would refrain from applying a negative binomial GLMM for the moment.



