##Here we are going to take the data we simulated, "find" the number of POPs
## and setup the data in a way to feed the model


##Load in the data we simulated
sim = readRDS("./sim.rds")
pop = sim$pop
sampledA = sim$sampledA
sampledJ = sim$sampledJ

 names(sampledA) = paste0(names(sampledA),"_A")
  names(sampledJ) = paste0(names(sampledJ),"_J")



##Check which sampled adults are mothers of juveniles, these are POPs
## This way is really really inefficient but ok for a small example and is easy to understand
library(reshape)
##we compare every row of the A df with every row of the J df
POPdf = expand.grid.df(sampledA,sampledJ)

##Remove rows where the juvenile is born before the adult (since they can't be the parent)
POPdf = POPdf[POPdf$byear_J > POPdf$byear_A,]
##Remove rows where the adult was sampled before the juvenile was born
POPdf = POPdf[POPdf$cyear_A >= POPdf$byear_J,]

##Find the parent offspring pairs
POPdf$is_POP = POPdf$id_A == POPdf$mother_J

library(dplyr)
##Since all parents have the same average fecundity, we only need to consider  
POPsum = POPdf |>
  group_by(byear_J) |>
  summarise(n_POP = sum(is_POP),n_comp=n())



##Finding HSPs within the juvenile samples

sampledJ2 = sampledJ
names(sampledJ2) = paste0(names(sampledJ2),"2")

HSPdf = expand.grid.df(sampledJ,sampledJ2)
##Remove rows where it's the same individual
HSPdf = HSPdf[HSPdf$id_J != HSPdf$id_J2,]
##Remove rows with J2 being older than J
HSPdf = HSPdf[HSPdf$byear_J < HSPdf$byear_J2,]


##Find the half sibling pairs
HSPdf$is_HSP = HSPdf$mother_J == HSPdf$mother_J2

HSPsum = HSPdf |>
  group_by(byear_J,byear_J2) |>
  summarise(n_HSP = sum(is_HSP),n_comp=n())

##save our data
mod_dat = list(POPs = POPsum,HSPs = HSPsum)
saveRDS(mod_dat,file="mod_dat.rds")
