##set a seed so it's reproducible
set.seed(42)


##Setup demographic parameters
beta = 2.05
s = 0.37

##Number of years to simulate
Y = 30


##Create Initial population
pop = data.frame(id=1:100,
		 mother= NA,
		 byear=-1,
		 cyear=0,
		 mature=TRUE,
		 alive=TRUE)

##Where to keep our sampled individuals
sampAdults = list()
sampJuvs = list()

##What years to sample
samp_years = 21:30

## How many to sample each year
n_s_adults = 40
n_s_juvs = 40


##function to create new rows of the pop data.frame for new offspring
## for each year of our simulation
make_kids <- function(ids,n_kids,cyear,max_id){
  ##Make the new mother column 
  new_m = rep(ids,n_kids)
  byear = cyear
  mature = FALSE
  alive = TRUE
  n_df = data.frame(id=seq(from=max_id+1,length.out=length(new_m)),
		    mother=new_m,
		    byear=byear,
		    cyear=cyear,
		    mature=mature,
		    alive=alive)
  n_df
}

##For comparison keep a record of how many adults are in the population each year
adult_pop = numeric(Y)

for(y in 1:Y){

  print(y)

  ##Find those that live to year y
  c_pop = pop[pop$alive == TRUE & pop$cyear == y-1,]

  ##determine which survive
  survived = rbinom(nrow(c_pop),1,s)
  c_pop$alive = ifelse(survived == 1,TRUE,FALSE)
  c_pop = c_pop[c_pop$alive == TRUE,]

  ##If they survived then they made it a year so they're mature
  c_pop$mature = TRUE

  ## Up the year
  c_pop$cyear = y

  ##Make the new kids
  n_births = rpois(nrow(c_pop),beta)
  new_kids = make_kids(c_pop$id,n_births,y,max(pop$id))

  ##Sample adults and juveniles
  if(y %in% samp_years){
    sampA = sample(1:nrow(c_pop),n_s_adults)
    sampJ = sample(1:nrow(new_kids),n_s_juvs)
    ##Make sampled individuals not alive anymore
    c_pop$alive[sampA] = FALSE
    new_kids$alive[sampJ] = FALSE
    ##Store them
    sampAdults[[y]] = c_pop[sampA,]
    sampJuvs[[y]] = new_kids[sampJ,]
  }else{
    sampAdults[[y]] = NA
    sampJuvs[[y]] = NA
  }

  print(paste("Cur adult pop. size:",nrow(c_pop)))
  adult_pop[y] = nrow(c_pop)

  ## combine kids with adults  
  c_pop = rbind(c_pop,new_kids)
  ## put them back in the population
  pop = rbind(pop,c_pop)
}

##Get the samples and combine them into one data.frame
sampledA = sampAdults[samp_years]
sampledA = do.call(rbind,sampledA)
sampledJ = sampJuvs[samp_years]
sampledJ = do.call(rbind,sampledJ)

##Save our simulated population and samples
sim = list(pop=pop,sampledA = sampledA,sampledJ = sampledJ,adult_pop=adult_pop)
saveRDS(sim,file="./sim.rds")
