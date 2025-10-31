# Code based on:
# https://prioritizr.net/articles/prioritizr.html

# load packages
library(prioritizrdata)
library(prioritizr)
library(sf)
library(terra)
library(vegan)
library(cluster)

# set seed for reproducibility
set.seed(500)

# load planning unit data
tas_pu <- get_tas_pu()

# load feature data
tas_features <- get_tas_features()

# print planning unit data
print(tas_pu)

# set costs for existing protected areas to zero
tas_pu$cost <- tas_pu$cost * !tas_pu$locked_in

# plot map of planning unit costs
plot(st_as_sf(tas_pu[, "cost"]), main = "Planning unit costs")

# plot map of planning unit coverage by protected areas
plot(st_as_sf(tas_pu[, "locked_in"]), main = "Protected area coverage")

# print planning unit data
print(tas_features)

# plot map of the first four vegetation classes
plot(tas_features[[1:4]])

# build problem
p1 <-
  problem(tas_pu, tas_features, cost_column = "cost") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.005) %>%
  add_relative_targets(0.17) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

# print problem
print(p1)

# Installing a solver 
# Supported solvers:
#   https://prioritizr.net/reference/solvers.html

# If you have access to Gurobi, that's the recommended solver:
# https://prioritizr.net/reference/add_gurobi_solver.html

# If you want to use an open source solver CBC is the recommendation:
# https://prioritizr.net/reference/add_cbc_solver.html

# if that doesn't work Symphony is an alternative:
# https://prioritizr.net/reference/add_rsymphony_solver.html

# solve problem
s1 <- solve(p1)

# plot map of prioritization
plot(
  st_as_sf(s1[, "solution_1"]), main = "Prioritization",
  pal = c("grey90", "darkgreen")
)


# create column with existing protected areas
tas_pu$pa <- round(tas_pu$locked_in)

# calculate feature representation statistics based on existing protected areas
tc_pa <- eval_target_coverage_summary(p1, tas_pu[, "pa"])
print(tc_pa)

# calculate  feature representation statistics based on the prioritization
tc_s1 <- eval_target_coverage_summary(p1, s1[, "solution_1"])
print(tc_s1)

# explore representation by existing protected areas
## calculate number of features adequately represented by existing protected
## areas
sum(tc_pa$met)

## summarize representation (values show percent coverage)
summary(tc_pa$relative_held * 100)

## visualize representation  (values show percent coverage)
hist(tc_pa$relative_held * 100,
     main = "Feature representation by existing protected areas",
     xlim = c(0, 100),
     xlab = "Percent coverage of features (%)")

# explore representation by prioritization
## summarize representation (values show percent coverage)
summary(tc_s1$relative_held * 100)


## calculate number of features adequately represented by the prioritization
sum(tc_s1$met)

## visualize representation  (values show percent coverage)
hist(
  tc_s1$relative_held * 100,
  main = "Feature representation by prioritization",
  xlim = c(0, 100),
  xlab = "Percent coverage of features (%)"
)


# calculate relative importance
imp_s1 <- eval_rank_importance(p1, s1["solution_1"], n = 10)
print(imp_s1)

# manually set locked in planning units to -1 to help with visualization,
# this way we can easily see the importance scores for the priority areas
imp_s1$rs[tas_pu$locked_in] <- -1


# plot map of importance scores
plot(st_as_sf(imp_s1[, "rs"]), main = "Overall importance")


# create new problem with a portfolio added to it
p2 <-
  p1 %>%
  add_gap_portfolio(number_solutions = 1000, pool_gap = 0.2)

# print problem
print(p2)

# generate prioritizations
prt <- solve(p2)

print(prt)


# extract solutions
prt_results <- sf::st_drop_geometry(prt)
prt_results <- prt_results[, startsWith(names(prt_results), "solution_")]

# calculate pair-wise distances between different prioritizations for analysis
prt_dists <- vegan::vegdist(t(prt_results), method = "jaccard", binary = TRUE)

# run cluster analysis
prt_clust <- hclust(as.dist(prt_dists), method = "average")

# visualize clusters
opar <- par()
par(oma = c(0, 0, 0, 0), mar= c(0, 4.1, 1.5, 2.1))
plot(
  prt_clust, labels = FALSE, sub = NA, xlab = "",
  main = "Different prioritizations in portfolio"
)
suppressWarnings(par(opar))

# run k-medoids analysis
prt_med <- pam(prt_dists, k = 6)

# extract names of prioritizations that are most central for each group.
prt_med_names <- prt_med$medoids
print(prt_med_names)


# create a copy of prt and set values for locked in planning units to -1
# so we can easily visualize differences between prioritizations
prt2 <- prt[, prt_med_names]
prt2[which(tas_pu$locked_in > 0.5), prt_med_names] <- -1

# plot a map showing main different prioritizations
# dark grey: locked in planning units
# grey: planning units not selected
# green: selected planning units
plot(st_as_sf(prt2), pal = c("grey60", "grey90", "darkgreen"))





