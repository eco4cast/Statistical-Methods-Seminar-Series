Make sure you have the packages reshape, dplyr and RTMB installed.
If not they can be install using the following R commands:

install.packages("reshape")
install.packages("dplyr")
install.packages("RTMB")

The file "sim.R" contains the individual based simulation used in this example.
It is used to generate the "sim.rds" file that is included and contains the
results of the simulation.

The file "makedata.R" "finds" the kin pairs from the simulation and sets up
the data in a way that can be used by the model. The result of this is
"mod_dat.rds".

Finally "model.R" is our simplified CKMR model along with code to run it
and plot results.
