**NOTE:**

Fit_models_in_secr_beardata.R will run all of the same models covered in the webinar, 
including a null model, a model with sex specific detection parameters, a behavioral model, 
and a model with forest cover as a covariate for density.

To run the script, you will need the secr and terra packages.
You will also need data from the main Data folder and the "beartraps.csv" file, which
has the trap coordinates in UTMs, which is required for secr. 

Please note, the script provided here is just for simple demonstration to match the models being fit in the Statistical Methods Seminar Series.  This is not meant to be a tutorial on the secr package, but if you are having trouble running Nimble, this will be a way to fit the same models for comparison.
Also, the estimated density from secr is in units of # of animals/ha.
&#x20;

