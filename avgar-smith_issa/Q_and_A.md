# Follow-up Q & A

## How does this approach (mostly) resolve statistical issues of temporal autocorrelation?

In unconditional HSA, used points are considered independent samples from the animal’s utilization distribution. However, positional data are collected at high frequency and are hence characterized by strong temporal autocorrelation – the spatial position of the animal at one time is much more likely to be proximate to the spatial position of the animal a short time ago, or after, than to the spatial position of the animal a long time ago, or after. This means that the effective sample size could be much smaller than the number of points (as these points are not actually independent) and hence standard errors may be grossly underestimated. In conditional HSA, such as iSSA, the inference is based on within-cluster contrasts, and hence whether used points (belonging to different clusters) are independent or not has negligible consequences.

## Can you the difference between the tentative and updated distributions?

The tentative step-length and turn-angle distributions reflect our initial guesses as to the ‘true’ step-length and turn-angle distributions (the distributions that, according to the assumption of the underlying model, our animal truly samples from when deciding where to go next). As long as habitat selection plays a part in the animal’s decision-making process, the directly observed distributions of step lengths and turn angles always deviate from these true distributions (but we still recommend basing our initial guess on these observations). This is perhaps best explained by imagining an animal moving across a heterogeneous landscape, being more likely to move into (‘select’) preferred habitats as long as those are available to it. Assuming that these preferred habitats are as available as any other habitat across the landscape, we would expect the animal spends more time in these preferred habitats. Not only that, whenever it is within a preferred habitat, moving outside of it is less appealing and so the animal will do less of such movements. The result is that our animal is moving less frequently and for shorter distances than we might expect based on its true step-length distribution alone. In fact, the only way the true step-length and turn-angle distribution could be fully manifested (without distortion) in the observed movement pattern, is when movement occurs across a homogeneous landscape or in the absence of any habit selection. This is the reason we refer to these ‘true’ distributions as defining the ‘selection-free’ movement kernel. iSSA allows the user to estimate the parameters of these true distributions without ever observing them directly.       
## What does it mean if the parameters of the selection-free movement kernel (the adjusted shape/scale/concentration of the step-length or turn-angle distributions) come out negative?

If the von Mises concentration parameter is negative, the adjusted turn angle distribution is centered at 𝜋 (180°) rather than 0 (negative directional autocorrelation; the animal is more likely to turn back). If any of the parameters of the step-length distribution is negative, the model is ill-fitted. Try a different step-length distribution. Include an interaction between step-length and turn-angle (e.g., ln[step length]  × cos[turn angle]). Include additional step-length interactions. Remove ‘non-movement’ steps (e.g., based on a step-length threshold). Resample the data to a coarser temporal resolution (longer step duration)

## If the agent in the model always has to be stepping, how do you approach resting times or states? Do you use a different distribution during resting times with much smaller step sizes? Are non-movement steps those steps where the animal is sleeping/chilling out? Or would foraging movement within 5 m areas considered non-movement? Could a sleeping animal be mistaken for foraging movement if there is some error in positional accuracy? Can I (and should I) combine iSSA with behavioral-state assignment or path segmentation analysis (e.g., a Hidden Markov Model)?

iSSA assumes step-length distribution is continuous. As such, it cannot and should not be used if the data is ‘0-inflated’. The typical signature of GPS error around a stationary position is a half-normal step length distribution (with very small standard deviation; 3-30 meters), and a π-centered turn-angle distribution (indicating back and forth positional jumps). It is possible to fit iSSF to these distributions, but it is entirely not clear what might one learn from that. If the data contains a mix of ‘non-movement’ and movement ‘steps’, it is recommended to either resample to a coarser temporal grain until ‘non-movement’ steps become infrequent, or leave the ‘non-movement’ steps out of the iSSA (instead, one can examine the probability of being stationary as a function of environmental covariates in a separate analysis using, e.g., a logistic regression). It is up to the user, based on the biology of the system, the temporal grain of the data, and the accuracy of the GPS, to decide what should be considered part of the animal’s movement process, and what should be considered ‘stationary’. It is generally a good idea, if the data likely contain a mixture of distinct behavioral states (e.g., resting, foraging, and dispersing) to first segment the data into those states (using, e.g., HMM) and then fit a separate iSSF to each one.

## When we update our movement distribution or turning angle distribution we don't bring any error ... is it possible to update our distribution and show the uncertainty? This is particularly interesting for interactions as you did for day/night.

It is possible, and recommended, to incorporate error estimates into our selection-free movement kernel, particularly if one infers on the resulting expectations (e.g., the mean step-length under certain conditions). In some cases, it is straightforward to calculate the 95% confidence bounds analytically based on the iSSF parameter’s standard error. If you are trying to update a distribution with a single parameter -- *e.g.* if you are using the coefficient for `cos_ta_` to update the von Mises concentration parameter or if you are using the coefficient for `sl_` to update the exponential rate parameter -- then you can simply use the upper and lower bounds of the 95% confidence interval for that coefficient to update the tentative distribution to the upper and lower bounds of the selection-free movement distribution.

 It becomes trickier when you are trying to use two parameters to update the distribution -- *e.g.* if you are using the coefficients for `log_sl_` and `sl_` to update the gamma distribution. This is because the coefficients for `log_sl_` and `sl_` will *covary*, so combining the upper and lower bounds for each coefficient will tend to overestimate the uncertainty in the selection-free movement distribution. In that case, one might choose to use parametric bootstrapping (relying on the model’s variance-covariance matrix), non-parametric bootstrapping (sampling clusters with replacement), or empirical distribution across different data instances (e.g., individuals, study area, or years) or tentative parameter values. *I.e.*, you can use any method that allows you to incorporate the covariance into the prediction of the parameters. 
 
## How can I combine data from multiple individuals?

If all you want is ‘population-level’ inference, include the same number of clusters from each individual in a single model. 

Mixed-effects approach: tricky in conditional logistic regression, use a Poisson regression instead (see [Muff et al. 2020](https://doi.org/10.1111/1365-2656.13087)). 

"Two-step" approach: first conduct independent iSSA for each individual, then summarize across  individual coefficient estimates using weighted bootstrapping [e.g., using AIC weights], or linear models with inverse-variance weights.

## Do you think that iSSF can be used in the case of geese? They have very small step lengths when they walk and when they change habitat, they fly.

Yes, iSSA can be used with geese (or other birds), but it is probably a good idea to fit separate models for walking and flying. If the animal could reasonably get to any point within the study area with a single step, reconsider using iSSA – an exponential HSA (where availability is assumed constant across space and time) might be a better choice.   

## Is there more levels in Time of Day than just day vs night?

Time of day may be represented as a categorical variable with 2 (‘day’/’night’), 3 (‘day’/’twilight’/’night’), or 4 (‘day’, ‘dusk’, ‘night’, ‘dawn’) levels. Time of day can also be represented as continuous variable/s (e.g., time to/from midnight and/or midday and or sunrise and or sunset, or cosine and sine transforms of time/24 hours). Either way, time will always have the same value across all steps in a cluster, and hence can only be inferred on as an interaction with step attributes that do vary within a cluster. 

The `time_of_day()` function in `amt` is a wrapper for `maptools` functions `maptools::sunriset()` and `maptools::crepuscule()`. See `?maptools::crepuscule` for more details on how that works.

## The extract_covariates is great. Is it possible to use this function with time-dependent covariates?

Yes, see `?amt::extract_covariates_var_time`. 

Note that your `Raster*` object will need to have a z-column (the time stamp). See `?raster::setZ`.

## A technical question, I tried saving my cleaned steps to use later in a different script, but didn't know how best to do it. I tried save_rds(steps) and then load_rds(steps), which didn't work.

I generally would recommend this in a workflow, and in fact, I do it myself. However, I use base R's `saveRDS()` and `readRDS()` without issue.

I am not familiar with `save_rds()` -- perhaps you are referring to `futures::save_rds()`? You may also be referring to `readr::write_rds()`, but in that case, I also have no problem reading and writing a `steps_xyt` object.

My general recommendation would be to make sure all of your packages are up-to-date and fall back to `base::saveRDS()` and `base::readRDS()` if you still have issues.

## Is there a way of incorporating random intercepts and slopes in your framework? Or would you always estimate models for each individual separately?
A: It is theoretically possible, yet technically challenging, to fit mixed-effects conditional logistic-regression models (e.g., [Duchesne et al. 2010](https://doi.org/10.1111/j.1365-2656.2010.01670.x), [Craiu et al. 2011](https://doi.org/10.1198/jcgs.2011.09189)). 

An alternative and arguably simpler approach, developed and demonstrated by [Muff et al. 2020](https://doi.org/10.1111/1365-2656.13087), is to cast the iSSA as a Poisson regression with cluster-specific random effects.   

## Is there a way to focus on the final destination of the track instead of running iSSA all along the trajectory? i.e. focusing on last locations that could be territory or settlement habitat selection? it might be a question of changing of movement (time or spatial) scale?

There are two ways to incorporate a bias towards (or away from) a point attractor (e.g. a den), even if that attractor itself is moving (e.g. another animal): 

- Given a hypothetical preferred direction, each step can be characterized by its angular deviation from this preferred direction, and the coefficient associated with the cosine of this angular deviation is an estimator of the concentration parameter of the corresponding von Mises distribution (the larger it is, the stronger the bias; [Duchesne et al. 2015](https://doi.org/10.1371/journal.pone.0122947)). If the attractor itself moves, the preferred direction and the angular deviations form it must be updated on step-by-step basis.

- Include a raster of ‘distance to attractor’ as one of the habitat covariates (affecting the selection of the step’s endpoint). Standard transformations of this covariate can be used to represent various biological processes, including ‘home ranging’ (weaken attraction when near the attractor; `x^2`), and spatially decaying attraction (so selection for proximity or avoidance of the target is stronger near the target; `ln(x)` or  `x^(-1)`).  

## I usually scale continuous habitat covariates. Would I scale movement parameters also (ta_ and sl_)?

Unless that is the only way to avoid model failure, it is highly recommended not to scale movement covariates, as that would make the interpretation of the parameters very challenging. 

## What range of tolerance would be appropriate for resample_track( ). Could you tell us a bit more about the impact irregularity in sampling frequency will have on the outcome and how to choose the tolerance?

As a rule of thumb, do not use temporal tolerance exceeding 20% of the common step duration in the data. So, for data that were supposed to have a fix-rate of τ minutes, but in reality ended up having some fixes taken more frequently and others taken less frequently, tolerance < 0.2τ. It is the user’s responsibility to make sure, before fitting the model, that all steps included in the model fit could be reasonably considered as having the same step duration (there will be no warning generated if you try to fit a model with very different step durations, but you should not trust the results). The issue of accounting for step-duration within the iSSA is, as far as we know, not resolved yet, mainly because it is unclear how should each of the kernels scale with time (but see [Munden et al. 2021](https://doi.org/10.1111/2041-210X.13574) for some exciting progress).     

## If you are interested in looking at the effect of a specific categorical habitat variable, should you limit steps included in the analysis only to matched case-control clusters that have at least one location within the habitat?

Not necessarily, you can still make inference even if not all clusters include variation in your focal covariate (whether it is categorical or continuous). We do recommend however having a sufficient number of clusters in the data to enable reliable inference (at the very least, more than 10 per degree of freedom).

Note that [Fortin et al. 2005](https://doi.org/10.1890/04-0953) used 200 available steps for each observed step to capture a rare categorical habitat variable.

## Is there a way to spatially constrain where available steps may be generated (i.e., to prevent available steps ending in a building when modeling urban species)?

Yes, ‘thinning’ - sample available steps as recommended, but reject an available step and sample again if it ends in a habitat the species cannot possibly occur in.  

## Fieberg et al 2021 recommends using infinite weighting when using RSF and SSF? but it looks like it is decreasing deviance explained of the models, apparently based on how weight is implemented in the deviance formula? 

From the Zoom chat:

> I recently had a conversation with John Fieberg and Stepahanie Muff regarding their 2020 (i think) manuscript infinite weighting is only done for RSF NOT for SSF.  The way it is written doesn't make that distinction clear. We never really got to the bottom of why but this is the quote from John: "Just a quick follow up - I added weights to the SSF using the fisher data and then tried different values for the fixed variance.  This led to really erratic behavior and increasing SEs. So, this: 1) reinforces that you should *not* using weights when fitting SSFs; and 2) likely explains the strange results you were getting Cory.  Hopefully, removing the weights argument solves the problem."

As far as we know, weights should not be used in iSSA, only in exponential HSA.

## Do you foresee any issues in amt with the emerging "terra" package for raster data?

It's hard to say, but `terra` is designed to be very similar to `raster`, so hopefully not. At some point there will likely be a full switch from `raster` to `terra` in `amt`, but I don't think Johannes Signer (`amt`'s author/creator) has a timeline for that yet. You can stay up-to-date with the latest updates on [GitHub](https://github.com/jmsigner/amt).

## Why is it necessary to include both log step length and step length in the model for the purpose of adjusting it ? and Why are these models biased?

The underlying mechanistic movement model is a Biased Correlated Random Walk, because it incorporates bias towards or away from spatial positions based on the local environmental conditions. A step-selection analysis that does not include any step-length covariates leads to biased inference about habitat selection strength (demonstrated first by [Forester, Im, and Rathous 2009](https://doi.org/10.1890/08-0874.1)). It is not necessary to include both step length (`sl_`) and the natural logarithm of step length (`log_sl_`) in the model. Doing so is only adequate if available steps were sampled from a Gamma distribution, and the user wants to adjust both the shape and scale parameters of this distribution. Otherwise, including only `sl_` as a covariate assumes the Gamma tentative shape parameter is the true ‘selection-free’ value, while including only `log_sl_` assumes the Gamma tentative scale parameter is the true ‘selection-free’ value. 

## do people usually/ever check some kind of predictive performance holding out x-number of steps (e.g., at the end of the individual's trajectory) fitting the model to the not-held-out data and predicting the holdout step trajectory to see how closely the predicted end-steps match the actual end-steps

Not often, mostly because this is fairly technically difficult, and theoretically questionable, to do using `survival::predict.coxph()`, but also because, while parameter estimate accuracy should increase with the number of available steps per cluster, correctly pointing out the used step is less likely the more available step there are (regardless of how close the model is to the truth). It is possible to calculate a likelihood-based pseudo R^2^ or concordance criteria (the probability that a used step is ranked higher than an available step; within-cluster ROC AUC) for out-of-sample data, but this is again rather technically difficult. We recommend Used-Habitat Calibration Plots ([Fieberg et al. 2018](https://doi.org/10.1111/ecog.03123)), or simulation-based validation, both soon to be implemented in amt.     

## Is it possible with amt to simulate tracks from iSSA models based on empirical data? and or based on specific distribution parameters chosen by the user?

This should be possible soon using the iSSF simulator in `amt`. For an idea of how this works, see the script `extra/sim_tracks.R` in this repository.

## Is it straightforward to add covariates collected along the expected path rather than just at the start and the end locations of the step? (when time intervals are long)
A: Yes (see function `amt::extract_covariates_along`), but this is not necessarily recommended. Keep in mind that the underlying mechanistic model is a positional jump process, where the ‘jump’ implies the animal does not interact with the environment along the step. Another way to think about it is that the straight-line step is a useful heuristic construct, but real animal do not travel along straight lines and turn at regular time intervals.  

## Is it possible to implement non-linear relationship with iSSA? or is possible to run conditional logistic regression within Generalized Additive Models framework?

If by "non-linear relationship"" you mean including variable transformations (*e.g.*, `x^2`, `ln(x)`, `x^(-1)`, `cos(x)`, *etc.*) in the model formula, then obviously yes. The basic link function is always the natural logarithm.

## Is there a way to incorporate vertical movement (e.g., dive depth) into the iSSA?

We see no reason why not, but… It could be incorporated as the depth at the step’s endpoint (either in absolute terms, or in relation to the depth at the start point), which would make it part of the ‘movement-free habitat-selection kernel’. This should be straightforward. It could (and ideally should) also be incorporated into the ‘selection-free movement kernel’ which will then become 3D, meaning that turn angles should be sampled from the spherical (rather than a circular) von Mises–Fisher distribution, and step lengths are calculated in 3D (x,y,z).   

## For attributes can you use fuzzy habitat attributes or some composite around the end of the step?

Yes

## Can iSSA be run with parallel processing? Can it be done in such a way that it will combine the results in order?

In principle, there should be no problem maximizing the conditional likelihood across multiple computational threads in parallel; the order of the clusters is of no significance. That said, `clogit()` calls `survival::coxph()` to maximize the likelihood, and it is hard to imagine a situation where parallel processing will be necessary, as `coxph()` is extremely efficient.

