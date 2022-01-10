# Webinar Q & A

## From Poll EV

### Fitting a time-series is great, but what do the basis coefficients tell us about the data? Specifically in the temp example, what do the fitted coefficients tell us about the change in temperature over time?

In the HadCRUT example, the fitted coefficients themselves don't tell us how temperature changes over time, but when we use these to weight the basis functions to give the fitted spline, we have

1. an estimate of how the mean (or expected) temperature anomaly changes over time,
2. a fitted model with which we can assess over aspects of the changing temperature.

As an example of 2., we can compute the derivative of the fitted smooth trend to get rates of change in temperature at any point in the series. This would be done in practice with the method of finite differences using components of the fitted mdoel, and for individual smooths in a `gam()` can be returned using the `derivatives()` function from my {gratia} package. If you want derivatives of the response on the response scale then you'll need to wait until I code that specific option for {gratia}, which should hopefully be available in version 0.7.0.

I also have several blog posts on this topic [here](https://fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/), [here](https://fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/), and again [here](https://fromthebottomoftheheap.net/2017/03/21/simultaneous-intervals-for-derivatives-of-smooths/). Note the warning at the top of the second post(!) and also note that `fderiv()` has been superceeded by `derivatives()` now.

See also a more general question along these line and response below in the *From chat* section.

### The first example examined temperature over time. Is there a way to address and evaluate temporal autocorrelation in GAMs?

Yes, for Gaussian responses mainly. `gamm()` allows one to specify the `correlation` argument of `nlme::lme()` which allows for several types of ARMA stochastic trend components in the model residuals. See blog posts [here](https://fromthebottomoftheheap.net/2011/07/21/smoothing-temporally-correlated-data/) and [here](https://fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/) for examples of doing this.

An alternative is the found in the `mgcv::bam()` function which allows estimation of mdoels with AR(1) processes in the residuals for models with `family = gaussian()` only.

For non-Gaussian responses/models, this is more difficult as these typically don't have separable mean and covariance terms. If you can stomach the PQL approach to fitting GLMMs then you can use `gamm()` with `family` set to one of the families understood by `glm()` whereupon fitting will now take place using `MASS::glmmPQL()` which iteratively calls `nlme::lme()` to fit the model. Hopefully you can see that iterating `gamm()` <-> `glmmPQL()` <-> `nlme:lme()` is a lot more complex an estimation procedure and potentially less tolerable of problems with data and model issues. You should also read Bolker et al (2009) before you embark on using PQL to estimate a GAMM/GLMM.

Better alternatives would require fitting in a more flexible language where model specification. {brms} for example would allow for ARMA terms in the mean part of a model with a non-Gaussian response. And there are other options as well, although the downside here is that one needs to specify the model by hand and this typically requires deeper knowledge than that required for fitting mdoel with R's formula interface.

More generally, as I have summarised (Simpson 2018), we can run into problems with models like this where there is an identifiability problem in that the model can't distinguish wiggly trends and high/strong autocorrelation. Inevitably something has to give and the model will typically throw everything in on *either* a wiggly trend and no autocorrelation *or* a very smooth (often linear) trend with very strong autocorrelation (for an AR(1), rho > 0.9 for example). This is really only reflecting that strong autocorrelation is really just a very wiggly trend, mathematically speaking. So when modelling time series data, we might simply make `k` as large as we can an accept the wiggly trend.

If we are talking about situations where the response may have been collected over time and thus show autocorrelation but we are estimating smooth function of covariates that aren't time, then we do need to be mindful of autocorrelation as I discuss [here](https://fromthebottomoftheheap.net/2011/07/21/smoothing-temporally-correlated-data/).

Bolker, B.M., Brooks, M.E., Clark, C.J., Geange, S.W., Poulsen, J.R., Stevens, M.H.H., White, J.-S.S., 2009. Generalized linear mixed models: a practical guide for ecology and evolution. Trends Ecol. Evol. 24, 127–135. https://doi.org/10.1016/j.tree.2008.10.008

Simpson, G.L., 2018. Modelling Palaeoecological Time Series Using Generalised Additive Models. Frontiers in Ecology and Evolution 6, 149. https://doi.org/10.3389/fevo.2018.00149

### How do you assess statistically significant predictors in a GAM model, especially when the sample size is very large (with large n, most predictors end up being significant)?

This is no different to how you'd do it other statistical models, although I would recommend to stop worrying about *p* values and instead focus on effect sizes and the shape and magnitudes of estimated smooth functions. Information criteria are no solution as they will typically use less stringent *p* value thresholds than those a human might as thus make the problem worse.

If one wants *p* values then `summary()` will provide them, but if your *n* is so large as to render even small effects "statistically significant" then one should probably just ignore the *p* value and focus on the effect sizes of the estimated parameter or smooth functions.

### Biologically speaking, how do I interpret GAM results when the model is considered significant?

You interpret them as you would for a (G)LM or (G)LMM, except that in the GAM case we look at the estimated function, not a single coefficient. One should look at the estimated function to see how it increases or decreases over the range of the function; does it reach an assymptote where there is no further change in the effect of the function on the response (it is flat), or perhaps the curve is unimodal, indicating the response is maximised around some value of the covariate.

The biological interpretation is up to you, the biologist. Just focus on the estimated functions and look where they are not flat, which would indicate that the response isn't changing as you vary the covariate.

An underused approach to help you understand what the model has encoded, is to predict from your model for specific data scenarios or slices through your data. This involves generating new data at which to predict values of the response at certain combinations of covariates. Typically one would vary one or two covariates only, while holding the other covariates constant at typical or chosen values.

### Is there a minimum sample size to be able to fit GAMs?

Yes, but it's not that much higher than the minimum sample size you'd need to fit a GLM. A more pertinent question is do you have enough unique values of the covariates you want to use represent as smooth functions to be able to estimate them? Practical problems arise when you have too few unique values; in which case you need to reduce `k` to be equal to or less than the number of unique values for a specific covariate.

Where you want to include multiple smooth terms in a model, you are going to need at least as many observations as you have basis functions over all your smooths, plus one for the intercept, plus one for good luck. However, trying to estimate as many parameters as you have data is asking for trouble and unlikely to work unless the true functions are very simple (in which case you could have gotten away with smaller values for your `k`s) or they have huge effects on the response (in which case one could reasonably as whether you needed to fit the model at all). More generally then, you are going to multiple observations per basis function used across your multiple smooths in order to get reliable estimates.

Other than that, I can't really reply with much more than "How long is a piece of string?" :-)

### Do GAMS work on timeseries with large datagaps? Will those gaps be interpolated linearly?

Yes they can work with large data gaps and the model will estimate smooth functions that interpolate over the gaps. I answered this question in the webinar so I go into more detail there, but I also mentioned [this blog post](https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/) that I wrote exploring the issue a little (but mainly focusing on extrapolation). The smooths don't necessarily interpolate linearly; they'll smooth across the gap in such a way as to minimise the penalty. This need not imply linear interpolation. But for the default setting, the interpolation will be approximately linear.

### What exactly are the assumptions for GAMs that you are looking for with the diagnostic plots? Are they the same as those for linear regression models?

These are really the same as for other GLMs or GLMMs. We're really checking the distributional assumptions of the models we have fitted, and looking for problems with the fits. Model checking probably deserves it's own webinar so I won't go into much detail here, other than to say I did answer this a bit in the webinar and I would also suggest that you look at the DHARMa package and rootograms as two other diagnostic tools for models like GAMs.

### How would you fit multiple seasonal components using GAM ?

Depends on what you mean by "seasonal"? If you meant multiple cyclical components, such as variation within a day, within a week, within a year, then I would model each as a cyclic smooth using `s(x, bs = "cc")` and be sure to set appropriate knot values for each covariate using the `knots` argument.

### When comparing the performance of different GAMs, which metric reported in mgcv (adjusted R-squared or deviance explained) is the most appropriate to use (and in what instances)?

If you mean comparing simpler models with more complex ones fitted to the same data, then if I were going to do that I would probably use `AIC()` which in newer versions of {mgcv} adjusted for us having estimated (selected) the smoothness parameters. Often however, where I have a set of covariates that represent the known or hypothesised relationships I wish to test, I would just fit the full model and see what the model comes up with (while appropriately setting `k` for each term, etc). Failing that, I would use do the same but add `select = TRUE` to my models, to do feature selection using extra penalisation, which can shrink terms out of the model entirely (a bit like the LASSO penalty can).

Adjusted R^2 and deviance explained tell you something about individual models but they aren't a good way to compare models unless you really want to compare them on one or both of those very specific metrics and I'm not sure when and where that would be a sensible thing to do --- but I am open to arguments in their favour.

### Are there ways of dealing with datasets in which the collection interval was not constant? For example, annual for part of the record but every 5 years for the rest of the record.

If we're talking about in a similar setting to the arthropod example, then GAMs can work well in those settings. Although I simplified the example to remove the forest plots as there were sampled at lower frequency than the the grassland plots, I could have modelled everything at the same time. If your data start out with one frequency and then at change to a different frequency, this is also OK.

The main issue to be aware of is that where you have fewer observations the effects (trends, say) that you can estimate will typically be less complex as there is less data to inform more complex trends. Also, the estimates in the less-frequently sampled periods or types of sites (forests versus gransslands) will be less precise (more uncertain) all else equal than in the more frequently samples periods or types of site. This will limit the shapes of the functions you will be able to estimate too.

A practical point to be aware of is that you will likely want to specify the knot locations for cubic splines or other knot-based splines in the case where you start collecting data at one frequency and then at some later time change to a different frequency. You'll want to avoid having basis functions not supported by any data, so you might want to override the default knot placement (which spreads them regularly and evenly over the range of the covariate), and put the knots where you have data. You don't need to do this of course, but it might help you estimate more complex functions where you have the data to estimate them rather than wasting basis functions where you have fewer data.

### How might GAMs be useful in near-term iterative forecasting when they rely on existing data so heavily? They have vast capabilities for fitting existing data, but what about prediction? 

Forecasting from GAMs where one has used smooths for the long term trend is going to be problematic for all but short extrapolation periods. I have covered some of the issues in [this blog post](https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/) The issue goes away for example if the trend is not modelled with a smooth, but via say lag-1 (AR(1)) terms say where we include the previous value of the response as a covariate. But this necessarily limits GAMs to situations where you have evenly-spaced data (without getting more complex).

I mentioned in the webinar that Gaussian Processes might be a better option where extrapoltion is required. I'm not as familiar with that literature, but the preprint paper that inspired the blog post referenced in the paragraph above is one example that might get you into the broader literature.

If by prediction you mean observations not included in the training data, but where values lie within the observed ranges of the training data, then this isn't really a problem for GAMs - in fact predictions for data hypotheticals is a good way to understand what you model is telling you about your data and system. Where it might be problematic is where you are asking the model to predict beyond the range of the combinations of the covariates observed in the training data. An example here would bea 2D smoother, say in space where you aren't extrapolating in either spatial covariate alone, but the combination of spatial covariates you want to predict at lies outside the domain of support of the observed data. Then we're back into an extrapolation situation as mentioned above.

### Since we are getting several coefficient estimated and smooth term in the summary output, how can you write the results for gam when working on a manuscript?

The output from `summary()` gives a single test for each smooth in the model against the null hypothesis that the indicated smooth is a flat constant function with value 0. The information in the `summary()` output could be turned into a table for the smooth terms and another for the parametric terms.

I would also suggest including the partial effects of the smooths, as produced by `plot.gam()` or `draw()`, but these are of lesser importance in some cases, especially that of more complex models with "main effect" smooths and "tensor interation smooths" of the form `y ~ s(x1) + s(x2) + s(x3) + ti(x1, x2) + ti(x1, x3)` where it isn't as easy to visualise the overall effect of `x1` and `x2` on `y` as this is spread other four separate smooth functions in the model. In this case, it is better to create new data evenly over the space of `x1` and `x2` while keeping `x3` constant at some typical or reference value of interest, and then predict from the model at those new data values/combinations and visualise the fitted surface in your favourite plotting software.

### Do you think GAMs (mathematical) complexity has the potential to hinder comparisons with other studies and/or the development of our understanding of the processes underlying trends; in a way that linear models, which are more readily described mathematically, do not?

No I do not. It might complicate comparison with models fitted using linear functions only, but in many instances linear effects are often used without any justification too. And if the estmiated effects you see in your data / system are robust then you're not loosing anything by fitting them with a GAM and comparing smooth functions with linear effects estimated elsewhere can be informative, but you have to assume that the linear effects aren't linear simply because the previous analyst didn't know about GAMs or choose to fit them instead.

I went off on a bit of rant in the webinar about this, so go listen to that if you're interested. I talked about the diversity setting in particular, where we have moved beyond asking whether or not there is a trend in biodiversity time series (although fights still happen about that too) to asking more nuanced questions such as what is the rate of change etc. And if the trend is non-linear, summarising it by a linear function might answer the "is there a trend" question (it also just as well might not give the right answer, depending on the trend), but it will get the "what is the rate of change?" answer wrong a lot of the time. For an example of that, compare the linear trend fitted to the HadCRUT temperature data with the estimate smooth function (or even just the data) and you'll see that it gets this important estimate wrong through almost the entire time series.

I do accept that non-linear trends make explaining change like this more complex, but that's a problem on us to solve through better communication. We don't solve that by fitting demonstrably poor models just because they are easier to interpret.

## From chat

### REML smoothness selection & collinearity

> I notice you are using REML. To what extent are we fitting these splines to the "residual" after taking out fixed effects - i.e., it functions like a random effect. I ask as I’m curious about how GAMs handle collinearity with other fixed linear predictors. In mixed models, we assume that the random effects don’t correlate with fixed predictors, versus in, say, multiple linear regression where our estimates of coefficients control for the effects of other covariates - and it’s OK as long as collinearity isn’t too too high.

REML smoothness selection methods cast the GAM estimation problem as a mixed effects / hierarchical model in which the wiggly parts of the spline basis are incorporated into the random effects model matrix **Z**, while the perfectly smooth parts of the basis, the functions in the so-called null space of the penalty, are incoporated into the fixed effects model matrix **X**. The smoothness parameter(s) are proportional to the inverse of the variance parameters for the radom effects; hence the wiggliness of the smooth(s) is estimated by estimating the variance parameters of the random effects. Once the basis expansions have been performed and the relevant model matrices have been formed, fitting can be done in standard LMM/GLMM software as if the model were a true mixed effects / hierarchial model. In practice however, {mgcv} employs some sophisticated estimation algorithms to fit the GAM efficiently. However, where we want GAMMs with many random effects terms or many levels (subjects) in each, the algorithms for fitting GAMs become inefficient as they do not currently exploit the sparsity of the random effects model matrix. In those cases you can have {mgcv} or {gamm4} do the neccessary conversion to the GLMM form and have the model fitted using {nlme} or {lme4} respectively.

If your concern relates to the theory of GLMMs and parameter estimates, then we might be less concerned about this in the GAM setting as we are only using the mixed effects representation as a nice computational trick to fit the model and estimate the functions. However, ...

...In GAMs we are concerned with multicollinearity as in other model settings, but we also have the additional problem of non-linear correlations --- something called concurvity, "the presence of covariates which are themselves well modelled as smooth functions of other covariates" (Wood 2008). The methods in {mgcv} have largely been developed to work in the presence of significant concurvity: Wood (2008) discusses a PIRLS algorithm that is robust to concurvity problems which is the basis for the GCV/AIC-based smoothness slection in {mgcv} and Wood (2011) extends this work to the REML and ML cases of smoothness selection for GAMs as we in applied setting prefer these over GCV/AIC based approaches.

So yes, while we do still have to be careful when fitting in the presence of significant concurvity, the methods in {mgcv} have been developed with this in mind and have been shown to work in such cases, for example in the two papers of Wood (2008, 2011). You can check the concurvity of model covariates using `mgcv::concurvity()`.

If one is fitting with `gamm()` or `gamm4()`, I am less certain about the performance of those algorithms as they were not designed with GAMs and smoothness selection in mind.

Wood, S.N., 2008. Fast stable direct fitting and smoothness selection for generalized additive models. J. R. Stat. Soc. Series B Stat. Methodol. 70, 495–518. https://doi.org/10.1111/j.1467-9868.2007.00646.x

Wood, S.N., 2011. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. J. R. Stat. Soc. Series B Stat. Methodol. 73, 3–36. https://doi.org/10.1111/j.1467-9868.2010.00749.x

### Is there a reason you switch between the tweedie and negative binomial families?

If you meant while fitting `m2` where I used `family = tw()`, that was a typo arising from initially preparing examples for arthropod biomass. I have corrected this in the file `arthropod.R` and I'll rerun this and update it again with comments if it would change anything I said.

If you meant while fitting with `family = twlss()`, then this was intentional; the point here was to try to allow different plots to have different variances (dispersions) as we cannot currently model the extra dispersion parameter of the negative binomial distribution in {mgcv}.

### How can I present the real response values on the y axis (not the effect) using gratia? In plot.gam I found “shift” to extract the model intercept

You can use the `constant` and `fun` arguments to `draw.gam()`, with the former corresponding to the `shift` argument of `plot.gam()` for example. The inverse of the link function (or some intermediary, such as using `fun = exp` in a binomial GAM to go to the odds scale) can be passed to `fun` to apply a transformation function.

That said, I find it more useful to use `predict.gam()` and include or exclude terms using the `terms` or `exclude` arguments or to just predict for data slices wherein one (or two) variables are varied across their range while all other predictors are held constant at some typical values.

### Are GAMs useful only to fit any response variable against time or where time is involved?

No, GAMs are useful where one might want to use a linear model, a GLM, or a (G)LMM. As I showed towards the end of the arthropod example, we can include smooth functions of any type of covariates into our models. One typically can use lower values of `k` (the upper limit on the size of the basis) in those situations, unless you expect very complex and non-linear relationships between a covariate and the response.

### How should one interprete these basis functions? Do they represent algebraic functions?

For most of the spline types in {mgcv} there isn't a straightforward interpretation of the basis functions themselves. And I wouldn't worry about that; they are a means to an end - a way to flexibly represent the unknown function *f*() in our models.

### What if missing some days of the year or some months?

I will assume this is in relation to the brief answer I gave on modelling seasonal variation. Missing days or months is fine, although in the latter, it will make it more difficult to use a spline if you say only collected samples quarterly in the same four months each year of sampling or at each location, as you'll need to set `k` quite low in that instance (`k = 4` is the largest value you could use in that case, but I'd have to check). That said, if you have large gaps with no data to inform the seasonal smooth you should expect wider uncertainty bands and thence reduced power to detect the seasonal effect if present.

### How much does non constant variance matter with these GAMs and is there a way to account for this in mcgv?

As with GLMs, deviation from the implied mean-variance relationship --- constant variance in the `gaussian()` case, variance = mean in the `poisson()` case etc --- is a concern for inference in GAMs. As I mentioned in the webinar, GAMs are, once we have created the basis expansion, just fancy GLMs, exactly so if we ignore the smoothness selection, so the same considerations apply.

A way to handle departures from the fixed mean-variance relationship implied by the `family` used when fitting is to fit a distributional model (these are also called GAMLSS --- GAMs for Location, Scale, Shape --- in the GAM setting), where we have linear predictors for all parameters of the conditional distribution of the response. In the case of `family = gaussian()` for example {mgcv} has the `gaulss()` family which allows one to model both the mean and the variance of response, which is a way to model the non-constant variance in the Gaussian setting. {mgcv} has added quite a few LSS families now, but they don't cover all the standard distributions in {mgcv}. The {gamlss} and {bamlss} packages have a wider range of LSS families available, and the {brms} package allows modelling of distributional parameters too using {mgcv}-like syntax if you can't find something appropriate within {mgcv} itself.