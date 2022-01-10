# Webinar Q & A

## From Poll EV

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