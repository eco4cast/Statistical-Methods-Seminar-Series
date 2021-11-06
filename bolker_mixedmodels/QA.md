---
title: 'Questions and answers'
---
 
1 .  What are your thoughts on the pros and cons of fitting these types of models in a Bayesian vs MLE framework? 

 **Answer**: 
 If difficulty were no object (in terms of (1) complexity of implementing the model; (2) computational intensity of fitting the model; (3) troubleshooting; (4) thoughtfully choosing priors; (5) convincing reviewers etc.), I would probably default to Bayesian approaches most of the time. Bayesian statisticians would say that (3) and (4) are things you ought to be doing anyway for scientific reasons (the requirement for troubleshooting often, although not always, comes from models that are subtly poorly posed or unidentifiable). They would say that (1), (2), and (5) are becoming progressively less problematic. 

2 .  What are the best practices for when you have covariates measured at different scales? e.g., if you had some predictors measured at the biome level and some measured at the ecoregion level? 

 **Answer**: 
 For the most part, modern mixed-model frameworks just take care of this for you. Unlike the classical ANOVA framework, you don’t have to memorize a lot of formulas for crossed, nested, split-plot, split-split-plot, Latin squares, etc.. On the other hand, you do have to think about which random effects terms to include (i.e., which terms vary within groups at which levels and thus can be included in the maximal model?) 

3 .  Is it valid to use AICc (or other IC) to choose between allowing variation of intercepts only (1|g) vs also allowing slopes to vary (f|g)? 

 **Answer**: 
 Maybe; see [here](http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#can-i-use-aic-for-mixed-models-how-do-i-count-the-number-of-degrees-of-freedom-for-a-random-effect for the discussion of some of the issues (counting parameters; boundary effects). I don’t know that anyone has really looked at how reliable AICc is in)  the mixed-model context (lots of people have incorrectly assessed AIC based on how reliably it selects the “correct” model; the right way to assess AIC is to whether it correctly chooses the model with the minimum Kullback-Leibler distance – see Richards (2005, Ecology) for an example, in a case where the true model is maximally complex). Because of boundary effects AIC is likely to be conservative (i.e. telling you a random effect is not useful for prediction when it really is); AICc and BIC will both be  even more conservative. So I would probably stick with AIC, and maybe even use a penalty term less than 2 (i.e. $k<2$ in $-2 \log(L) + k p$).  

4 .  Hi Ben, should we always try to present mixed models with e.g. marginal and conditional R2 values? What are the arguments against this (your 'double bend' roadsigns)? Are there any useful alternatives that can be presented to approximate goodness of fit for mixed models? Thanks! 

 **Answer**: 
 I think marginal and conditional R^2 values (aka the Schielzeth/Nakagawa/Johnson framework) is probably pretty good. I’ve mostly seen it give weird answers in cases where the data were small and noisy otherwise. When using R^2 values for mixed models and especially for GLMMs it’s worth noting that the possibilities have proliferated because, beyond the simple linear model framework i.e. (G)LM(M)s, it’s impossible to specify a single metric that shares all of the nice properties of the simple R^2. My caution in talk was specifically about the partial R^2 values, which I took from a series of papers by Edwards, Jaeger et al (DOIs: 10.1002/sim.3429, 10.1080/02664763.2016.1193725). I had to choose between a Kenward-Roger and a standardized generalized variance (SGV) approach; K-R seemed more principled, but in some cases I got termwise partial R^2 values that were larger than the computed R^2 for the full model. I spoke to Byron Jaeger about this, he said he wasn’t sure what was going on (but had seen similar phenomena), it remains to be explored. See also http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#how-do-i-compute-a-coefficient-of-determination-r2-or-an-analogue-for-glmms 

5 .  when to choose a variation over the intercept or a variation over the slopes in the random part? 

 **Answer**: 
 Is this in the context of model simplification? In general if the design allows it and you have sufficient data to estimate random intercepts and slopes, you should estimate both. Dropping the random slopes seems like a reasonable possible avenue of simplification. Barr (2013) talks about the alternative possibility of dropping the intercepts as well, but in most cases this doesn’t make much sense to me (to be honest I don’t remember their justification very well); you’d have to have a situation where it was plausible that the groups had varying slopes but all had the same intercept … 

6 .  is there any specific reason for choosing log scaling over other scaling methods such as Z-score scaling for NPP data? 

 **Answer**: 
 In this analysis we actually did both (we log scaled, then Z-scaled). The two types of scaling do quite different things. Log-scaling means that parameters measure changes on a proportional/multiplicative rather than an additive scale; it changes the predicted values, the fit of the model, and the distribution of the residuals as well as the units of the parameters.  Z-scaling has no effect on the overall fit of the model (predictions, log-likelihood, AIC, etc.); it is made up of centering, which changes the interpretation of intercepts and of main effects in a model with interactions; and scaling, which changes the units of the parameters (making them unitless and hence comparable). See Schielzeth 2010 for more. 

7 .  every time I apply the mixed-effects model  I receive the isSingular comment, how to cope with a singular fit ? 

 **Answer**: 
 You probably have a small, and possibly noisy, data set, and not very many groups. This is discussed [here]( https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-models-random-effect-variances-estimated-as-zero-or-correlations-estimated-as---1) and at `?lme4::isSingular`. You can leave the model as is; drop or simplify the random effect term that is causing the problem (this may mean that you no longer have any random effects at all, i.e. you now have a (G)LM instead of a (G)LMM); or convert the random effect to a fixed effect (the latter is probably the most conservative). 

8 .  Can you expand on the notion that Mixed Effect Models should not be used to test differences between groups? 

 **Answer**: 
 From a theoretical point of view, random effect deviations from the population-average parameter aren’t estimates in the formal sense (they don’t have a true value, sampling distributions, etc.); they are random variables. All of the frequentist statistical theory for working with parameters goes out the window. Another way of looking at this is that not being able to test difference 

9 .  If you wanted to examine variation among biomes on slopes (i.e, effect of biome on slope), how would you proceed given that you said hypothesis testing is not valid for random effects? 

 **Answer**: 
 I would just examine it, most typically with a “caterpillar plot” – a plot of the conditional modes, with ± 2 conditional-SD intervals, ordered by their values (this looks a little like a caterpillar). I can make perfectly reasonable statements about which biomes are (e.g.) very sensitive to fire, which are not very sensitive, whether the apparent signs differ for some biomes (i.e. richness increases with increasing fire in some biomes, decreases in others), and use the intervals to assess my certainty. I just can’t get p-values. (Note that while I can easily get intervals for the deviations from the population mean slope, getting intervals for the biome-level slopes [i.e. pop mean + deviation] presents some challenges: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#confidence-intervals-on-conditional-meansblupsrandom-effects ) 

10 .  "log response & predictor means coefficients interpreted as elasticities" Could you expand on just what this means? 

 **Answer**: 
 When predictors or response variables are logged it means that we are looking at effects of proportional changes, e.g. a change of 0.01 in log-NPP means (approximately) a 1% change in NPP (delta-NPP = 0.01*NPP). An elasticity represents a “proportional change caused by a proportional change”, e.g. an elasticity of 0.5 means “a 1% change in NPP causes a 0.5% change in richness”, which is exactly what we get if we do a log-log regression. 

11 .  I have a question. Why are you log transforming the variables? 

 **Answer**: 
 Log transforms have many nice properties: they often make a regression curve more linear; they make it impossible to get negative predicted values of the richness; they often improve the fit of the data to the assumptions (i.e., they take care of certain forms of heteroscedasticity, when the standard deviation is approximately proportional to the mean). The parameters have reasonably nice interpretations in terms of elasticities (see previous answer). If log(y) = a + b*log(x), then the x-y relationship is a power law (y = exp(a)*x^b), which also has some nice biological interpretations. 

12 .  Should you always take coefficients from the full model, or use model averaging via AIC, or from the best model estimated by model selection? 

 **Answer**: 
 You should never do inference (i.e. compute Cis/p-values) on coefficients estimated from a model that was selected in a data-driven way. If you want to do standard statistical analyses (p-values/Cis/etc.), use the full model (with non-data-driven, a priori choices to reduce the model to a sensible size: see F. Harrell Regression Modeling Strategies). If you are only interested in prediction, and don’t mind confidence intervals that are probably too optimistic/narrow, then model averaging is OK, although regularization methods (ridge/lasso/etc.) are better (more efficient with better theoretical grounding). 

13 .  When using mixed models as a tool for "(2) estimating the variability among groups", what statistics should you report? Relative measures such as variance partition coefficient become hard to interpret for complicated RE structures. 

 **Answer**: 
 You can still report the variances, and the relative magnitudes of the different variances, even if partitioning gets dicey. 

14 .  Does it make sense if we regress the condition modes of random terms (estimation of each group in `(1 | g)`) with another variable about ` g`? For example, (1|species) then regress with species traits? 

 **Answer**: 
 As far as I understand the question, this is a bad idea. See Hadfield 2010 “The misuse of BLUP in ecology and evolution”. I **think** that the best way to do this is to put the species traits into the model as fixed effects (`lme4` doesn’t care that some covariates vary at the individual level and others vary at the species level); then the species effects will be decomposed into a component that’s explainable by the covariates and a residual component that’s not ... 

15 .  Is it valid to use means separation using `emmeans` with `glmer` in `lme4`? 

 **Answer**: 
 I’m not sure what “means separation” means. Emmeans generally works fine for GLMMs; the only tricky bit for GLMMs with non-identity links is that it can generate biased predictions of the mean if you’re not careful (see the vignettes/documentation) 

16 .  why when including random slopes the CI is bigger than just intercept 

 **Answer**: 
 It’s often the case that adding predictors to the model will make the Cis of individual parameters larger (this is the basis of computations like the Variance Inflation Factor). The more stuff is in the model, the less certainly you can attribute the patterns to the effect of any one particular predictor. 

17 .  I didn't follow the model selection procedure you used. is there a function like `dredge` that allows you to compare all the possible models for fixed and random effects? 

 **Answer**: 
 The `dredge` function in the `MuMIn` package only (AFAIK) does selection on the fixed effects.  I would actually want the opposite (i.e. only selecting on the random effects). I explicitly do not want to do selection on the fixed effects, which is the focal component of my model (i.e. the part I’m really interested in); drawing conclusions after model selection messes up inference. There might be a function somewhere that can simultaneously do selection on both fixed and random components, but I wouldn’t recommend using it. 

18 .  What are your thoughts on using variance weighting (`varIdent` and other related functions?) 

 **Answer**: 
 It’s very useful if there’s heteroscedasticity of the level of the residuals that you can’t fix by transforming the response variable. It’s unfortunate that we still haven’t implemented this functionality in `lme4`. `glmmTMB` can do some, but not all of these models (it will allow the variance to depend on covariates, but not on the fitted value of the mean) 

19 .  What are your thoughts on including random effects using `s(variable, bs="re")`  in`mgcv`? 

 **Answer**: 
 Random effects via `mgcv` are great.  They’re less efficient for really large random effects (`mgcv` doesn’t use a sparse model matrix at this step), and I don’t know how flexible they are for random-slope models. `gamm4` is clunky, but uses sparse matrices for random effects and allows all of the `mgcv`-style `s()` terms in the model. 

20 .  When fitting polynomial terms in `lme4`, is it meaningful to include polynomial terms in the random terms too? e.g., `(1 + f^2 | g)` 

 **Answer**: 
 Yes. It hurts my brain a little bit to think of these models, but if you have measured multiple values of x in each level then there’s no reason to think that your polynomial effect (or spline, or whatever) doesn’t vary among groups, and to try to estimate that variation 

21 .  sorry, what was the point of using `gamm4`? 

 **Answer**: 
 To allow us to use a spherical spline term from mgcv() to incorporate spatial autocorrelation. 

22 .  the goal was not to explore the random factors in this model but can we use this way and set a model with random factors where our goal is to interpret the random factors or is not recomended? 

 **Answer**: 
 I don’t quite understand the question. The model I presented would be fine for exploring variation among biomes etc., but (1) I would probably need more data (2) I’d have to be more careful about doing model selection on the RE terms (so as not to mess up my inference). 

23 .  For what diagnostics should studentized/ standardized residuals vs. normal residuals be used? 

 **Answer**: 
 Not really a GLMM question … “standardized” = “studentized” residuals scale the  residuals by their leverage (potential influence). In general you should use standardized/studentized residuals when you’re considering heteroscedasticity (e.g. a scale-location plot) or influence. Looking at `?plot.lm` in base R, standardized residuals are used for the scale-location plot (assessing heteroscedasticity), residual-leverage (assessing influence), and Q-Q (assessing Normality). I may have been a little bit sloppy about this in my examples. [This StackOverflow Q](https://stackoverflow.com/questions/31976284/how-can-i-extract-studentized-residuals-from-mixed-model-lmer) shows how to use `residuals()` + `hatvalues()` to compute studentized resids … 

24 .  What do you think about use variables that violates the assumptions? 

 **Answer**: 
 It depends on how bad the violations are, and what kind of violations. For example, people worry a lot about non-Normality of the residuals/conditional distribution of the data, but that’s one of the least important problems (most linear models, LMMs, etc. are quite robust to non-Normality) … note that I didn’t even show a Q-Q plot in the presentation … I would generally try to address the violations if I can (by transformation, using robust models, adding model components that address the violations etc., but it’s important to think about the consequences of assumption-violation (bias? Type I error/undercoverage of confidence intervals?) and their likely magnitudes  

25 .  I would like to know what´s your opinion of using weights (VarIdent, etc...) vs log-transforming variables. 

 **Answer**: 
 Weights (or using a GLM, which includes a specification for the mean-variance relationship) are more flexible than transformation. Transformation affects lots of components of the model simultaneously: linearity, meaning of interactions, variance-mean relationship, distribution of the residuals.  However, when it works (as in this example) log-transformation is great because (1) it’s very easy and (2) it does address all of the components listed at once. 

26 .  the spatial autocorrelation can be estimated through distance decay relationship too, right? 

 **Answer**: 
 Not sure I understand the question. Instead of drawing a picture I could have computed a spatial variogram or correlogram of the residuals (“distance decay relationship”); I didn’t because it’s a little bit fussy to compute/draw variograms when I have to worry about points on a sphere (great circle distances etc.). I could have used a different model to account for spatial autocorrelation – Gaussian process, Markov random field, etc. (some of this is hinted at in the ‘extras’ notes) 

27 .  Is there any way to make the mixed models accept negative variances, instead of forcing negative vari to zero and generate the Singularity warning? 

 **Answer**: 
 No, not really. In the population genetics literature there is some stuff about how to deal with negative variances. The way forward with this (although it’s not always simple) is to use a compound symmetric structure that allows negative within-group correlations. See e.g. Molenberghs, Geert, and Geert Verbeke. “A Note on a Hierarchical Interpretation for Negative Variance Components.” Statistical Modelling 11, no. 5 (2011): 389–408. https://doi.org/10.1177/1471082X1001100501.
 

