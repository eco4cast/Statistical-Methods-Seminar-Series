# brain dump for EFI talk

- 15 minutes for overview (subject-matter + essential mixed model stuff)
- 30 minutes for code (with interspersed additional theory)
- extras (spatial autocorr etc.)

# stats intro

- models
   - what is a random effect?
    (table from Fox chapter)
	
- describe tools (lme4; mixture of lattice/ggplot graphs)
- scope
   - skimming very quickly over (G)LM-specific stuff
   - mostly LMMs, come back to GLMMs at the end
   - mostly lme4, describe wider landscape at end

# science intro

- describe data and questions
- pix from mmd project
- univariate ggplot graphs?

## need for mixed models: geographic variation

- three-level structure (biome/realm/interaction); random slopes
- Yackulic and Uriarte/Gelman graphs
- log scale, then center (Schielzeth); don't scale, for now

## simplest model

- start with simplest case (1-level, intercept-only)
- basic diagnostics

## three-level model

- discuss nesting/crossing

## Nested vs crossed designs

**Nested**: sub-unit IDs only measured within a single larger unit.
e.g.: Plot1 in Block1 independent of Plot1 in Block2

![](pix/CV_nested.png)

**Crossed**: sub-unit IDs can be measured in multiple larger units.
e.g. year, site

![](pix/CV_crossed.png)

**Unique coding**: removes ambiguity

![](pix/CV_unique.png)

Robert Long, [Cross Validated](https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified)

## maximal model

- idea
- why it usually doesn't work
- (confounding with residual variance)

## model simplification

- avoid singularity
- convergence warnings

## AIC table/strategy

- `for` loop over table

## diagnostics

basic
DHARMa

## spatial correlation

- diagnosis
- choices 


## display/description

coefficient plots
predictions
partial residuals
R^2 values

# extras

- more on regularization
- more on model simplification (compound symmetry, factor-analytic)
- more complex structures (AR etc.)
- more on autocorrelation (INLA, gamm4, brms; soap-film, MRF, ?)
- more on available packages (Google sheet)
