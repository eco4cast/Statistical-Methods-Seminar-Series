# brain dump for EFI talk

## stats intro

- describe tools (lme4; mixture of lattice/ggplot graphs)
- scope
   - skimming very quickly over (G)LM-specific stuff
   - mostly LMMs, come back to GLMMs at the end
   - mostly lme4, describe wider landscape at end

## science intro

- describe data and questions
- basic graphs (univariate ggplot)

## need for mixed models: geographic variation

- three-level structure (biome/realm/interaction); random slopes
- Yackulic and Uriarte/Gelman graphs
- log scale, then center (Schielzeth); don't scale, for now

## simplest model

- start with simplest case (1-level, intercept-only)
- basic diagnostics

## three-level model

- mention nesting/crossing

## maximal model

- idea
- why it usually doesn't work

## model simplification

- avoid singularity
- convergence warnings

## AIC table/strategy

## diagnostics

basic
DHARMa

## spatial correlation

diagnosis
choices (INLA, gamm4, brms; soap-film, MRF, ?)


## display/description

coefficient plots
predictions
partial residuals
R^2 values

