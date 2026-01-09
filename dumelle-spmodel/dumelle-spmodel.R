## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: setup
#| include: false

# load background packages
library(tidyverse) # for general data wrangling
library(broom) # for tidy() of lm() objects
library(kableExtra) # for tables
library(sf) # for spatial data operations
library(spmodel) # for spatial modeling
library(tigris) # for southwestern state boundaries


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false

# simulate data
set.seed(1)
n <- 20
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- rnorm(n)
dat <- data.frame(x1, x2, y)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lmod <- lm(formula = y ~ x1 + x2, data = dat)
summary(lmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

# set reproducible seed
set.seed(29)

# simulate data 
n <- 200
xc <- runif(n)
yc <- runif(n)
dat <- tibble(xc, yc)
params <- spcov_params(spcov_type = "exponential", de = 2, ie = 0, range = 1)
dat$x1 <- sprnorm(params, data = dat, xcoord = xc, ycoord = yc)
dat$x2 <- rnorm(n, mean = 0, sd = 0.1)
dat$y <- sprnorm(params, mean = 1 * dat$x2 , data = dat, xcoord = xc, ycoord = yc)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-realized_y
#| fig-align: center
#| fig-cap: "The spatial distribution of the response variable, y."
ggplot(dat, aes(x = xc, y = yc, color = y)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

lmod <- lm(y ~ x1 + x2, data = dat)
tidy_lmod <- tidy(lmod, conf.int = TRUE)

spmod <- splm(y ~ x1 + x2, data = dat, spcov_type = "exponential",
              xcoord = xc, ycoord = yc)
tidy_spmod <- tidy(spmod, conf.int = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

b1 <- tibble(
  method = c("Nonspatial", "Spatial"),
  est = c(tidy_lmod$estimate[2], tidy_spmod$estimate[2]),
  se = c(tidy_lmod$std.error[2], tidy_spmod$std.error[2]),
  p.value = c(tidy_lmod$p.value[2], tidy_spmod$p.value[2]),
  conf.low = c(tidy_lmod$conf.low[2], tidy_spmod$conf.low[2]),
  conf.high = c(tidy_lmod$conf.high[2], tidy_spmod$conf.high[2])
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: tab-fe_tidy1
#| tbl-cap: "Parameter inference for $\\beta_1$ in the nonspatial and spatial linear models."

b1 |>
  kbl(digits = 2) |>
  kable_classic(full_width = FALSE) |>
  column_spec(4, color = c("red", "green")) |>
  row_spec(seq(1, NROW(b1)), extra_css = "border-bottom-style: none;")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

b2 <- tibble(
  method = c("Nonspatial", "Spatial"),
  est = c(tidy_lmod$estimate[3], tidy_spmod$estimate[3]),
  se = c(tidy_lmod$std.error[3], tidy_spmod$std.error[3]),
  p.value = c(tidy_lmod$p.value[3], tidy_spmod$p.value[3]),
  conf.low = c(tidy_lmod$conf.low[3], tidy_spmod$conf.low[3]),
  conf.high = c(tidy_lmod$conf.high[3], tidy_spmod$conf.high[3])
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: tab-fe_tidy2
#| tbl-cap: "Parameter inference for $\\beta_2$ in the nonspatial and spatial linear models."

b2 |>
  kbl(digits = 2) |>
  kable_classic(full_width = FALSE) |>
  column_spec(4, color = c("red", "green")) |>
  row_spec(seq(1, NROW(b2)), extra_css = "border-bottom-style: none;")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

lmod <- update(spmod, spcov_type = "none")
lmod_loocv <- loocv(lmod, cv_predict = TRUE)
spmod_loocv <- loocv(spmod, cv_predict = TRUE)
dat$Nonspatial <- lmod_loocv$cv_predict
dat$Spatial <- spmod_loocv$cv_predict
long_dat <- dat |>
  pivot_longer(c(Nonspatial, Spatial), names_to = "approach", values_to = "pred")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-loocv
#| fig-align: center
#| fig-cap: "Observations vs leave-one-out cross validation predictions for the nonspatial and spatial linear models."

ggplot(long_dat, aes(x = pred, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "orange", linewidth = 1.5) +
  facet_wrap(~ approach) +
  labs(x = "Leave-One-Out Prediction")  +
  geom_text(data = data.frame(
    pred = -2.1,
    y = 1.8,
    approach = c("Nonspatial", "Spatial"),
    label = c("RMSPE = 0.86", "RMSPE = 0.31")
  ), aes(label = label)) +
  geom_text(data = data.frame(
    pred = -2.05,
    y = 1.5,
    approach = c("Nonspatial", "Spatial"),
    label = c("R2 = 0.36", "R2 = 0.92")
  ), aes(label = label)) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-pos_cov
#| fig-align: center
#| fig-cap: "Positive covariance between two random variables, x and y."

set.seed(0)
x <- rnorm(50)
y <- 2 * x + rnorm(50)
dat <- tibble(x, y)
ggplot(dat, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-neg_cov
#| fig-align: center
#| fig-cap: "Negative covariance between two random variables, x and y."

set.seed(1)
x <- rnorm(50)
y <- -2 * x + rnorm(50)
dat <- tibble(x, y)
ggplot(dat, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-no_cov
#| fig-align: center
#| fig-cap: "No covariance between two random variables, x and y."

set.seed(2)
x <- rnorm(50)
y <- rnorm(50)
dat <- tibble(x, y)
ggplot(dat, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

set.seed(1)
n <- 200
xc <- runif(n)
yc <- runif(n)
dat <- tibble(xc, yc)

p1 <- spcov_params(spcov_type = "exponential", de = 1, ie = 0, range = 1)
dat$y11 <- sprnorm(p1, data = dat, xcoord = xc, ycoord = yc)
dat$y12 <- sprnorm(p1, data = dat, xcoord = xc, ycoord = yc)
dat$y13 <- sprnorm(p1, data = dat, xcoord = xc, ycoord = yc)

p2 <- spcov_params(spcov_type = "exponential", de = 0.01, ie = 0.99, range = 1)
dat$y21 <- sprnorm(p2, data = dat, xcoord = xc, ycoord = yc)
dat$y22 <- sprnorm(p2, data = dat, xcoord = xc, ycoord = yc)
dat$y23 <- sprnorm(p2, data = dat, xcoord = xc, ycoord = yc)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-high_spcov1
#| fig-align: center
#| fig-cap: "A random variable with strong spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y11)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-high_spcov2
#| fig-align: center
#| fig-cap: "A random variable with strong spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y12)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-high_spcov3
#| fig-align: center
#| fig-cap: "A random variable with strong spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y13)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-low_spcov1
#| fig-align: center
#| fig-cap: "A random variable with weak spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y21)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-low_spcov2
#| fig-align: center
#| fig-cap: "A random variable with weak spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y22)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-low_spcov3
#| fig-align: center
#| fig-cap: "A random variable with weak spatial covariance."

ggplot(dat, aes(x = xc, y = yc, color = y23)) +
  geom_point() +
  scale_color_viridis_c(name = "y") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

h <- seq(1e-4, 4, length.out = 1000)
sill <- 3
spcov1 <- tibble(h = h, cov = sill * exp(-h^2), Example = 1) # gaussian
spcov2 <- tibble(h = h, cov = sill * exp(-h), Example = 2) # exponential
spcov3 <- tibble(h = h,
                 cov = sill * 0.72 * (1 - 1.5 * h/0.47 + 0.5 * (h/0.47)^3) * (h <= 0.47),
                 Example = 3) # spherical
spcov4 <- tibble(h = h, cov = sill * 0.4 * exp(-h/3), Example = 4) # exponential
spcov5 <- tibble(h = h, cov = sill * 0.25 * exp(-h^2/3), Example = 5) # gaussian
spcov6 <- tibble(h = h,
                 cov = sill * 0.15 * (1 - 1.5 * h/1.6 + 0.5 * (h/1.6)^3) * (h <= 1.6),
                 Example = 6) # spherical
spcovs <- bind_rows(spcov1, spcov2, spcov3, spcov4, spcov5, spcov6) |>
  mutate(Example = factor(Example))
spcov_total <- data.frame(h = 0, cov = sill)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-ex_spcov
#| fig-align: center
#| fig-cap: "Examples of several different spatial covariance functions."

ggplot() +
  geom_point(data = spcov_total, aes(x = h, y = cov), size = 4) +
  geom_line(data = spcovs, aes(x = h, y = cov, color = Example, linetype = Example), linewidth = 2) +
  scale_color_viridis_d(direction = -1, begin = 0.2) +
  labs(x = "Distance", y = "Spatial Covariance") +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lake


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false
states <- states(cb = TRUE) |>
  filter(STUSPS %in% c("AZ", "CO", "NV", "UT")) |>
  st_transform(crs = 5070)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-logcond
#| fig-align: center
#| fig-cap: "The (natural) logarithm of conductivity at lakes in the southwestern United States."

ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = lake, aes(color = log_cond, shape = state)) + 
  scale_color_viridis_c() +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spmod <- splm(
  formula = log_cond ~ elev + origin,
  data = lake,
  spcov_type = "exponential"
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
tidy(spmod)
tidy(spmod, conf.int = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
glance(spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
augment(spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-plot_stdresid_fitted
#| fig-align: center
#| fig-cap: "Standardized residuals vs fitted values. There is no obvious pattern between the standardized residuals and fitted values."

plot(spmod, which = 1)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-logcond_explanatory_spcov
#| fig-align: center
#| fig-cap: "The fitted spatial covariance of log conductivity (after accounting for the explanatory variables)."

plot(spmod, which = 7)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lake_preds


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-logcond_preds_spcov
#| fig-align: center
#| fig-cap: "The (natural) logarithm of conductivity at lakes in the southwestern United States alongside log conductivity predictions at new lakes (black triangles)."

ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = lake, aes(color = log_cond)) + 
  scale_color_viridis_c() +
  geom_sf(data = lake_preds, pch = 17, size = 6) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
predict(spmod, newdata = lake_preds)
augment(spmod, newdata = lake_preds)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-logcond_preds_spcov2
#| fig-align: center
#| fig-cap: "The (natural) logarithm of conductivity at lakes in the southwestern United States alongside log conductivity predictions (triangles)."

aug_spmod <- augment(spmod, newdata = lake_preds)
ggplot() +
  geom_sf(data = states, fill = "white") +
  geom_sf(data = lake, aes(color = log_cond)) + 
  geom_sf(data = aug_spmod, aes(color = .fitted), pch = 17, size = 6) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
predict(spmod, newdata = lake_preds, interval = "prediction")
augment(spmod, newdata = lake_preds, interval = "prediction")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lmod <- splm(log_cond ~ elev + origin, data = lake, spcov_type = "none")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
glances(lmod, spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
loocv(lmod)
loocv(spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
pseudoR2(lmod)
pseudoR2(spmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spmod_list <- splm(
  formula = log_cond ~ elev + origin,
  data = lake,
  spcov_type = c("exponential", "gaussian", "none")
)
glances(spmod_list)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
citation(package = "spmodel")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1)
n <- 2500
xc <- runif(n)
yc <- runif(n)
x <- runif(n)
dat <- data.frame(xc, yc, x)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
params <- spcov_params(spcov_type = "exponential", de = 1, ie = 0.2, range = 1)
mu <- 1 + x
dat$y <- sprnorm(params, mean = mu, data = dat, xcoord = xc, ycoord = yc)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_train <- 2000
train_dat <- dat[seq(1, n_train), ]
test_dat <- dat[-seq(1, n_train), ]


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-train_y
#| fig-align: center
#| fig-cap: "The simulated response variable in the training data."

ggplot(data = train_dat, aes(x = xc, y = yc, color = y)) +
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-train_x
#| fig-align: center
#| fig-cap: "The simulated explanatory variable in the training data."

ggplot(data = train_dat, aes(x = xc, y = yc, color = x)) +
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-train_y_x
#| fig-align: center
#| fig-cap: "The simulated response and explanatory variables in the training data."

ggplot(data = train_dat, aes(x = x, y = y)) +
  geom_point() +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
coords <- data.frame(train_dat$xc, train_dat$yc)
dists <- coords |>
  dist() |>
  as.matrix()
kmeans_output <- kmeans(dists, centers = 4)
train_dat$kmeans <- kmeans_output$cluster |>
  as.factor()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
index <- rep(1:4, length.out = n_train)
train_dat$random <- index |>
  sample() |>
  as.factor()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-spin_kmeans
#| fig-align: center
#| fig-cap: "Spatial indexing groups chosen via kmeans clustering."

ggplot(data = train_dat, aes(x = xc, y = yc, color = kmeans, shape = kmeans)) +
  geom_point() +
  scale_color_viridis_d() + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-spin_random
#| fig-align: center
#| fig-cap: "Spatial indexing groups chosen via random assignment."

ggplot(data = train_dat, aes(x = xc, y = yc, color = random, shape = random)) +
  geom_point() +
  scale_color_viridis_d() + 
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| eval: false

# local_list <- list(
#   method = "kmeans",
#   size = 100,
#   var_adjust = "theoretical",
#   parallel = FALSE
# )
# splm(..., local = local_list)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spmod_direct_start <- Sys.time()
spmod_direct <- splm(
  formula = y ~ x,
  data = train_dat,
  spcov_type = "exponential",
  xcoord = xc,
  ycoord = yc
)
spmod_direct_time <- Sys.time() - spmod_direct_start


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spmod_bigdata_start <- Sys.time()
spmod_bigdata <- splm(
  formula = y ~ x,
  data = train_dat,
  spcov_type = "exponential",
  xcoord = xc,
  ycoord = yc,
  local = TRUE
)
spmod_bigdata_time <- Sys.time() - spmod_bigdata_start


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spmod_direct_time
spmod_bigdata_time


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
tidy(spmod_direct) |>
  filter(term == "x")

tidy(spmod_bigdata) |>
  filter(term == "x")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
loocv(spmod_direct)
loocv(spmod_bigdata)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
coef(spmod_direct, type = "spcov")[1:3]
coef(spmod_bigdata, type = "spcov")[1:3]


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-direct_fit
#| fig-align: center
#| fig-cap: "Spatial covariance from the direct model."

plot(spmod_direct, which = 7)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-spin_fit
#| fig-align: center
#| fig-cap: "Spatial covariance from the spatial indexing (SPIN) model."

plot(spmod_bigdata, which = 7)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| eval: false

# local_list <- list(
#   method = "covariance",
#   size = 100,
#   parallel = FALSE
# )
# predict(object, ..., local = local_list)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
preds_direct_start <- Sys.time()
preds_direct <- predict(spmod_direct, newdata = test_dat, se.fit = TRUE)
preds_direct_time <- Sys.time() - preds_direct_start


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
preds_bigdata_start <- Sys.time()
preds_bigdata <- predict(spmod_bigdata, newdata = test_dat, se.fit = TRUE)
preds_bigdata_time <- Sys.time() - preds_bigdata_start


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
preds_direct_time
preds_bigdata_time


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-test_y
#| fig-align: center
#| fig-cap: "The simulated response variable in the test data."

ggplot(test_dat, aes(x = xc, y = yc, color = y)) +
  geom_point() +
  scale_color_viridis_c(limits = c(-3.25, 2.5)) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

test_dat$direct <- preds_direct$fit
test_dat$bigdata <- preds_bigdata$fit
test_dat_long <- test_dat |>
  pivot_longer(cols = c(direct, bigdata), names_to = "Approach", values_to = "Pred") |>
  mutate(Approach = factor(Approach, levels = c("direct", "bigdata")))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-test_y_preds
#| fig-align: center
#| fig-cap: "Predictions for the direct and big data (LNBH) approaches."

ggplot(test_dat_long, aes(x = xc, y = yc, color = Pred)) +
  geom_point() +
  scale_color_viridis_c(limits = c(-3.25, 2.5)) +
  facet_wrap(~ Approach) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-test_y_se
#| fig-align: center
#| fig-cap: "Standard Errors for the direct and big data (LNBH) approaches."

test_dat$direct <- preds_direct$se.fit
test_dat$bigdata <- preds_bigdata$se.fit
test_dat_long <- test_dat |>
  pivot_longer(cols = c(direct, bigdata), names_to = "Approach", values_to = "StdErr") |>
  mutate(Approach = factor(Approach, levels = c("direct", "bigdata")))
ggplot(test_dat_long, aes(x = xc, y = yc, color = StdErr)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~ Approach) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
rmspe_direct <- sqrt(mean((test_dat$y - test_dat$direct)^2))
rmspe_direct
rmspe_bigdata <- sqrt(mean((test_dat$y - test_dat$bigdata)^2))
rmspe_bigdata


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
moose


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-moose_pres
#| fig-align: center
#| fig-cap: "Moose presence in the Togiak Region."

ggplot(moose, aes(color = presence)) +
  geom_sf() +
  scale_color_viridis_d() +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
spgmod <- spglm(
  formula = presence ~ elev * strat,
  family = binomial,
  data = moose,
  spcov_type = "spherical"
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
tidy(spgmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
glance(spgmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
augment(spgmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
AUROC(spgmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
loocv(spgmod)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
moose_preds


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-moose_pres_pred
#| fig-align: center
#| fig-cap: "Moose observed data and prediction locations."

ggplot() +
  geom_sf(data = moose) +
  geom_sf(data = moose_preds, color = "orange", shape = 17) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
preds <- predict(spgmod, newdata = moose_preds, type = "response")
head(preds)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
augment(spgmod, newdata = moose_preds, interval = "prediction")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| include: false

moose$.fitted <- fitted(spgmod)
moose_preds$.fitted <- preds


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: fig-moose_fit_pred
#| fig-align: center
#| fig-cap: "Fitted values and predictions for moose presence probability."

ggplot() +
  geom_sf(data = moose, aes(color = .fitted)) +
  geom_sf(data = moose_preds, aes(color = .fitted), shape = 17) +
  scale_color_viridis_c(limits = c(0, 1)) +
  theme_bw(base_size = 14)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#| echo: false
#| label: tab-glm
#| tbl-cap: "Generalized linear model response variable types, families, and link functions supported by `spmodel`."

dat <- tibble(
  type = c("binary", "count", "count", "proportion", "skewed", "skewed"),
  family = c("binomial", "poisson", "nbinomial", "beta", "Gamma", "inverse.gaussian"),
  link = c("logit", "log", "log", "logit", "log", "log")
)
dat |>
  kbl() |>
  kable_classic()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

