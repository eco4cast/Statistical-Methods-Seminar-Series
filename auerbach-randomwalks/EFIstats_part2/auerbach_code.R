## Code for "Modeling biological processes as stopped random walks"
## Presented at EFI and Statistical Ecology Section Webinar
## Date: December 2, 2024 
## Author: Jonathan Auerbach
## Contact: jauerba@gmu.edu
## Website: jauerbach.github.io


############
# 0. Setup #
############

library("tidyverse")
library("rstan")
options(mc.cores = parallel::detectCores())


############################################
# CLT for Stopped Random Walks: Simulation #
############################################

set.seed(1)

simulation_data <-
  tibble(temp = rnorm(121, 10, 5),
         date = seq(as.Date("2024-01-01"), as.Date("2024-04-30"), "day"))


bloom_date <-
  simulation_data %>%
  filter(date >= as.Date("2024-02-01")) %>%
  mutate(status = cumsum(temp) > 300) %>%
  filter(status == TRUE) %>%
  summarize(bloom_date = first(date)) %>%
  pull(bloom_date)

simulation_data %>% 
  ggplot() + 
  aes(date, 
      cumsum(ifelse((date >= as.Date("2024-02-01")) &
                      (date <= as.Date(bloom_date)),
                    temp, 
                    0))) + 
  geom_line() +
  theme_bw() +
  labs(x = "", y = "cumulative temperature (°C)") +
  geom_vline(xintercept = as.Date("2024-02-01"), linetype = 2) +
  geom_vline(xintercept = as.Date(bloom_date), linetype = 2) +
  annotate("label", 
           y = c(50, 250), 
           x = c(as.Date("2024-02-01"), as.Date(bloom_date)),
           label = c("last\nfrost", "first\nbloom"))

simulation_data <- 
  function(mean_temp)
  tibble(temp = rnorm(121, mean_temp, 5),
         date = seq(as.Date("2024-01-01"),
                    as.Date("2024-04-30"),
                    "day")) %>%
  filter(date >= as.Date("2024-02-01")) %>% 
  mutate(status = cumsum(temp) > 300) %>%
  filter(status == TRUE) %>%
  summarize(bloom_date = first(date)) %>%
  mutate(mean_temp = mean_temp)

lapply(rep(c(5, 10, 15, 20, 25), 20), 
       simulation_data) %>%
  Reduce(bind_rows, x = .) %>%
  mutate(doy = as.numeric(format(bloom_date, "%j"))) %>%
  lm(formula = doy ~ I(1/mean_temp),
     weight = mean_temp^3, 
     data = .) %>%
  broom::tidy()


####################################
# Application 1: Experimental Data #
####################################

charrier11 <- 
  read_csv("data/dwsl.csv") %>%
  filter(datasetID == "charrier11", 
         study == "exp2")

figure_charrier <-
  charrier11 %>%
  ggplot() +
  aes(forceday, response.time) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(linewidth = 0.5, colour = "darkgray")) +
  labs(x = "Temperature (°C)",
       y =  "Mean Time until Budburst (days)") +
  ylim(0, 250)

figure_charrier +
  geom_smooth(aes(weight = forceday^3), 
              method = "lm", formula = y ~ I(1/x), color = "red") +
  geom_smooth(method = "glm", formula = y ~ I(1/x),
              method.args = list(family = inverse.gaussian(link = "identity")))


stan_mod1 <- stan_model("stan_code/stan_mod1.stan")
fit1 <- 
  sampling(stan_mod1,
           data = list(n = nrow(charrier11),
                       m_lower = 5,
                       m_upper = 25,
                       mu = charrier11$forceday,
                       y = charrier11$response.time),
           control = list(adapt_delta = 0.95),
           refresh = 0)

figure_charrier +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey", alpha = .5,
              data = 
                tibble(forceday = 5:25,
                       response.time =  colMeans(extract(fit1, "y_pred")[[1]]),
                       ymin = apply(extract(fit1, "y_pred")[[1]], 2, quantile, prob = .025),
                       ymax = apply(extract(fit1, "y_pred")[[1]], 2, quantile, prob = .975))) +
  geom_line(data = 
              tibble(forceday = 5:25,
                     response.time =  colMeans(extract(fit1, "y_pred")[[1]]),),
            color = "orange", linewidth = 1)


#####################################
# Application 2: Observational Data #
#####################################

marsham <- 
  read_csv("data/Marsham-Combes_UK.csv")  %>%
  select(year, response.time = oak) %>%
  na.omit()

temp <- 
  read_csv("data/meantemp_daily_totals.csv") %>%
  mutate(date  = seq(as.Date("1772-01-01"), as.Date("2024-10-07"), "day"),
         year  = as.numeric(format(date, "%Y")),
         month = as.numeric(format(date, "%m")),
         doy   = as.numeric(format(date, "%j")))

marsham %<>%
  left_join(temp %>% 
              filter(month >= 2,
                     month <= 4) %>% 
              group_by(year) %>% 
              summarize(spring.temp = sum(pmax(0, Value)))) %>%
  na.omit()

figure_marsham <-
  marsham %>%
  ggplot() +
  aes(spring.temp, response.time) +
  geom_point(alpha = .5) +
  labs(x = "Cumulative Daily Temperature (°C, 2/1 to 4/30)",
       y = "Time until Budburst (number of days after January 1)") +
  theme_bw() +
  ylim(70, 180)

figure_marsham +
  geom_smooth(aes(weight = spring.temp^3), 
              method = "lm", formula = y ~ I(1/x), color = "red") +
  geom_smooth(method = "glm", formula = y ~ I(1/x),
              method.args = list(family = inverse.gaussian(link = "identity")))

figure_marsham +
  geom_smooth(aes(weight = spring.temp^3), 
              method = "lm", formula = y ~ I(1/x), color = "red") +
  geom_smooth(method = "glm", formula = y ~ I(1/x),
              method.args = list(family = inverse.gaussian(link = "identity")))


########################
# Non-asymptotic Model #
########################

stan_mod3 <- stan_model("stan_code/stan_mod3.stan")

fit3 <- 
  sampling(stan_mod3,
           data = list(n = nrow(charrier11),
                       m_lower = 5,
                       m_upper = 25,
                       mu = charrier11$forceday,
                       y = charrier11$response.time),
           control = list(adapt_delta = 0.95),
           refresh = 0,
           init = rep(list(list(sigma = 1e4)), 4))

figure_charrier +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey", alpha = .5,
              data = 
                tibble(forceday = 5:25,
                       response.time =  colMeans(extract(fit3, "y_pred")[[1]]),
                       ymin = apply(extract(fit3, "y_pred")[[1]], 2, quantile, prob = .025),
                       ymax = apply(extract(fit3, "y_pred")[[1]], 2, quantile, prob = .975))) +
  geom_line(data = 
              tibble(forceday = 5:25,
                     response.time =  colMeans(extract(fit3, "y_pred")[[1]]),),
            color = "orange", linewidth = 1)

stan_mod4 <- stan_model("stan_code/stan_mod4.stan")

temp_matrix <-
  temp %>%
  filter(year %in% marsham$year) %>%
  transmute(doy, year, temp = pmax(0, Value)) %>%
  pivot_wider(names_from = doy, 
              values_from = temp,
              values_fill = 0) %>%
  select(-year) %>%
  as.matrix() 

fit4 <- 
  sampling(stan_mod4, 
           data = list(n = nrow(marsham),
                       m_lower = 300,
                       m_upper = 800,
                       temp = temp_matrix,
                       y = marsham$response.time),
           control = list(adapt_delta = 0.99),
           refresh = 0,
           init = rep(list(list(sigma = 1e6)), 4))

figure_marsham +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey", alpha = .5,
              data = 
                tibble(spring.temp = 300:800,
                       response.time =  colMeans(extract(fit4, "y_pred")[[1]]),
                       ymin = apply(extract(fit4, "y_pred")[[1]], 2, quantile, prob = .025),
                       ymax = apply(extract(fit4, "y_pred")[[1]], 2, quantile, prob = .975))) +
  geom_line(data = 
              tibble(spring.temp = 300:800,
                     response.time =  colMeans(extract(fit4, "y_pred")[[1]])),
            color = "orange", linewidth = 1)


##############
# Comparison #
##############

tibble(data = c(rep("Charrier", 2), rep("Marsham", 2)), 
       model = rep(c("asymptotic", "non-asymptotic"), 2),
       gamma = c(mean(extract(fit1, "gamma")[[1]]),
                 mean(extract(fit3, "gamma")[[1]]),
                 mean(extract(fit2, "gamma")[[1]]) / 90,
                 mean(extract(fit4, "gamma")[[1]])),
       lower = c(quantile(extract(fit1, "gamma")[[1]], p = .025),
                 quantile(extract(fit3, "gamma")[[1]], p = .025),
                 quantile(extract(fit2, "gamma")[[1]], p = .025) / 90,
                 quantile(extract(fit4, "gamma")[[1]], p = .025)),
       upper = c(quantile(extract(fit1, "gamma")[[1]], p = .975),
                 quantile(extract(fit3, "gamma")[[1]], p = .975),
                 quantile(extract(fit2, "gamma")[[1]], p = .975) / 90,
                 quantile(extract(fit4, "gamma")[[1]], p = .975)))
