# Poisson point process for presence only data featuring Joshua Trees

# packages
pkgs <- c("here", "readr", "janitor", "mgcv", "ppgam", "gratia", "dplyr",
          "ggplot2", "mvnfast", "patchwork", "tidyr")
vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

# data - from Godsoe et al 2009 New Phytologist
# https://doi.org/10.1111/j.1469-8137.2009.02942.x
# download from dryad https://doi.org/10.5061/dryad.6s67t
joshua <- read_tsv(here("data/godsoe-et-al/joshua-tree.txt"),
                   col_types = "ddc")

joshua %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point() +
  coord_quickmap()

# fit the Poisson point process
# about 5 seconds
m <- ppgam(~ s(longitude, latitude, k = 100), data = joshua)

summary(m)
draw(m, rug = FALSE)

m_ds <- ppgam(~ s(longitude, latitude, k = 100, bs = "ds", m = c(2,0.5)),
              data = joshua)

summary(m_ds)
draw(m_ds)

# validation
lon_valid <- range(joshua$longitude) + 1e-6 * c(-1, 1)
lat_valid <- range(joshua$latitude) + 1e-6 * c(-1, 1)

mids <- function(x) x[-1] - .5 * diff(x)

nlon_valid <- 50
nlat_valid <- 50

lon_seq <- seq(lon_valid[1], lon_valid[2], l = nlon_valid + 1)
lat_seq <- seq(lat_valid[1], lat_valid[2], l = nlat_valid + 1)
lon_cut <- cut(joshua$longitude, lon_seq)
lat_cut <- cut(joshua$latitude, lat_seq)
n_valid <- nlon_valid * nlat_valid # total number of integration points

lon_brks_valid <- seq(lon_valid[1], lon_valid[2], l = nlon_valid + 1)
lon_nodes_valid <- mids(lon_brks_valid)
lat_brks_valid <- seq(lat_valid[1], lat_valid[2], l = nlat_valid + 1)
lat_nodes_valid <- mids(lat_brks_valid)

quad_t_valid <- data.frame(longitude = cut(joshua$longitude, lon_brks_valid),
                           latitude = cut(joshua$latitude, lat_brks_valid))
tab_t_valid <- table(quad_t_valid)

mult <- nlon_valid * nlat_valid
nodes_wts_valid <- rep(1, nlon_valid) %o% rep(1, nlat_valid)
nodes_wts_valid <- nodes_wts_valid / nlon_valid / nlat_valid
df_valid <- expand.grid(longitude = lon_nodes_valid,
                        latitude = lat_nodes_valid) # integration nodes
Xp_valid <- predict(m_ds, df_valid, type = "lpmatrix")
E <- array(c(Xp_valid %*% coef(m_ds)), dim = c(nlon_valid, nlat_valid))
E <- exp(E) * nodes_wts_valid

(mu1 <- sum(E))
(mu2 <- sum(tab_t_valid))
mu1 / mu2

O <- apply(tab_t_valid, 1:2, sum)

n_samp <- 1e3
mu <- coef(m_ds)
n_par <- length(mu)
set.seed(42)
samps <- t(rmvn(n_samp, mu, vcov(m_ds), ncores = 4))
m2 <- c(Xp_valid %*% samps)
m2 <- matrix(m2, ncol = n_samp)
m2 <- exp(m2) * c(nodes_wts_valid)
sd_valid <- apply(m2, 1, sd)

fv <- tibble(observed = as.vector(O),
             expected = as.vector(E),
             sd = sd_valid)
fv <- bind_cols(as_tibble(df_valid), fv)

oe <- fv %>%
  select(!sd) %>%
  pivot_longer(observed:expected,
               names_to = "type",
               values_to = "count")

p_oe <- ggplot(oe, aes(x = longitude, y = latitude, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  coord_quickmap() +
  facet_wrap(~ type)

p_oe

p_e <- ggplot(fv, aes(x = longitude, y = latitude, fill = expected)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  coord_quickmap()

p_o <- ggplot(fv, aes(x = longitude, y = latitude, fill = observed)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  coord_quickmap()

p_sd <- ggplot(fv, aes(x = longitude, y = latitude, fill = sd)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  coord_quickmap()

p_o

p_e

p_sd

p_o + p_e
