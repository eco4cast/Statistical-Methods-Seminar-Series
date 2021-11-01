dwplot_ordered <- function(model, ..., dodge = 0.5) {
  if (!require("broom.mixed")) stop("please install broom.mixed")
  multi_model <- is.list(model)
  tt <- if (!multi_model) tidy(model, ...) else purrr::map_dfr(model, tidy, .id = "model", ...)
  tt <- (tt
    %>% filter(term != "(Intercept)")
    %>% mutate(across(term, fct_reorder, estimate))
  )
  gg <- (ggplot(tt, aes(y=term, x = estimate,
                       xmin = estimate - 2*std.error,
                       xmax = estimate + 2*std.error)))
  if (multi_model) {
    gg <- gg + aes(colour = model)
    pd <- position_dodge(width = dodge)
  } else pd <- position_identity()
  gg1 <- (gg
    + geom_pointrange(position = pd)
    + geom_vline(xintercept = 0, lty = 2)
  )
  return(gg1)
}

has_warning <- function(x) {
  m <- x@optinfo$conv$lme4$messages
  length(m)>0 && any(grepl("failed", m))
}

CRAN_pkgs <- c("lme4", "gamm4", "DHARMa", "broom.mixed",
               "gridExtra", "tidyverse", "remotes", "tidyverse", "car",
               "dotwhisker", "emmeans")
GH_pkgs <- c("bbolker/r2glmm", "hohenstein/remef")
install_all_pkgs <- function() {
  i1 <- installed.packages()
  to_install <- setdiff(CRAN_pkgs, rownames(i1))
  if (length(to_install) > 0) {
    install.packages(to_install)
  }
  sapply(GH_pkgs, remotes::install_github)
}
