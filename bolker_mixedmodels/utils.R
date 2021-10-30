dwplot_ordered <- function(model, ...) {
  tt <- (broom.mixed::tidy(model, ...)
    %>% filter(term != "(Intercept)")
    %>% mutate(across(term, fct_reorder, estimate))
    %>% ggplot(aes(y=term, x = estimate,
                   xmin = estimate - 2*std.error,
                   xmax = estimate + 2*std.error))
    + geom_pointrange()
    + geom_vline(xintercept = 0, lty = 2)
  )
  return(tt)
}
