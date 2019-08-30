# some convenience plot defaults

scale_x_count <- function(name="Count", ...) scale_x_continuous(name, ...)

scale_x_cases <- function(name="Cases", ...) scale_x_continuous(name, ...)

facet_labels <- labeller(
  exp_sd = function(s) 
)
