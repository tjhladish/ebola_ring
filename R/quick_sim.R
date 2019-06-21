# non-network assessment

#' modelling concept:
#'  - contacts based on negative binomial model
#'  - R0 sets probability that contacts will be exposed - solve for p given R0 = <contacts * p transmission>
#'  - only cases come from (1) patient 0, (2) their contact exposures
#'  - only negatives come from (2) patient zero contact non-exposures, and contacts of contacts
#'  
simulator <- function (R0, negbinom_pars, vac_eff, vac_coverage, i) {
  # preallocate return structure
  res <- c(vax_pos = 0, vax_neg = 0, unvax_pos = 0, unvax_neg = 0)
  set.seed(i)
  
  covered <- runif(1) < vac_coverage
  resisted <- covered & (runif(1) < vac_eff)
  # if patient 0 resists exposure, no info
  if (covered & resisted) return(res)
  
  if (covered) res[1] <- 1
  
  p0_exps <- with(negbinom_pars, rnbinom(1, r, p))
  
  return(res)
}