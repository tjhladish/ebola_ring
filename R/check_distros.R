require(jsonlite)
require(data.table)
require(cowplot)

.args <- c("timedistros.json")
.args <- commandArgs(trailingOnly = T)

# standard is wikipedia k, theta formulation
# so mean = k*theta
# sd^2 = k theta^2 = mn * theta => theta = sd^2/mn
# and k => (mn/sd)^2
# cpp gamma_distribution uses alpha = k, beta = theta
# rgamma uses shape = k, scale = theta

gamma_alpha <- function(mn, sd) (mn/sd)^2
gamma_beta <- function(mn, sd) (sd^2)/mn
adgamma <- function(mn, sd) function(n) rgamma(n, shape = gamma_alpha(mn,sd), scale = gamma_beta(mn,sd))

xformdistros <- lapply(read_json(.args[1], simplifyVector = T), function(xs) {
  list(shape=gamma_alpha(xs[1],xs[2]), scale=gamma_beta(xs[1], xs[2]))
})


dt <- rbindlist(mapply(function(pars, nm, samples) {
  data.table(
    tm = do.call(rgamma, c(pars, n=samples)),
    event=nm, sample = 1:samples
  )
}, pars=xformdistros, nm=names(xformdistros), MoreArgs = list(samples=1000), SIMPLIFY = F))

histop <- ggplot(dt) + aes(x=tm) + facet_grid(event ~ ., scales="free_y") + geom_histogram()

# what happens
tst <- dt[event == "hospital"][dt[event == "recover"], on=.(sample)][,
  .(tm=pmin(tm, i.tm), outcome=c("hospital","recover")[(tm<i.tm)+1])
]

histop2 <- ggplot(tst) + aes(x=tm, fill=outcome) + geom_histogram()

# what quantitle of the exposure distro the min event time corresponds to
# this is the probability that transmission will occur along any given connection
tst[, qtle := pgamma(
  tm, shape = xformdistros$exposure$shape, scale=xformdistros$exposure$scale
)]

tst[, n0 := sample(10, .N, rep = T)]
tst[, R0 := n0*qtle ]

histof <- ggplot(tst) + aes(x=R0) + geom_histogram()
tst[, mean(R0)]
