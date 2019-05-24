require(jsonlite)
require(data.table)
require(cowplot)

.args <- c("timedistros.json", "../tests/networks")
.args <- commandArgs(trailingOnly = T)

netpath <- .args[2]
nets <- list.files(netpath, pattern=".csv", full.names = T)

require(igraph)

p0k <- rle(sort(sapply(nets, function(pth) {
  dt <- fread(pth)
  if (dim(dt)[2] == 2) {
    g <- igraph::graph_from_data_frame(dt, directed = F)
    unname(degree(g, V(g)[1]))
  } else 0
}, USE.NAMES = F)))

gamma_alpha <- function(mn, sd) (mn/sd)^2
gamma_beta <- function(mn, sd) (sd^2)/mn
adgamma <- function(mn, sd) function(n) rgamma(n, shape = gamma_alpha(mn,sd), scale = gamma_beta(mn,sd))

# standard is wikipedia k, theta formulation
# so mean = k*theta
# sd^2 = k theta^2 = mn * theta => theta = sd^2/mn
# and k => (mn/sd)^2
# cpp gamma_distribution uses alpha = k, beta = theta
# rgamma uses shape = k, scale = theta

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

tst[, n0 := sample(p0k$values, .N, rep = T, prob = p0k$lengths)]
tst[, Robs := rbinom(.N,n0,qtle) ]

histof <- ggplot(tst) + aes(x=Robs) + geom_histogram()
tst[, .(med=median(Robs),mn=mean(Robs))]
