suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
  require(ggplot2)
  require(cowplot)
})

.args <- c("digestlogs.rds", "../tests/abc_ebola_km.sqlite", "iddconf_efficacy.png")
.args <- commandArgs(trailingOnly = T)

drv <- dbDriver("SQLite")
db <- dbConnect(drv, .args[2], flags=SQLITE_RO)
pars <- data.table(dbGetQuery(
    db,
    "SELECT P.*, count_pre_vax FROM par P JOIN met USING(serial) JOIN job USING(serial) WHERE status=='D';"
))
dbDisconnect(db)
## these are the serials corresponding to demo network
# baseline.sers <- pars[net_rep == 1 & bveff == 0 & exp_sd == 4, serial]
baseline.sers <- pars[
  (back_vac_mech == 0) & (exp_sd == 4) & (bveff == 0.0)
][, .(
  net_rep, epi_rep, serial, coverage, bveff,
  # indices=c("multi","one")[(count_pre_vax == 1)+1],
  netsrc = factor(
    c("Realistic", "De-Clustered", "Randomized")[netsrc+1],
    levels = c("De-Clustered", "Randomized", "Realistic"),
    ordered = T
  )) ]
baseline.sers

intervention.sers <- pars[
  (back_vac_mech == 0) & (exp_sd == 4) & (bveff != 0.0)
][, .(
  net_rep, epi_rep, serial, coverage, bveff,
  #indices=c("multi","one")[(count_pre_vax == 1)+1],
  netsrc = factor(
  c("Realistic", "De-Clustered", "Randomized")[netsrc+1],
  levels = c("De-Clustered", "Randomized", "Realistic"),
  ordered = T
)) ]
intervention.sers

#compref <- intervention.sers[baseline.sers[,.SD,.SDcols=-c("bveff")], on=.(net_rep, epi_rep, coverage, netsrc)]

rm(pars)

dt <- readRDS(.args[1])

digested <- dt[, .(
  cases = sum(!is.infinite(infection_time))
), by=serial]
digested

rm(dt)

int.slice <- digested[intervention.sers, on = .(serial=serial), .(net_rep, epi_rep, coverage, bveff, netsrc, cases)]
int.slice[is.na(cases), cases := 0L ]

base.slice <- digested[baseline.sers, on = .(serial=serial), .(net_rep, epi_rep, coverage, bveff=0, netsrc, cases)]
base.slice[is.na(cases), cases := 0L ]

comp <- int.slice[base.slice[,.SD,.SDcols=-c("bveff")], on=.(net_rep, epi_rep, coverage, netsrc)]

comp[, ind.eff := ifelse(cases == i.cases, 0, (i.cases-cases)/i.cases) ]

# ggplot(melt(
#   intervention.dt,
#   id.vars = c("serial", "coverage", "bveff", "indices", "netsrc"),
#   measure.vars = c("cases","pop"))[variable == "pop"]
# ) +
#   aes(value, color=factor(coverage), linetype=factor(bveff), group=interaction(coverage, bveff)) + facet_grid(netsrc ~ indices) +
#   geom_density()
# 
# ggplot(intervention.dt[, .N, by=.(bveff, coverage, netsrc)]) +
#   aes(bveff, N) + facet_grid(netsrc ~ coverage) + geom_line() +
#   scale_x_continuous(breaks=seq(1,3)/4)

# ggplot(comp[!(cases == 0 & i.cases == 0)]) + aes(bveff, ind.eff, color=netsrc) + facet_grid(netsrc ~ coverage) + geom_boxplot()

eff.dt <- comp[,{
  uncensored <- sum(!((i.cases == 0) & (cases == 0)))
  eff.est <- (sum(i.cases)-sum(cases))/sum(i.cases)
  ci <- as.numeric(binom.test(round(uncensored*c(eff.est,1-eff.est)), n=uncensored, p=eff.est)$conf.int)
  .(eff = eff.est, lo=ci[1], hi=ci[2])
}, by=.(coverage, bveff, netsrc)]

p1 <- ggplot(eff.dt) + aes(bveff, eff, color=netsrc) +
  facet_grid(. ~ coverage, labeller = labeller(
    coverage=function(l) paste0(c("\n","Intervention Coverage\n","\n"),scales::percent(as.numeric(l), accuracy = 1)))
  ) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=netsrc, color = NULL), alpha=0.2) +
  geom_line() +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  scale_y_continuous("Intervention Effectiveness\n(Fraction of Cases Prevented)", expand=c(0,0)) +
  scale_x_continuous("Intervention Efficacy", expand=c(0,0), labels = scales::percent) +
  scale_color_discrete("Network Structure", labels=c(
    `De-clustetered`="Realistic Population Only",
    `Reclustered`="Realistic Degree Only",
    `Realistic`="Realistic Clustering"
  ), aesthetics = c("color","fill")) +
  theme_minimal() + theme(
    panel.spacing = unit(24,"pt"),
    legend.position = c(.1, .9),
    legend.justification = c(0, 1)
  )

save_plot(tail(.args,1), p1, nrow=3, ncol=3)
