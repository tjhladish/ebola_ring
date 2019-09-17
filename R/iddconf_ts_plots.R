suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
  require(ggplot2)
  require(cowplot)
})

.args <- c("timedistros.json", "timedistros.png")
.args <- commandArgs(trailingOnly = T)

pars <- read_json(.args[1])

gamma_alpha <- function(mn, sd) (mn/sd)^2
gamma_beta <- function(mn, sd) (sd^2)/mn

rd <- function(mn, sd) {
  shape = gamma_alpha(mn, sd)
  scale = gamma_beta(mn, sd)
  list(r=function(n) rgamma(n, shape=shape, scale=scale), d=function(x) dgamma(x, shape=shape, scale=scale))
}

rdhosp <- rd(pars$hospital[[1]], pars$hospital[[2]])
rdinf <- rd(pars$incubate[[1]], pars$incubate[[2]])
rdrem <- rd(pars$recover[[1]], pars$recover[[2]])
rdexp <- rd(pars$exposure[[1]], pars$exposure[[2]]+1)

funref <- list(hospitalisation=rdhosp, removal=rdrem, onset=rdinf, exposure=rdexp)

refgrid <- data.table(expand.grid(time=seq(0,20,by=0.01), distro=c("hospitalisation","removal","onset","exposure")))
refgrid[, p := {
  FUN <- funref[[as.character(distro)]]
  FUN$d(time)
},by=distro]

refgrid[, cump := cumsum(p*0.01), by=distro]
iqs <- refgrid[, .(lo=time[which.min(cump < 0.25)], hi=time[which.max(cump > 0.75)]), by=distro]

arearef <- refgrid[iqs, on=.(distro)][between(time, lo, hi)]

alldistros <- ggplot(refgrid[distro!="onset"]) + aes(time, p, color=distro) +
  geom_area(aes(fill=distro), arearef[distro!="onset"], position = "identity", alpha=0.5, linetype="dashed", show.legend = F) +
  geom_line() +
  scale_color_manual(
    "Waiting Time Distributions",
    labels=c(hospitalisation="Onset to Hospital", removal="Onset to Other Removal", exposure="Onset to Infectious Contact"),
    breaks=c("hospitalisation","removal","exposure"),
    values=c(hospitalisation="dodgerblue", removal="forestgreen", exposure="firebrick"),
    aesthetics = c("color","fill")
  ) +
  scale_x_continuous("Days after Initiating Event", expand=c(0,0)) +
  scale_y_continuous("Density", expand=c(0,0), labels = scales::percent) +
  theme_minimal_grid() + theme(
    axis.line.y.left = element_line(),
    axis.line.x.bottom = element_line(),
    legend.position = c(.97,1),
    legend.justification = c(1,1),
    plot.margin = margin(10,15,10,10)
  )

save_plot(tail(.args, 1), alldistros, base_asp = 1.75, base_width = 10)
