suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
  require(ggplot2)
  require(cowplot)
})

.args <- c("digestlogs.rds", "../tests/abc_ebola_km.sqlite", "km_plot.rda")
.args <- commandArgs(trailingOnly = T)

dt <- readRDS(.args[1])
drv <- dbDriver("SQLite")
db <- dbConnect(drv, .args[2], flags=SQLITE_RO)
pars <- data.table(dbGetQuery(
    db,
    "SELECT P.* FROM par P JOIN job USING(SERIAL) WHERE status=='D';"
))
dbDisconnect(db)
## these are the serials corresponding to demo network
# baseline.sers <- pars[net_rep == 1 & bveff == 0 & exp_sd == 4, serial]
intervention.sers <- pars[
  (back_vac_mech == 0) & (exp_sd == 4) & (coverage == 0.75) & (bveff != 0.0)
][, .(serial, bveff, netsrc = factor(
  c("Realistic", "De-Clustered", "Randomized")[netsrc+1],
  levels = c("De-Clustered", "Randomized", "Realistic"),
  ordered = T
)) ]

intervention.dt <- dt[
  intervention.sers, on=.(serial), c("bveff","netsrc") := .(bveff, netsrc)][!is.na(bveff)]
ref.pop <- intervention.dt[infection_time > trace_time, .N, by=.(has_background, bveff, netsrc)]
# baseline.dt[,sum(!is.infinite(infection_time))>1,by=serial][V1==TRUE, .N]
# intervention.dt[,sum(!is.infinite(infection_time))>1,by=serial][V1==TRUE, .N]

cases <- intervention.dt[!is.infinite(infection_time - trace_time)][infection_time > trace_time]
#cases[, length(unique(serial))]
#tlim <- 25
#filt <- expression(V1 < tlim)
km <- cases[,
  infection_time - trace_time, by=.(serial, has_background, bveff, netsrc)
][V1 < 25, .(time=sort(V1)), by=.(serial, has_background, bveff, netsrc)]
plot.dt <- km[, .(time=sort(time), cases=1:.N), keyby=.(bveff, has_background, netsrc)]
plot.dt[ref.pop, arper := cases/N*100, on=.(bveff, has_background, netsrc)]

calc <- plot.dt[, arper[which.min(time < 10)], by=.(bveff, netsrc, has_background)]

non <- calc[has_background==0, .(ar=V1), keyby=.(bveff, netsrc)]
vac <- calc[has_background==1, .(ar=V1), keyby=.(bveff, netsrc)]

veff <- vac[non, .(Veff = 1 - ar/i.ar, bveff, netsrc)]

ulim <- 1.0

extension.dt <- plot.dt[,.(time=25, cases=cases[.N], arper=arper[.N]),by=.(bveff, has_background, netsrc)]

p1 <- ggplot(rbind(plot.dt, extension.dt)) + facet_grid(netsrc ~ .) +
  aes(time, arper, color=factor(has_background), alpha=bveff, group=interaction(bveff, has_background)) +
  annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, fill="grey65", alpha=0.3) +
  geom_step(size=1.5) +
  geom_text(
    aes(25-2, label=sprintf("%02i%% VE", bveff*100), alpha=NULL),
    data=function(dt) dt[
      has_background == 1,
      .(arper = arper[which.min(time < (25-1))]-0.05),
      by=.(bveff, has_background)
  ], show.legend = F) +
  coord_cartesian(x=c(0,25), y=c(0,ulim)) +
  scale_y_continuous("Prevalence of Confirmed EVD (%)", expand=c(0,0), breaks=seq(0,ulim,by=0.1), position = "right") +
  scale_x_continuous("Time after contact tracing (days)", expand=c(0,0)) +
  scale_alpha_continuous(guide="none", range = c(0.4, 1)) +
  scale_color_manual(
    "Vaccine Status",
    values = c(`0`="black", `1`="forestgreen"),
    labels=c(`0`="Unvaccinated", `1`="Vaccinated"),
    guide = guide_legend(override.aes = list(size=2))
  ) +
  theme(
    legend.position = c(0.25, .9),
    panel.background = element_blank(),
    legend.key = element_blank(),
    panel.grid = element_line(color="grey85"),
    rect = element_rect(fill="transparent"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing.y = unit(36, "pt"),
    axis.line.y.right = element_line(),
    axis.line.x.bottom = element_line(),
    axis.title.y.right = element_text(margin = margin(l=10))
  )

saveRDS(p1,tail(.args,1))
