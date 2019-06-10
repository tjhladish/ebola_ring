pkgs <- lapply(c("data.table","RSQLite", "jsonlite", "cowplot"), require, character.only=T)

.args <- c("../tests/test_abc_example.sqlite", paste0("../tests/work",1:6,".log"), "test_eff_calc.png")
.args <- commandArgs(trailingOnly = T)

abcref <- jsonlite::read_json("../tests/abc_ebola_sim.json")

extractcolname <- function(el) if(is.null(el$short_name)) el$name else el$short_name

parcols <- c("serial", sapply(abcref$parameters, extractcolname))
metcols <- sapply(abcref$metrics, extractcolname)

selcols <- c(parcols, metcols)

drv <- dbDriver("SQLite")
db <- dbConnect(drv, .args[1], flags=SQLITE_RO)
pars <- data.table(
  dbGetQuery(db,
    sprintf("SELECT %s FROM par JOIN met USING(serial) JOIN job USING(SERIAL) WHERE status=='D';", paste0(selcols, collapse = ", "))
  )
)
dbDisconnect(db)
pars[, replicate := net_rep*10 + epi_rep ]

thing <- pars[,{
  res <- hist(realized_coverage, breaks = seq(0,1,by=0.2), plot=F)
  .(ulim=tail(res$breaks, -1), count=res$counts)
}, keyby=.(bveff, trace_prob, exp_sd, vaccine_delay)]

thing[, share := count/sum(count), by=key(thing) ]

ggplot(thing) + aes(
  x=bveff, y=share, fill=as.character(ulim)
) + facet_grid(trace_prob + vaccine_delay ~ exp_sd, labeller = labeller(
  exp_sd=function(sd) c(`3`="R0 = 1",`4`="R0 = 1.5",`5`="R0 = 2"),
  trace_prob=function(p) {sprintf("contact tracing\np = %i%%", as.numeric(p)*100)},
  vaccine_delay=function(d) {sprintf("notification\nto vaccination\n%s days", d)}
)) + geom_col() +
  theme_minimal() +
  scale_y_continuous("Share of clusters observed by background coverage") +
  scale_x_continuous("Vaccine efficacy", breaks=seq(0,1,by=0.2)) +
  scale_fill_brewer(
    "Vaccine\nCoverage\nCategory",
    labels=function(b) sprintf("%i-%i%%",as.integer(100*(as.numeric(b)-0.2)), as.integer(100*as.numeric(b))),
    palette="Greens"
  )

coverage_distro <- ggplot(pars[bveff %in% seq(0,1,by=0.2)]) +
  aes(realized_coverage, stat(density), color=factor(exp_sd), linetype=factor(trace_prob)) +
  facet_grid(. ~ bveff, labeller = labeller(.cols = function(b) {
    res <- paste0(c("\n","\n","Vaccine Efficacy\n","\n","\n"), sprintf("%i%%", as.numeric(b)*100))
    res
  })) +
  geom_freqpoly(bins=100, alpha=0.5) + theme_minimal() +
  scale_color_manual(expression("approx. "*R[0]), labels=c(`3`=1,`4`=1.5,`5`=2), values=c(`3`="black",`4`="dodgerblue",`5`="firebrick")) +
  scale_x_continuous("Coverage in Area of Observed Case") +
  scale_y_continuous("Relative Likelihood") +
  scale_linetype("Contact\nTracing\nProbability") +
  coord_cartesian(xlim=c(0.1,0.9)) + theme(
    legend.position = c(0.01, .99), legend.justification = c(0, 1),
    legend.direction = "horizontal"
  )

save_plot("coverage_bias.png", coverage_distro, ncol = 5, base_width = 1.5)

# logs <- grep("log$", .args, value=T)
# 
# aggsims <- rbindlist(lapply(
#   logs,
#   function(log, ...) {
#     fread(log, ...)[level == 1 | !is.nan(onset), {
#       rvtm <- if(!all(is.nan(ring))) min(ring, na.rm = T) else 0
#       .SD[, .N, keyby=.(background, onset=round(onset-rvtm))]
#     }, by=serial]
#   },
#   col.names = c("serial","id","level","background","trace","ring","onset")
# ))
# 
# ref.dt <- aggsims[pars[,.(epi_rep,bveff,trace_prob,exp_sd,realized_coverage,count_pre_vax),by=serial], on=.(serial)]

id.vars <- c("serial","bveff","trace_prob","exp_sd")

testpos <- melt.data.table(
  pars[,.(vaccine=vaccine_pos_0, none=unvax_pos_0, coverage=realized_coverage), by=id.vars],
  id.vars = c(id.vars,"coverage"), variable.name = "background", value.name = "N"
)[, outcome := "positive" ]
testneg <- melt.data.table(
  pars[,.(vaccine=vaccine_neg_0, none=unvax_neg_0, coverage=realized_coverage), by=id.vars],
  id.vars = c(id.vars,"coverage"), variable.name = "background", value.name = "N"
)[, outcome := "negative" ]

# testpos <- ref.dt[
#   !is.nan(onset) & onset <= 0,
#   .(N=sum(N)),
#   by = .(serial, background=c("none","vaccine")[background+1], bveff, trace_prob, exp_sd)
# ][, outcome := "positive"] 
# testneg <- ref.dt[is.nan(onset),N,by=.(serial, background=c("none","vaccine")[background+1], bveff, trace_prob, exp_sd)][, outcome := "negative"]

results <- dcast(
  rbind(testpos, testneg),
  serial + bveff + trace_prob + exp_sd + coverage ~ background + outcome,
  value.var = "N", fill = 0
)[!(vaccine_negative+vaccine_positive+none_negative+none_positive == 0)]

# ggplot(
#   melt(results[exp_sd==3 & bveff==0.3,.SD,.SDcols=-c(2,4)], id.vars=c("serial","coverage","trace_prob"))
# ) + facet_grid(variable ~ trace_prob, scales = "free") + aes(x=value) + geom_histogram(bins=100) +
#   coord_cartesian(xlim=c(0,15))

samp <- results[exp_sd==3, {
  set.seed(1L)
  randomize <- sample(.N, replace = F)
  res <- .SD[randomize]
  alt <- lapply(res, cumsum)
  names(alt) <- paste0("c.", names(alt))
  c(res, alt, .(clusters=1:.N))
}, by=.(bveff, trace_prob, exp_sd), .SDcols=-c("coverage", "serial")]

slice <- samp[,c(.SD,.(neg_rate=seq(.2,.8,by=.2))),by=.(bveff,trace_prob,exp_sd,clusters)]

slice[, RR := (
    c.vaccine_positive/(c.vaccine_positive+round(c.vaccine_negative*neg_rate))
  )/(
    c.none_positive/(c.none_positive+round(c.none_negative*neg_rate))
  )
]

ggplot(slice[between(clusters,10,100) & trace_prob %in% c(.5,.7,.9, 1)]) + aes(clusters/10, 1-RR, color="simulated", alpha=neg_rate, group=neg_rate) + 
  facet_grid(trace_prob ~ bveff, labeller = labeller(.cols = function(b) {
    sprintf("%i%%", as.numeric(b)*100)
  })) +
  geom_hline(aes(yintercept=1-RR, color="asymptote"), slice[clusters > 2000 & trace_prob %in% c(.5,.7,.9, 1), .(RR=mean(RR)), by=.(bveff, trace_prob, neg_rate)]) +
  geom_hline(aes(yintercept=bveff, color="actual")) +
  geom_line() + theme_minimal() + labs(
    y=expression('Measured '*v[eff]*' = 1-RR'),
    x="10s of clusters observed",
    color="Efficacy\nMeasure"
  ) + ggtitle("Simulated Vaccine Efficacy") +
  scale_alpha(breaks=c(.2,.4,.6,.8)) + scale_x_continuous(breaks=seq(2,10,by=2)) +
  scale_color_manual(
    labels=c(actual="Reference", simulated="Measured", asymptote="Measured\nafter 2k clusters"),
    values=c(actual="black", simulated="firebrick", asymptote="dodgerblue")
  ) + coord_cartesian(ylim=c(-0.5,1)) + theme(
    panel.grid = element_line(),
    panel.grid.minor.x = element_blank()
  )



samp <- rbindlist(lapply(1:200, function(sid, negattackrate, maxsamples=100) results[, {
  ord <- order(sample(.N, maxsamples))
  # A = vac, +; B = no, +
  # C = vac, -; D = no, -
  cum.A <- cumsum(vaccine_positive[ord])
  cum.B <- cumsum(none_positive[ord])
  vacneg <- rbinom(maxsamples, vaccine_negative[ord], negattackrate)
  notneg <- rbinom(maxsamples, none_negative[ord], negattackrate)
  cum.C <- cumsum(vacneg)
  cum.D <- cumsum(notneg)
  RR <- (cum.A/(cum.A+cum.C))/(cum.B/(cum.B+cum.D))
  z <- 1.96*sqrt( (cum.C/cum.A)/(cum.A+cum.C) + (cum.D/cum.B)/(cum.B+cum.D) )
  .(clusters=1:maxsamples, RR=RR, z=z)
}, by=.(bveff, trace_prob, exp_sd)][, sample_id := sid ][, negattackrate := negattackrate ], negattackrate = 0.1))

samp[, veff := 1-RR ][, veffl := 1-exp(log(RR)+z) ][, veffh := 1-exp(log(RR)-z) ][, TP0 := (bveff > 0) & (veffl>0) ][, TP50 := (bveff > 0.5) & (veffl > 0.5)]

pwr <- samp[,.(
    power0=sum(TP0, na.rm = T)/sum(!is.na(TP0)),
    power50=sum(TP50, na.rm = T)/sum(!is.na(TP50))
  ),
  keyby=.(bveff, trace_prob, exp_sd, clusters)
]

ggplot(pwr[bveff !=0 & bveff!=1][]) + aes(x=clusters, color=factor(exp_sd)) +
  facet_grid(trace_prob ~ bveff) +
  geom_line(aes(y=power0, linetype="vs 0")) +
  geom_line(aes(y=power50, linetype="vs 0.5")) +
  scale_linetype("Comparison...") +
  scale_color_manual("R0 roughly", labels=c(`3`=1,`4`=1.5,`5`=2), values=c(`3`="black",`4`="dodgerblue",`5`="firebrick")) +
  scale_y_continuous("Power") +
  scale_x_continuous("# of Rings") +
  theme_minimal()

q.dt <- samp[,{
  qs <- quantile(veff, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T)
  names(qs) <- c("lo95","lo50","med","hi50","hi95")
  as.list(qs)
}, keyby=.(bveff, trace_prob, exp_sd, clusters)]

nrs <- samp[,length(unique(trace_prob))]
ncs <- samp[,length(unique(bveff))] - 1

qp <- ggplot(q.dt[bveff != 1.0]) + aes(
  x=log10(clusters), fill=factor(exp_sd),
  group=interaction(trace_prob, bveff, exp_sd)) +
  facet_grid(trace_prob ~ bveff) +
  geom_hline(aes(yintercept = bveff, linetype="reference"), color="black") +
  geom_ribbon(aes(alpha="95", ymin=lo95, ymax=hi95)) +
  geom_ribbon(aes(alpha="50", ymin=lo50, ymax=hi50)) +
  geom_line(aes(y=med, color=factor(exp_sd), linetype="observed")) +
  scale_y_continuous("Estimated Veff") +
  scale_x_continuous("log10(# of Rings)", breaks = log10(c(10,18,32,56,100)), labels = function(b) round(10^b)) +
  scale_color_discrete("Transmission\nScenario") +
  scale_alpha_manual("Interval", values=c("95"=0.1,"50"=0.2)) +
  coord_cartesian(ylim=c(0,1), xlim=c(1,2)) +
  theme_minimal()

save_plot(tail(.args, 1), qp, nrow = nrs, ncol = ncs)

p <- ggplot(samp) + aes(
  x=log10(clusters), y=veff, color=factor(exp_sd),
  group=interaction(trace_prob, bveff, exp_sd, sample_id)) +
  facet_grid(trace_prob ~ bveff) +
  geom_hline(aes(yintercept = bveff, linetype="reference"), color="black") +
  geom_line(aes(linetype="observed"), alpha=0.2) +
  scale_y_continuous("Estimated Veff") +
  scale_x_continuous("log10(# of Rings)") +
  scale_color_discrete("Transmission\nScenario") +
  coord_cartesian(ylim=c(-0.5,1), xlim=c(1,3)) +
  theme_minimal()

save_plot(tail(.args, 1), p, nrow = nrs, ncol = ncs)

plot.dt <- ref.dt[!is.nan(onset) & !is.na(background), {
  bnds <- range(onset, na.rm = T)
  refind <- bnds[1]-1
  ind <- onset - refind
  if (any(is.infinite(bnds))) browser()
  full <- seq(bnds[1], bnds[2])
  cases <- rep(0, length.out=length(full))
  cases[ind] <- N
  .(onset=full, cases=cases, cc=cumsum(cases))
}, by=.(serial, bveff, trace_prob, exp_sd, background=c("none","prophylactic")[background+1])]

p <- ggplot(plot.dt[trace_prob==0.75 & bveff == 0.5]) + aes(x=onset, y=cc, color=factor(background), group=serial) +
  facet_grid(bveff ~ trace_prob) +
  geom_step(alpha=0.1) + geom_rug(sides="r", alpha=0.1) +
  scale_color_manual("Background\nVaccination",values=c(none="firebrick", prophylactic="dodgerblue")) +
  theme_minimal()

hist.dt <- plot.dt[onset <= 0,.(peak=max(cc)),by=.(serial, bveff, trace_prob, exp_sd, background)]

p <- ggplot(hist.dt[trace_prob == 0.5 & peak != 0]) + aes(x=log(peak), fill=factor(background)) +
  facet_grid(trace_prob ~ bveff + exp_sd) +
  geom_histogram() +
  scale_y_log10() +
  scale_fill_manual("Background\nVaccination",values=c(none="firebrick", prophylactic="dodgerblue")) +
  theme_minimal()
