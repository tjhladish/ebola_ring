pkgs <- lapply(c("data.table","RSQLite", "jsonlite", "cowplot"), require, character.only=T)

.args <- c("../../tests/test_abc_example.sqlite", "../../tests/abc_ebola_sim.json", "stem-.png")
.args <- commandArgs(trailingOnly = T)

abcref <- jsonlite::read_json(.args[2])

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

coveragebinwidth <- 0.2

binnedcoverage <- pars[,{
  res <- hist(realized_coverage, breaks = seq(0,1,by=coveragebinwidth), plot=F)
  .(ulim=tail(res$breaks, -1), count=res$counts)
}, keyby=.(bveff, trace_prob, exp_sd, vaccine_delay)]

binnedcoverage[, share := count/sum(count), by=key(thing) ]

binnedplotter <- function(dt) ggplot(dt) + aes(
  x=bveff, y=share, fill=as.character(ulim)
) + facet_grid(trace_prob + vaccine_delay ~ exp_sd, labeller = labeller(
  exp_sd=function(sd) c(`3`="R0 = 1",`4`="R0 = 1.5",`5`="R0 = 2")[as.character(sd)],
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

nr <- binnedcoverage[,.(ind=1),by=.(trace_prob,vaccine_delay)][,sum(ind)]
nc <- binnedcoverage[,.(ind=1),by=.(exp_sd)][,sum(ind)]

fullp <- binnedplotter(binnedcoverage)
save_plot("coverage_bin_full.png", fullp, ncol = nc, nrow = nr, base_width = 2.5, base_height = 1.5)


zoomp <- binnedplotter(binnedcoverage[exp_sd==3 & trace_prob==0.5 & vaccine_delay == 2])
save_plot("coverage_bin_zoom.png", zoomp, ncol = 1, nrow = 1, base_width = 3.5*nc, base_height = 1*nr)

coverage_distro <- ggplot(pars[bveff %in% seq(0,1,by=0.2)]) +
  aes(realized_coverage, stat(density), color=factor(exp_sd), linetype=factor(trace_prob)) +
  facet_grid(. ~ bveff, labeller = labeller(.cols = function(b) {
    res <- paste0(c("\n","\n","Vaccine Efficacy\n","\n","\n"), sprintf("%i%%", as.numeric(b)*100))
    res
  })) +
  # geom_freqpoly(bins=100, alpha=0.5) +
  stat_bin(geom="point", bins=100, alpha=0.5) +
  theme_minimal() +
  scale_color_manual(expression("approx. "*R[0]), labels=c(`3`=1,`4`=1.5,`5`=2), values=c(`3`="black",`4`="dodgerblue",`5`="firebrick")) +
  scale_x_continuous("Coverage in Area of Observed Case") +
  scale_y_continuous("Relative Likelihood") +
  scale_linetype("Contact\nTracing\nProbability") +
  coord_cartesian(xlim=c(0.1,0.9)) + theme(
    legend.position = c(0.01, .99), legend.justification = c(0, 1),
    legend.direction = "horizontal"
  )

save_plot("coverage_bias.png", coverage_distro, ncol = 5, base_width = 1.5)

