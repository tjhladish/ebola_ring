suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})
# digest a bunch of samples

pkgs <- lapply(c("data.table","ggplot2","cowplot"), require, character.only=T)

.args <- c(list.files(path="simple", pattern="\\d+\\.rds", full.names = T), "aggregate_simple.rds")
.args <- c(list.files(path="stratified", pattern="\\d+\\.rds", full.names = T), "aggregate_stratified.rds")
.args <- c(list.files(path="coverage_lvls", pattern="\\d+\\.rds", full.names = T), "aggregate_covlevels.rds")
.args <- commandArgs(trailingOnly = TRUE)

srcs <- head(.args, -1)
dms <- c("sample_id", key(readRDS(srcs[1])))

samples <- rbindlist(lapply(srcs, function(fn) {
  res <- readRDS(fn)[, sample_id := as.integer(gsub("^.+/(\\d+)\\.rds$","\\1",fn)) ]
  setkeyv(res, dms)
  # how to handle case where vaccine_pos == 0?
  # seems like lim vaccine_pos == 0 => varhat = 1/unvax pos + 1/unvax tot
  # TODO re-evaluate this as permutation test for small samples, then switch to normal approx?
  # res[, varhat := {
  #   ref <- cumsum(vaccine_pos)
  #   alt <- 1/cumsum(unvax_pos) + 1/cumsum(unvax_pos + unvax_neg)
  #   alt[ref != 0] <- alt[ref != 0] + 1/ref[ref != 0] + (1/cumsum(vaccine_pos + vaccine_neg))[ref != 0]
  #   sqrt(alt)
  # }, by = c(head(dms, -1))]
  res
}))

RRcalc <- function(vac_pos, vac_neg, non_pos, non_neg, recruit_fraction, B, gencov) {
  if (length(B)==1) B <- rep(B, length(vac_pos))
  ## browser()
  vac_neg_recruits <- rbinom(length(vac_neg), vac_neg, recruit_fraction)
  non_neg_recruits <- rbinom(length(non_neg), non_neg, recruit_fraction)
  back_non_recruits <- rbinom(length(B), B, 1-gencov)
  back_vac_recruits <- rbinom(length(B), B, gencov)
  return(list((
    cumsum(vac_pos)/cumsum(vac_neg_recruits + back_vac_recruits)
  )/(
    cumsum(non_pos)/cumsum(non_neg_recruits + back_non_recruits)
  ), vac_neg_recruits, non_neg_recruits, back_vac_recruits, back_non_recruits))
}

recruit_fraction <- 0.5
B <- 10
gencov <- 0.5
samples[, c("RR", "vac_neg_recruits", "non_neg_recruits", "back_vac_recruits", "back_non_recruits") := {
  RRcalc(vaccine_pos, vaccine_neg, unvax_pos, unvax_neg, recruit_fraction, B, gencov)
}, by=c(setdiff(dms,"cluster"))]
samples[, Veff := 1 - RR ]
plot.dt <- samples[,.(Veff.ave=mean(Veff), Veff.med=median(Veff)), by=c(dms[-1])]
ref.dt <- samples[cluster == 500,.(Veff.ave=mean(Veff)), by=c(dms[-1])]

filt <- expression(round(bveff %% 0.3, 1)==0)

ggplot(plot.dt[eval(filt)]) + aes(x=cluster, y=Veff.ave, color=factor(back_vac_mech), linetype=index, group=interaction(back_vac_mech, index)) + facet_grid(
  exp_sd ~ bveff
) + geom_line() +
  geom_hline(aes(yintercept=bveff), color = "black", alpha = 0.5) +
  geom_hline(aes(yintercept=Veff.ave, color=factor(back_vac_mech), linetype=index), ref.dt[eval(filt)], alpha = 0.5) +
  coord_cartesian(xlim = c(10,100), ylim=c(0,1)) +
  scale_color_discrete("Vaccine Mechanism", labels = c(`0`="Leaky", `1`="All-or-Nothing"))

# zeroslice <- samples[bveff == 0]
# 
# thingcmp <- samples[bveff != 0]
# 
# thing <- zeroslice[,.(tot = vaccine_pos + unvax_pos), by=.(exp_sd, trace_prob)]
# 
# thing[, mean(tot), by=.(exp_sd, trace_prob)]
# 
# ggplot(thing) + aes(x=tot) + facet_grid(trace_prob ~ exp_sd, labeller = facet_labels) + stat_bin(aes(y=log10(..count..))) 

zhat <- 1.96 # 95% interval z score

# # how to handle case where vaccine_pos == 0?
# # seems like lim vaccine_pos == 0 => varhat = 1/unvax pos + 1/unvax tot
# # TODO re-evaluate this as permutation test for small samples, then switch to normal approx?
# samples[, varhat := {
#   ref <- cumsum(vaccine_pos)
#   alt <- 1/cumsum(unvax_pos) + 1/cumsum(unvax_pos + unvax_neg)
#   alt[ref != 0] <- alt[ref != 0] + 1/ref[ref != 0] + (1/cumsum(vaccine_pos + vaccine_neg))[ref != 0]
#   sqrt(alt)
# }, by = c(head(dms, -1))]

samples[, Veff := 1 - RR ]
# samples[is.infinite(varhat) & RR == 0, Veff.lo := -Inf ]
# samples[!(is.infinite(varhat) & RR == 0), Veff.lo := 1-exp(log(RR)+zhat*varhat) ]
samples[, Veff.lo := 1-exp(log(RR)+zhat*varhat) ]
samples[, Veff.hi := 1-exp(log(RR)-zhat*varhat) ]

plot.dt <- samples[,.(Veff.ave=mean(Veff), Veff.med=median(Veff)), by=c(dms[-1])]

ggplot(plot.dt[round(bveff %% 0.3, 1)==0]) + aes(x=cluster, y=Veff.ave, color=factor(back_vac_mech), linetype=index, group=interaction(back_vac_mech, index)) + facet_grid(
  exp_sd ~ bveff
) + geom_line() + coord_cartesian(xlim = c(10,100), ylim=c(0,1))

samples[, pass0 := Veff.lo > 0 ]

power <- samples[, .(power=sum(pass0)/.N), by=c(grep("sample_id", dms, invert = T, value=T))]
power[is.na(power), power := 0]

required_samples <- power[, {
  fnd <- any(power > 0.8)
  .(cluster=if (fnd) as.numeric(which.max(power > 0.8)) else Inf)
}, by=c(grep("(sample_id|clusters)", dms, invert = T, value=T))]

plotter <- function(dt) ggplot(dt[!is.infinite(cluster)]) +
  aes(bveff, cluster, color=factor(back_vac_mech), linetype=factor(use_bias)) +
  facet_grid(exp_sd ~ trace_prob, labeller = label_bquote(
    cols = .(sprintf("%i%% Tracing", trace_prob*100)),
    rows = R[0]%~~%.(exp_sd/2-0.5)
  )) + geom_line() +
  coord_cartesian(ylim=c(0,100), xlim=c(0,1)) +
  scale_color_manual("Vaccine\nMechanism", labels=c(`0`="Leaky",`1`="All-or-Nothing"), values=c(`1`="firebrick",`0`="dodgerblue")) +
  scale_linetype_manual("Coverage Sampling", labels=c(`0`="Random",`1`="Adjusted"), values=c(`0`="dashed",`1`="solid")) +
  scale_x_continuous("Simulated Vaccine Efficacy") +
  scale_y_continuous("Clusters Required for >80% Power", breaks=seq(0,100,by=20)) +
  theme_minimal() + theme(
    panel.spacing.y = unit(1, "line"),
    panel.spacing.x = unit(1, "line")
  )

prow <- plotter(required_samples[exp_sd == 3])
save_plot("sample_size_row.png", prow, ncol = 4, nrow = 1, base_height = 3, base_aspect_ratio = 1)
pfull <- plotter(required_samples)
save_plot("sample_size_full.png", pfull, ncol = 4, nrow = 3, base_height = 3, base_aspect_ratio = 0.8)
