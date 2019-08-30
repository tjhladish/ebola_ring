pkgs <- lapply(c("data.table","cowplot"), require, character.only=T)

.args <- c(list.files(path="simple", pattern="\\d+\\.rds", full.names = T), "thing.png")
.args <- commandArgs(trailingOnly = T)

srcs <- head(.args, -1)
samples <- rbindlist(lapply(srcs, function(fn) {
  readRDS(fn)[, fn := fn ]
}))
samples[, sample_id := as.integer(gsub("^.+/(\\d+)\\.rds$","\\1",fn)) ]

setkeyv(samples, c("bveff", "trace_prob", "exp_sd", "vaccine_delay", "back_vac_mech", "use_bias", "window", "clusters"))

subview <- samples[trace_prob == .6 & bveff != 1 ]

bnd <- subview[,{
 res <- quantile(1-RR, probs = (1:3)/4, na.rm = T)
 names(res) <- c("lo.Veff","md.Veff","hi.Veff")
 as.list(res)
}, by=key(samples)]
 
asymp <- bnd[, .(Veff=mean(md.Veff[clusters > .9*max(clusters)])), by=c(setdiff(key(samples),"clusters")) ]

thing <- samples[, sum(vaccine_neg + unvax_neg) / sum(vaccine_pos + unvax_pos), by=c(setdiff(key(samples),"clusters"))]

ggplot(thing) + aes(V1) + facet_grid(trace_prob ~ exp_sd) + geom_histogram()

fig1filt <- expression(bveff %in% c(0.2,0.4,0.6,0.8) & exp_sd == 3 & use_bias == 1)

coreplotter <- function(dt, asymp.dt, filt) ggplot(dt[clusters <= 100][eval(filt)]) +
  aes(x=clusters, color = interaction(back_vac_mech, use_bias)) +
  facet_grid(exp_sd ~ bveff, labeller = labeller(.cols = function(b) {
    sprintf("%i%%", as.numeric(b)*100)
  }, .rows = function(r) rep("", length(r))
  )) +
  geom_hline(aes(yintercept=bveff, color="reference"), linetype="dashed") +
  geom_hline(aes(yintercept=0, color="reference"), linetype="dashed") +
  geom_ribbon(aes(ymin=lo.Veff, ymax=hi.Veff, fill=interaction(back_vac_mech, use_bias), color=NULL), alpha=0.3) +
  geom_hline(aes(yintercept=Veff, color = interaction(back_vac_mech, use_bias)), alpha=0.5, asymp.dt[eval(filt)]) +
  geom_text(aes(x=80, y=Veff, label="Asymptotic\nEstimate"), asymp.dt[eval(filt)], alpha=0.5, show.legend = F) +
  geom_line(aes(y=md.Veff), dt[between(clusters,10,100)][eval(filt)]) + theme_minimal() + labs(
    y=expression('Estimated '*v[eff]*' = 1-RR'),
    x="# of clusters observed"
  ) + ggtitle("Simulated Vaccine Efficacy") +
  scale_alpha(breaks=c(.2,.4,.6,.8)) + scale_x_continuous(breaks=seq(10,100,by=10)) +
  scale_color_manual("Mechanism x Coverage",
                     labels=c(reference="(reference)", `0.0`="Leaky,\nRandom Coverage",`0.1`="Leaky,\nSampled Coverage",`1.0`="All-or-nothing,\nRandom Coverage",`1.1`="All-or-nothing,\nSampled Coverage"),
                     values=c(reference="black", `0.0`="dodgerblue", `1.1`="firebrick", `0.1`="blue", `1.0`="red"),
                     aesthetics = c("color", "fill")
  ) +
  coord_cartesian(ylim=c(-0.5,1)) + theme(
    panel.grid = element_line(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  ) + scale_size(guide = "none")

effplot1 <- coreplotter(bnd, asymp, fig1filt)
effplot1
save_plot("example_efficacy.png", effplot1, ncol = 4, nrow = 1, base_height = 5, base_aspect_ratio = 0.6)


sifig1filt <- expression(bveff %in% c(0.2,0.4,0.6,0.8))

effplot2 <- coreplotter(bnd, asymp, sifig1filt)
effplot2
save_plot("scan_efficacy.png", effplot2, ncol = 4, nrow = 3, base_height = 3, base_aspect_ratio = 0.8)
# 
# save_plot("nonleaky_sim_vax_eff_lines.png", effplot1, ncol = 4, nrow = 4, base_height = 4)
# save_plot("nonleaky_sim_vax_eff_ribbon.png", effplot2, ncol = 9, nrow = 4, base_height = 4)
# 
# samp <- rbindlist(lapply(1:200, function(sid, negattackrate, maxsamples=100) results[, {
#   ord <- order(sample(.N, maxsamples))
#   # A = vac, +; B = no, +
#   # C = vac, -; D = no, -
#   cum.A <- cumsum(vaccine_positive[ord])
#   cum.B <- cumsum(none_positive[ord])
#   vacneg <- rbinom(maxsamples, vaccine_negative[ord], negattackrate)
#   notneg <- rbinom(maxsamples, none_negative[ord], negattackrate)
#   cum.C <- cumsum(vacneg)
#   cum.D <- cumsum(notneg)
#   RR <- (cum.A/(cum.A+cum.C))/(cum.B/(cum.B+cum.D))
#   z <- 1.96*sqrt( (cum.C/cum.A)/(cum.A+cum.C) + (cum.D/cum.B)/(cum.B+cum.D) )
#   .(clusters=1:maxsamples, RR=RR, z=z)
# }, by=.(bveff, trace_prob, exp_sd)][, sample_id := sid ][, negattackrate := negattackrate ], negattackrate = 0.1))
# 
# samp[, veff := 1-RR ][, veffl := 1-exp(log(RR)+z) ][, veffh := 1-exp(log(RR)-z) ][, TP0 := (bveff > 0) & (veffl>0) ][, TP50 := (bveff > 0.5) & (veffl > 0.5)]
# 
# pwr <- samp[,.(
#     power0=sum(TP0, na.rm = T)/sum(!is.na(TP0)),
#     power50=sum(TP50, na.rm = T)/sum(!is.na(TP50))
#   ),
#   keyby=.(bveff, trace_prob, exp_sd, clusters)
# ]
# 
# ggplot(pwr[bveff !=0 & bveff!=1][]) + aes(x=clusters, color=factor(exp_sd)) +
#   facet_grid(trace_prob ~ bveff) +
#   geom_line(aes(y=power0, linetype="vs 0")) +
#   geom_line(aes(y=power50, linetype="vs 0.5")) +
#   scale_linetype("Comparison...") +
#   scale_color_manual("R0 roughly", labels=c(`3`=1,`4`=1.5,`5`=2), values=c(`3`="black",`4`="dodgerblue",`5`="firebrick")) +
#   scale_y_continuous("Power") +
#   scale_x_continuous("# of Rings") +
#   theme_minimal()
# 
# q.dt <- samp[,{
#   qs <- quantile(veff, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = T)
#   names(qs) <- c("lo95","lo50","med","hi50","hi95")
#   as.list(qs)
# }, keyby=.(bveff, trace_prob, exp_sd, clusters)]
# 
# nrs <- samp[,length(unique(trace_prob))]
# ncs <- samp[,length(unique(bveff))] - 1
# 
# qp <- ggplot(q.dt[bveff != 1.0]) + aes(
#   x=log10(clusters), fill=factor(exp_sd),
#   group=interaction(trace_prob, bveff, exp_sd)) +
#   facet_grid(trace_prob ~ bveff) +
#   geom_hline(aes(yintercept = bveff, linetype="reference"), color="black") +
#   geom_ribbon(aes(alpha="95", ymin=lo95, ymax=hi95)) +
#   geom_ribbon(aes(alpha="50", ymin=lo50, ymax=hi50)) +
#   geom_line(aes(y=med, color=factor(exp_sd), linetype="observed")) +
#   scale_y_continuous("Estimated Veff") +
#   scale_x_continuous("log10(# of Rings)", breaks = log10(c(10,18,32,56,100)), labels = function(b) round(10^b)) +
#   scale_color_discrete("Transmission\nScenario") +
#   scale_alpha_manual("Interval", values=c("95"=0.1,"50"=0.2)) +
#   coord_cartesian(ylim=c(0,1), xlim=c(1,2)) +
#   theme_minimal()
# 
# save_plot(tail(.args, 1), qp, nrow = nrs, ncol = ncs)
# 
# p <- ggplot(samp) + aes(
#   x=log10(clusters), y=veff, color=factor(exp_sd),
#   group=interaction(trace_prob, bveff, exp_sd, sample_id)) +
#   facet_grid(trace_prob ~ bveff) +
#   geom_hline(aes(yintercept = bveff, linetype="reference"), color="black") +
#   geom_line(aes(linetype="observed"), alpha=0.2) +
#   scale_y_continuous("Estimated Veff") +
#   scale_x_continuous("log10(# of Rings)") +
#   scale_color_discrete("Transmission\nScenario") +
#   coord_cartesian(ylim=c(-0.5,1), xlim=c(1,3)) +
#   theme_minimal()
# 
# save_plot(tail(.args, 1), p, nrow = nrs, ncol = ncs)
# 
# plot.dt <- ref.dt[!is.nan(onset) & !is.na(background), {
#   bnds <- range(onset, na.rm = T)
#   refind <- bnds[1]-1
#   ind <- onset - refind
#   if (any(is.infinite(bnds))) browser()
#   full <- seq(bnds[1], bnds[2])
#   cases <- rep(0, length.out=length(full))
#   cases[ind] <- N
#   .(onset=full, cases=cases, cc=cumsum(cases))
# }, by=.(serial, bveff, trace_prob, exp_sd, background=c("none","prophylactic")[background+1])]
# 
# p <- ggplot(plot.dt[trace_prob==0.75 & bveff == 0.5]) + aes(x=onset, y=cc, color=factor(background), group=serial) +
#   facet_grid(bveff ~ trace_prob) +
#   geom_step(alpha=0.1) + geom_rug(sides="r", alpha=0.1) +
#   scale_color_manual("Background\nVaccination",values=c(none="firebrick", prophylactic="dodgerblue")) +
#   theme_minimal()
# 
# hist.dt <- plot.dt[onset <= 0,.(peak=max(cc)),by=.(serial, bveff, trace_prob, exp_sd, background)]
# 
# p <- ggplot(hist.dt[trace_prob == 0.5 & peak != 0]) + aes(x=log(peak), fill=factor(background)) +
#   facet_grid(trace_prob ~ bveff + exp_sd) +
#   geom_histogram() +
#   scale_y_log10() +
#   scale_fill_manual("Background\nVaccination",values=c(none="firebrick", prophylactic="dodgerblue")) +
#   theme_minimal()
