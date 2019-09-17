pkgs <- lapply(c("data.table"), require, character.only=T)

.args <- c("digested_nobias.rds", "simple_sum.png")
.args <- commandArgs(trailingOnly = T)

recast <- readRDS(.args[1])

tarfile <- tail(.args, 1)

res <- recast[window == 0, {
  cases <- vaccine_pos+unvax_pos
  qs <- quantile(cases, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(qs) <- c("lo","lo.q","md","hi.q","hi")
  c(as.list(qs), mn=mean(cases))
}, by=.(exp_sd, bveff, back_vac_mech)]

ggplot(res) + aes(bveff, color=factor(back_vac_mech)) + facet_grid(. ~ exp_sd, labeller = label_bquote(cols=R[0]%~~%.(exp_sd/2-.5))) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=factor(back_vac_mech), color=NULL), alpha = 0.3) +
  geom_ribbon(aes(ymin=lo.q, ymax=hi.q, fill=factor(back_vac_mech), color=NULL), alpha = 0.3) +
  geom_line(aes(y=mn)) +
  scale_x_continuous("Vaccine Efficacy") +
  scale_y_continuous("Number of Cases per Ring") +
  scale_color_manual("Vaccine\nMechanism",
    values=c(`0`="black",`1`="dodgerblue"), aesthetics = c("color","fill")
  ) + theme_minimal()
