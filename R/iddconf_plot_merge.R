suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(cowplot)
})


.args <- c("net_comp.rda","km_plot.rda","poster.png")
.args <- commandArgs(trailingOnly = T)

pleft <- readRDS(.args[1])
pright <- readRDS(.args[2])

mergep <- plot_grid(pleft, pright, align = "h")

save_plot(tail(.args,1), mergep, base_width=10, nrow=3, base_asp = 0.8, bg="transparent")
