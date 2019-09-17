suppressPackageStartupMessages({
  require(data.table)
  require(igraph)
  require(ggplot2)
  require(cowplot)
})

readepifirenet <- function(pth) {
  dtorig <- fread(pth)
  # igraph convention indexes from 1
  dtorig[, c("V1","V2") := .(V1+1L, V2+1L) ]
  return(dtorig)
}

extractL2 <- function(dt) {
  ig <- graph_from_data_frame(dt)
  V(ig)$L <- distances(ig, 1, to = V(ig))
  return(induced_subgraph(ig, V(ig)[L<3]))
}

mkdtsfromig <- function(ig, name, amp=5) {
  E(ig)$weight <- 1
  L0 <- V(ig)[L==0]
  L1 <- V(ig)[L==1]
  E(ig)[L0 %--% L1]$weight <- amp
  E(ig)[L1 %--% L1]$weight <- 2*amp
  L2 <- V(ig)[L==2]
  E(ig)[L2 %--% L2]$weight <- 1.5*amp
  
  pos.dt <- data.table(layout_with_fr(ig, niter=1000, start.temp = vcount(ig)), v=1:gorder(ig), L=V(ig)$L, name = name)
  pos.dt[, V1 := V1 - mean(V1) ]
  pos.dt[, V2 := V2 - mean(V2) ]
  maxr <- pos.dt[, max(sqrt(V1^2+V2^2))]
  pos.dt[, V1 := V1/maxr ]
  pos.dt[, V2 := V2/maxr ]
  el.dt <- data.table(as_edgelist(ig, names = FALSE), eid = 1:ecount(ig), name = name)
  el.dt[pos.dt, c("x","y","fromL") := .(i.V1, i.V2, L), on=.(V1=v)]
  el.dt[pos.dt, c("xend","yend","toL") := .(i.V1, i.V2, L), on=.(V2=v)]
  list(pos.dt, el.dt)
}

plotter <- function(igref, igdec, igrec) {
  processed <- mapply(mkdtsfromig, ig=list(igref,igdec,igrec), name=c("ref","dec","rec"), SIMPLIFY = F)
  reduced <- Reduce(function(l, r) {
    list(rbind(l[[1]], r[[1]]), rbind(l[[2]], r[[2]]))
  }, processed)
  pos.dt <- reduced[[1]]
  el.dt <- reduced[[2]]
  
  ggplot(pos.dt) + aes(color=factor(L), size=factor(L)) +
    facet_wrap(
      ~name, nrow =3,
      labeller = labeller(
        name=c(ref="Realistic Clustering", rec="Realistic Degree Only", dec="Realistic Population Size Only")
      )
    ) +
    geom_segment(aes(x, y, xend=xend, yend=yend, color=factor(fromL), size=NULL), el.dt, alpha=0.4, size=0.5) +
    geom_point(aes(V1, V2), data=function(dt) dt[order(-L)]) +
    # geom_text(aes(0,1,
    #   label=unique(c(ref="Realistic Clustering", rec="Realistic Degree", dec="Realistic Number")[name]),
    #   size=NULL, color=NULL
    # )) +
    scale_color_manual(values=c(`0`="firebrick", `1`="dodgerblue", `2`="darkgrey")) +
    scale_size_manual(values=c(`0`=4, `1`=3, `2`=2)) +
    guides(color="none", size="none") +
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) +
    theme(
      rect = element_rect(fill="transparent"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face="bold", size = rel(2))
    )
}

.args <- c(
  sprintf("../tests/%s/4%03i_lvl4.csv",c("networks","declustered","reclustered"),42), # 1, 42
  "net_comp.rda"
)
.args <- commandArgs(trailingOnly = T)

igs <- lapply(head(.args, 3), function(fn) extractL2(readepifirenet(fn)))

origp <- plotter(igs[[1]], igs[[2]], igs[[3]])

saveRDS(origp, tail(.args, 1))

# save_plot(tail(.args,1), origp, base_width=5, nrow=3, base_asp = 0.9, bg="transparent")
