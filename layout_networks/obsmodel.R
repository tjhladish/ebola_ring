## observation model

require(data.table)

args <- c("~/Downloads/sample_networks","obsplot-85.png")
args <- commandArgs(trailingOnly = TRUE)

set.seed(0L)

pth <- args[1]
tar <- args[2]
whichid <- as.integer(gsub(".+-(\\d+)\\..+","\\1",tar))
sides <- as.integer(unlist(strsplit(gsub("^.+[^\\d]+(\\d+)x(\\d+)\\..+$", "\\1 \\2", tar)," ")))

layouts = list.files(paste0(pth,"/layouts"), pattern = "\\.layout$", full.names = T)
networks = list.files(paste0(pth,"/posterior_networks"), pattern = "\\.csv$", full.names = T)
labels = fread(list.files(pth, pattern = "\\.out", full.names = T, include.dirs = F)[1], col.names = c("id","vertex","L"))

net.dt <- rbindlist(lapply(networks, function(fn) {
  id <- as.integer(gsub("^.+/([^/]+)_.+\\.csv$","\\1",fn))
  res <- fread(fn, col.names = c("from","to"))[, id := id ]
  res
}))

layout.dt <- rbindlist(lapply(layouts, function(fn) {
  id <- as.integer(gsub("^.+/([^/]+)_.+\\.layout$","\\1",fn))
  res <- fread(fn, col.names = c("vertex","x","y"))[, id := id ]
  res
}))[labels, on=c("id", "vertex")]

## get the networks by layer
ref <- net.dt[layout.dt, .(from, to, from.L=L, id), on=.(id, from=vertex)][!is.na(to)][
  layout.dt, .(from, to, from.L, to.L=L, id), on=.(id, to=vertex)
][!is.na(from)][id == whichid]
## would like to assert about labeling: all L1 labels < all L2 labels, however:
# ref[,max(from),keyby=.(id, L)]

## if we assume all nodes surveyed, according to those probabilities, then the following:

# ref[, forward_outcome := factor(ifelse(
#   runif(.N) < name_p, ifelse(runif(.N) < misname_p,"MISSED","ID"),
#   "NONID"
# ), levels = c("MISSED","NONID","ID"), ordered = T) ]
# 
# ref[, reverse_outcome := factor(ifelse(
#   runif(.N) < name_p, ifelse(runif(.N) < misname_p,"MISSED","ID"),
#   "NONID"
# ), levels = c("MISSED","NONID","ID"), ordered = T) ]

sampler <- function(reportp, namep, ref) {
  subref <- copy(ref[from.L < 2 & to.L < 2])[, c("report","name") := .(reportp, namep) ]
  
  subref[, forward_outcome := factor(ifelse(
    runif(.N) < reportp, ifelse(runif(.N) < namep, "ID", "MISSED"),
    "NONID"
  ), levels = c("MISSED","NONID","ID"), ordered = T) ]
  
  subref[, reverse_outcome := factor(ifelse(
    runif(.N) < reportp, ifelse(runif(.N) < namep, "ID", "MISSED"),
    "NONID"
  ), levels = c("MISSED","NONID","ID"), ordered = T) ]
  
  # assert: patient zero is always id == 0
  new_finds <- subref[from == 0 & forward_outcome == "ID", .(tar=to, src=from), by=id]
  found_edges <- copy(new_finds)
  
  while(new_finds[,.N]) {
    # found_nodes <- rbind(new_finds, found_nodes)
    new_edges <- rbind(
      subref[new_finds, on=.(id, from=tar)][forward_outcome == "ID", .(tar=to, src=from), by=id],
      subref[new_finds, on=.(id, to=tar)][reverse_outcome == "ID", .(tar=from, src=to), by=id]
    )
    new_finds <- new_edges[!found_edges, on=.(id, tar)][,.(src=src[1]),by=.(id,tar)]
    found_edges <- rbind(new_edges, found_edges)
  }
  
  subref[, found := FALSE ]
  subref[found_edges, found := TRUE, on=.(id, from==src, to==tar)]
  subref[found_edges, found := TRUE, on=.(id, from==tar, to==src)]
  subref
}

grd <- expand.grid(reportp = c(seq(.1,.9,by=.2),1), namep = c(seq(.1,.9,by=.2),1))

dt <- rbindlist(
  mapply(
    sampler, reportp = grd$reportp, namep = grd$namep, MoreArgs = list(ref=ref), SIMPLIFY = F
  )
)[,.(from, to, from.L, to.L, found),by=.(id, report, name)]

require(ggplot2)
require(geomnet)

aug <- dt[,.(missing_from=setdiff(to, from)),by=.(id, report, name)][
  dt, .(id, report, name, from=missing_from, to=NA, from.L=to.L, to.L=NA, found=found),
  on=.(id, report, name, missing_from==to)
]

plot.dt <- rbind(dt, aug)

colors = c(`0`='firebrick',`1`='dodgerblue',`2`='grey65')

ggsave(tar,
ggplot(
  plot.dt[found==TRUE],
  aes(
    from_id=as.character(from),
    to_id=as.character(to),
    color=as.factor(from.L),
    linetype = found
  )
) +
  facet_grid(report~name, scales = "free") +
  geom_net(
    # layout.alg="fruchtermanreingold",
    directed = FALSE, labelon = FALSE, size = 1,
    ecolour = "lightgray", ealpha = 0.5,
    # linewidth = 0.2,
    fiteach = TRUE
  ) + theme_net() +
  theme(strip.background = element_blank()) +
  scale_color_manual(guide="none", values = colors) +
  scale_linetype_manual(values = c(`TRUE`="solid", `FALSE`="dotted")),
  width = 6, height = 6, units = "in"
)
# ggplot(
#   plot.dt[id == 19156 & found == TRUE],
#   aes(
#     from_id=as.character(from),
#     to_id=as.character(to),
#     color=as.factor(from.L),
#     linetype = found
#   )
# ) +
#   facet_grid(report~name, scales = "free") +
#   geom_net(
#     # layout.alg="fruchtermanreingold",
#     directed = FALSE, labelon = FALSE, size = 1,
#     ecolour = "lightgray", ealpha = 0.5,
#     # linewidth = 0.2,
#     fiteach = TRUE
#   ) + theme_net() +
#   theme(strip.background = element_blank()) +
#   scale_color_manual(guide="none", values = colors) +
#   scale_linetype_manual(values = c(`TRUE`="solid", `FALSE`="dotted"))