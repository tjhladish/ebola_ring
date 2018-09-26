require(data.table)
require(ggplot2)

args <- c("~/Downloads/sample_networks","adsfasdfasdf4x4.png")
args <- commandArgs(trailingOnly = TRUE)

set.seed(0L)

pth <- args[1]
tar <- args[2]
sides <- as.integer(gsub("^.+(\\d+)\\..+$", "\\1", tar))

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

positions.dt <- {
  tmp <- net.dt[layout.dt,
    .(id, from, to, from.x=x, from.y=y), on=.(from=vertex, id=id)
  ][!is.na(to)]
  tmp[layout.dt,
    .(id, from, to, from.x, from.y, to.x=x, to.y=y), on=.(to=vertex, id=id)
  ][!is.na(from), .SD, keyby=.(id, from, to)]
}

picks <- sample(unique(net.dt$id), sides^2)

colors = c(`0`='firebrick',`1`='dodgerblue',`2`='grey65')
vertices = layout.dt[id %in% picks]
edges = positions.dt[id %in% picks]
    
# ggplot(vertices) +
#   facet_wrap(~id, scales = "free") +
#   geom_segment(
#     aes(x=from.x, y=from.y, xend=to.x, yend=to.y), data=edges,
#     color="lightgrey", alpha=0.5
#   ) +
#   aes(x=x, y=y, color=factor(L)) +
#   geom_point() +
#   theme_minimal() + theme(
#     axis.text = element_blank(), axis.title = element_blank(),
#     panel.grid = element_blank()
#   ) + scale_color_manual(guide="none", values = colors)

require(geomnet)

ggsave(tar,
  ggplot(
    edges[vertices, .(from, to, id, L), on=.(id, from=vertex)],
    aes(from_id=as.character(from), to_id=as.character(to), color=as.factor(L))
  ) +
    facet_wrap(~id, scales = "free") +
    geom_net(
      layout.alg="fruchtermanreingold",
      directed = FALSE, labelon = FALSE, size = 1,
      ecolour = "lightgray", ealpha = 0.5,
      linewidth = 0.2, fiteach = TRUE
    ) + theme_net() +
    theme(strip.background = element_blank()) +
    scale_color_manual(guide="none", values = colors),
  width = 6, height = 6, units = "in"
)