# compute assorted statistics for all the contact networks

suppressPackageStartupMessages({
  require(igraph)
  require(data.table)
})

.args <- c("../tests/networks", "network_properties.rds")
.args <- commandArgs(trailingOnly = T)

# list all the network files
network_files <- list.files(path=.args[1], pattern = "*.csv", full.names = T)

# for each network file, do some analysis, then bind result up into a single data.table
results <- rbindlist(lapply(network_files, function(file_name) {
  dt <- fread(file_name)
  # defaults for singleton networks
  trns <- 0
  if (any(grepl("V2", names(dt)))) { # for non-singleton networks...
    # make an igraph object
    ig <- graph_from_data_frame(fread(file_name)[, .(A=V1+1, B=V2+1) ], directed = F)
    
    # n.b., this is in the order defined by V(ig), so no concerns about the vertex ids
    # being "out of order" from distances
    V(ig)$L <- distances(ig, to = 1)
    trns <- transitivity(induced_subgraph(ig, V(ig)[L <= 2]))
    
    # TODO: other measurements...
  } else {
    warning(sprintf("%s is singleton.", file_name))
  }
  
  data.table(fn = file_name, transitivity = trns)
}))

results[, netid := as.integer(gsub("^.+/4(\\d{3})_.+$","\\1", fn)) ]
results$fn <- NULL

results[, partition := {
  brks <- quantile(transitivity, probs = (1:2)/3)
  cut(transitivity, breaks=c(0, brks, 1), labels=c("lo","md","hi"), right=FALSE, include.lowest = T)
}]

# results[,.N,by=partition] # confirmed divided into roughly equal partitions

saveRDS(results, tail(.args, 1))