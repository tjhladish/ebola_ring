require(data.table)
require(statnet)
require(RColorBrewer)
require(colorRamps)
require(network)
require(ggplot2)
require(igraph)

# [heavily] modified from a script by Mario Jendrossek


# set location
if(Sys.info()[["user"]] == "Rosalind") {
  network.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/networks/" # these ones had simulations
  network.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/posterior_network_csvs/" # these ones were the networks only
} else {
  # add your own location
  output.location <- "~/TomOrCarlsFolders"
}

# which networks to plot
networks <- list.files(network.location)
networks <- gsub(".csv", "", networks)
networks <- networks[grep("_1", networks)] # needed when comparing clustered or unclustered networks
networks <- gsub("_1", "", networks)
networks <- networks[order(as.numeric(networks))]
keeps  <- seq(1, 200, 2) # needed when the values increase in mulitples because of simulations
networks <- networks[keeps]

# select a network
k <- 1
par(mfcol=c(5,5))
events$set <- floor(events$serial/1000)*1000 # needed when there are sets of networks

# plot the networks
for(k in 71:95) {
  net.num <- networks[k]
  print(net.num)
  this.net <- read.csv(paste0(network.location, net.num, "_1.csv")) %>% data.table
  # match to add level to it
  node.chars <- events[events$serial==net.num, list(node_id, level)]
  node.chars <- cbind(node_name=node.chars$node_name, node_id=node.chars$node_id, level=node.chars$level)
  
  net <- graph.data.frame(this.net, node.chars, directed=F)
  
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  #plot(net, edge.arrow.size=.4, vertex.label=NA, 
  #     vertex.color=new.pal[as.numeric(V(net)$level)+1], vertex.frame.color=new.pal[as.numeric(V(net)$level)+1],
  #     vertex.size=5, layout=layout.fruchterman.reingold)
  V(net)$level <- as.numeric(V(net)$level)
  
  L2s <- V(net)[V(net)$level=="2"]
  
  small.net <- delete_vertices(net, L2s)
  plot(small.net, edge.arrow.size=.4, vertex.label=NA, 
       vertex.color=new.pal[as.numeric(V(small.net)$level)+1], vertex.frame.color=new.pal[as.numeric(V(small.net)$level)+1],
       vertex.size=5, layout=layout.fruchterman.reingold) 
}

this.net <- read.csv(paste0(network.location, 145000, "_1.csv")) %>% data.table
this.net2 <- read.csv(paste0(network.location, 147000, "_1.csv")) %>% data.table


