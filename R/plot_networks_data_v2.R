require(data.table)
require(statnet)
require(RColorBrewer)
require(colorRamps)
require(network)
require(ggplot2)
require(igraph)


# set location
if(Sys.info()[["user"]] == "Rosalind") {
  network.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/networks/"
  network.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/posterior_network_csvs/"
} else {
  # add your own location
  output.location <- "~/TomOrCarlsFolders"
}

# which networks to plot
networks <- list.files(network.location)
networks <- gsub(".csv", "", networks)
networks <- networks[order(as.numeric(networks))]

# select a network
k <- 1
par(mfcol=c(5,5))
levels.file <- read.table(paste0(network.location, "../node_levels.out"), header=F, sep=" ") %>% data.table
colnames(levels.file) <- c("network", "node", "level")
levels.file[network==457,]

for(k in 1:25) {
  net.num <- networks[k]
  print(net.num)
  this.net <- read.csv(paste0(network.location, net.num, ".csv")) %>% data.table
  # match to add level to it
  node.chars <- levels.file[levels.file$network==net.num, list(node, level)]
  node.chars <- cbind(node_name=node.chars$node, level=node.chars$level)
  
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
       vertex.size=8, layout=layout.fruchterman.reingold) 
}

this.net <- read.csv(paste0(network.location, 145000, "_1.csv")) %>% data.table
this.net2 <- read.csv(paste0(network.location, 147000, "_1.csv")) %>% data.table

