## read in files
library(data.table)
library(magrittr)

# set location
if(Sys.info()[["user"]] == "Rosalind") {
  output.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/"
  output.location <- "~/Sync/LSHTM/Collaboration/Ring Simulator/Outputs/v2/" # we need a system for each iteration of outputs
} else {
  # add your own location
  output.location <- "~/TomOrCarlsFolders"
}

# this is the column labels for the event.time, which doesnt have a header column 
event.time.labels <- read.table(paste0(output.location, "event_time.columns"))[1,] %>% unlist %>% unname %>% as.character
# serial: is how you match output to the database (0-399,999)
# and (using arithmetic) match up with the network files
# unique serial for each simulation
# replicate: for each set of parametesr we did 1000 replicates, with different seeds
# node_id: sequential integers, starting at 0. node_id 0 == L0 (patient zero)
# level: L0, L1, L2
# inf_by: who the node was infected by
# NA = -1

# event_time.log == big file
# paired with event_time.columns
# serials are used for matching
events <- fread(paste0(output.location, "event_times-75-90.out")) %>% data.table
colnames(events) <- event.time.labels

# network files
# the _0 means there is no clustering
# the _1 means there is no clustering, but they are paired
# 0 goes with 1000
# it goes up by 1000 because there are 1000 sims per network

# get list of network files 
network.files <- list.files(paste0(output.location, "networks/")) #%>% gsub(".csv", "", .)

# read in an example net
this.net <- read.csv(paste0(output.location, "networks/", network.files[1]))
