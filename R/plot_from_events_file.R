# process events file
these.serial <- c(0:99)
#this.replicate <- 
some.events <- events[serial %in% these.serial] #replicate==this.replicate & 

# plot a separated replicate
plot.ts <- function(this.ts, this.replicate) {
  plot(this.ts, xlim=c(0, max(as.numeric(names(this.ts)))), axes=F)
  axis(1)
  axis(2)
  legend("topleft", legend=this.replicate, bty="n", lwd=1, col="white")
}

# make plots for each replicate
for(this.replicate in 0:99) {
  one.epi <- some.events[replicate==this.replicate]
  e.ts <- one.epi$e_time[one.epi$e_time != -1] %>% floor %>% table 
  i.ts <- one.epi$i_time[one.epi$i_time != -1] %>% floor %>% table
  d.ts <- one.epi$d_time[one.epi$d_time != -1] %>% floor %>% table
  plot.ts(i.ts, this.replicate)
  }

# # make a table of number infected in each simulation
# num.infected <- data.table(ID=0:99, infected=NA) #num infected in ring in total
# num.infected.index <- data.table(ID=0:99, infected=NA) # num infected by index case
# num.inf.by.L1 <- data.table(ID=NA, L1=NA, infected=NA) #num infected by L1s
# for(this.replicate in 0:99) {
#   one.epi <- some.events[replicate==this.replicate]
#   # infections in the ring in total
#   num.infected$infected[this.replicate-1] <- one.epi[, list(sum(one.epi$inf_by > -1))] %>% as.numeric
#   # infections caused by indexes
#   num.infected.index$infected[this.replicate-1] <- one.epi[, list(sum(one.epi$inf_by == 0))] %>% as.numeric
#   
#   # infections caused by L1s
#   L1s <- one.epi$node_id[one.epi$level==1]
#   L1s.inf <- one.epi$node_id[one.epi$level==1 & one.epi$inf_by > -1]
#   num.L1s.inf <- length(L1s.inf)
#   if(num.L1s.inf > 1) {
#     add.values <- c(ID=this.replicate, L1s.inf=one.epi[, list(sum(one.epi$inf_by %in% L1s))] %>% as.numeric
#     
#     num.inf.by.L1 <- rbind(num.inf.by.L1, ) 
#     
#   } else {
#   }
#   }

# histograms of offspring distributions 
par(las=1)
hist(num.infected$infected, col="grey65", border=NA, 
     breaks=seq(0, max(num.infected$infected), 1), main="total infected in ring")
hist(num.infected.index$infected, col="grey45", border=NA, 
     breaks=seq(0, max(num.infected$infected), 1), main="Offspring dist index case")
abline(v=mean(num.infected.index$infected), col="firebrick")
legend("topright", legend=paste0("mean: ", mean(num.infected.index$infected)), bty="n", col="white")
hist(num.inf.by.L1$infected, col="grey45", border=NA, 
     breaks=seq(0, max(num.infected$infected), 1), main="Offspring dist from L1s")
abline(v=mean(num.inf.by.L1$infected), col="firebrick")
legend("topright", legend=paste0("mean: ", mean(num.inf.by.L1$infected)), bty="n", col="white")
