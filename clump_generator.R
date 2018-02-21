png("pop_dist.png", width=2000, height=2000, res=200)

grid_size = 5
hhs = 10**seq(log10(10),log10(500),length.out = grid_size)
ss  = 10**seq(log10(0.01), log10(0.25), length.out = grid_size)
par(mfrow=c(grid_size,grid_size))
par(oma=c(4,4,0,0))
par(mar=c(1,1,0,0))

for (r in 1:grid_size) {
    for (c in 1:grid_size) {
        N = 500
        hh = hhs[r]
        hh_size = N/hh
        s = ss[c]
        
        x = runif(hh); y = runif(hh)
        
        y2 = rnorm(hh_size*hh, y, s)
        x2 = rnorm(hh_size*hh, x, s)

        plot(x, y, pch=20, cex=0.75, xlab='', ylab='', axes=F, col=1)
        box()
        points(x2, y2, pch=20, cex=0.75)
    }
}

mtext(signif(ss,3),1,0,outer=T,at=seq(0.1,0.9,0.2))
mtext(signif(hhs,3),2,0,outer=T,at=seq(0.9,0.1,-0.2))
mtext("SD",1,outer=T,at=0.5,line=2)
mtext("# households",2,outer=T,at=0.5,line=2)
dev.off()