
layouts = Sys.glob('./*gauss.layout')
networks = Sys.glob('../posterior_networks/*gauss.csv')

levels = read.table('gauss_levels.out')
colors = c('firebrick','dodgerblue','grey65')
par(mfrow=c(3,5))
par(mar=c(0,0,0,0))
for (i in 1:15) {
#for (i in 1:length(networks)) {
#i = 14
    d = read.table(layouts[i])
    n = read.table(networks[i])
    net_id = regmatches(networks[i], regexpr('\\d+', networks[i]))    
    names(d) = c('id1', 'x1', 'y1')
    d = d[order(d$id1),]
    plot(d$x1, d$y1, pch=20, cex=3, col=colors[levels$V3[levels$V1==net_id]+1])
    names(n) = c('id1', 'id2')
    layout = merge(n, d)
    names(d) = c('id2', 'x2', 'y2')
    layout = merge(layout, d)
    segments(layout$x1, layout$y1, layout$x2, layout$y2, lwd=2, col='#00000044')
}
