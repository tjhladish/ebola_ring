png('whole_pop.png', width=1500, height=1500, res=200)
par(mar=c(3,3,1,1))
pop = read.table('coords', header=T, sep=',')
net = read.table('ebola_ring.csv', header=F, sep=',')

plot(pop$x, pop$y, pch=20, col='#bbbbbb', xlab='', ylab='')

names(pop) = c('V1', 'x1', 'y1')
segs = merge(net, pop)
names(pop) = c('V2', 'x2', 'y2')
segs = merge(segs, pop)[c('x1', 'y1', 'x2', 'y2')]
segments(segs$x1, segs$y1, segs$x2, segs$y2, col='#555555')

all_coords = rbind(as.matrix(segs[1:2]), as.matrix(segs[3:4]))
all_coords = unique(all_coords)
points(all_coords, pch=20, col='#0000bb')
dev.off()

png('ring_only.png', width=1500, height=1500, res=200)
par(mar=c(3,3,1,1))
plot(all_coords, xlab='', ylab='', type='n')
segments(segs$x1, segs$y1, segs$x2, segs$y2, col='#555555bb', lwd=2)
points(all_coords, pch=20, col='#0000bbbb', cex=3)
dev.off()
