require("RSQLite")
require(beanplot)
drv = dbDriver("SQLite")
db = dbConnect(drv, "ebola_rings_popfit-hh_cliques.sqlite")
#db = dbConnect(drv, "ebola_rings_popfit-test.sqlite")
abc = dbGetQuery(db, 'select J.*, P.*, M.* from job J, par P, met M where J.serial = P.serial and J.serial = M.serial')
dbDisconnect(db)

extra_serials = grepl('serial', names(abc))
extra_serials[1] = FALSE
abc = abc[,!extra_serials]
abc = subset(abc, select=-c(startTime, duration, attempts, seed))
abc$post_bool = abc$posterior >= 0

par_cols = 6:9
met_cols = c(10:29)
obs_met = c(8, 14, 23, 17.37, 47.5, 64.5, 81.25, 72.966667, 2.652256, 4.840336, 9.087912, 10.03201, 0.02122255, 0.1005773, 0.2675585, 0.0110294, 0.12861, 0.23708025, 0, 0, 0);


#obs_met = c(8.75, 20, 25, 18.75, 40, 85.75, 0, 0, 2, 0.27403825, 0.373397, 0.65661775, 0.02122255, 0.1005773, 0.2675585, 0.0110294, 0.12861, 0.23708025, 0, 0, 0, -0.2950645, 0.261722, 0.81093, 0.599706, 0.692489, 0.74585775)


incomplete_sets = unique(abc$smcSet[abc$status!='D']) # 6
all_sets = unique(abc$smcSet)
complete_sets = setdiff(all_sets, incomplete_sets)
last_complete_set = max(complete_sets)

pdf('marginal_pars_v5.pdf', width=8.5, height=11)
par(mfrow=c(4,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
for (col in par_cols) {
    colname = names(abc)[col];
    cat(paste(colname, '\n'))
    beanplot( abc[,colname] ~ abc$post_bool*abc$smcSet, what=c(1,1,1,0), col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1,1, grey(0.7))), main='', ylab=colname, side='both', bw='nrd0')
}
dev.off()

#ylims = list(c(0,800), c(0,5), c(0,200), c(0,350), c(0,600), c(10,30000), c(10e-1, 10e4), c(-2, 6), c(0, 1), c(0, 1), c(0, 1), c(-6, 4), c(-0.1, 0.1))

# Plot metrics - tends to be trickier, because distributions can be weird (e.g, highly skewed or long-tailed)
pdf('marginal_mets_v5.pdf', width=5, height=20)
par(mfrow=c(17,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
#alpha = 0.025 # plot middle 95% of distributions
alpha = 0.0 # plot middle 95% of distributions
for (col in met_cols) {
    met_idx = col - max(par_cols)
    colname = names(abc)[col]

    # calculate reasonable plot limits
    obs_val = obs_met[met_idx]
    val_lims = quantile(abc[,colname], na.rm=T, probs=c(alpha, 1-alpha))
    val_min = min(val_lims[1], obs_val) 
    val_max = max(val_lims[2], obs_val) 

    # filter out NAs
    complete = complete.cases(abc[,colname]) & abc[,colname] > val_lims[1] & abc[,colname] < val_lims[2]
    cat(paste0(colname, ' ', val_lims[1], ', ', val_lims[2], '\n'))

    # plot the stuff
    beanplot( abc[complete, colname] ~ abc$post_bool[complete] * abc$smcSet[complete], what=c(0,1,1,0), 
              col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1, 1, grey(0.7))), 
              #main='', ylab=colname, side='both', ylim=unlist(ylims[met_idx]) )
              main='', ylab=colname, side='both', ylim=c(val_min, val_max), na.rm=T, log='' , bw='nrd0')
    abline(h=mean(abc[abc$smcSet==last_complete_set, colname], na.rm=T), lty=3)
    abline(h=obs_val, col=2, lwd=1)
}

dev.off()
