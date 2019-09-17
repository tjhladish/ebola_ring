suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
  require(jsonlite)
})

.args <- c("../tests/abc_ebola_sim_nonleaky.json", "../tests/test_abc_example_nonleaky.sqlite", "digested_stratify.rds")
.args <- commandArgs(trailingOnly = T)

abcref <- read_json(.args[1])

extractcolname <- function(el) if(is.null(el$short_name)) el$name else el$short_name

parcols <- c("serial", sapply(abcref$parameters, extractcolname))
metcols <- sapply(abcref$metrics, extractcolname)

selcols <- c(parcols, metcols)

drv <- dbDriver("SQLite")
db <- dbConnect(drv, .args[2], flags=SQLITE_RO)
pars <- data.table(
  dbGetQuery(
    db,
    sprintf("SELECT %s FROM par JOIN met USING(serial) JOIN job USING(SERIAL) WHERE status=='D';", paste0(selcols, collapse = ", "))
  )
)[count_pre_vax != 0]
dbDisconnect(db)

pars[, index := c("multi","one")[(count_pre_vax==1)+1] ]

pars$count_pre_vax <- NULL

id.vars <- c(parcols, "realized_coverage", "index")

mlt <- melt.data.table(pars, id.vars = id.vars, value.name = "count")

mlt[, group   := factor(gsub("^(vaccine|unvax)_.*$", "\\1", variable)) ]
mlt[, outcome := factor(gsub("^.*_(pos|neg)_.*$", "\\1", variable)) ]
mlt[, window  := as.integer(gsub("^.*_(\\d+)$", "\\1", variable)) ]

mlt$variable <- NULL

recast <- mlt[, dcast.data.table(.SD, ... ~ group + outcome, value.var = "count")][unvax_pos + vaccine_pos != 0]
#' key(recast)
#' require(ggplot2)
#' ggplot(pars) + aes(x=coverage, fill=factor(back_vac_mech)) + facet_grid(
#'   (count_pre_vax == 0) + back_vac_mech ~ bveff, scales = "free_y"
#' ) + geom_histogram()

saveRDS(recast, tail(.args, 1))
