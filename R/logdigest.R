suppressPackageStartupMessages({
  require(data.table)
})

.args <- c("../tests", "digestlogs.rds")
.args <- commandArgs(trailingOnly = T)

logfiles <- list.files(.args[1], pattern = "\\.log$", full.names = T)

res <- rbindlist(lapply(logfiles, function(fn) {
  dt <- fread(fn)
  colnames(dt) <- c("serial","node_id", "L", "has_background", "trace_time", "rv_time", "infection_time")
  dt
}))

saveRDS(res, tail(.args, 1))
