require(data.table)

.args <- c("../tests/test_result.csv")
.args <- commandArgs(trailingOnly = T)

dt <- fread(.args[1])

# so we can appropriately baseline time...
tracetm <- dt[,min(trace, na.rm = T)]
rvtm <- dt[, min(rv, na.rm = T)]

dt[, .N, by=.(background, onset=round(onset-rvtm)) ]
