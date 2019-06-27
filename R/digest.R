pkgs <- lapply(c("data.table","RSQLite"), require, character.only=T)

.args <- c("../tests/test_abc_example_nobias.sqlite", "../tests/abc_ebola_sim_nobias.json", "digest_nobias.rds")
.args <- commandArgs(trailingOnly = T)

sqldb <- .args[1]
abcref <- jsonlite::read_json(.args[2])
tarfile <- .args[3]

extractcolname <- function(el) if(is.null(el$short_name)) el$name else el$short_name

parcols <- c("serial", sapply(abcref$parameters, extractcolname))
metcols <- sapply(abcref$metrics, extractcolname)

selcols <- c(parcols, metcols)

drv <- dbDriver("SQLite")
db <- dbConnect(drv, .args[1], flags=SQLITE_RO)
pars <- data.table(
  dbGetQuery(db,
    sprintf("SELECT %s FROM par JOIN met USING(serial) JOIN job USING(SERIAL) WHERE status=='D';", paste0(selcols, collapse = ", "))
  )
)
dbDisconnect(db)

pars[, unvax_correct_prob := 1 ]
pars[, coverage := realized_coverage ]
pars[back_vac_mech == 1, coverage := realized_coverage/bveff ]
pars[back_vac_mech == 1, unvax_correct_prob := (1-coverage)/(1-realized_coverage) ]

pars$count_pre_vax <- NULL

id.vars <- c("serial", "net_rep", "epi_rep", "bveff", "trace_prob", "exp_sd", "vaccine_delay", "back_vac_mech",
             "use_bias", "realized_coverage", "unvax_correct_prob", "coverage")

mlt <- melt.data.table(pars, id.vars = id.vars)

mlt[, group   := factor(gsub("^(vaccine|unvax)_.*$", "\\1", variable)) ]
mlt[, outcome := factor(gsub("^.*_(pos|neg)_.*$", "\\1", variable)) ]
mlt[, window  := as.integer(gsub("^.*_(\\d+)$", "\\1", variable)) ]

mlt$variable <- NULL

# fm <- as.formula(
#   sprintf("%s ~ %s",paste(c(id.vars,"coverage"), collapse=" + "),paste(c("background", "outcome"), collapse=" + "))
# )
# 
# results <- dcast(
#   rbind(testpos, testneg),
#   fm,
#   value.var = "N", fill = 0
# )

saveRDS(mlt, tarfile)
