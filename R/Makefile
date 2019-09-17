
network_properties.rds: analyze_nets.R ../tests/networks
	Rscript $^ $@

digested_%.rds: digest.R ../tests/test_abc_example_%.sqlite ../tests/abc_ebola_sim_%.json
	Rscript $^ $@

digested_stratify.rds: multi_vs_single_index.R ../tests/abc_ebola_sim_nonleaky.json ../tests/test_abc_example_nonleaky.sqlite
	Rscript $^ $@

digested_stratify_nobias.rds: coverage_varying.R ../tests/abc_ebola_sim_nobias.json ../tests/test_abc_example_nobias.sqlite
	Rscript $^ $@


SAMPS := $(addprefix tertiles/,$(addsuffix .rds,$(shell seq -w 00 40)))

SIMPSAMPS := $(addprefix simple/,$(addsuffix .rds,$(shell seq -w 00 80)))

STRATSAMPS := $(addprefix stratified/,$(addsuffix .rds,$(shell seq -w 00 80)))

COVSAMPS := $(addprefix coverage_lvls/,$(addsuffix .rds,$(shell seq -w 00 80)))

$(SAMPS): analyze_with_epi_metrics.R digested_nobias.rds
	Rscript $^ $@

$(SIMPSAMPS): analyze_simple.R digested_nobias.rds
	Rscript $^ $@

$(STRATSAMPS): analyze_simple.R digested_stratify.rds | stratified
	Rscript $^ $@

$(COVSAMPS): analyze_covlvls.R digested_stratify_nobias.rds | coverage_lvls
	Rscript $^ $@

stratified coverage_lvls:
	mkdir -p $@

samples: $(SAMPS)

simpsamples: $(SIMPSAMPS)

stratsamples: $(STRATSAMPS)

covsamples: $(COVSAMPS)

digestlogs.rds: logdigest.R ../tests
	Rscript $^ $@

km_plot.rds: iddconf_km_plots.R digestlogs.rds ../tests/abc_ebola_km.sqlite
	Rscript $^ $@

net_comp.rds: iddconf_net_plots.R $(patsubst %,../tests/%/4042_lvl4.csv,networks declustered reclustered)
	Rscript $^ $@


iddconf_poster_big.png: iddconf_plot_merge.R net_comp.rds km_plot.rds
	Rscript $^ $@

iddconf_tdistros.png: iddconf_ts_plots.R timedistros.json
	Rscript $^ $@

iddconfigs: iddconf_poster_big.png iddconf_tdistros.png