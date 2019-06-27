
network_properties.rds: analyze_nets.R ../tests/networks
	Rscript $^ $@

digested_%.rds: digest.R ../tests/test_abc_example_%.sqlite ../tests/abc_ebola_sim_%.json
	Rscript $^ $@