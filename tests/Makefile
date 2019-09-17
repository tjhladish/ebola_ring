-include ../local.mk

FLAGS ?= -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic -fmax-errors=5
#
# WORKSPACE ?= $(HOME)/work
#
ABC_PATH = $(WORKSPACE)/AbcSmc
# # override these to use system copies, rather than AbcSmc copies
# GSL_PATH ?= $(ABC_PATH)/gsl_local

SQL_PATH ?= $(ABC_PATH)/sqdb
EPI_PATH = $(WORKSPACE)/EpiFire/src
ABC_LIB = -L$(ABC_PATH) -labc -ljsoncpp -lsqdb -Wl,--no-as-needed -ldl $(ABC_PATH)/sqlite3.o
GSL_LIB = -lm -lgsl -lgslcblas -pthread

#INCLUDE = -I$(ABC_PATH) -I$(GSL_PATH)/include/ -I$(EPI_PATH)
INCLUDE = -I$(EPI_PATH)
LINK = -L$(EPI_PATH)
#
# LDFLAGS=  $(EPI_PATH)/*.o

$(EPI_PATH)/libsim.a: $(EPI_PATH)
	$(MAKE) -C $^

epifire: $(EPI_PATH)/libsim.a

event_template.o core_event_queue.o: ../NetworkSimplate.h

id_community.o: ../IDCommunity.h

test_es.o: ../IDCommunity.h ../EbolaSim.h ../NetworkSimplate.h

%.o: %.cpp | epifire
	g++ $(FLAGS) $< -o $@ $(INCLUDE) $(LINK) -lsim $(GSL_LIB)

id_community.dot: id_community.o test_net.csv
	./$^ $@

test_es.o: test_es.cpp ../EbolaSim.h
	g++ $(FLAGS) $< -o $@ $(INCLUDE) $(LINK) -lsim $(GSL_LIB)

test_es_abc: test_es_abc.cpp ../EbolaSim.h ../IDCommunity.h
	g++ $(FLAGS) $< -o $@ $(INCLUDE) -I$(SQL_PATH) -I$(ABC_PATH) $(ABC_LIB) $(LINK) -lsim $(GSL_LIB)

%.sqlite: test_es_abc %.json
	./$^ --process

test_abc_example.sqlite: test_es_abc abc_ebola_sim.json
	make cleanwork
	./$^ --process

test_abc_example_nonleaky.sqlite: test_es_abc abc_ebola_sim_nonleaky.json
	./$^ --process

test_abc_example_nobias.sqlite: test_es_abc abc_ebola_sim_nobias.json
	./$^ --process

SRCDIR := networks
SRCS := $(wildcard $(SRCDIR)/*_lvl4.csv)
TARDIR := declustered
TARS := $(subst $(SRCDIR),$(TARDIR),$(SRCS))
TAR2DIR := reclustered
TARS2 := $(subst $(SRCDIR),$(TAR2DIR),$(SRCS))

$(TARDIR) $(TAR2DIR):
	mkdir -p $@

$(TARDIR)/%.csv $(TAR2DIR)/%.csv: ~/workspaces/network_prune/decluster $(SRCDIR)/%.csv | $(TARDIR) $(TAR2DIR)
	$^ $(TARDIR)/$*.csv $(TAR2DIR)/$*.csv

declusterall: $(TARS)

FORCE:

TARGET := abc_ebola_km

THREADS = 6
SEQ = $(strip $(shell seq $(THREADS)))
LOGS = $(patsubst %,work%.log,$(SEQ))

NPER = $(shell [ -f $(TARGET).sqlite ] && sqlite3 $(TARGET).sqlite 'select count(*)/$(THREADS) from par;')

work%.log: test_es_abc $(TARGET).json FORCE | $(TARGET).sqlite
	./$(filter-out FORCE,$^) --simulate -n $(NPER) >> $@

network_filenames.txt: $(wildcard networks/*.csv) | networks
	cd $|; ls *.csv > ../$@

dowork: $(LOGS)
	sqlite3 $(TARGET).sqlite 'select count(*) from job where status="Q";'
	wc -l *.log

cleanwork:
	rm -f $(TARGET).sqlite
	rm -f $(LOGS)

test_result.csv: test_es.o test_net.csv
	./$^ > $@
