-include local.mk

FLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic

WORKSPACE ?= $(HOME)/work

ABC_PATH = $(WORKSPACE)/AbcSmc
# override these to use system copies, rather than AbcSmc copies
GSL_PATH ?= $(ABC_PATH)/gsl_local
SQL_PATH ?= $(ABC_PATH)/sqdb

EPI_PATH = $(WORKSPACE)/EpiFire/src
ABC_LIB = -L$(ABC_PATH) -labc -ljsoncpp -Wl,--no-as-needed -ldl -lsqdb $(ABC_PATH)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -pthread -Wl,--no-as-needed -ldl

INCLUDE = -I$(ABC_PATH) -I$(GSL_PATH)/include/ -I$(EPI_PATH)

LDFLAGS=  $(EPI_PATH)/*.o

default: ebola_sim

## EpiFire dependency

EPIFIREREPO = https://github.com/tjhladish/EpiFire.git

$(EPI_PATH):
	cd $(WORKSPACE) && git clone $(EPIFIREREPO)

$(EPI_PATH)/libsim.a: $(EPI_PATH)
	$(MAKE) -C $^

epifire: $(EPI_PATH)/libsim.a

## AbcSmc dependency

ABCSMCREPO = https://github.com/tjhladish/AbcSmc.git

$(ABC_PATH):
	cd $(WORKSPACE) && git clone $(ABCSMCREPO)

$(ABC_PATH)/libabc.a: $(ABC_PATH)
	$(MAKE) -C $^

libabc: $(ABC_PATH)/libabc.a

## Simulation proper
## convention for targets:
## resulting_binary: source_with_main.cpp assorted other dependencies
##   commands

ebola_net: main_net_gen.cpp Gaussian_Ring_Generator.h epifire libabc
	g++ $(FLAGS) $< -o $@ $(INCLUDE) -I$(SQL_PATH) $(GSL_LIB) $(ABC_LIB) $(LDFLAGS)

degvtime: degree_vs_time.cpp epifire
	g++ $(FLAGS) $< -o $@ $(INCLUDE) $(LDFLAGS)

ebola_sim: main_transmission_sim.cpp Ring_Generator.h Event_Driven_Ebola_Sim.h epifire libabc
	g++ $(FLAGS) $< -o $@ $(INCLUDE) -I$(SQL_PATH) $(GSL_LIB) $(ABC_LIB) $(LDFLAGS)

SIMS = ebola_net degvtime ebola_sim

clean:
	rm $(SIMS)
