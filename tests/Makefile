-include ../local.mk

FLAGS ?= -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
#
# WORKSPACE ?= $(HOME)/work
#
# ABC_PATH = $(WORKSPACE)/AbcSmc
# # override these to use system copies, rather than AbcSmc copies
# GSL_PATH ?= $(ABC_PATH)/gsl_local
# SQL_PATH ?= $(ABC_PATH)/sqdb

EPI_PATH = $(WORKSPACE)/EpiFire/src
# ABC_LIB = -L$(ABC_PATH) -labc -ljsoncpp -lsqdb $(ABC_PATH)/sqlite3.o
GSL_LIB = -lm -lgsl -lgslcblas -pthread -Wl,--no-as-needed -ldl

#INCLUDE = -I$(ABC_PATH) -I$(GSL_PATH)/include/ -I$(EPI_PATH)
INCLUDE = -I$(EPI_PATH)
LINK = -L$(EPI_PATH)
#
# LDFLAGS=  $(EPI_PATH)/*.o

$(EPI_PATH)/libsim.a: $(EPI_PATH)
	$(MAKE) -C $^

epifire: $(EPI_PATH)/libsim.a

test_id_comm.o: test_id_comm.cpp epifire
	g++ $(FLAGS) $< -o $@ $(INCLUDE) $(LINK) -lsim $(GSL_LIB)

test_id_comm.dot: test_id_comm.o test_net.csv
	./$^ $@

test_priorityq_comp.o: test_priorityq_comp.cpp
	g++ $(FLAGS) $< -o $@

test_eventdriven.o: test_eventdriven.cpp
	g++ $(FLAGS) $< -o $@

test_classmethodmap.o: test_classmethodmap.cpp
	g++ $(FLAGS) $< -o $@

test_rnghandling.o: test_rnghandling.cpp
	g++ $(FLAGS) $< -o $@

test_es.o: test_es.cpp ../EbolaSim.h
	g++ $(FLAGS) $< -o $@ $(INCLUDE) $(LINK) -lsim $(GSL_LIB)

.PHONY: run_test_es

run_test_es: test_es.o test_net.csv
	./$^ $@
