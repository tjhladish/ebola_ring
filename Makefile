FLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic

ABC_PATH = $(HOME)/work/AbcSmc
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
EPI_PATH = $(HOME)/work/EpiFire/src
SQLDIR = $(ABC_PATH)/sqdb
ABC_LIB = -L$(ABC_PATH) -labc -ljsoncpp -lsqdb $(ABC_PATH)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

INCLUDE = -I$(ABC_PATH) -I$(GSL_PATH)/include/ -I$(EPI_PATH)

LDFLAGS=  $(EPI_PATH)/*.o

default: ebola_net

epifire: 
	$(MAKE) -C $(EPI_PATH) -f Makefile

libabc:
	$(MAKE) -C $(ABC_PATH) -f Makefile

ebola_net: epifire main.cpp
	g++ $(FLAGS) main.cpp -o ebola $(INCLUDE) $(GSL_LIB) $(LDFLAGS)
