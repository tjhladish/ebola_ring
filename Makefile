FLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic

ABC_PATH = $(HOME)/work/AbcSmc
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
EPI_PATH = $(HOME)/work/EpiFire/src
SQL_PATH = $(ABC_PATH)/sqdb
ABC_LIB = -L$(ABC_PATH) -labc -ljsoncpp -lsqdb $(ABC_PATH)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -pthread -Wl,--no-as-needed -ldl

INCLUDE = -I$(ABC_PATH) -I$(GSL_PATH)/include/ -I$(EPI_PATH)

LDFLAGS=  $(EPI_PATH)/*.o

default: ebola_net

epifire: 
	$(MAKE) -C $(EPI_PATH) -f Makefile

libabc:
	$(MAKE) -C $(ABC_PATH) -f Makefile

ebola_net: epifire main.cpp
	g++ $(FLAGS) main.cpp -o ebola $(INCLUDE) $(GSL_LIB) $(LDFLAGS)

ebola_abc: epifire libabc main_abc.cpp
	#g++ $(FLAGS) -Wno-ignored-attributes -Wno-misleading-indentation main_abc.cpp -o ebola_abc $(INCLUDE) -I$(SQL_PATH) $(GSL_LIB) $(ABC_LIB) $(LDFLAGS)
	g++ $(FLAGS) main_abc.cpp -o ebola_abc $(INCLUDE) -I$(SQL_PATH) $(GSL_LIB) $(ABC_LIB) $(LDFLAGS)

degvtime: epifire degree_vs_time.cpp
	g++ $(FLAGS) degree_vs_time.cpp -o degvtime $(INCLUDE) $(LDFLAGS)

ebola_sim: epifire Ring_Generator.h Event_Driven_Ebola_Sim.h main_transmission_sim.cpp
	g++ $(FLAGS) main_transmission_sim.cpp -o etsim $(INCLUDE) $(LDFLAGS)
