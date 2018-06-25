#include "Ring_Generator.h"
#include "Event_Driven_Ebola_Sim.h"

int main() { 
    Node* p_zero;
    Network* net = generate_ebola_network(p_zero);
    //Network net = Network("ebola_pop", Network::Undirected);
    //net.populate(10000);
    //net.fast_random_graph(8);

    Vaccine vaccine;
    vaccine.efficacy = {1.0, 0.0};   // single dose vaccine
    vaccine.efficacy = {0.8, 0.9}; // two dose vaccine
    vaccine.coverage = 0.5;
    vaccine.isLeaky  = false;

    for(int i=0; i<1; i++ ) {
        Event_Driven_Ebola_Sim sim(net, vaccine);
        sim.expose(p_zero);
        sim.run_simulation();
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
