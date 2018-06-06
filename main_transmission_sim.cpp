#include "Ring_Generator.h"
#include "Event_Driven_Polio_Sim.h"

int main() { 
    Node* p_zero;
    Network* net = generate_ebola_network(p_zero);
    //Network net = Network("ebola_pop", Network::Undirected);
    //net.populate(10000);
    //net.fast_random_graph(8);

    double lambda  = 0.1; // dummy vals
    double epsilon = 0.1;
    double nu      = 0.1;

    for(int i=0; i<1; i++ ) {
        Event_Driven_Polio_Sim sim(net, lambda, epsilon, nu);
        //sim.rand_infect(1);
        sim.infect(p_zero);
        sim.run_simulation(1000);
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
