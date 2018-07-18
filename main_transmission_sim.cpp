#include "Ring_Generator.h"
#include "Event_Driven_Ebola_Sim.h"

int main() { 
    const unsigned int seed = 0;
    Network* net = generate_ebola_network(seed); // omit seed argument for seed based on current time
    Node* p_zero = net->get_nodes()[0];          // not elegant, but works for now

    Vaccine vaccine;
    //vaccine.efficacy = {1.0, 0.0};   // single dose vaccine
    vaccine.efficacy = {0.8, 0.9}; // two dose vaccine
    vaccine.coverage = {0.5, 1.0};
    vaccine.isLeaky  = false;

    for(int i=0; i<1; i++ ) {
        Event_Driven_Ebola_Sim sim(net, vaccine);
        sim.expose(p_zero);
        sim.run_simulation();
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
