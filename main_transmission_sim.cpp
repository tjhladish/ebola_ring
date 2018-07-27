#include "Ring_Generator.h"
#include "Event_Driven_Ebola_Sim.h"

/*
enum EventType { StoE_EVENT,
                 EtoI_EVENT,
                 ItoR_EVENT,
                 ItoH_EVENT,
                 ItoH_PZERO_EVENT,
                 ItoD_EVENT,
                 V1_EVENT,
                 V2_EVENT,
                 NUM_OF_EVENT_TYPES }; // must be last*/

// mean  = alpha/beta
// sd    = sqrt(alpha)/beta
// alpha = (mean/sd)**2
// beta  = mean/alpha --> mean/(sd**2 * beta**2) --> (mean/sd**2)**(1/3)

inline double gamma_alpha(double mean, double sd) {return pow(mean/sd, 2);}
inline double gamma_beta(double mean, double sd) {return pow(mean/pow(sd,2), 1.0/3);}

map<EventType, function<double(mt19937&)> > initialize_event_generator() {
    map<EventType, function<double(mt19937&)> > event_generator = {
        {StoE_EVENT,      gamma_distribution<double>(gamma_alpha(2,1), gamma_beta(2,1))},   // 2,1 are totally arbitary--prob needs to be fit
        {EtoI_EVENT,      gamma_distribution<double>(gamma_alpha(6,2), gamma_beta(6,2))},   // 6,2
        {ItoR_EVENT,      gamma_distribution<double>(gamma_alpha(9,4), gamma_beta(9,4))},   // 9,4
        {ItoH_EVENT,      gamma_distribution<double>(gamma_alpha(2,1), gamma_beta(2,1))},   // 2,1
        {ItoH_PZERO_EVENT,gamma_distribution<double>(gamma_alpha(5,3), gamma_beta(5,3))},   // 5,3
        {ItoD_EVENT,      gamma_distribution<double>(gamma_alpha(8,4), gamma_beta(8,4))},   // 8,4 arbitary numbers
        {V1_EVENT,        gamma_distribution<double>(gamma_alpha(5,3), gamma_beta(5,3))},   // 5,3
        {V2_EVENT,        uniform_real_distribution<double>(28.0,28.0)},                    // 5,3
    };
    return event_generator;
}

int main() { 
    const unsigned int seed = 0;
    Network* net = generate_ebola_network(seed); // omit seed argument for seed based on current time
    Node* p_zero = net->get_nodes()[0];          // not elegant, but works for now

    Vaccine vaccine;
    //vaccine.efficacy = {1.0, 0.0};   // single dose vaccine
    vaccine.efficacy = {0.8, 0.9}; // two dose vaccine
    vaccine.coverage = {0.5, 1.0};
    vaccine.timeToProtection = 7;
    vaccine.isLeaky  = false;

    map<EventType, function<double(mt19937&)> > event_generator = initialize_event_generator();

    for(int i=0; i<1; i++ ) {
        Event_Driven_Ebola_Sim sim(net, event_generator, vaccine);
        sim.expose(p_zero);
        sim.run_simulation();
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
