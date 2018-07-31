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
random_device true_rng;

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

void initialize_parameters(vector<double> &abc_pars, NetParameters &netpar, SimParameters &simpar) {
    // NB: There are two different random seeds, one for network generation and one for transmission
    // modeling.  This is so that they can separately be varied or fixed.  Variable seeds can be
    // chosen in a number of ways; two that both work on the dev machine are as follows:
    //
    // seed = chrono::system_clock::now().time_since_epoch().count(); // set seed using system clock
    // seed = true_rng(); // use a true (not pseudo) random number generator to seed
    //
    // It seems plausible that one or the other of these might not work, given alternate sys architecture

    // Network construction parameters
    netpar.N = 1e3;
    netpar.clusters = 200;
    netpar.mean_deg = 16.0;
    netpar.cluster_kernel_sd = 0.01;
    netpar.wiring_kernel_sd = 0.094;
    netpar.seed = 0;

    // Transmission model parameters
    Vaccine vac;
    //vac.efficacy = {1.0, 0.0};   // single dose vaccine
    vac.efficacy = {0.8, 0.9}; // two dose vaccine
    vac.coverage = {0.5, 1.0};
    vac.timeToProtection = 7;
    vac.isLeaky  = true;

    simpar.vaccine = vac;
    simpar.seed = true_rng();
    simpar.event_generator = initialize_event_generator();
    simpar.prob_quarantine = 0.33;     // start 1/3, max 2/3
    simpar.prob_community_death = 0.8; // 80%
}

int main() { 
    // pull out the network parameterization
    // parameterize quarantine & death probs
    vector<double> abc_pars;
    NetParameters netpar = {};
    SimParameters simpar = {};
    initialize_parameters(abc_pars, netpar, simpar);

    Network* net = generate_ebola_network(netpar); // omit seed argument for seed based on current time
    Node* p_zero = net->get_nodes()[0];          // not elegant, but works for now

    simpar.network    = net;
    simpar.index_case = p_zero;

    for(int i=0; i<1; i++ ) {
        Event_Driven_Ebola_Sim sim(simpar);
        //Event_Driven_Ebola_Sim sim(net, event_generator, vaccine);
        sim.expose(p_zero);
        sim.run_simulation();
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
