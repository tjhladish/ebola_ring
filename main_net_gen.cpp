#include "Ring_Generator.h"
#include "Event_Driven_Ebola_Sim.h"
#include "AbcSmc.h"

const gsl_rng* GSL_RNG = gsl_rng_alloc (gsl_rng_taus2); // RNG for AbcSmc

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

void initialize_parameters(vector<double> &abc_args, NetParameters &netpar, SimParameters &simpar) {
    // NB: There are two different random seeds, one for network generation and one for transmission
    // modeling.  This is so that they can separately be varied or fixed.  Variable seeds can be
    // chosen in a number of ways; two that both work on the dev machine are as follows:
    //
    // seed = chrono::system_clock::now().time_since_epoch().count(); // set seed using system clock
    // seed = true_rng(); // use a true (not pseudo) random number generator to seed
    //
    // It seems plausible that one or the other of these might not work, given alternate sys architecture

    // Network construction parameters
    netpar.desired_levels = 3; // index case is a level
    netpar.N = 1e4;
    netpar.clusters = 2000;
    netpar.mean_deg = 16.0;
    netpar.cluster_kernel_sd = abc_args[0]; //0.01;
    netpar.wiring_kernel_sd  = abc_args[1]; //0.094;
    netpar.seed = abc_args[0];

    // Transmission model parameters
    Vaccine vac;
    //vac.efficacy = {1.0, 0.0};   // single dose vaccine
    vac.efficacy = {0.8, 0.9}; // two dose vaccine
    vac.coverage = {0.5, 1.0};
    vac.timeToProtection = 7;
    vac.isLeaky  = true;

    simpar.vaccine = vac;
    simpar.seed = 0;//true_rng();
    simpar.event_generator = initialize_event_generator();
    // TODO - finish implementing:
    //simpar.prob_quarantine = 0.33;     // start 1/3, max 2/3
    //simpar.prob_community_death = 0.8; // 80%
}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* mp) {
    // pull out the network parameterization
    // parameterize quarantine & death probs
    //vector<double> abc_pars = {4.0};
    //vector<double> abc_pars = {(double) atoi(argv[1])};
    NetParameters netpar = {};
    SimParameters simpar = {};
    initialize_parameters(args, netpar, simpar);
    netpar.seed = rng_seed;
    simpar.seed = rng_seed;

    Network* net = generate_ebola_network(netpar); // omit seed argument for seed based on current time
    Node* p_zero = net->get_nodes()[0];          // not elegant, but works for now

    simpar.network    = net;
    simpar.index_case = p_zero;

//    net->write_edgelist("net_w_clustering.csv", Network::NodeIDs);
//    remove_clustering(net, rng);
//    net->validate();
//    net->write_edgelist("net_wo_clustering.csv", Network::NodeIDs);
    //net->dumper();

//    cerr << "Total size after pruning: " << net->size() << endl;
//    cerr << "Transitivity clustering coefficient after pruning: " << net->transitivity() << endl;

/*    for(int i=0; i<1; i++ ) {
        Event_Driven_Ebola_Sim sim(simpar);
        sim.expose(p_zero);
        sim.run_simulation();
        //cout << sim.current_epidemic_size() << endl;
    }*/
    vector<double> metrics = {(double) p_zero->deg(), (double) net->size()};

    return metrics;
}




void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --process --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(GSL_RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(GSL_RNG);
    }

    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
