#include "Gaussian_Ring_Generator.h"
#include "AbcSmc.h"

const gsl_rng* GSL_RNG = gsl_rng_alloc (gsl_rng_taus2); // RNG for AbcSmc
random_device true_rng;

void initialize_parameters(vector<double> &abc_args, NetParameters &netpar) {
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
    netpar.mean_deg = 16.0;
    // N = 10x larger than the net size given simple branching process
    // (using sum of first n terms of geometric series)
    netpar.N = (int) (1*(1.0-pow(netpar.mean_deg, netpar.desired_levels)) / (1.0 - netpar.mean_deg));//1e4;
    cerr << netpar.N << endl;
    //netpar.clusters = (int) (netpar.N / 1.0);
    netpar.clusters = (int) (netpar.N / 10.0);
    netpar.density = abc_args[0];
    netpar.cluster_kernel_sd = abc_args[1]; //0.01;
    netpar.seed = abc_args[0];
}

vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par* /*mp*/) {
    // pull out the network parameterization
    // parameterize quarantine & death probs
    //vector<double> abc_pars = {4.0};
    //vector<double> abc_pars = {(double) atoi(argv[1])};
    NetParameters netpar = {};
    initialize_parameters(args, netpar);
    const int REPS = 1;  //100;
    vector<vector<double>> level_sizes(2, vector<double>(REPS,0.0));

    for (unsigned int rep = 0; rep < REPS; ++rep) {
        netpar.seed = rng_seed + rep;

        map<Node*, int> level_of;
        Network* net = generate_ebola_network(netpar, level_of); // omit seed argument for seed based on current time
        //Node* p_zero = net->get_nodes()[0];          // not elegant, but works for now

        net->write_edgelist("./revpar_testnets/net_" + ABC::toString(serial) + ".csv", Network::NodeIDs);
    //    remove_clustering(net, rng);
    //    net->validate();
    //    net->write_edgelist("net_wo_clustering.csv", Network::NodeIDs);
        //net->dumper();

    //    cerr << "Total size after pruning: " << net->size() << endl;
    //    cerr << "Transitivity clustering coefficient after pruning: " << net->transitivity() << endl;

        // not fitting to transitivity, but want to know it for a bit of analysis

        vector<Node*> inner_nodes;
        for (Node* n: net->get_nodes()) {
            const int lev = level_of[n];
            //cout << serial << " " << n->get_id() << " " << lev << endl;
            switch (lev) {
                case 1:
                    ++level_sizes[0][rep];
                    break;
                case 2:
                    ++level_sizes[1][rep];
                    break;
                default:
                    break;
            }
            if (lev < 2) inner_nodes.push_back(n);
        }

        //double trans = net->transitivity();
        //trans = isfinite(trans) ? trans : -99999.9;

        //double inner_trans = net->transitivity(inner_nodes);
        //inner_trans = isfinite(inner_trans) ? inner_trans : -99999.9;
        //net->validate();
        //net->write_edgelist(ABC::toString(serial) + "_gauss.csv", Network::NodeIDs);
        //cerr << net->size() << " " << inner_nodes.size() << " " << trans << " " << inner_trans << endl;
        //vector<double> metrics = {(double) p_zero->deg(), (double) net->size(), trans, inner_trans};
        delete net;
    }
    vector<double> metrics(4, 0.0);
    //vector<double> metrics = {(double) p_zero->deg(), (double) net->size()};
    /*vector<double> metrics = {mean(  level_sizes[0] ),
                              stdev( level_sizes[0] ),
                              mean(  level_sizes[1] ),
                              stdev( level_sizes[1] )};
    */
    //cerr << metrics[0] << " " << metrics[1] << " " << metrics[2] << " " << metrics[3] << endl;

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
