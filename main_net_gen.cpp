#include "Gaussian_Ring_Generator.h"
#include "AbcSmc.h"

const vector<double>
  hh_nbinom = {0.0092685312656875431386, 0.029773374281162414551, 0.056065554030988888623, 0.080734397804623894701,
               0.098371758555788293976, 0.10678632744148358613, 0.10645775412627803136, 0.099391766929326080815,
               0.088076396540448381733, 0.074797185985119141183, 0.061310677989033279811, 0.048774144951275709425,
               0.037818721623758405626, 0.028679570197046585361, 0.021330036383913390796, 0.015593897368670185319,
               0.011227606105442555079, 0.0079741405172409872415, 0.0055941662705567551939, 0.0038810377818729062047,
               0.0026653773320431668574, 0.0018136281802210226389, 0.0012236282994945771839, 0.00081913525226365942102,
               0.00054440989073523289848, 0.00035941103432662023249, 0.00023580766559843656968, 0.00015381915417497982928,
               9.9796514972427065098e-05, 6.4420371044269903967e-05, 4.1387610689364768039e-05, 2.6471639035460803064e-05,
               1.6860397785662649417e-05, 1.0696330679128077722e-05, 6.7604681866996148138e-06, 4.2577577221552508418e-06,
               2.672561770214363352e-06, 1.6722124539786208069e-06, 1.0431356069191315491e-06, 6.4884269230378000514e-07,
               4.0248211313366749343e-07, 2.4901032200179913214e-07, 1.53674941578253143e-07, 9.4613180489395049441e-08,
               5.8117634784534615102e-08, 3.562163953562825933e-08, 2.1787576716641175238e-08, 1.3299336517575369752e-08,
               8.1023650168612471966e-09, 4.9270519826238454307e-09, 2.9907963542524132368e-09, 1.8123413925858968321e-09,
               1.0964129229466391519e-09, 6.6223977070169719924e-10, 3.9938152325394750004e-10, 2.4050029182128833673e-10,
               1.4461731833473386342e-10, 8.6840650426752217931e-11, 5.2076748669708836478e-11, 3.1189015995980902304e-11,
               1.8655829875749597709e-11, 1.1145505626652290984e-11, 6.6507970052043342322e-12, 3.964167357827279896e-12,
               2.3602042576602429935e-12, 1.4037210079641657934e-12, 8.3398697088555893618e-13, 4.9499185358094166446e-13,
               2.9350105200564112427e-13, 1.7386256298073286386e-13, 1.0289606804275167619e-13, 6.0841341381378175556e-14,
               3.5943192446844916651e-14, 2.1215952245547885132e-14, 1.2512559943902100415e-14, 7.3735553244804332372e-15,
               4.3417404307385869419e-15, 2.5545568872005925849e-15, 1.501898061019234089e-15, 8.8236145481653870023e-16,
               5.1801404793536392352e-16, 3.0390157478874766475e-16, 1.7816728722324022027e-16, 1.0438390388419359135e-16,
               6.1116202186261735985e-17, 3.5760444964297906486e-17, 2.0911224060697337661e-17, 1.222058590199670071e-17,
               7.1375058359075278974e-18, 4.1662799838708078652e-18, 2.4305436459751034724e-18, 1.4171528193134194482e-18,
               8.2583514126144413715e-19, 4.8099261627018396098e-19, 2.7999910598196868199e-19, 1.6291186849748790051e-19,
               9.4739517372385194929e-20, 5.5067626211946162105e-20, 3.1992821206463204844e-20, 1.8578208901990910102e-20,
               1.0783364083912647632e-20 };

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
    netpar.desired_levels = 4; // index case is a level
    netpar.mean_deg = 16.0;
    // N = 10x larger than the net size given simple branching process
    // (using sum of first n terms of geometric series)
    const double safety_factor = 2; // was 10
    double desired_N = safety_factor*(1.0-pow(netpar.mean_deg, netpar.desired_levels)) / (1.0 - netpar.mean_deg);//1e4;
    netpar.clusters = 0;
    double avg_hh_size = 0;
    for (unsigned int i=0; i<hh_nbinom.size(); ++i) avg_hh_size += i*hh_nbinom[i];

    netpar.clusters = (int) (desired_N / avg_hh_size);
    cerr << "Desired N based on # levels, expected degree & safety factor: " << desired_N << endl;
    cerr << "Expected household size based on HH size dist: " << avg_hh_size << endl;
    cerr << "Resulting # households: " << netpar.clusters << endl;

    netpar.hh_dist = discrete_distribution<int>(hh_nbinom.begin(), hh_nbinom.end());

    netpar.between_cluster_sd = abc_args[0];
    netpar.within_cluster_sd  = abc_args[1]; //0.01;
    netpar.seed = abc_args[0];
}


struct InterviewProbabilities {
    InterviewProbabilities() : fixed(1.0) {};
    InterviewProbabilities(double f) : fixed(f) {};
    double fixed;
};

enum InterviewState {UNDETECTED, DETECTED, INTERVIEWED, NUM_OF_INTERVIEW_STATES};

Network* interview_network(Network* net, const map<Node*, int> &level_of, const InterviewProbabilities &ip, const unsigned long int seed, map<Node*, int> &ilevel_of, vector<set<Node*, NodePtrComp> > &ilevels
) {
    mt19937 rng(seed);
    uniform_real_distribution<double> runif(0.0, 1.0);

    Network* inet = net->duplicate();
    set<Node*, NodePtrComp> interviewed_nodes;
    for (Node* node: inet->get_nodes()) {
        const int lvl = level_of.at(net->get_node(node->get_id()));
        ilevel_of[node] = lvl;
        ilevels[lvl].insert(node);
        switch (lvl) {
          case 0:
            interviewed_nodes.insert(node);
            break;
          case 1:
          case 2:
            if (runif(rng) < ip.fixed) interviewed_nodes.insert(node);
            break;
          default:
            break;
        }
    }

    set<Node*, NodePtrComp> detected_nodes(interviewed_nodes);
    for (Node* node: interviewed_nodes) {
        for (Node* neighbor: node->get_neighbors()) detected_nodes.insert(neighbor);
    }

    vector<Node*> inodes = inet->get_nodes(); // still have everyone at this point
    for (Node* inode: inodes) { // label nodes with interviewed/detected status
        const int id = inode->get_id();
        int state = interviewed_nodes.count(inode) != 0 ? INTERVIEWED :
                    detected_nodes.count(inode)    != 0 ? DETECTED : UNDETECTED;
        net->get_node(id)->set_state(state);
    }

    vector<Node*> undetected_nodes;
    set_difference(inodes.begin(), inodes.end(), 
                   detected_nodes.begin(), detected_nodes.end(), 
                   inserter(undetected_nodes, undetected_nodes.begin()));
    for (Node* unode: undetected_nodes) inet->delete_node(unode);
    return inet;
}

struct TrialRawMetrics {
    vector<int> l1_size;
    vector<int> l2_size;
    vector<double> l1_l2_ratio;
};

struct InterviewRawMetrics {
    vector<int> l1_size;
    vector<int> l2_size;
    vector<int> l3_size;
    vector<double> l1_l2_ratio;
    vector<double> l111_trans; 
    vector<double> l112_trans;
    vector<double> l122_trans;
    vector<double> l1_log_component_to_size_ratio;
    vector<double> mean_path_diameter_ratio;
};

enum TransitivityType {L111, L112, L122, L222, NUM_OF_TRANSITIVITY_TYPES}; // 0, 1, 2, or 3 nodes in level 2

vector<double> special_transitivity (Network* net, map<Node*, int> level_of) {
    vector<Node*> node_set = net->get_nodes();
    vector<int> triangles(NUM_OF_TRANSITIVITY_TYPES, 0);
    vector<int> tripples(NUM_OF_TRANSITIVITY_TYPES, 0);

    for (Node* a: node_set) {
        if (level_of[a] < 1 or level_of[a] > 2) continue;

        for (Node* b: a->get_neighbors()) {
            if (level_of[a] < 1 or level_of[a] > 2) continue;

            for (Node* c: b->get_neighbors()) {
                if (level_of[a] < 1 or level_of[a] > 2) continue;
                if ( c == a ) continue;
                const int l2_ct = (level_of[a] == 2) + (level_of[b] == 2) + (level_of[c] == 2);
                if ( c->is_neighbor(a) ) triangles[l2_ct]++;
                tripples[l2_ct]++;
            }
        }
    }

    vector<double> trans(NUM_OF_TRANSITIVITY_TYPES);
    for (int i = 0; i < NUM_OF_TRANSITIVITY_TYPES; ++i) trans[i] = (double) triangles[i] / tripples[i];
    return trans;
}


pair<double, double> calc_mean_and_max(vector< vector<double> > distance_matrix) {

    double total = 0;
    int ct = 0;
    double max = 0;
    for( unsigned int i=0; i < distance_matrix.size(); i++ ) {
        for( unsigned int j=0; j < distance_matrix[i].size(); j++ ) {
            if (i != j) {  // don't consider distance from nodes to themselves
                const double d = distance_matrix[i][j];
                total += d;
                max = d > max ? d : max;
                ct++;
            }
        }
    }
    const double mean = total / ct;
    return make_pair(mean, max);
}



void raw_metrics(Network* net, vector<set<Node*, NodePtrComp> > levels, map<Node*, int> level_of, TrialRawMetrics& trm, InterviewRawMetrics& irm, bool do_interview, const unsigned long int rng_seed) {
    //net->write_edgelist("./revpar_testnets/net_" + ABC::toString(serial) + ".csv", Network::NodeIDs);
    trm.l1_size.push_back(levels[1].size());
    trm.l2_size.push_back(levels[2].size());
    trm.l1_l2_ratio.push_back((double) levels[1].size() / levels[2].size());

   if (do_interview) { 
        InterviewProbabilities ip(0.5);
        map<Node*, int> ilevel_of;
        vector<set<Node*, NodePtrComp> > ilevels(levels.size());

        Network* inet = interview_network(net, level_of, ip, rng_seed+1, ilevel_of, ilevels);
        vector<double> s_trans = special_transitivity(net, level_of);

        // Determine # components (after removing L0) and # of interviewed L1 nodes
        Network* tmp_net = net->duplicate();
        tmp_net->delete_node(tmp_net->get_node(0)); // delete the index case
        const int comp_ct = tmp_net->get_components().size();
        int l1i_size = 0;
        for (Node* n: levels[1]) l1i_size += n->get_state() == INTERVIEWED;

        vector< vector<double> > distance_matrix;
        vector<Node*> inet_nodes = inet->get_nodes();
        inet->calculate_distances(inet_nodes, distance_matrix);
        pair<double, double> mean_and_max = calc_mean_and_max(distance_matrix);

        irm.l1_size.push_back(ilevels[1].size());
        irm.l2_size.push_back(ilevels[2].size());
        irm.l3_size.push_back(ilevels[3].size());
        irm.l1_l2_ratio.push_back((double) ilevels[1].size() / ilevels[2].size());
        irm.l111_trans.push_back(s_trans[L111]);
        irm.l112_trans.push_back(s_trans[L112]);
        irm.l122_trans.push_back(s_trans[L122]);
        irm.l1_log_component_to_size_ratio.push_back(log((double) comp_ct/l1i_size));
        irm.mean_path_diameter_ratio.push_back(mean_and_max.first / mean_and_max.second);
        //cout << "Interviewed network size: " << inet->size() << endl;
    }

        //double trans = net->transitivity();
        //trans = isfinite(trans) ? trans : -99999.9;

        //double inner_trans = net->transitivity(inner_nodes);
        //inner_trans = isfinite(inner_trans) ? inner_trans : -99999.9;
        //net->validate();
        //net->write_edgelist(ABC::toString(serial) + "_lvl4.csv", Network::NodeIDs);
        //cerr << net->size() << " " << inner_nodes.size() << " " << trans << " " << inner_trans << endl;
        //vector<double> metrics = {(double) p_zero->deg(), (double) net->size(), trans, inner_trans};
}


vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int /*serial*/, const ABC::MPI_par* /*mp*/) {
    // pull out the network parameterization
    // parameterize quarantine & death probs
    //vector<double> abc_pars = {4.0};
    //vector<double> abc_pars = {(double) atoi(argv[1])};
    NetParameters netpar = {};
    initialize_parameters(args, netpar);
    const int trial_networks = 2;  //100;
    const int interviewed_networks = 1;
    assert(trial_networks >= interviewed_networks);

    TrialRawMetrics trm;
    InterviewRawMetrics irm;

    //vector<vector<double>> level_sizes(2, vector<double>(REPS,0.0));

    for (unsigned int rep = 0; rep < trial_networks; ++rep) {
        netpar.seed = rng_seed + rep;

        vector<set<Node*, NodePtrComp> > levels(netpar.desired_levels, set<Node*, NodePtrComp>());
        map<Node*, int> level_of;

        Network* net = generate_ebola_network(netpar, levels, level_of); // omit seed argument for seed based on current time

        const bool do_interview = rep < interviewed_networks;
        raw_metrics(net, levels, level_of, trm, irm, do_interview, rng_seed);
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
