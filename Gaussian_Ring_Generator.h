#ifndef GAUS_RING_GEN_H
#define GAUS_RING_GEN_H

#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <Network.h>
#include <set>
#include <numeric>
#include <unistd.h>

using namespace std;

struct Coord {
    Coord() {
        x = 0;
        y = 0;
    }

    Coord(double _x, double _y): x(_x), y(_y) {}

    double x;
    double y;
};


std::ostream& operator<<(std::ostream &out, const Coord &c) {
    return out << c.x << "," << c.y << endl;
}


struct NetParameters {
    NetParameters() {
        desired_levels = 3; // index case counts as a level
        expected_N = 1e3;
        clusters = 100;
        hh_dist = discrete_distribution<int>();
        mean_deg = 16.0;
        between_cluster_sd = 0.01;
        within_cluster_sd = 0.01;
        wiring_kernel_sd = 1.0;
        pzero_total_weight = mean_deg;
        seed = 0;
    }

    NetParameters(int n, int c, discrete_distribution<int> hhd, double md, double btn_sd, double wtn_sd, unsigned int s) {
        expected_N = n;
        clusters = c;
        hh_dist = hhd;
        mean_deg = md;
        between_cluster_sd = btn_sd;
        within_cluster_sd = wtn_sd;
        seed = s;
    }

    int desired_levels;
    int expected_N;
    int clusters;
    discrete_distribution<int> hh_dist;
    double mean_deg;
    double between_cluster_sd; // a scaling par
    double within_cluster_sd;
    double wiring_kernel_sd;
    double pzero_total_weight;
    unsigned int seed;
};


vector<Coord> generate_spatial_distribution(int clusters, discrete_distribution<int> hh_dist, double between_cluster_sd, double within_cluster_sd, mt19937& rng) {
    vector<Coord> cluster_coords(clusters);
    vector<int> cluster_sizes(clusters);
    vector<Coord> node_coords;
    normal_distribution<double> rnorm_loc(0.0, sqrt((double) clusters)*between_cluster_sd);

    int N = 0;
    for (unsigned int i = 0; i < cluster_sizes.size(); ++i) {
        cluster_coords[i] = Coord(rnorm_loc(rng), rnorm_loc(rng));
        cluster_sizes[i] = hh_dist(rng);
        N += cluster_sizes[i];
    }

    normal_distribution<double> cluster_noise(0.0, within_cluster_sd);

    double radius = numeric_limits<double>::infinity();
    unsigned int test_pzero = 0;
    for (unsigned int i = 0; i < cluster_sizes.size(); ++i) { // for each remaining person
        for (int j = 0; j < cluster_sizes[i]; ++j) {
            const double x = cluster_coords[i].x + cluster_noise(rng);
            const double y = cluster_coords[i].y + cluster_noise(rng);
            node_coords.push_back(Coord(x, y));
            const double dist = sqrt(x*x + y*y); // distance of this node from the origin
            if (dist < radius) {
                radius = dist;
                test_pzero = node_coords.size() - 1;
            }
//cout << i << " " << x << " " << y << endl;
        }
    }
    // NB: node_coords is no longer organized by household membership.  This isn't something we've been using, but may become of interest.
    Coord tmp = node_coords[0];
    node_coords[0] = node_coords[test_pzero];
    node_coords[test_pzero] = tmp;

    return node_coords;
}


long double normal_weight(long double x, long double mu, long double var) { return exp(-pow(x-mu,2) / (2.0*var)); }


double kahan_sum(vector<double> nums) {
    double sum = 0.0;
    double c = 0.0;
    for (auto num: nums) {
        double y = num - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

double calc_weights(const vector<Coord> &coords, const double wiring_kernel_sd) {
    const unsigned int N = coords.size();
    vector<double> pzero_basic_weights(N);

    const double wiring_kernel_var = pow(wiring_kernel_sd, 2);
    const double gaussian_threshold = 10*sqrt(wiring_kernel_var); // distances greater than this are equally likely
    double min_wt = normal_weight(gaussian_threshold, 0.0, wiring_kernel_var);

    const int i = 0;
    const double x1 = coords[i].x;
    const double y1 = coords[i].y;
    for (unsigned int j = i+1; j < N; ++j) {
        const double x2 = coords[j].x;
        const double y2 = coords[j].y;
        const double distance = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
        pzero_basic_weights[j] = distance > gaussian_threshold ? min_wt : normal_weight(distance, 0.0, wiring_kernel_var);
    }

    sort(pzero_basic_weights.begin(), pzero_basic_weights.end());

    return kahan_sum(pzero_basic_weights);
}


vector<vector<double>> recalc_weights(const vector<Node*> &nodes, const vector<Coord> &coords, const double wiring_kernel_sd, const double pzero_total_weight, const double mean_deg) {

    cerr << "pzero_total_weight: " << pzero_total_weight << endl;
    cerr << "net size: " << nodes.size() << endl;
    const unsigned int N = nodes.size();
    assert(pzero_total_weight > 0);
//    assert(pzero_total_weight < N); // this is now sometimes failing, with the change in algorithm.  I *think* that's okay
    assert(coords.size() == N);

    const double wiring_kernel_var = pow(wiring_kernel_sd, 2);
    const double gaussian_threshold = 10*sqrt(wiring_kernel_var); // distances greater than this are equally likely
    double min_wt = normal_weight(gaussian_threshold, 0.0, wiring_kernel_var);
    vector<vector<double>> wiring_probs(N, vector<double>(N, 0.0));

    const double weight_coef = pzero_total_weight < mean_deg ?
                               (mean_deg - pzero_total_weight) / (N - pzero_total_weight) :
                               mean_deg/pzero_total_weight;

    for (unsigned int i = 0; i < N; ++i) {
        const double x1 = coords[i].x;
        const double y1 = coords[i].y;
        for (unsigned int j = i+1; j < N; ++j) {
            const double x2 = coords[j].x;
            const double y2 = coords[j].y;
            const double distance = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
            const double basic_weight = distance > gaussian_threshold ? min_wt : normal_weight(distance, 0.0, wiring_kernel_var);
            wiring_probs[i][j] = pzero_total_weight < mean_deg ?
                                 basic_weight + weight_coef*(1.0 - basic_weight) :
                                 basic_weight * weight_coef;
            wiring_probs[j][i] = wiring_probs[i][j];
        }
    }

    return wiring_probs;
}


struct NodePtrComp { bool operator()(const Node* A, const Node* B) const { return A->get_id() < B->get_id(); } };


bool is_node_in_level(Node* const n, const set<const Node*, NodePtrComp> &level) { return level.count(n) > 0; }


Network* generate_ebola_network(const NetParameters &par, vector<Coord> &coords, vector<set<const Node*, NodePtrComp>> &levels, map<const Node*, int> &level_of) {
    const double mean_deg = par.mean_deg;
    const double wiring_kernel_sd = 1.0;
    const unsigned int seed = par.seed;
    discrete_distribution<int> hh_dist = par.hh_dist;
    const double pzero_total_weight = par.pzero_total_weight;

    mt19937 rng(seed);
    const int N = coords.size();

    uniform_real_distribution<double> runif(0.0, 1.0);

    Network* ebola_ring = new Network(Network::Undirected);
    ebola_ring->populate(N);

    vector<Node*> nodes = ebola_ring->get_nodes();
    const int p_zero_idx = 0; // important--we use index to identify p_zero for triggering mass vaccination of group
    Node* p_zero = nodes[p_zero_idx];
    assert(p_zero_idx == nodes[p_zero_idx]->get_id()); // only needs to be true for p_zero
    levels[0].insert(p_zero);
    level_of[p_zero] = 0;

    const vector<vector<double>> wiring_probs = recalc_weights(nodes, coords, wiring_kernel_sd, pzero_total_weight, mean_deg);

    set<const Node*, NodePtrComp> zero_weight_nodes = {p_zero}; // no incoming edges to p_zero

    // this yields an expected degree for each node
    for (unsigned int level_idx = 0; level_idx < levels.size(); ++level_idx) { // for each level, inner to outer
        for (const auto& inner_level_node: levels[level_idx]) {                 // look at each node in level
            const int self_id = inner_level_node->get_id();
            // calculate weights for whether that node should be connected to others
            // connections from outer levels back to inner should be disallowed, as well as within-level
            // connections that have already been considered (no double jeapardy)
            // so . . . how to determine the zero_weight_nodes membership?
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                Node* n = nodes[i];
                if (level_idx == levels.size()-1 and not is_node_in_level(n, levels[level_idx])) {
                    // if we're looking at the outermost level, we only consider within-level connections
                    continue;
                } else if (is_node_in_level(n, levels[level_idx]) and self_id <= n->get_id()) {
                    // each potential connection should only be considered once--make sure we don't ask whether
                    // 'A' should be connected to 'B' *and* whether 'B' should be connected to 'A'
                    // Also, disallow connections to self
                    continue;
                } else if ((runif(rng) < wiring_probs[self_id][i]) and (zero_weight_nodes.count(n) == 0)) {
                    Node* unconst_node = ebola_ring->get_node(inner_level_node->get_id());
                    unconst_node->connect_to(n);
                    // if other node not in this level, put it in the next level
                    if (not is_node_in_level(n, levels[level_idx])) {
                        levels[level_idx+1].insert(n);
                        level_of[n] = level_idx+1;
                    }
                }
            }
        }
        // Now that we're done with that level, add nodes to zero_weight set so they won't be considered further.
        // I think this needs to be done as a separate loop to allow within-level edges.
        for (const auto& inner_level_node: levels[level_idx]) {
            zero_weight_nodes.insert(inner_level_node);
        }
    }

    for (unsigned int i = 1; i < nodes.size(); ++i) { // always leave p_zero
        if (nodes[i]->deg() == 0) ebola_ring->delete_node(nodes[i]);
    }

    ebola_ring->reset_node_ids();
    /*
    ebola_ring->write_edgelist("ebola_ring.csv", Network::NodeIDs);

    cout << "idx,x,y\n";
    for (unsigned int i = 0; i < coords.size(); ++i) {
        cout << i << "," << coords[i].first << "," << coords[i].second << endl;
    }*/

    return ebola_ring;
}

void remove_clustering(Network* net, mt19937& rng) {
    vector<Node*> nodes = net->get_nodes();
    const Node* p_zero = nodes[0];
    vector<double> level_num_vec = p_zero->min_paths(nodes);
    map<Node*, int> level_of;
    for (unsigned int i = 0; i < nodes.size(); ++i) level_of[nodes[i]] = (int) level_num_vec[i];
    set<Edge*> within_level_edges_to_delete;

    for (Node* n: nodes) {
        const int l_n = level_of[n];
//        cerr << n->get_id() << ", " << l_n << endl;
        vector<Node*> inner_nodes; // nodes that are one level closer to index

        for (Edge* e: n->get_edges_out()) { // let's take a look at n's neighbors
            Node* m = e->get_end();         // some neighbor m
            const int l_m = level_of[m]; // m's level is l_m
            if (l_n == l_m) {
                // remove edges between nodes in same level
                // NB: edges in EpiFire are inherently directional;
                // need to delete edge and its complement
                within_level_edges_to_delete.insert(e);
                within_level_edges_to_delete.insert(e->get_complement());
            } else if (l_n == l_m + 1) {
                inner_nodes.push_back(m);
            } else {
                assert(l_n == l_m - 1); // only acceptable alternative is one level farther
            }
        }
        if (inner_nodes.size() > 1) {
            // Node has multiple connections to inner level
            shuffle(inner_nodes.begin(), inner_nodes.end(), rng);
            for (unsigned int j = 1; j < inner_nodes.size(); ++j) {
                n->disconnect_from(inner_nodes[j]);
            }
        }
    }
    for (Edge* e: within_level_edges_to_delete) e->delete_edge();
}

#endif
