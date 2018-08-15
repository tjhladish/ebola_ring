#ifndef RING_GEN_H
#define RING_GEN_H

#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <Network.h>
#include <set>
#include <numeric>

using namespace std;

struct NetParameters {
    NetParameters() {
        desired_levels = 3; // index case counts as a level
        N = 1e3;
        clusters = 200;
        mean_deg = 16.0;
        cluster_kernel_sd = 0.01;
        wiring_kernel_sd = 0.094;
        seed = 0;
    }

    NetParameters(int n, int c, double md, double ck_sd, double wk_sd, unsigned int s) {
        N = n;
        clusters = c;
        mean_deg = md;
        cluster_kernel_sd = ck_sd;
        wiring_kernel_sd = wk_sd;
        seed = s;
    }

    int desired_levels;
    int N;
    int clusters;
    double mean_deg;
    double cluster_kernel_sd;
    double wiring_kernel_sd;
    unsigned int seed;
};


vector<pair<double, double>> generate_spatial_distribution(int N, int clusters, double cluster_kernel_sd, mt19937& rng) {
    assert(N >= clusters);
    vector<pair<double, double>> coords(N);
    coords[0]  = make_pair(0.0, 0.0);  // for convenience, always put p_zero at the origin
    uniform_real_distribution<double> runif_loc(-1.0, 1.0);

    for (int i = 1; i < clusters; ++i) { // p_zero anchors one cluster; where are the others?
        const double x = runif_loc(rng);
        const double y = runif_loc(rng);
        coords[i] = make_pair(x, y);
    }

    uniform_int_distribution<int> rand_clust(0,clusters-1); // [a,b] inclusive
    normal_distribution<double> cluster_noise(0.0,cluster_kernel_sd);
    for (int i = clusters; i < N; ++i) { // for each remaining person
        const int cluster_idx = rand_clust(rng);
        const double x = coords[cluster_idx].first + cluster_noise(rng);
        const double y = coords[cluster_idx].second + cluster_noise(rng);
        coords[i] = make_pair(x, y);
    }

    return coords;
}


//vector<double> calc_weights(const vector<pair<double, double>> &coords, int reference_node, double wiring_kernel_sd) {
vector<double> calc_weights(const int reference_node_idx, const vector<Node*> &nodes, const vector<pair<double, double>> &coords, double wiring_kernel_sd, set<Node*> zero_weight_nodes) {
    assert(reference_node_idx < (signed) coords.size());
    assert(nodes.size() == coords.size());
    vector<double> weights(coords.size(), 0.0); // important to initialize vals to 0
    const double wiring_kernel_var = pow(wiring_kernel_sd, 2);

    const double x1 = coords[reference_node_idx].first;
    const double y1 = coords[reference_node_idx].second;
    for (unsigned int i = 0; i < weights.size(); ++i) {
        if (zero_weight_nodes.count(nodes[i]) > 0) {
            // never connect back to p_zero, to self (sometimes that's the same thing), or to an inner ring (already wired)
            // leave weight at 0
            continue;
        }
        const double x2 = coords[i].first;
        const double y2 = coords[i].second;
        const double distance = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
        weights[i] = normal_pdf(distance, 0.0, wiring_kernel_var);
    }

    return weights;
}


struct NodePtrComp { bool operator()(const Node* A, const Node* B) const { return A->get_id() < B->get_id(); } };


bool is_node_in_level(Node* const n, const set<Node*, NodePtrComp> &level) { return level.count(n) > 0; }


Network* generate_ebola_network(const NetParameters &par, map<Node*, int> &level_of) {
    const int N = par.N;
    const int clusters = par.clusters;
    const double mean_deg = par.mean_deg;
    const double cluster_kernel_sd = par.cluster_kernel_sd;
    const double wiring_kernel_sd = par.wiring_kernel_sd;
    const unsigned int seed = par.seed;
    const unsigned int desired_levels = par.desired_levels;
    vector<set<Node*, NodePtrComp> > levels(desired_levels, set<Node*, NodePtrComp>());

    mt19937 rng(seed);
    uniform_real_distribution<double> runif(0.0, 1.0);

    Network* ebola_ring = new Network("ebola_ring", Network::Undirected);
    ebola_ring->populate(N);

    vector<Node*> nodes = ebola_ring->get_nodes();
    const int p_zero_idx = 0; // important--we use index to identify p_zero for triggering mass vaccination of group
    Node* p_zero = nodes[p_zero_idx];
    assert(p_zero_idx == nodes[p_zero_idx]->get_id()); // only needs to be true for p_zero
    levels[0].insert(p_zero);
    level_of[p_zero] = 0;

    // locate nodes, get weights for wiring p_zero
    vector<pair<double, double>> coords = generate_spatial_distribution(N, clusters, cluster_kernel_sd, rng);
    set<Node*> zero_weight_nodes = {p_zero}; // no incoming edges to p_zero

    // this yields an expected degree for each node
    for (unsigned int level_idx = 0; level_idx < levels.size(); ++level_idx) { // for each level, inner to outer
    //cerr << "level idx, size: " << level_idx << ", " << levels[level_idx].size() << endl;
        for (const auto& inner_level_node: levels[level_idx]) {                 // look at each node in level
            // calculate weights for whether that node should be connected to others
            // connections from outer levels back to inner should be disallowed, as well as within-level
            // connections that have already been considered (no double jeapardy)
            // so . . . how to determine the zero_weight_nodes membership?
            // TODO - calc_weights does not consider that nodes with lower id will not get wired to nodes with
            // higher id.  I *think* that's okay, because the weight matrix is symmetrical and that potential edge
            // gets considered separately
            const vector<double> weights = calc_weights(inner_level_node->get_id(), nodes, coords, wiring_kernel_sd, zero_weight_nodes);
            const double remaining_expected_degree = level_idx == 0 ? mean_deg : mean_deg - 1;
            const double weight_coef = remaining_expected_degree / accumulate(weights.begin(), weights.end(), 0.0);
            assert(remaining_expected_degree > 0);
            assert(weight_coef > 0);
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                Node* n = nodes[i];
                if (level_idx == levels.size()-1 and not is_node_in_level(n, levels[level_idx])) {
                    // if we're looking at the outermost level, we only consider within-level connections
                    continue;
                } else if (is_node_in_level(n, levels[level_idx]) and inner_level_node->get_id() < n->get_id()) {
                    // each potential connection should only be considered once--make sure we don't ask whether
                    // 'A' should be connected to 'B' *and* whether 'B' should be connected to 'A'
                    continue;
                } else if (runif(rng) < weight_coef*weights[i]) {
                    inner_level_node->connect_to(n);
                    // if other node not in this level, put it in the next level
                    if (not is_node_in_level(n, levels[level_idx])) {
                        levels[level_idx+1].insert(n);
                        level_of[n] = level_idx+1;
                        //cerr << "linking " << level_idx << " to " << level_idx+1 << endl;
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

//    cerr << "p-zero degree: " << p_zero->deg() << endl;

	for (unsigned int i = 1; i < nodes.size(); ++i) { // always leave p_zero
	    if (nodes[i]->deg() == 0) ebola_ring->delete_node(nodes[i]);
    }
    
    ebola_ring->reset_node_ids();
//    cerr << "Total size: " << ebola_ring->size() << endl;
//    cerr << "Level 1 size: " << levels[1].size() << endl;
//    cerr << "Level 2 size: " << levels[2].size() << endl;
//    cerr << "Transitivity clustering coefficient: " << ebola_ring->transitivity() << endl;
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
