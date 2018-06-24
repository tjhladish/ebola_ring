#ifndef RING_GEN_H
#define RING_GEN_H

#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <Network.h>
#include <unordered_set>
#include <numeric>

using namespace std;


vector<pair<double, double>> generate_spatial_distribution(int N, int clusters, double cluster_kernel_sd, default_random_engine& rng) {
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


vector<double> calc_weights(const vector<pair<double, double>> &coords, int reference_node, double wiring_kernel_sd) {
    assert(reference_node < (signed) coords.size());
    vector<double> weights(coords.size(), 0.0); // important to initialize vals to 0
    const double wiring_kernel_var = pow(wiring_kernel_sd, 2);

    const double x1 = coords[reference_node].first;
    const double y1 = coords[reference_node].second;
    for (unsigned int i = 0; i < weights.size(); ++i) {
        if (i == 0 or (signed) i == reference_node) {
            // never connect back to p_zero, or to self (sometimes that's the same thing)
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


Network* generate_ebola_network(Node* p_zero) {
    const int N = 1e3;
    const int clusters = 200;
    const double mean_deg = 16.0;
    const double cluster_kernel_sd = 0.01;
    const double wiring_kernel_sd = 0.094;

    const unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine rng(seed);
    uniform_real_distribution<double> runif(0.0, 1.0);

    Network* ebola_ring = new Network("ebola_ring", Network::Undirected);
    ebola_ring->populate(N);

    vector<Node*> nodes = ebola_ring->get_nodes();
    const int p_zero_idx = 0; // important--we use index to identify p_zero for triggering mass vaccination of group
    p_zero = nodes[p_zero_idx];

    // locate nodes, get weights for wiring p_zero
    vector<pair<double, double>> coords = generate_spatial_distribution(N, clusters, cluster_kernel_sd, rng);
    vector<double> weights = calc_weights(coords, p_zero_idx, wiring_kernel_sd);
    double weight_coef = mean_deg / accumulate(weights.begin(), weights.end(), 0.0);

    unordered_set<Node*> ring_1;

    // this yields an expected degree for p_zero.
    // we may want to change to a while loop that samples
    // until exactly the desired degree is achieved
    for (unsigned int i = 1; i < nodes.size(); ++i) {
        if (runif(rng) < weight_coef*weights[i]) {
            p_zero->connect_to(nodes[i]);
            ring_1.insert(nodes[i]);
        }
    }
    cerr << "p-zero degree: " << p_zero->deg() << " " << ring_1.size() << endl;

    // now define the 2-step neighborhood, disallowing
    // connections back to p_zero

    unordered_set<Node*> ring_2;

    for (const auto& ring1_node: ring_1) {
        // Currently, this is a 2-ring generator, as only wiring back to p_zero
        // is disallowed.  To make it an n-ring generator, wiring back to any
        // inner loop needs to be disallowed, by setting all of those weights
        // to 0 (to get the re-weighting correct) and possibly also only iterating
        // over only non-inner-ring nodes when wiring (not necessary, may or may
        // not be faster)
        weights = calc_weights(coords, ring1_node->get_id(), wiring_kernel_sd);
        weight_coef = mean_deg / accumulate(weights.begin(), weights.end(), 0.0);

        // this yields an expected degree of mean_deg
        for (unsigned int i = 1; i < nodes.size(); ++i) {
            if (runif(rng) < weight_coef*weights[i]) {
                ring1_node->connect_to(nodes[i]);
                ring_2.insert(nodes[i]);
            }
        }
    }

	for (unsigned int i = 0; i < nodes.size(); ++i) {
	    if (nodes[i]->deg() == 0) ebola_ring->delete_node(nodes[i]);
    }

    cerr << "Total size: " << ebola_ring->size() << endl;
    cerr << "Ring 1 size: " << ring_1.size() << endl;
    cerr << "Ring 2 size: " << ring_2.size() << endl;
    cerr << "Transitivity clustering coefficient: " << ebola_ring->transitivity() << endl;
    /*
    ebola_ring->write_edgelist("ebola_ring.csv", Network::NodeIDs);

    cout << "idx,x,y\n";
    for (unsigned int i = 0; i < coords.size(); ++i) {
        cout << i << "," << coords[i].first << "," << coords[i].second << endl;
    }*/
    
    return ebola_ring;
}

#endif
