#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <Network.h>
#include <unordered_set>

using namespace std;


int main() {
    const unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine rng(seed);
    uniform_real_distribution<double> runif(0.0, 1.0);
    uniform_real_distribution<double> runif_loc(-1.0, 1.0);
    //normal_distribution<double> rnorm(0.0, 0.1); 
    const int N = 1e3;
    const int mean_deg = 10;
    const double kernel_sd = 0.2;

    Network ebola_ring("asdf", Network::Undirected);
    ebola_ring.populate(N);

    vector<Node*> nodes = ebola_ring.get_nodes();
    Node* p_zero = nodes[0];

    vector<pair<double, double>> coords(nodes.size());
    coords[0]  = make_pair(0.0, 0.0);  // for convenience, put p_zero at the origin
   
    vector<double> weights(nodes.size(), 0.0);
    weights[0] = 0.0;                  // p-zero shouldn't get connected to itself
    
    double total_weight = 0.0;
    for (unsigned int i = 1; i < nodes.size(); ++i) {
        assert(i = nodes[i]->get_id()); // this allows us to use index in nodes and node id interchangeably later
        const double x = runif_loc(rng);
        const double y = runif_loc(rng);
        coords[i] = make_pair(x, y);
        const double distance = sqrt(x*x + y*y); // get away with pythagorean, because ref is origin
        weights[i] = normal_pdf(distance, 0.0, pow(kernel_sd, 2));
        total_weight += weights[i];
    }

    double weight_coef = (double) mean_deg / total_weight;
    vector<Node*> ring_1;

    // this yields an expected degree for p_zero.
    // we may want to change to a while loop that samples
    // until exactly the desired degree is achieved
    for (unsigned int i = 1; i < nodes.size(); ++i) {
        if (runif(rng) < weight_coef*weights[i]) {
            p_zero->connect_to(nodes[i]);
            ring_1.push_back(nodes[i]);
        }
    }
    cerr << "p-zero degree: " << p_zero->deg() << " " << ring_1.size() << endl;

    // now define the 2-step neighborhood, disallowing
    // connections back to p_zero

    unordered_set<Node*> ring_2;

    for (unsigned int j = 0; j < ring_1.size(); ++j) {
        Node* ring1_node = ring_1[j];
        const double x1 = coords[ring1_node->get_id()].first;
        const double y1 = coords[ring1_node->get_id()].second;
        weights.clear();
        weights.resize(nodes.size(), 0.0);
        
        double total_weight = 0.0;
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            if (i == 0 or i == j) {
                // ring 1 nodes are already connected to p_zero
                // node can't be connected to itself
                continue;
            }
            const double x2 = coords[ring1_node->get_id()].first;
            const double y2 = coords[ring1_node->get_id()].second;
            const double distance = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
            weights[i] = normal_pdf(distance, 0.0, pow(kernel_sd, 2));
            total_weight += weights[i];
        }

        double weight_coef = (double) mean_deg / total_weight;

        // this yields an expected degree of mean_deg
        for (unsigned int i = 1; i < nodes.size(); ++i) {
            if (runif(rng) < weight_coef*weights[i]) {
                ring1_node->connect_to(nodes[i]);
                ring_2.insert(nodes[i]);
            }
        }
    }

	for (unsigned int i = 0; i < nodes.size(); ++i) {
	    if (nodes[i]->deg() == 0) ebola_ring.delete_node(nodes[i]);
    }

    cerr << "Total size: " << ebola_ring.size() << endl;
    cerr << "Ring 1 size: " << ring_1.size() << endl;
    cerr << "Ring 2 size: " << ring_2.size() << endl;
    cerr << "Transitivity clustering coefficient: " << ebola_ring.transitivity() << endl;
    ebola_ring.write_edgelist("ebola_ring", Network::NodeIDs);
}
