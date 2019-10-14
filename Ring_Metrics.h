#ifndef RINGMETRICS_H
#define RINGMETRICS_H

#include <vector>
#include <map>
#include <set>
#include <random>
#include <Network.h>
#include "Gaussian_Ring_Generator.h"

struct TrialRawMetrics {
    vector<double> l1_size;
    vector<double> l2_size;
    vector<double> l1_l2_ratio;

    void dumper(ostream &os) {
        assert(l1_size.size() == l2_size.size());
        assert(l1_size.size() == l1_l2_ratio.size());

        for (unsigned int i = 0; i < l1_size.size(); ++i) {
            os << l1_size[i] << "\t" 
               << l2_size[i] << "\t" 
               << l1_l2_ratio[i] << endl; 
        }
    }
};

struct InterviewRawMetrics {
    vector<double> l1_size;
    vector<double> l2_size;
    vector<double> l3_size;
    vector<double> l1_l2_ratio;
    vector<double> l111_trans; 
    vector<double> l112_trans;
    vector<double> l122_trans;
    vector<double> l1_log_component_to_size_ratio;
    vector<double> mean_path_diameter_ratio;       

    void dumper(ostream &os) {
        assert(l1_size.size() == l2_size.size());
        assert(l1_size.size() == l3_size.size());
        assert(l1_size.size() == l1_l2_ratio.size());
        assert(l1_size.size() == l111_trans.size());
        assert(l1_size.size() == l112_trans.size());
        assert(l1_size.size() == l122_trans.size());
        assert(l1_size.size() == l1_log_component_to_size_ratio.size());
        assert(l1_size.size() == mean_path_diameter_ratio.size());

        for (unsigned int i = 0; i < l1_size.size(); ++i) {
            os << l1_size[i] << "\t" 
               << l2_size[i] << "\t" 
               << l3_size[i] << "\t" 
               << l1_l2_ratio[i] << "\t"
               << l111_trans[i]  << "\t"
               << l112_trans[i] << "\t"
               << l122_trans[i] << "\t"
               << l1_log_component_to_size_ratio[i] << "\t"
               << mean_path_diameter_ratio[i] << endl;
        }
    }
};

struct InterviewProbabilities {
    InterviewProbabilities() : fixed(1.0) {};
    InterviewProbabilities(double f) : fixed(f) {};
    double fixed;
};

enum InterviewState {UNDETECTED, DETECTED, INTERVIEWED, NUM_OF_INTERVIEW_STATES};

Network* interview_network(Network* net, const map<const Node*, int> &level_of, const InterviewProbabilities &ip, const unsigned long int seed) {
    mt19937 rng(seed);
    uniform_real_distribution<double> runif(0.0, 1.0);

    Network* inet = net->duplicate();
    set<Node*, NodePtrComp> interviewed_nodes;
    for (Node* node: inet->get_nodes()) {
        const int lvl = level_of.at(net->get_node(node->get_id()));
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
        inet->get_node(id)->set_state(state);
    }

    vector<Node*> undetected_nodes;
    set_difference(inodes.begin(), inodes.end(), 
                   detected_nodes.begin(), detected_nodes.end(), 
                   inserter(undetected_nodes, undetected_nodes.begin()));
    for (Node* unode: undetected_nodes) inet->delete_node(unode);
    cerr << "interviewed ct: " << interviewed_nodes.size() << endl;
    cerr << "detected ct: " << detected_nodes.size() << endl;
    cerr << "undetected ct: " << undetected_nodes.size() << endl;
    cerr << "Trial size, interviewed size: " << net->size() << ", " << inet->size() << endl;
    return inet;
}

enum TransitivityType {L111, L112, L122, L222, NUM_OF_TRANSITIVITY_TYPES}; // 0, 1, 2, or 3 nodes in level 2

vector<double> special_transitivity (Network* net, map<const Node*, int> level_of) {
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


void determine_level(const DistanceMatrix& dm, map<const Node*, int> &ilevel_of, vector<set<const Node*, NodePtrComp> > &ilevels) {
    for (pair<const Node*, double> nd: dm) {
        const Node* n = nd.first;
        const int level = nd.second;
        ilevel_of[n] = level;
        if ((unsigned) level >= ilevels.size()) ilevels.resize(level+1);
        ilevels[level].insert(n);
    }
}


pair<double, double> calc_mean_and_max(PairwiseDistanceMatrix& pdm) {
    double total = 0;
    int ct = 0;
    double max = 0;
    for (pair<const Node*, DistanceMatrix> pair_ndm: pdm) {
        const Node* n = pair_ndm.first;
        for (pair<const Node*, double> pair_nd: pair_ndm.second) {
            const Node* m = pair_nd.first;
            if (n != m) {  // don't consider distance from nodes to themselves
                const double d =  pair_nd.second;
                total += d;
                max = d > max ? d : max;
                ct++;
            }
        }
    }

    const double mean = total / ct;
    return make_pair(mean, max);
}


void raw_trial_metrics(vector<set<const Node*, NodePtrComp> > levels, TrialRawMetrics& trm) {
    trm.l1_size.push_back(levels[1].size());
    trm.l2_size.push_back(levels[2].size());
    trm.l1_l2_ratio.push_back((double) levels[1].size() / levels[2].size());
}


void raw_interview_metrics(Network* inet, vector<set<const Node*, NodePtrComp> > ilevels, map<const Node*, int> ilevel_of, PairwiseDistanceMatrix& pdm, InterviewRawMetrics& irm) {

    vector<double> s_trans = special_transitivity(inet, ilevel_of);

    // Determine # components (after removing L0) and # of interviewed L1 nodes
    Network* tmp_net = inet->duplicate();
    tmp_net->delete_node(tmp_net->get_node(0)); // delete the index case
    const int comp_ct = tmp_net->get_components().size();
    int l1i_size = 0;
    if (ilevels.size() > 1) {
        for (const Node* n: ilevels[1]) l1i_size += n->get_state() == INTERVIEWED;
    }

    if (pdm.size() == 0) pdm = inet->calculate_distances_map();

    pair<double, double> mean_and_max = calc_mean_and_max(pdm);
    
    irm.l1_size.push_back(ilevels[1].size());
    irm.l2_size.push_back(ilevels[2].size());
    irm.l3_size.push_back(ilevels[3].size());
    irm.l1_l2_ratio.push_back((double) ilevels[1].size() / ilevels[2].size());
    irm.l111_trans.push_back(s_trans[L111]);
    irm.l112_trans.push_back(s_trans[L112]);
    irm.l122_trans.push_back(s_trans[L122]);
    irm.l1_log_component_to_size_ratio.push_back(log((double) comp_ct/l1i_size));
    irm.mean_path_diameter_ratio.push_back(mean_and_max.first / mean_and_max.second);
}


unsigned int get_pzero_idx(Network* net) {
    unsigned int pzero_idx;
    vector<Node*> net_nodes = net->get_nodes();

    for (pzero_idx = 0; pzero_idx < net_nodes.size(); ++pzero_idx) {
        if (net_nodes[pzero_idx]->get_id() == 0) break;
    }
    assert(net_nodes[pzero_idx]->get_id() == 0);

    return pzero_idx;
}


void raw_metrics(Network* net, vector<set<const Node*, NodePtrComp> > levels, map<const Node*, int> level_of, TrialRawMetrics& trm, InterviewRawMetrics& irm, bool do_interview, InterviewProbabilities& ip, const unsigned long int rng_seed) {
    raw_trial_metrics(levels, trm);

    if (do_interview) { 
        Network* inet = interview_network(net, level_of, ip, rng_seed+1);

        PairwiseDistanceMatrix pdm = inet->calculate_distances_map();

        // build level look-up data structures
        map<const Node*, int> ilevel_of;
        vector<set<const Node*, NodePtrComp> > ilevels(4);

        unsigned int pzero_idx = get_pzero_idx(inet);
        vector<Node*> inet_nodes = inet->get_nodes();

        DistanceMatrix dm = pdm[inet_nodes[pzero_idx]];
        determine_level(dm, ilevel_of, ilevels);
        raw_interview_metrics(inet, ilevels, ilevel_of, pdm, irm);
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

#endif
