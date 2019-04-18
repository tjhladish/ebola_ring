#ifndef IDCOMMUNITY_H
#define IDCOMMUNITY_H

#include <vector>
#include <iostream>
#include "Network.h"

using namespace std;

// for anything that should be counts
typedef unsigned int capita;

// the mutually exclusive disease states
enum DiseaseState {
  SUSCEPTIBLE, // expected to be first for data log
  EXPOSED,
  INFECTIOUS,
  HOSPITALIZED,
  RECOVERED,
  DEAD,
  NUM_OF_STATE_TYPES }; // must be last

class IDCommunity {

  public:
    IDCommunity(Network* n) : network(n), state_counts(NUM_OF_STATE_TYPES, 0) { reset(); }

    Network* network;           // population
    vector<capita> state_counts;   // S, E, I, R, etc. counts

    capita current_epidemic_size() {
      return state_counts[EXPOSED] + state_counts[INFECTIOUS];
    }

    capita final_size() {
      return state_counts[EXPOSED] + state_counts[INFECTIOUS] + state_counts[HOSPITALIZED] + state_counts[RECOVERED] + state_counts[DEAD];
    }

    void reset() {
      for (auto node: network->get_nodes()) node->set_state(SUSCEPTIBLE);
      state_counts.clear();
      state_counts.resize(NUM_OF_STATE_TYPES, 0);
      state_counts[SUSCEPTIBLE] = network->size();
    }

    void update_state(Node* node, DiseaseState new_state) {
      auto old_state = static_cast<DiseaseState>(node->get_state());
      node->set_state(new_state);
      state_counts[old_state]--; state_counts[new_state]++;
    }

    void log(ostream o, double Now) {
      o << "(SvV | EI | HRD) " << Now << " :\t"
        << state_counts[SUSCEPTIBLE] << "\t"
        << state_counts[VAC_PARTIAL] << "\t"
        << state_counts[VAC_FULL] << "\t|\t"
        << state_counts[EXPOSED] << "\t"
        << state_counts[INFECTIOUS] << "\t|\t"
        << state_counts[HOSPITALIZED] << "\t"
        << state_counts[RECOVERED] << "\t"
        << state_counts[DEAD] << endl;
    }

};

#endif
