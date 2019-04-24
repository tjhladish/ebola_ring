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
  RECOVERED,
  DEAD,
  N_STATES
}; // must be last

enum ControlCondition {
  NONE,
  TRACED,
  VACCINATED,
  HOSPITALIZED,
  N_CONDITIONS
};

class IDCommunity {

  public:
    IDCommunity(Network* n) :
    network(n), state_counts(N_STATES, 0), control_counts(N_CONDITIONS, 0),
    control_condition(n->size(), NONE) {
        reset();
    }

    Network* network;           // population
    vector<capita> state_counts;   // S, E, I, R, etc. counts
    vector<capita> control_counts;
    // accounting for lack of multi-state on Node
    vector<ControlCondition> control_condition;

    capita size() { return network->size(); }

    capita current_epidemic_size() {
      return state_counts[EXPOSED] + state_counts[INFECTIOUS];
    }

    capita final_size() {
      return state_counts[EXPOSED] + state_counts[INFECTIOUS] + state_counts[HOSPITALIZED] + state_counts[RECOVERED] + state_counts[DEAD];
    }

    void reset() {
      for (auto node: network->get_nodes()) node->set_state(SUSCEPTIBLE);
      state_counts.clear();
      state_counts.resize(N_STATES, 0);
      state_counts[SUSCEPTIBLE] = network->size();

      control_counts.clear();
      control_counts.resize(N_CONDITIONS, 0);
      control_counts[NONE] = network->size();

      control_condition.clear();
      control_condition.resize(network->size(), NONE);
    }

    void update_state(Node* node, DiseaseState new_state) {
      auto old_state = static_cast<DiseaseState>(node->get_state());
      node->set_state(new_state);
      state_counts[old_state]--; state_counts[new_state]++;
    }

    void update_condition(Node* node, ControlCondition new_state) {
      auto old_state = control_condition[node->get_id()];
      control_condition[node->get_id()] = new_state;
      control_counts[old_state]--; control_counts[new_state]++;
    }

    ControlCondition get_condition(Node* node) {
      return control_condition[node->get_id()];
    }

    void log(ostream o, double Now) {
      o << "(S | EI | RD | UTVH ) " << Now << " :\t"
        << state_counts[SUSCEPTIBLE] << "\t"
        << state_counts[EXPOSED] << "\t"
        << state_counts[INFECTIOUS] << "\t|\t"
        << state_counts[RECOVERED] << "\t"
        << state_counts[DEAD] << "\t|\t"
        << control_counts[NONE] << "\t"
        << control_counts[TRACED] << "\t"
        << control_counts[VACCINATED] << "\t"
        << control_counts[HOSPITALIZED] << endl;
    }

};

#endif
