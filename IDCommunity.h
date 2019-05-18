#ifndef IDCOMMUNITY_H
#define IDCOMMUNITY_H

#include <vector>
#include <iostream>
#include "Network.h"

using namespace std;

// the mutually exclusive disease states
enum DiseaseState {
  SUSCEPTIBLE, // expected to be first for data log
  EXPOSED,
  INFECTIOUS,
  RECOVERED,
//  DEAD,
  N_STATES
}; // must be last

std::ostream& operator<<(std::ostream &out, DiseaseState &ds) {
  string res;
  switch (ds) {
    case SUSCEPTIBLE: res = "SUSCEPTIBLE"; break;
    case EXPOSED    : res = "EXPOSED"; break;
    case INFECTIOUS : res = "INFECTIOUS"; break;
    case RECOVERED  : res = "RECOVERED"; break;
//    case DEAD       : res = "DEAD"; break;
    case N_STATES   : res = "N_STATES"; break;
  }
  return out << res;
}

enum ControlCondition {
  NONE,
  TRACED,
  VACCINATED,
  HOSPITALIZED,
  N_CONDITIONS
};

std::ostream& operator<<(std::ostream &out, ControlCondition &cc) {
  string res;
  switch (cc) {
    case NONE         : res = "NONE"; break;
    case TRACED       : res = "TRACED"; break;
    case VACCINATED   : res = "VACCINATED"; break;
    case HOSPITALIZED : res = "HOSPITALIZED"; break;
    case N_CONDITIONS : res = "N_CONDITIONS"; break;
  }
  return out << res;
}

// for anything that should be counts
typedef unsigned int capita;

class IDCommunity {

  static const capita zero = 0;

  public:
    IDCommunity(Network* n) :
    network(n),
    state_counts(N_STATES, zero),
//    control_counts(N_CONDITIONS, 0),
    traced(n->size(), -1),
    quarantined(n->size(), false),
    reactive_vaccine(n->size(), -1.0),
    prophylactic_vaccine(n->size(), -1.0)
    {
        reset();
    }

    Network* network;           // population
    vector<capita> state_counts;   // S, E, I, R, etc. counts
    // vector<capita> control_counts;
    // accounting for lack of multi-state on Node
    bool anytrace = false;
    bool isTraced() { return anytrace; }

    vector<int> traced;
    int level(Node* n) { return level(n->get_id()); }
    int level(int id) { return traced[id]; }


    vector<bool> quarantined;

    vector<double> reactive_vaccine;
    vector<double> prophylactic_vaccine;

    capita size() { return network->size(); }

    capita current_epidemic_size() {
      return state_counts[EXPOSED] + state_counts[INFECTIOUS];
    }

    capita final_size() {
      return current_epidemic_size() + state_counts[HOSPITALIZED] + state_counts[RECOVERED];// + state_counts[DEAD];
    }

    template<typename T>
    void reset(vector<T>& tar, T def) {
      int sz = tar.size();
      tar.clear(); tar.resize(sz, def);
    }

    void reset() {
      for (auto node: network->get_nodes()) node->set_state(SUSCEPTIBLE);

      // control_counts.clear();
      // control_counts.resize(N_CONDITIONS, 0);
      // control_counts[NONE] = network->size();

      reset(traced, -1);
      reset(quarantined, false);
      reset(reactive_vaccine, -1.0);
      reset(prophylactic_vaccine, -1.0);
      reset(state_counts, zero);
      state_counts[SUSCEPTIBLE] = network->size();

    }

    void update_state(int id, DiseaseState new_state) { update_state(network->get_node(id), new_state); }

    void update_state(Node* node, DiseaseState new_state) {
      auto old_state = static_cast<DiseaseState>(node->get_state());
      node->set_state(new_state);
      state_counts[old_state]--; state_counts[new_state]++;
    }

    static const double found_cost;
    static const double missed_cost;
    static const double unknown_cost;

    // found = 0
    // missed = Inf
    // un-assessed
    void trace(Edge* e, bool success = true) {
      anytrace = true;
      assert(e->get_cost() == unknown_cost);
      e->set_cost(success ? found_cost : missed_cost);
      if (success) {
        cout << "tracing edge " << e->get_id() << "..." << endl;
        // ? traced[e->get_end()->get_id()] = true;
        e->get_complement()->set_cost(found_cost);
      }
    }

    void quarantine(int n) { quarantined[n] = true; }
    void quarantine(Node* n) { quarantine(n->get_id()); }

    bool isQuarantined(int n) { return quarantined[n]; }
    bool isQuarantined(Node* n) { return isQuarantined(n->get_id()); }

    // void update_condition(Node* node, ControlCondition new_state) {
    //   auto old_state = control_condition[node->get_id()];
    //   control_condition[node->get_id()] = new_state;
    //   control_counts[old_state]--; control_counts[new_state]++;
    // }
    //
    // ControlCondition get_condition(Node* node) {
    //   return control_condition[node->get_id()];
    // }

    static const string header;// = "(S | EI | R )";

};

const string IDCommunity::header = "S\t|\tE\tI\t|\tR";
const double IDCommunity::found_cost = 0.0;
const double IDCommunity::missed_cost = std::numeric_limits<double>::infinity();
const double IDCommunity::unknown_cost = 1.0;

std::ostream& operator<<(std::ostream &out, IDCommunity &cc) {
  auto sc = cc.state_counts;
  return out
    << sc[SUSCEPTIBLE] <<
    "\t|\t"
    << sc[EXPOSED] << "\t" << sc[INFECTIOUS] <<
    "\t|\t"
    << sc[RECOVERED];
}

#endif
