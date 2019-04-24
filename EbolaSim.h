#ifndef EBOLA_SIM_H
#define EBOLA_SIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <random>
#include <iomanip>
#include <map>
#include <functional>
#include <unistd.h>
#include "IDCommunity.h"
#include "NetworkSimplate.h"
#include "Utility.h"

using namespace std;
using namespace std::placeholders; // gives _1, _2, etc

//int seed = 0;             // to use a static seed
//mt19937 gen(seed);

// uniform_real_distribution<double> runif(0.0, 1.0);
// mt19937 rng;                      // random number generator

enum EventType {
  EXPOSE,
  INCUBATE,
  HOSPITAL,
  RECOVER,
  DIE,
  VACCINATE,
  TRACE,
  N_EVENTS // must be last
};

// want vaccine object to answer question:
// is an individual protected?
// what is potentially needed to answer - when vaccinated, what time now, what kind of vaccine

// struct SimParameters {
//     SimParameters() {
//         network = nullptr;
//         index_case = nullptr;
//         seed = 0;
//         prob_quarantine = 0.0;
//         prob_community_death = 0.0;
//     }
//
//     Network* network;
//     Node* index_case;
map<EventType, function<double(mt19937&)> > event_generator;
//     Vaccine vaccine;
//     unsigned int seed;
//     double prob_quarantine;
//     double prob_community_death;
//     int control_radius = 2;
// };

//template class Event<EventType>;
using EboEvent = Event<EventType>;

class EbolaSim : public EventDrivenSim<EboEvent> {
  public:

    EbolaSim(Network* network) :
    community(IDCommunity(network)),
    index_case(network->get_node(0)),
    disease_log_data(network->size(), vector<double>(N_STATES, std::numeric_limits<double>::quiet_NaN())) {}

    mt19937 localrng; // random number generator
    int day;
    IDCommunity community;
    Node* index_case;
    vector< vector<double> > disease_log_data;
    int control_radius = 2;

    virtual void reset() {
      disease_log_data.clear();
      disease_log_data.resize(community.size(), vector<double>(N_STATES, std::numeric_limits<double>::quiet_NaN()));
      community.reset();
      EventDrivenSim<EboEvent>::reset();
    }

    vector<EboEvent> defaultEvents() {
      return { EboEvent(0, EXPOSE, index_case) };
    }

    map<EventType, function<double(mt19937&)> > event_generator = {
      { EXPOSE,   [](mt19937& rng) { return 5.0; }},
      { INCUBATE, [](mt19937& rng) { return 5.0; }},
      { DIE, [](mt19937& rng) { return 10.0; }},
      { RECOVER, [](mt19937& rng) { return 5.0; }},
      { HOSPITAL, [](mt19937& rng) { return 2.0; }}
    };

    bool canTransmit(Node* source) {
      return !source ? true :
        source->get_state() == INFECTIOUS && community.get_condition(source) != HOSPITALIZED;
    }

    bool exposure(EboEvent event) {
      if (event.node->get_state() == SUSCEPTIBLE && canTransmit(event.source)) {

        auto newEvent = event;
        newEvent.type = INCUBATE;
        newEvent.source = nullptr;
        newEvent.time += event_generator[newEvent.type](localrng);

        add_event(newEvent);

        community.update_state(event.node, EXPOSED);

        return true;
      } else return false;
    }

    bool incubation(EboEvent event) {
      if (event.node->get_state() == EXPOSED) {

        auto newEvent = event;

        double tDeath    = event_generator[DIE](localrng),
               tRecover  = event_generator[RECOVER](localrng),
               tHospital = event_generator[HOSPITAL](localrng);

        newEvent.type = tDeath < tRecover ? DIE : RECOVER;
        newEvent.time += min(tDeath, tRecover);

        add_event(newEvent);

        if (tHospital < min(tDeath, tRecover)) {
          auto hospEvent = event;
          hospEvent.type = HOSPITAL;
          hospEvent.time += tHospital;
          add_event(hospEvent);
        }

        community.update_state(event.node, INFECTIOUS);

        return true;
      } else return false;
    }

    bool notYetTraced = true;

    bool hospital(EboEvent event) {
      if (event.node->get_state() == INFECTIOUS) {
        community.update_condition(event.node, HOSPITALIZED);
        if (notYetTraced) {
          notYetTraced = false;
          // TODO trace events to neighbors
        }
        return true;
      } else return false;
    }

    bool recovery(EboEvent event) {
      if (event.node->get_state() == INFECTIOUS) {
        community.update_state(event.node, RECOVERED);
        return true;
      } else return false;
    }

    bool death(EboEvent event) {
      if (event.node->get_state() == INFECTIOUS) {
        community.update_state(event.node, DEAD);
        return true;
      } else return false;
    }

    bool vaccination(EboEvent event) {
      return false;
    }

    bool tracing(EboEvent event) {
      return false;
    }

    bool eventerror(EboEvent event) {
      cerr << "ERROR: Unsupported event type: " << event.type << endl;
      exit(-2);
      return false;
    }

    map< EventType, function<bool(EboEvent&)> > process_steps = {
      { EXPOSE,    bind(&EbolaSim::exposure,   this, _1) },
      { INCUBATE,  bind(&EbolaSim::incubation, this, _1) },
      { HOSPITAL,  bind(&EbolaSim::hospital,   this, _1) },
      { RECOVER,   bind(&EbolaSim::recovery,   this, _1) },
      { DIE,       bind(&EbolaSim::death,      this, _1) },
      { VACCINATE, bind(&EbolaSim::vaccination,  this, _1) },
      { TRACE,     bind(&EbolaSim::tracing,      this, _1) },
      { N_EVENTS,  bind(&EbolaSim::eventerror, this, _1) },
    };

    map< EventType, DiseaseState > disease_outcome = {
      { EXPOSE,   EXPOSED      },
      { INCUBATE, INFECTIOUS   },
      { HOSPITAL, N_STATES },
      { RECOVER,  RECOVERED },
      { DIE, DEAD },
      { VACCINATE, N_STATES },
      { TRACE, N_STATES },
      { N_EVENTS, N_STATES }
    };

    map< EventType, ControlCondition > control_outcome = {
      { EXPOSE,   N_CONDITIONS },
      { INCUBATE, N_CONDITIONS   },
      { HOSPITAL, HOSPITALIZED },
      { RECOVER,  N_CONDITIONS },
      { DIE, N_CONDITIONS },
      { VACCINATE, VACCINATED },
      { TRACE, TRACED },
      { N_EVENTS, N_CONDITIONS }
    };

    void process(EboEvent event) {
      bool success = process_steps[event.type](event);
      DiseaseState updateTo = disease_outcome[event.type];
      //DiseaseState updateTo = disease_outcome[event.type];
      if (success) logupdate(event, updateTo);
      if (updateTo == EXPOSED) {
        // update infector log
      }
      // don't need to update TRACED log? since that's already a thing that can be printed from
    }

    void logupdate(EboEvent event, DiseaseState ds) {
      if (ds != N_STATES) {
        cout << "DiseaseState: " << ds << endl;
        cout << "Event: " << event.type << endl;
        cout << "Target Node: " << event.node->get_id() << endl;
        disease_log_data[event.node->get_id()][ds] = event.time;
      }
    }

};
#endif
