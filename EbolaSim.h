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
//     map<EventType, function<double(mt19937&)> > event_generator;
//     Vaccine vaccine;
//     unsigned int seed;
//     double prob_quarantine;
//     double prob_community_death;
//     int control_radius = 2;
// };

//template class Event<EventType>;
using EboEvent = Event<EventType>;

struct SimPars {
  Network* net;
  map<EventType, function<double(mt19937&)>> event_time_distribution;
  int simseed;
};

class EbolaSim : public EventDrivenSim<EboEvent> {
  public:

    EbolaSim(SimPars& pars) :
    // assorted simple constructions
      community(IDCommunity(pars.net)),
      index_case(pars.net->get_node(0)),
      disease_log_data(
        pars.net->size(),
        vector<double>(N_STATES, numeric_limits<double>::quiet_NaN())
      ),
      rngseed(pars.simseed)
    // construction with a bit more complexity
    {
        // concept: input parameters provides the distributions
        // but internally, the simulation code only invokes the generators via ()
        // with all generators locked to the same rng
        for_each(
          pars.event_time_distribution.begin(), pars.event_time_distribution.end(),
          [&](pair<const EventType, function<double(mt19937&)>>& p) {
            event_time_generator[p.first] = bind(p.second, std::ref(localrng));
          }
        );
        derivedReset();
    }

    int rngseed;
    mt19937 localrng; // random number generator
    map<EventType, function<double()>> event_time_generator;
    double runif() { return uniform_real_distribution<double>(0,1)(localrng); }

    int day;
    IDCommunity community;
    Node* index_case;
    vector< vector<double> > disease_log_data;
    int control_radius = 2;

    virtual void reset() {
      localrng.seed(rngseed);
      disease_log_data.clear();
      disease_log_data.resize(community.size(), vector<double>(N_STATES, std::numeric_limits<double>::quiet_NaN()));
      community.reset();
      EventDrivenSim<EboEvent>::reset(); // also reset eventQ
    }

    vector<EboEvent> defaultEvents() {
      return { EboEvent(0, EXPOSE, index_case) };
    }

    bool canTransmit(Node* source) {
      return !source ? true : // if source is null, means external introduction
        // otherwise, can transmit if source is infectious & not hospitalized
        source->get_state() == INFECTIOUS && community.get_condition(source) != HOSPITALIZED;
    }

    bool exposure(EboEvent event) {
      if (event.node->get_state() == SUSCEPTIBLE && canTransmit(event.source)) {

        auto newEvent = event; // copy old event
        newEvent.type = INCUBATE; // update event type
        newEvent.source = nullptr; // remove source
        newEvent.time += event_time_generator[newEvent.type](); // draw new time

        add_event(newEvent); // add the event to the queue

        community.update_state(event.node, EXPOSED); // update community accounting

        return true; // success
      } else return false; // no exposure
    }

    bool incubation(EboEvent event) {
      if (event.node->get_state() == EXPOSED) {

        auto newEvent = event; // definitely having at least one new event, so copy previous

        // draw times for possible outcomes
        double tDeath    = event_time_generator[DIE](),
               tRecover  = event_time_generator[RECOVER](),
               tHospital = event_time_generator[HOSPITAL]();
        // which disease outcome?
        newEvent.type = tDeath < tRecover ? DIE : RECOVER;
        newEvent.time += min(tDeath, tRecover);

        add_event(newEvent);

        // hospitalized first?
        if (tHospital < min(tDeath, tRecover)) {
          auto hospEvent = event;
          hospEvent.type = HOSPITAL;
          hospEvent.time += tHospital;
          add_event(hospEvent);
        }

        community.update_state(event.node, INFECTIOUS);

        double tInfectiousInCommunity = min(tHospital, tDeath, tRecover);
        auto refExpEvent = event;
        refExpEvent.source = event.node;
        refExpEvent.type = EXPOSE;

        for (auto neighbor : event.node->get_neighbors()) {
          double expTime = event_time_generator[EXPOSE]();
          if (expTime < tInfectiousInCommunity) {
            auto expEvent = refExpEvent;
            expEvent.node = neighbor;
            expEvent.time += expTime;
            add_event(expEvent);
          }
        }

        return true;
      } else return false;
    }

    bool notYetTraced = true;

    bool hospital(EboEvent event) {
      if (event.node->get_state() == INFECTIOUS) {
        if (notYetTraced) {
          notYetTraced = false;
          auto refTraceEvent = event;
          refTraceEvent.type = TRACE;
          refTraceEvent.source = event.node;
          for (auto neighbor : event.node->get_neighbors()) {
            // TODO trace events to neighbors - at this point, no one traced, so only need to check obs prob
            // but still need to record edges as checked vs not
          }
        }
        community.update_condition(event.node, HOSPITALIZED);
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
        // TODO if tracing hasn't happened yet, draw time-to-detection event
        community.update_state(event.node, DEAD);
        return true;
      } else return false;
    }

    bool vaccination(EboEvent event) {
      // TODO vaccination tracking
      return false;
    }

    vector<int> level;
    map<Node*,set<Node*>> traces;

    bool tracing(EboEvent event) {
      int myLevel = level[event.source->get_id()]+1;

      if (myLevel <= control_radius) {
        if (level[event.node->get_id()] < 0) { // previously undetected node
          myLevel = event.node->get_state() == INFECTIOUS ? 0 : myLevel;
          if (myLevel and runif() < pCoverage) {
            auto vacEvent = event;
            vacEvent.type = VACCINATE;
            vacEvent.time += event_time_generator[VACCINATE]();
            add_event(vacEvent);
          }
          set<Node*> knownContacts = *(traces.emplace(event.node, set<Node*>()).first);
          for (auto neighbor : event.node->get_neighbors()) if (runif() < pTrace) knownContacts.add(neighbor);
          auto refTrace = event;
          refTrace.source = event.node;
          for (auto contact : knownContacts) {
            auto traceEvent = refTrace;
            traceEvent.node = contact;
            traceEvent.time += event_time_generator[TRACE]();
            // add event for everyone; alternative is get lowest of levels first
            // then only add trace for relevant ones
            add_event(traceEvent);
            if (level[contact->get_id] >= 0 and level[contact->get_id]+1 < myLevel) myLevel = level[contact->get_id]+1;
          }
        } else if (myLevel < level[event.node->get_id()]) { // need to lower my level, resend
          set<Node*> knownContacts = traces[event.node]; // contacts should already be known
          for (auto contact : knownContacts) {
            auto traceEvent = refTrace;
            traceEvent.node = contact;
            traceEvent.time += event_time_generator[TRACE]();
            // add event for everyone; alternative is get lowest of levels first
            // then only add trace for relevant ones
            add_event(traceEvent);
            if (level[contact->get_id()] >= 0 and level[contact->get_id()]+1 < myLevel) myLevel = level[contact->get_id()]+1;
          }
        }
        level[event.node->get_id()] = myLevel;
        return true;
      } else return false;
    }

    bool eventerror(EboEvent event) {
      cerr << "ERROR: Unsupported event type: " << event.type << endl;
      exit(-2);
      return false;
    }

    // TODO re-write next in terms of this + transform
    // map<EventType, decltype(&EbolaSim::exposure)> links {
    //   { EXPOSE, &EbolaSim::exposure },
    //   { EXPOSE, &EbolaSim::exposure }
    // };

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
