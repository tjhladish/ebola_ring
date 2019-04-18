#ifndef ED_EBOLA_SIM_H
#define ED_EBOLA_SIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include "Utility.h"
#include "Network.h"
#include <random>
#include <iomanip>
#include <map>
#include <functional>
#include <unistd.h>

using namespace std;

//int seed = 0;             // to use a static seed
//mt19937 gen(seed);

uniform_real_distribution<double> runif(0.0, 1.0);
mt19937 rng;                      // random number generator

enum StateType { SUSCEPTIBLE,          // expected to be first for data log
                 VAC_PARTIAL,
                 VAC_FULL,
                 EXPOSED,
                 INFECTIOUS,
                 HOSPITALIZED,
                 RECOVERED,
                 DEAD,
                 NUM_OF_STATE_TYPES }; // must be last

enum EventType { StoE_EVENT,
                 EtoI_EVENT,
                 ItoR_EVENT,
                 ItoH_EVENT,
                 ItoH_PZERO_EVENT,
                 ItoD_EVENT,
                 V1_EVENT,
                 V2_EVENT,
                 NUM_OF_EVENT_TYPES }; // must be last


class Event {
  public:
    double time;
    EventType type;
    Node* node;
    Node* source; // for exposure events
    Event(const Event& o) {  time=o.time; type=o.type; node=o.node; source=o.source; }
    Event(double t, EventType e, Node* n, Node* s = nullptr) { time=t; type=e; node=n; source=s; }
    Event& operator=(const Event& o) { time=o.time; type=o.type; node=o.node; source=o.source; return *this; }
};

class compTime {
  public:
    bool operator() (const Event* lhs, const Event* rhs) const {
        return (lhs->time>rhs->time);
    }

    bool operator() (const Event& lhs, const Event& rhs) const {
        return (lhs.time>rhs.time);
    }
};

class Vaccine {
  public:
    Vaccine () : efficacy({1.0,1.0}), coverage({1.0,1.0}), isLeaky(true), timeToProtection(numeric_limits<double>::max()) {}
    vector<double> efficacy; // doses 1 and 2
    vector<double> coverage; // same
    bool isLeaky;
    double timeToProtection;
};

struct SimParameters {
    SimParameters() {
        network = nullptr;
        index_case = nullptr;
        seed = 0;
        prob_quarantine = 0.0;
        prob_community_death = 0.0;
    }

    Network* network;
    Node* index_case;
    map<EventType, function<double(mt19937&)> > event_generator;
    Vaccine vaccine;
    unsigned int seed;
    double prob_quarantine;
    double prob_community_death;
    };


class Event_Driven_Ebola_Sim {
  public:
                                // constructor
    Event_Driven_Ebola_Sim (SimParameters& par) :
      vaccine(par.vaccine),
      log_data(par.network->size(), vector<double>(NUM_OF_STATE_TYPES, -1))
    {
        rng.seed(par.seed);
        community = IDCommunity(par.network);
        event_generator = par.event_generator;
        prob_quarantine = par.prob_quarantine;
        prob_community_death = par.prob_community_death;
        // node_vac_immunity.clear();
        // for (Node* n: network->get_nodes()) node_vac_immunity[n->get_id()] = make_pair(0.0, 0.0);
        // vaccine_doses_used = {0, 0}; // doses 1 and 2
        reset();
    }

    int day;
    vector< vector<double> > log_data;
    IDCommunity community;
    priority_queue<Event, vector<Event>, compTime > EventQ;

    void reset() {
        day = 0;
        EventQ = priority_queue<Event, vector<Event>, compTime >();
        for (auto el: log_data) { el.assign(NUM_OF_STATE_TYPES, -1); }
        community.reset();
    }

    // TODO implement vaccination model(s)
    // rough idea: vaccine object(s) that hold universal information about vaccine
    // and nodes hold state information about their particular vaccine events (e.g., when, what kind(s), assignment for non-leaky vax)
    bool protected(Node* node) {
      return false;
    }

    void ascertain(Node* node) {
      // TODO: mark node as ascertained
      // if currently infectious, move immediately to hospitalized
      // and prune future infection events from queue
    }

    bool expose(Node* node, const int source_id, double eventtime) {
        //cerr << "node: " << node->get_id() << endl;
        auto state = static_cast<StateType>(node->get_state());
        if (state == NUM_OF_STATE_TYPES) {
          cerr << "ERROR: Exposure of node in unsupported state.  Node state is " << state << endl;
          exit(-1);
        }

        bool didTransmissionOccur = (static_cast<DiseaseState>(node->get_state())==SUSCEPTIBLE) and (not protected(node));

        // add infectiousness event
        if (didTransmissionOccur) {
          add_event(eventtime + time_to_event(EtoI_EVENT), EtoI_EVENT, node);
          community.update_state(node, EXPOSED);
        }

        return didTransmissionOccur;
    }

    void onset(Node* node) {
      // is the node ascertained yet
      // if yes - distribute recovery or hospital death (even if death would happen sooner than hospitalization) + hospitalization (if it occurs sooner than death / recovery)
      // if not, will the node be ascertained without any intervention?
      // draw time to death, time to recovery, time to ascertainment
      // if death < recovery => death instead of recovery
      // if ascertainment < min(death, recovery) => hosp first, then next thing
      // whatever min disease state event:
      // draw n exposure times (n == num neighbors)
      // for the ones less than min disease state event, add them to eventQ
    }

    // little effect on dynamics - book-keeping
    void recover(Node* node) {
      community.update_state(node, RECOVERED);
    }

    // potential to trigger vaccination
    void hospitalize(Node* node) {
      // if already ascertained
      community.update_state(node, HOSPITALIZED);
    }

    // if (state == SUSCEPTIBLE or state == VAC_PARTIAL or state == VAC_FULL) {
    //     if (unprotected_by_vaccine) {
    //         double Ti = Now + time_to_event(EtoI_EVENT);      // time to become infectious
    //         add_event(Ti, EtoI_EVENT, node);
    //
    //         // TODO - use probabilities of quarantine, community death
    //         // TODO - figure out how to handle the special case of patient zero
    //         double Th = Ti;
    //         Th += isPatientZero(node) ? time_to_event(ItoH_PZERO_EVENT) : time_to_event(ItoH_EVENT);
    //         double Tr = Ti + time_to_event(ItoR_EVENT);
    //         double Td = Ti + time_to_event(ItoD_EVENT);
    //
    //         double Ti_end;
    //         if (Th < Tr and Th < Td) {                        // end of infectious period is whichever happens first
    //             Ti_end = Th;
    //             add_event(Th, ItoH_EVENT, node);
    //         } else if (Tr < Th and Tr < Td) {
    //             Ti_end = Tr;
    //             add_event(Tr, ItoR_EVENT, node);
    //         } else {
    //             Ti_end = Td;
    //             add_event(Td, ItoD_EVENT, node);
    //         }
    //         for (Node* neighbor: node->get_neighbors()) {      // density-dependent assumption! more neighbors --> more contact per unit time
    //             double Tc = Ti + time_to_event(StoE_EVENT);    // time to next contact--we'll worry about whether contact is susceptible at Tc
    //             while ( Tc < Ti_end ) {                        // does contact occur before recovery?
    //                 add_event(Tc, StoE_EVENT, neighbor, node); // potential transmission event
    //                 Tc += time_to_event(StoE_EVENT);           // if vaccine is leaky, or successful transmission is probabilistic, this is impt
    //             }
    //         }
    //         didTransmissionOccur = true;
    //         update_state(node, state, EXPOSED);
    //     }
    // }

    // void schedule_vaccinations() {
    //     cerr << "vaccination triggered\n";
    //     // schedule initial vaccinations for everyone in network who is eligible
    //     // first and second doses are scheduled if efficacy is > 0 (dose 2 is conditional on dose 1)
    //     if (vaccine.efficacy[0] > 0.0) {
    //         const double Tv1 = Now + time_to_event(V1_EVENT);
    //         add_event(Tv1, V1_EVENT);
    //         if (vaccine.efficacy[1] > 0.0) {
    //             const double Tv2 = Tv1 + time_to_event(V2_EVENT);
    //             add_event(Tv2, V2_EVENT);
    //         }
    //     }
    // }
    //
    // void vaccine_campaign(EventType et) {
    //     StateType eligible_group;
    //     StateType converts_to;
    //     if (et == V1_EVENT or et == V2_EVENT) {
    //         int dose_idx = 0;
    //         if (et == V1_EVENT) {
    //             assert(vaccine_doses_used[0] == 0); // V1_EVENT is only expected once
    //             eligible_group = SUSCEPTIBLE;
    //             converts_to    = VAC_PARTIAL;
    //         } else {
    //             assert(vaccine_doses_used[1] == 0); // V2_EVENT is only expected once
    //             eligible_group = VAC_PARTIAL;
    //             converts_to    = VAC_FULL;
    //             dose_idx = 1;
    //         }
    //
    //         for (Node* node: network->get_nodes()) {
    //             const StateType state = (StateType) node->get_state();
    //             if (state == eligible_group and runif(rng) < vaccine.coverage[dose_idx]) {
    //                 vaccine_doses_used[dose_idx]++;
    //                 log_data[node->get_id()][converts_to] = Now;
    //                 update_state(node, eligible_group, converts_to);
    //                 const double vac_eff = vaccine.efficacy[dose_idx];
    //                 if (vaccine.isLeaky) {
    //                     node_vac_immunity[node->get_id()] = make_pair(Now, vac_eff);
    //                 } else {
    //                     const double r = runif(rng);
    //                     node_vac_immunity[node->get_id()] = make_pair(Now, r < vac_eff ? 1.0 : 0.0);
    //                 }
    //             }
    //         }
    //     } else {
    //         cerr << "ERROR: Unsupported vaccination event type: " << et << endl;
    //         exit(-4);
    //     }
    // }

    vector< vector<double> > run_simulation(
      vector<Event> evts, // initial list of events - e.g., background vaccination, incipient infection
      const double timelimit = numeric_limits<double>::max()
    ) {
        for (auto e : evts) add_event(e);
        while (next_event(timelimit) < timelimit) continue;
        return log_data;
    }

    double next_event(const double maxtime) {
        if ( EventQ.empty() ) {
            return maxtime;
        } else {
            auto event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            switch(event.type) {
                case EtoI_EVENT: onset(event.node);                          break;
                case ItoR_EVENT: recover(event.node);                        break;
                case ItoH_EVENT: hospitalize(event.node);                    break;
                case HtoD_EVENT: hospitaldeath(event.node);                  break;
                case ItoD_EVENT: communitydeath(event.node);                 break;
                case StoE_EVENT: expose(event.node, event.source->get_id(), event.time); break;
                case V1_EVENT:
                case V2_EVENT:   vaccine_campaign(event.type);               break;
                default:
                    cerr << "ERROR: Unsupported event type: " << event.type << endl;
                    exit(-2);
                    break;          // superfluous after exit(), but included in case refactoring makes it necessary
            }
            log(event, event.time);
            return event.time;
        }
    }

    // template these?
    void add_event( Event et ) { EventQ.push(et); }
    void add_event( double time, EventType et, Node* node = nullptr, Node* source = nullptr) {
        add_event(Event(time, et, node, source));
    }

    void verbose(double Now) {
      if (static_cast<int>(Now) > day) {
        community.log(cerr, Now);
        day = static_cast<int>(Now);
      }
    }

    void log(EventType et, double Now) {
      verbose(Now);
      auto tar = log_data[et.node->get_id()];
      DiseaseState which = NUM_OF_STATE_TYPES;
      switch(et.type) {
          case EtoI_EVENT: which = INFECTIOUS; break;
          case ItoR_EVENT: which = RECOVERED; break;
          case ItoH_EVENT: which = HOSPITALIZED; break;
          case ItoD_EVENT: // intentional fall through
          case HtoD_EVENT: which = DEAD; break;
          case StoE_EVENT: which = EXPOSED;
            tar[0] = event.source->get_id();
            break;
          case V1_EVENT: // intentional fall-through to V2_EVENT
          case V2_EVENT: break;
          default:
              cerr << "ERROR: Unsupported event type: " << event.type << endl;
              exit(-2);
              break;          // superfluous after exit(), but included in case refactoring makes it necessary
      }
      if (which != NUM_OF_STATE_TYPES) tar[which] = Now;
    }

};
#endif
