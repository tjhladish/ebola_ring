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
#include "IDCommunity.h"

using namespace std;

//int seed = 0;             // to use a static seed
//mt19937 gen(seed);

uniform_real_distribution<double> runif(0.0, 1.0);
mt19937 rng;                      // random number generator

enum EventType {
  StoE_EVENT,
  EtoI_EVENT,
  ItoR_EVENT,
  ItoH_EVENT,
  ItoD_EVENT,
  V_EVENT,
  TRACE_EVENT,
  NUM_OF_EVENT_TYPES // must be last
};

// want vaccine object to answer question:
// is an individual protected?
// what is potentially needed to answer - when vaccinated, what time now, what kind of vaccine

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
    int control_radius = 2;
};

class Event_Driven_Ebola_Sim : public EventDrivenSim<Event<EventType>> {
  public:
    
    
    int day;
    vector< vector<double> > log_data;
    IDCommunity community;
    int control_radius;

    bool expose(Node* node, const int source_id) {
        //cerr << "node: " << node->get_id() << endl;
        auto state = static_cast<StateType>(node->get_state());
        if (state == NUM_OF_STATE_TYPES) {
          cerr << "ERROR: Exposure of node in unsupported state.  Node state is " << state << endl;
          exit(-1);
        }

    enum TraceStatus {
      FOUND, MISSED, NUM_TRACE_STATUS;
    }

    vector<vector<TraceStatus>> edgetrace;
    vector<int> nodetrace;

    bool trace(Node* node, Node* src, double time) {
        if (nodetrace[node->get_id()] != -1) {
          int lvl = ;
          // attempt to trace neighbors for new finds
          // if untraced, create trace events
          // if already traced, check level difference
          // TODO refactor to use edges from EpiFire
          for (auto tar : node->get_neighbors()) 
            edgetrace[node->get_id()][tar->get_id()] = pTrace > runif(1) ? FOUND : MISSED;
          // of the found edges, what is the min lvl?
          
          int ref = 2; // TODO use algo stuff to filter edgetrace[node->get_id()] by FOUND, then get min
          int lvl = min(nodetrace[src->get_id()]+1, ref);
          
          
          
          
        } else {
          
        }
      // check node as ascertained - if already ascertained at equal or lower level that source level + 1,
      // don't need to re-notify other neighbors.  *Might* to need to notify source. I.e., if L2 node connects to an L0 node
      // otherwise, proceed
      // if currently infectious, new hospitalized event
      // test neighbors (or reuse previous identifications)
      // register contact tracing events for found neighbors (ignoring source neighbor)
    }


        // TODO actually implement this
        bool unprotected_by_vaccine = not vaccine.protected(node->get_vaccination());
        // if (state == VAC_PARTIAL or state == VAC_FULL) {     // determine if node is protected by vaccine
        //     const double vac_time = node_vac_immunity[node->get_id()].first;
        //                                                      // current assumption is vaccine abruptly assumes full efficacy after timeToProtection
        //     const double vac_protection = (Now - vac_time > vaccine.timeToProtection) ? node_vac_immunity[node->get_id()].second : 0.0;
        //     unprotected_by_vaccine = runif(rng) > vac_protection;
        // }

        bool didTransmissionOccur = unprotected_by_vaccine and state == SUSCEPTIBLE;


        if (didTransmissionOccur) {

        }

        if (state == SUSCEPTIBLE or state == VAC_PARTIAL or state == VAC_FULL) {
            if (unprotected_by_vaccine) {
                double Ti = Now + time_to_event(EtoI_EVENT);      // time to become infectious
                add_event(Ti, EtoI_EVENT, node);

                // TODO - use probabilities of quarantine, community death
                // TODO - figure out how to handle the special case of patient zero
                double Th = Ti;
                Th += isPatientZero(node) ? time_to_event(ItoH_PZERO_EVENT) : time_to_event(ItoH_EVENT);
                double Tr = Ti + time_to_event(ItoR_EVENT);
                double Td = Ti + time_to_event(ItoD_EVENT);

                double Ti_end;
                if (Th < Tr and Th < Td) {                        // end of infectious period is whichever happens first
                    Ti_end = Th;
                    add_event(Th, ItoH_EVENT, node);
                } else if (Tr < Th and Tr < Td) {
                    Ti_end = Tr;
                    add_event(Tr, ItoR_EVENT, node);
                } else {
                    Ti_end = Td;
                    add_event(Td, ItoD_EVENT, node);
                }
                for (Node* neighbor: node->get_neighbors()) {      // density-dependent assumption! more neighbors --> more contact per unit time
                    double Tc = Ti + time_to_event(StoE_EVENT);    // time to next contact--we'll worry about whether contact is susceptible at Tc
                    while ( Tc < Ti_end ) {                        // does contact occur before recovery?
                        add_event(Tc, StoE_EVENT, neighbor, node); // potential transmission event
                        Tc += time_to_event(StoE_EVENT);           // if vaccine is leaky, or successful transmission is probabilistic, this is impt
                    }
                }
                log_data[node->get_id()][EXPOSED] = Now;
                log_data[node->get_id()][0] = source_id;
                update_state(node, state, EXPOSED);
            }
        } else if (state == EXPOSED or state == INFECTIOUS or state == RECOVERED or state == HOSPITALIZED or state == DEAD) {
            didTransmissionOccur = false;
        } else {
            cerr << "ERROR: Exposure of node in unsupported state.  Node state is " << state << endl;
            exit(-1);
        }
        return didTransmissionOccur;
    }

    void schedule_vaccinations() {
        cerr << "vaccination triggered\n";
        // schedule initial vaccinations for everyone in network who is eligible
        // first and second doses are scheduled if efficacy is > 0 (dose 2 is conditional on dose 1)
        if (vaccine.efficacy[0] > 0.0) {
            const double Tv1 = Now + time_to_event(V1_EVENT);
            add_event(Tv1, V1_EVENT);
            if (vaccine.efficacy[1] > 0.0) {
                const double Tv2 = Tv1 + time_to_event(V2_EVENT);
                add_event(Tv2, V2_EVENT);
            }
        }
    }

    void vaccine_campaign(EventType et) {
        StateType eligible_group;
        StateType converts_to;
        if (et == V1_EVENT or et == V2_EVENT) {
            int dose_idx = 0;
            if (et == V1_EVENT) {
                assert(vaccine_doses_used[0] == 0); // V1_EVENT is only expected once
                eligible_group = SUSCEPTIBLE;
                converts_to    = VAC_PARTIAL;
            } else {
                assert(vaccine_doses_used[1] == 0); // V2_EVENT is only expected once
                eligible_group = VAC_PARTIAL;
                converts_to    = VAC_FULL;
                dose_idx = 1;
            }

            for (Node* node: network->get_nodes()) {
                const StateType state = (StateType) node->get_state();
                if (state == eligible_group and runif(rng) < vaccine.coverage[dose_idx]) {
                    vaccine_doses_used[dose_idx]++;
                    log_data[node->get_id()][converts_to] = Now;
                    update_state(node, eligible_group, converts_to);
                    const double vac_eff = vaccine.efficacy[dose_idx];
                    if (vaccine.isLeaky) {
                        node_vac_immunity[node->get_id()] = make_pair(Now, vac_eff);
                    } else {
                        const double r = runif(rng);
                        node_vac_immunity[node->get_id()] = make_pair(Now, r < vac_eff ? 1.0 : 0.0);
                    }
                }
            }
        } else {
            cerr << "ERROR: Unsupported vaccination event type: " << et << endl;
            exit(-4);
        }
    }

    int next_event() {
        if ( EventQ.empty() ) {
            return 0;
        } else {
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            Node* node = event.node;
            switch(event.type) {
                case EtoI_EVENT:
                    log_data[node->get_id()][INFECTIOUS] = Now;
                    update_state(node, EXPOSED, INFECTIOUS);
                    break;
                case ItoR_EVENT:    // recovery event
                    log_data[node->get_id()][RECOVERED] = Now;
                    update_state(node, INFECTIOUS, RECOVERED);
                    break;
                case ItoH_EVENT:    // hospitalization event
                    log_data[node->get_id()][HOSPITALIZED] = Now;
                    update_state(node, INFECTIOUS, HOSPITALIZED);
                    if (isPatientZero(node)) { schedule_vaccinations(); }
                    break;
                case ItoD_EVENT:    // death event
                    log_data[node->get_id()][DEAD] = Now;
                    update_state(node, INFECTIOUS, DEAD);
                    if (isPatientZero(node)) { schedule_vaccinations(); }
                    break;
                case StoE_EVENT:    // a contact event--contacted node might not actually be susceptible
                    //if (expose(node)) { update_state(node, SUSCEPTIBLE, EXPOSED); }
                    expose(node, event.source->get_id());
                    break;
                case V1_EVENT:      // intentional fall-through to V2_EVENT
                case V2_EVENT:
                    vaccine_campaign(event.type);
                    break;
                default:
                    cerr << "ERROR: Unsupported event type: " << event.type << endl;
                    exit(-2);
                    break;          // superfluous after exit(), but included in case refactoring makes it necessary
            }
            return 1;
        }
    }

};
#endif
