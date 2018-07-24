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

using namespace std;

//int seed = 0;             // to use a static seed
//mt19937 gen(seed);

uniform_real_distribution<double> runif(0.0, 1.0);
// mean  = alpha/beta
// sd    = sqrt(alpha)/beta
// alpha = (mean/sd)**2
// beta  = mean/alpha --> mean/(sd**2 * beta**2) --> (mean/sd**2)**(1/3)

inline double gamma_alpha(double mean, double sd) {return pow(mean/sd, 2);}
inline double gamma_beta(double mean, double sd) {return pow(mean/pow(sd,2), 1.0/3);}

gamma_distribution<double> rSE(gamma_alpha(2,1),       gamma_beta(2,1));   // 2,1 are totally arbitary--prob needs to be fit
gamma_distribution<double> rEI(gamma_alpha(6,2),       gamma_beta(6,2));   // 6,2
gamma_distribution<double> rIR(gamma_alpha(9,4),       gamma_beta(9,4));   // 9,4
gamma_distribution<double> rIH(gamma_alpha(2,1),       gamma_beta(2,1));   // 2,1
gamma_distribution<double> rIH_pzero(gamma_alpha(5,3), gamma_beta(5,3));   // 5,3
gamma_distribution<double> rID(gamma_alpha(8,4),       gamma_beta(8,4));   // 8,4 arbitary numbers
gamma_distribution<double> rV1(gamma_alpha(5,3),       gamma_beta(5,3));   // 5,3
random_device rd;                       // generates a random real number for the seed
mt19937 rng;                      // random number generator


enum StateType { SUSCEPTIBLE,
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
    Event(const Event& o) {  time=o.time; type=o.type; node=o.node; }
    Event(double t, EventType e, Node* n) { time=t; type=e; node=n; }
    Event& operator=(const Event& o) { time=o.time; type=o.type; node=o.node; return *this; }
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

class Event_Driven_Ebola_Sim {
  public:
                                // constructor
    Event_Driven_Ebola_Sim ( Network* net, const Vaccine& _vac, unsigned int seed = rd()) : vaccine(_vac) {
        rng.seed(seed);
        network = net;
        node_vac_immunity.clear();
        node_vac_immunity.resize(network->size(), make_pair(0.0, 0.0));
        vaccine_doses_used = {0, 0}; // doses 1 and 2
        reset();
    }

    Network* network;           // population
    const Vaccine vaccine;
    vector<int> vaccine_doses_used;
    vector<pair<double, double> > node_vac_immunity; // time of vaccination, efficacy

    priority_queue<Event, vector<Event>, compTime > EventQ;
    vector<int> state_counts;   // S, E, I, R, etc. counts
    double Now;                 // Current "time" in simulation

    bool isPatientZero(Node* node) { return node->get_id() == 0; }

    double time_to_event (EventType et) {
        double delta_t = 0.0; // time in days
        switch (et) {
            case StoE_EVENT:  // per-neighbor time until contact
                delta_t = rSE(rng);
                break;
            case EtoI_EVENT:  // incubation period
                delta_t = rEI(rng);
                break;
            case ItoR_EVENT:  // infectious period
                delta_t = rIR(rng);
                break;
            case ItoH_EVENT:  // time until hospitalization
                delta_t = rIH(rng);
                break;
            case ItoH_PZERO_EVENT:  // time until hospitalization
                delta_t = rIH_pzero(rng);
                break;
            case ItoD_EVENT:  // time until hospitalization
                delta_t = rID(rng);
                break;
            case V1_EVENT  :  // time from index case is detected until mass vaccination
                delta_t = rV1(rng);
                break;
            case V2_EVENT  :  // time between dose 1 and dose 2
                delta_t = 28;
                break;
            default:
                cerr << "ERROR: Unsupported event type in time_to_event(): " << et << endl;
                exit(-3);
                break;
        }
        return delta_t;
    }

    void run_simulation(double duration = numeric_limits<double>::max()) {
        double start_time = Now;
        int day = (int) Now;
        while (next_event() and Now < start_time + duration) {
            if ((int) Now > day) {
                cout << "(SvV | EI | HRD) " << Now << " :\t"
                                  << state_counts[SUSCEPTIBLE] << "\t"
                                  << state_counts[VAC_PARTIAL] << "\t"
                                  << state_counts[VAC_FULL] << "\t|\t"
                                  << state_counts[EXPOSED] << "\t"
                                  << state_counts[INFECTIOUS] << "\t|\t"
                                  << state_counts[HOSPITALIZED] << "\t"
                                  << state_counts[RECOVERED] << "\t"
                                  << state_counts[DEAD] << endl;
                day = (int) Now;
            }

            continue;
        }
    }

    int current_epidemic_size() {
        return state_counts[EXPOSED] + state_counts[INFECTIOUS];
    }

    void reset() {
        Now = 0.0;

        for (Node* node: network->get_nodes()) node->set_state(SUSCEPTIBLE);

        state_counts.clear();
        state_counts.resize(NUM_OF_STATE_TYPES, 0);
        state_counts[SUSCEPTIBLE] = network->size();

        EventQ = priority_queue<Event, vector<Event>, compTime > ();
    }

    bool expose(Node* node) {
        //cerr << "node: " << node->get_id() << endl;
        bool didTransmissionOccur = false;
        StateType state = (StateType) node->get_state();
        bool unprotected_by_vaccine = true;
        if (state == VAC_PARTIAL or state == VAC_FULL) {     // determine if node is protected by vaccine
            const double vac_time = node_vac_immunity[node->get_id()].first;
                                                             // current assumption is vaccine abruptly assumes full efficacy after timeToProtection
            const double vac_protection = (vac_time - Now > vaccine.timeToProtection) ? node_vac_immunity[node->get_id()].second : 0.0;
            const double asdf = runif(rng);
            unprotected_by_vaccine = runif(rng) > vac_protection;
        }

        if (state == SUSCEPTIBLE or state == VAC_PARTIAL or state == VAC_FULL) {
            if (unprotected_by_vaccine) {
                double Ti = Now + time_to_event(EtoI_EVENT);      // time to become infectious
                add_event(Ti, EtoI_EVENT, node);

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

                for (Node* neighbor: node->get_neighbors()) {     // density-dependent assumption! more neighbors --> more contact per unit time
                    double Tc = Ti + time_to_event(StoE_EVENT);   // time to next contact--we'll worry about whether contact is susceptible at Tc
                    while ( Tc < Ti_end ) {                       // does contact occur before recovery?
                        add_event(Tc, StoE_EVENT, neighbor);      // potential transmission event
                        Tc += time_to_event(StoE_EVENT);          // if vaccine is leaky, or successful transmission is probabilistic, this is impt
                    }
                }
                didTransmissionOccur = true;
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

    void update_state(Node* node, StateType old_state, StateType new_state) {
        node->set_state(new_state);
        state_counts[old_state]--; state_counts[new_state]++;
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
                if (state == INFECTIOUS) { continue; } // all others will be potentially vaccinated

                if (runif(rng) < vaccine.coverage[dose_idx]) {
                    vaccine_doses_used[dose_idx]++;

                    if (state == eligible_group) {
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
                    update_state(node, EXPOSED, INFECTIOUS);
                    break;
                case ItoR_EVENT:    // recovery event
                    update_state(node, INFECTIOUS, RECOVERED);
                    break;
                case ItoH_EVENT:    // hospitalization event
                    update_state(node, INFECTIOUS, HOSPITALIZED);
                    if (isPatientZero(node)) { schedule_vaccinations(); }
                    break;
                case ItoD_EVENT:    // death event
                    update_state(node, INFECTIOUS, DEAD);
                    if (isPatientZero(node)) { schedule_vaccinations(); }
                    break;
                case StoE_EVENT:    // a contact event--contacted node might not actually be susceptible
                    //if (expose(node)) { update_state(node, SUSCEPTIBLE, EXPOSED); }
                    expose(node);
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

    void add_event( double time, EventType et, Node* node = nullptr) {
        EventQ.push( Event(time,et,node) );
        return;
    }
};
#endif
