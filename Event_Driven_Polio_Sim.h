#ifndef ED_POLIO_SIM_H
#define ED_POLIO_SIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include "Utility.h"
#include "Network.h"

using namespace std;

class Event {
    public:
        double time;
        char type;
        Node* node;
        Event(const Event& o) {  time=o.time; type=o.type; node=o.node; }
        Event(double t, char e, Node* n) { time=t; type=e; node=n; }
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


class Event_Driven_Polio_Sim {
    public:
        // S -> E -> I -> R -> S
        typedef enum {
            SUSCEPTIBLE, VAC_PARTIAL, VAC_FULL, EXPOSED, INFECTIOUS, RESISTANT, STATE_SIZE // STATE_SIZE must be last
        } stateType;
                                    // constructor
        Event_Driven_Polio_Sim ( Network* net, double _lambda, double _epsilon, double _nu):lambda(_lambda), epsilon(_epsilon), nu(_nu) {
            network = net;
            reset();
        }

        Network* network;           // population
        double lambda;              // param for exponential exposed duration
        double epsilon;             // param for exponential time to transmission
        double nu;                  // param for exponential time to recovery

                                    // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        vector<int> state_counts;   // S, E, I, R counts
        double Now;                 // Current "time" in simulation

        MTRand mtrand;              // RNG

        void run_simulation(double duration) {
            double start_time = Now;
            int day = (int) Now;
            while (next_event() and Now < start_time + duration) {
                if ((int) Now > day) {
                    cout << (int) Now << " : "  << state_counts[SUSCEPTIBLE] << "\t" 
                                          << state_counts[EXPOSED] << "\t" 
                                          << state_counts[INFECTIOUS] << "\t" 
                                          << state_counts[RESISTANT] << endl; 
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
           
            vector<Node*> nodes = network->get_nodes();
            for (int i = 0; i < network->size(); i++) nodes[i]->set_state(SUSCEPTIBLE);

            state_counts.clear();
            state_counts.resize(STATE_SIZE, 0);
            state_counts[SUSCEPTIBLE] = network->size();
            
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        // choose n nodes without replacement
        vector<Node*> rand_choose_nodes (int n) {
            assert(n > -1 and n <= network->size());
            vector<Node*> nodes = network->get_nodes();
            vector<Node*> sample(n);
            vector<int> sample_ids(n);
            rand_nchoosek(network->size(), sample_ids, &mtrand);
            Node* node;
            for (unsigned int i = 0; i < sample_ids.size(); i++) {
                node = nodes[ sample_ids[i] ];
                sample[i] = node;
            };
            return sample;
        }


        void rand_infect(int n) {   // randomly infect k people
            vector<Node*> sample = rand_choose_nodes(n);
            for (unsigned int i = 0; i < sample.size(); i++) {
                infect(sample[i]);
            }
            return;
        }

        void infect(Node* node) {
            /*assert(state_counts[SUSCEPTIBLE] > 0);
            node->set_state(EXPOSED);
            state_counts[SUSCEPTIBLE]--;                // decrement susceptible group
            state_counts[EXPOSED]++;                    // increment exposed group

            double Ti = rand_exp(mu, &mtrand) + Now;    // time to become infectious
            add_event(Ti, 'i', node);
                                          
            double Tr = rand_exp(gamma, &mtrand) + Ti;  // time to recovery
            double Tc = rand_exp(beta, &mtrand) + Ti;   // time to next contact
            while ( Tc < Tr ) {                         // does contact occur before recovery?
                add_event(Tc, 'c', node);               // potential transmission event
                Tc += rand_exp(beta, &mtrand);
            }
            add_event(Tr, 'r', node);
                                    
            double Ts = Tr + immunity_duration;         // time to become susceptible again
            add_event(Ts, 's', node);*/
            return;
        }

        int next_event() {
            /*if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            Node* node = event.node;
            if (event.type == 'i') { 
                node->set_state(INFECTIOUS); 
                state_counts[EXPOSED]--;      // decrement Infected class
                state_counts[INFECTIOUS]++;   // increment Recovered class
            } else if (event.type == 'r') {   // recovery event
                node->set_state(RESISTANT);
                state_counts[INFECTIOUS]--;   // decrement Infected class
                state_counts[RESISTANT]++;    // increment Recovered class
            } else if (event.type == 's') {   // loss of immunity event
                node->set_state(SUSCEPTIBLE);
                state_counts[RESISTANT]--;
                state_counts[SUSCEPTIBLE]++;
            } else if (event.type == 'c') {                          // event type must be 'c'
                                 
                vector<Node*> neighbors = node->get_neighbors();
                if (neighbors.size() > 0) {
                    int rand_idx = mtrand.randInt(neighbors.size() - 1); // randInt includes endpoints
                    Node* contact = neighbors[rand_idx];
                    if ( contact->get_state() == SUSCEPTIBLE ) infect(contact);
                }
            } else {
                cerr << "Unknown event type encountered in simulator: " << event.type << "\nQuitting.\n";
            }*/
            return 1;
        }

        void add_event( double time, char type, Node* node) {
            EventQ.push( Event(time,type,node) );
            return;
        }

};
#endif
