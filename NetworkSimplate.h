#ifndef NETWORK_SIMPLATE_H
#define NETWORK_SIMPLATE_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include <random>
#include <iomanip>
#include <map>
#include <functional>
#include <unistd.h>

using namespace std;

// template<class ET, class NT = Node>
// class Event
// {
//   public:
//     double time;
//     ET type;
//     NT* node, source;
// 
//     Event(const Event& o) : time(o.time), type(o.type), node(o.node), source(o.source) {}
//     Event(double t, ET e, NT* n, NT* s = nullptr) : time(t), type(e), node(n), source(s) {}
//     Event& operator=(const Event& o) {
//       time=o.time; type=o.type; node=o.node; source=o.source;
//       return *this;
//     }
// 
//     bool operator<(const Event& right) const { return (time < right.time); }
//     bool operator>(const Event& right) const { return (time > right.time); }
// 
// };

template<class ET>
class EventDrivenSim {
  public:
                                // constructor
    EventDrivenSim() {
        reset();
    }

    priority_queue<ET, vector<ET>, greater<ET> > EventQ;

    void reset() { EventQ = priority_queue<ET, vector<ET>, greater<ET> >(); }

    void run_simulation(
      vector<ET> initial_events, // initial list of events - e.g., background vaccination, incipient infection
      const double timelimit = numeric_limits<double>::max()
    ) {
        for (auto e : initial_events) add_event(e);
        while (verbose(next_event(timelimit)) < timelimit) continue;
    }

    virtual void process(ET event) {}

    double next_event(const double maxtime) {
        if ( EventQ.empty() ) {
            return maxtime;
        } else {
            auto event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q
            process(event);
            return event.time;
        }
    }

    // template these?
    void add_event( ET et ) { EventQ.push(et); }
    // hook to do something as time passes
    virtual double verbose(double Now) { return(Now); }

};
#endif
