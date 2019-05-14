#ifndef NETWORK_SIMPLATE_H
#define NETWORK_SIMPLATE_H

#include <stdlib.h>
#include <iostream>   // provide pieces for stringify-ing Event<...>
#include <vector>     // provides vector
#include <queue>      // provides priority_queue
#include <functional> // provides greater (used w/ priority_queue)

#include "Network.h"  // provides Node

using namespace std;

// EE := Event Enum
// NT := Node Type, defaults to standard EpiFire Node
template<class EE, class NT = Node>
struct Event {
  double time;
  EE which;
  NT* node;
  NT* source;

  // Event(const Event& o) { time=o.time; type=o.type; node = o.node; source=o.source; }
  Event(double t, EE e, NT* n, NT* s = nullptr) : time(t), which(e), node(n), source(s) {}
  // ~Event() {}

  // Event& operator=(const Event& o) {
  //   time=o.time; type=o.type; node=o.node; source=o.source;
  //   return *this;
  // }

  bool operator<(const Event& right) const { return (time < right.time); }
  bool operator>(const Event& right) const { return (time > right.time); }

};

template<class EE, class NT = Node>
ostream& operator<<(std::ostream &out, const Event<EE,NT> &e) {
  return out << "Event(" << e.time << ", " << e.which << ")";
}

// ET := Event Type;
// intent: an Event<SomeEventEnum[,OptionalNodeType]>
// requirement: something that defines >(ET)->bool and time/time()->double
template<class ET>
class EventDrivenSim {
  public:

    EventDrivenSim() { reset(); }

    double run(
      // initial list of events - e.g., background vaccination, incipient infection
      const vector<ET> initial_events = vector<ET>(),
      const double timelimit = numeric_limits<double>::max()
    ) {
        for (auto e : initial_events) add_event(e);
        double time = verbose(next_event(timelimit));
        while (time < timelimit) time = verbose(next_event(timelimit));
        return time;
    }

    double next_event(const double maxtime) {
        if ( EventQ.empty() ) {
            return maxtime;
        } else {
            auto event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q
            process(event);
            return event.time();
        }
    }

    // template these? i.e. might be able to use deduction on ET constructor args
    void add_event( ET et ) { EventQ.push(et); }

    // hook to squawk as time passes
    virtual void verbose(const double Now) {}
    // hook to handle events
    virtual void process(ET event) {}

    priority_queue<ET, vector<ET>, greater<ET> > EventQ;
    virtual void reset() { EventQ = priority_queue<ET, vector<ET>, greater<ET> >(); }

};
#endif
