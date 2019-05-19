#ifndef NETWORK_SIMPLATE_H
#define NETWORK_SIMPLATE_H

#include <stdlib.h>
#include <iostream>   // provide pieces for stringify-ing Event<...>
#include <vector>     // provides vector
#include <queue>      // provides priority_queue
#include <functional> // provides greater (used w/ priority_queue)

#include "Network.h"  // provides Node

using namespace std;
using namespace std::placeholders;

// EE := Event Enum
// NT := Node Type, defaults to standard EpiFire Node
template<class EE, class NT = Node>
struct Event {
  double tm;
  EE which;
  NT* node;
  NT* source;

  // Event(const Event& o) { tm=o.tm; which=o.which; node = o.node; source=o.source; }
  Event(double t, EE e, NT* n, NT* s = nullptr) : tm(t), which(e), node(n), source(s) {}
  // ~Event() {}

  // Event& operator=(const Event& o) {
  //   tm=o.tm; which=o.which; node=o.node; source=o.source;
  //   return *this;
  // }

  bool operator<(const Event& right) const { return (tm < right.tm); }
  bool operator>(const Event& right) const { return (tm > right.tm); }
  double time() const { return tm; }
  void time(double shift) { tm += shift; }

};

template<class EE>
ostream& operator<<(std::ostream &out, const Event<EE> &e) {
  return out << "Event(" << e.time() << ", " << e.which << ")";
}

template<class EE, class ET, class SM, class RT = bool>
map<EE,function<bool(ET&)>> mapper(vector<EE> es, vector<RT (SM::*)(ET&)> fs, SM* that) {
  map<EE,function<RT(ET&)>> res;
  auto et = es.begin(); auto ft = fs.begin();
  for (; et != es.end() and ft != fs.end(); et++, ft++) res[*et] = bind(*ft, that, _1);
  return res;
}

template<class EE, class ET, class SM, class RT = bool>
pair<EE, function<RT(ET&)>> fpair(EE e, RT (SM::*f)(ET&), SM* that) {
  return { e, bind(f, that, _1) };
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
        double time = next_event(timelimit);
        while (time < timelimit) time = next_event(timelimit);
        return time;
    }

    double next_event(const double maxtime) {
        if ( EventQ.empty() ) {
            return maxtime;
        } else {
            auto event = EventQ.top(); // get the element
            EventQ.pop();              // remove from Q
            process(event);            // processing event may add items to Q, so this must be first
            verbose(event.time());     // squawk as necessary
            return event.time();
        }
    }

    // template these? i.e. might be able to use deduction on ET constructor args
    void add_event( ET et ) { EventQ.push(et); }

    // hook to squawk as time passes
    virtual void verbose(const double /* Now */) {}
    // hook to handle events
    virtual void process(ET /* event */) {}

    priority_queue<ET, vector<ET>, greater<ET> > EventQ;
    virtual void reset() { EventQ = priority_queue<ET, vector<ET>, greater<ET> >(); }

};
#endif
