#ifndef EBOLA_SIM_H
#define EBOLA_SIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <random>
#include <map>
#include <functional>
#include <unistd.h>
#include "IDCommunity.h"
#include "NetworkSimplate.h"
#include "Utility.h"

using namespace std;
using namespace std::placeholders; // gives _1, _2, etc

enum EventType {
  EXPOSE,
  INCUBATE,
  HOSPITAL,
  RECOVER,
//  DIE,
  VACCINATE,
  TRACE,
  N_EVENTS // must be last
};

std::ostream& operator<<(std::ostream &out, EventType &e) {
  string res;
  switch (e) {
    case EXPOSE:    res = "EXPOSE";    break;
    case INCUBATE:  res = "INCUBATE";  break;
    case HOSPITAL:  res = "HOSPITAL";  break;
    case RECOVER:   res = "RECOVER";   break;
//    case DIE:       res = "DIE";       break;
    case VACCINATE: res = "VACCINATE"; break;
    case TRACE:     res = "TRACE";     break;
    case N_EVENTS:  res = "N_EVENTS";  break;
  }
  return out << res;
}

using EboEvent = Event<EventType>; // from NetworkSimplate, POD event container for EventType

struct SimPars {
  Network* net;
  map<EventType, function<double(mt19937&)>> event_time_distribution;
  double traceProbability;
  double backgroundEff;
  double backgroundCoverage;
  mt19937& sharedrng;
};

class EbolaSim : public EventDrivenSim<EboEvent> {
  public:

    void process(EboEvent event) {
      bool success = process_steps[event.which](event);
      DiseaseState DSupdate = disease_outcome[event.which];
      ControlCondition CCupdate = control_outcome[event.which];
      //DiseaseState updateTo = disease_outcome[event.which];
      if (success) {
        logupdate(event, DSupdate);
        logupdate(event, CCupdate);
        if (DSupdate == EXPOSED) {
          // update infector log
        }
        // don't need to update TRACED log? since that's already a thing that can be printed from
      }
    }

    map< EventType, function<bool(EboEvent&)> > process_steps = {
      fpair(EXPOSE,   &EbolaSim::exposure, this),
      fpair(INCUBATE, &EbolaSim::incubation, this),
//      fpair(DIE,      &EbolaSim::death, this),
      fpair(HOSPITAL, &EbolaSim::hospital, this),
      fpair(RECOVER,  &EbolaSim::recovery, this),
      fpair(TRACE,  &EbolaSim::tracing, this),
      fpair(VACCINATE,  &EbolaSim::vaccinate, this)
    };
    //   { EXPOSE, INCUBATE, HOSPITAL, RECOVER,
    //     DIE, VACCINATE, TRACE, N_EVENTS },
    //   { &EbolaSim::exposure, &EbolaSim::incubation, &EbolaSim::hospital, &EbolaSim::recovery,
    //     &EbolaSim::death, &EbolaSim::vaccination, &EbolaSim::tracing, &EbolaSim::eventerror },
    //   this
    // );

    map< EventType, DiseaseState > disease_outcome = {
      { EXPOSE,   EXPOSED      },
      { INCUBATE, INFECTIOUS   },
      { RECOVER,  RECOVERED },
//      { DIE, DEAD },
      { HOSPITAL, N_STATES },
      { VACCINATE, N_STATES },
      { TRACE, N_STATES },
    };

    map< EventType, ControlCondition > control_outcome = {
      { EXPOSE,   N_CONDITIONS },
      { INCUBATE, N_CONDITIONS   },
      { RECOVER,  N_CONDITIONS },
//      { DIE, N_CONDITIONS },
      { HOSPITAL, HOSPITALIZED },
      { VACCINATE, VACCINATED },
      { TRACE, TRACED },
    };

    EbolaSim(SimPars& pars) :
    // assorted simple constructions
      rng(pars.sharedrng),
      community(IDCommunity(pars.net, pars.backgroundCoverage, pars.sharedrng)),
      index_case(pars.net->get_node(0)),
      disease_log_data(
        pars.net->size(),
        vector<double>(N_STATES, numeric_limits<double>::quiet_NaN())
      ),
      traceProbability(pars.traceProbability),
      backgroundEff(pars.backgroundEff)
    // construction with a bit more complexity
    {
        // concept: input parameters provides the distributions
        // but internally, the simulation code only invokes the generators via ()
        // with all generators locked to the same rng
        for_each(
          pars.event_time_distribution.begin(), pars.event_time_distribution.end(),
          [&](pair<const EventType, function<double(mt19937&)>>& p) {
            event_time_generator[p.first] = bind(p.second, std::ref(rng));
          }
        );
        reset();
    }

    mt19937& rng; // random number generator
    map<EventType, function<double()>> event_time_generator;
    double runif() { return uniform_real_distribution<double>(0,1)(rng); }

    int day;
    IDCommunity community;
    Node* index_case;
    vector< vector<double> > disease_log_data;
    int control_radius = 2;
    double traceProbability, backgroundEff;

    virtual void reset() {
      disease_log_data.clear();
      disease_log_data.resize(community.size(), vector<double>(N_STATES, std::numeric_limits<double>::quiet_NaN()));
      community.reset();
      EventDrivenSim<EboEvent>::reset(); // also reset eventQ
    }

    vector<EboEvent> defaultEvents() {
      return { EboEvent(0, INCUBATE, index_case) };
    }

    bool canTransmit(Node* source) {
      return !source ? true : // if source is null, means external introduction
        // otherwise, can transmit if source is infectious & not hospitalized
        source->get_state() == INFECTIOUS && not community.isQuarantined(source);
    }

    double logistic(const double x, const double offset, const double stretch, const double k = 1.25, const double mid = 3.5) const {
      return(stretch*(1.0/(1.0+exp(-k*(x-mid)))-offset));
    }

    const double refoffset = logistic(0.0,0.0,1.0);
    const double refstretch = 1.0/(1.0-2.0*refoffset);

    double ringVaccineEff(double timeSinceRingVaccination) {
      if (timeSinceRingVaccination < 0.0) {
        return 0.0;
      } else {
        // see R/ring_vax_model_pars.R for value picking
        return min(logistic(timeSinceRingVaccination, refoffset, refstretch), 1.0);
      }
    }

    bool isSusceptible(Node* target, double when) {
      if (target->get_state() == SUSCEPTIBLE) {
        double beff = community.hasBackground(target) ? backgroundEff : 0.0;
        double reff = ringVaccineEff(when - community.ringVaccineTime(target));
        double efficacy = max(beff, reff);
        return efficacy < runif();
      } else return false;
    }

    bool exposure(EboEvent& event) {
      if (isSusceptible(event.node, event.time()) && canTransmit(event.source)) {

        auto newEvent = event; // copy old event
        newEvent.which = INCUBATE; // update event type
        newEvent.time(event_time_generator[newEvent.which]()); // draw new time

        add_event(newEvent); // add the event to the queue

        community.update_state(event.node, EXPOSED); // update community accounting

        return true; // success
      } else return false; // no exposure
    }

    bool incubation(EboEvent& event) {
      assert(event.node->get_state() == EXPOSED | !event.source);

      auto newEvent = event; // definitely having at least one new control_radievent, so copy previous

      // draw times for possible outcomes
      double // tDeath    = event_time_generator[DIE](),
             tRecover  = event_time_generator[RECOVER](),
             tHospital = community.isNodeTraced(event.node) ? 2.0 : event_time_generator[HOSPITAL]();
             // ANOTHER MAGIC NUMBER - reduced hospitalization time is 2.0
      // which outcome?
      newEvent.which = tHospital < tRecover ? HOSPITAL : RECOVER;
      newEvent.time(min(tHospital, tRecover));

      add_event(newEvent);
      community.update_state(event.node, INFECTIOUS);

      double tInfectiousInCommunity = newEvent.time() - event.time();
      auto refExpEvent = event;
      refExpEvent.source = event.node;
      refExpEvent.which = EXPOSE;

      for (auto neighbor : event.node->get_neighbors()) {
        if (neighbor->get_state() == SUSCEPTIBLE) {
          double expTime = event_time_generator[EXPOSE]();
          if (expTime < tInfectiousInCommunity) {
            auto expEvent = refExpEvent;
            expEvent.node = neighbor;
            expEvent.time(expTime);
            add_event(expEvent);
          }
        }
      }

      return true;
    }

    bool notYetTraced = true;

    bool hospital(EboEvent& event) {
      if (event.node->get_state() == INFECTIOUS) {
        if (not community.isTraced()) {
          auto traceEvent = event;
          traceEvent.which = TRACE;
          traceEvent.source = traceEvent.node;
          add_event(traceEvent);
        }
        community.quarantine(event.node);
        return true;
      } else return false;
    }

    bool recovery(EboEvent& event) {
      if (event.node->get_state() == INFECTIOUS) {
        community.update_state(event.node, RECOVERED);
        return true;
      } else return false;
    }

    void sendTraceEvents(EboEvent& event) {
      assert(event.which == TRACE);
      EboEvent refEvent = event;
      // swap the reference event to be *from* this node
      refEvent.source = event.node;
      // the proposed new level for *targets*
      int refLvl = community.get_level(event.node)+1;
      // check all outgoing edges
      for (auto edge : event.node->get_edges_out()) {
        // if that edge has been found
        if (community.isFound(edge)) {
          auto target = edge->get_end();
          if (community.get_level(target) > refLvl) {
            EboEvent sendEvent = refEvent;
            sendEvent.node = target;
            add_event(sendEvent);
          }
        }
      }
    }

    bool initialTrace(EboEvent& event) {
      int level = event.node->get_state() >= INFECTIOUS ? 0 : community.get_level(event.source) + 1;
      for (auto edge : event.node->get_edges_out()) {
        // draw trace success
        bool isFound = runif() < traceProbability;
        // record trace success or not
        community.trace(edge, isFound);
        // potentially update level
        if (isFound) {
          int otherLvl = community.get_level(edge->get_end()) + 1;
          if (otherLvl < level) level = otherLvl;
        }
      }
      // set level
      community.set_level(event.node, level);
      if (level <= control_radius) {
        if (level < control_radius) sendTraceEvents(event);
        if (level != 0) {
          auto vacEvent = event;
          vacEvent.which = VACCINATE;
          vacEvent.time(2.0); // MAGIC NUMBER; vaccinate 2 days later
          add_event(vacEvent);
        }
      }
      return true;
    }

    bool reTrace(EboEvent& event) {
      int proposedLevel = community.get_level(event.source)+1;
      if (proposedLevel < community.get_level(event.node)) {
        // adjust my level
        community.set_level(event.node, proposedLevel);
        // notify my identified edges
        sendTraceEvents(event);
        return true;
      } else return false;
    }

    bool notification(EboEvent& event) {
      assert(event.node == event.source);
      if (community.isTraced()) {
        community.set_level(event.node, IDCommunity::NOTIFY_LEVEL);
        // for now: only one round of tracing
        return true;
      } else {
        community.setTraced();
        return initialTrace(event);
      }
    }

    bool tracing(EboEvent& event) {
      if (event.node == event.source) return notification(event);
      return community.isNodeTraced(event.node) ? reTrace(event) : initialTrace(event);
    }

    bool vaccinate(EboEvent& event) {
      // hasn't been previously vaccinated...
      assert(community.ringVaccineTime(event.node) == std::numeric_limits<double>::infinity());
      community.set_ringVaccineTime(event.node, event.time());
      return true;
    }


/*    bool death(EboEvent& event) {
      if (event.node->get_state() == INFECTIOUS) {
        // TODO if tracing hasn't happened yet, draw time-to-detection event
        community.update_state(event.node, DEAD);
        return true;
      } else return false;
    }

    bool eventerror(EboEvent event) {
      cerr << "ERROR: Unsupported event type: " << event.which << endl;
      exit(-2);
      return false;
    }
*/

    static const string loghead;

    void logupdate(EboEvent event, DiseaseState ds) {
      if (ds != N_STATES) {
        cout << "Event(" <<
        event.which << " on " << event.node->get_id() <<
        " @ " << event.time();
        if (event.source) cout << " from " << event.source->get_id();
        cout <<  ")" << endl;
        disease_log_data[event.node->get_id()][ds] = event.time();
      }
    }

    void logupdate(EboEvent event, ControlCondition cc) {
      if (cc != N_CONDITIONS) {
        cout << "Event(" <<
        event.which << " on " << event.node->get_id() <<
        " @ " << event.time();
        if (event.source) cout << " from " << event.source->get_id();
        cout <<  ")" << endl;
        disease_log_data[event.node->get_id()][cc] = event.time();
      }
    }

};

const string EbolaSim::loghead = "Event(WHICH on Target @ Time[ from Source])";

#endif
