
enum EventType {
  StoE_EVENT,
  EtoI_EVENT,
  ItoR_EVENT,
  ItoH_EVENT,
  ItoH_PZERO_EVENT,
  ItoD_EVENT,
  V_EVENT,
  TRACING_EVENT,
  NUM_OF_EVENT_TYPES
}; // must be last

class Event {
public:
  double time;
  EventType type;
  Node* node;
  Node* source; // for exposure events

  Event(double t, EventType e, Node* n, Node* s = nullptr) : time(t), type(e), node(n), source(s) {}
  Event(const Event& o) : time(o.time), type(o.type), node(o.node), source(o.source) {}
  Event& operator=(const Event& o) { time=o.time; type=o.type; node=o.node; source=o.source; return *this; }

  bool operator<(const Event& right) const {
    return (time < right.time);
  }
  bool operator>(const Event& right) const {
    return (time > right.time);
  }
};

class Event_Driven_Ebola_Sim : public NetworkSimplate<Event> {
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
        reset();
    }

    int day;
    vector< vector<double> > log_data;
    IDCommunity community;

    virtual void reset() {
      // derived class reset
      day = 0;
      for (auto el: log_data) { el.assign(NUM_OF_STATE_TYPES, -1); }
      community.reset();
      // base reset
      NetworkSimplate<Event>::reset();
    }

    // TODO implement vaccination model(s)
    // rough idea: vaccine object(s) that hold universal information about vaccine
    // and nodes hold state information about their particular vaccine events (e.g., when, what kind(s), assignment for non-leaky vax)
    bool vaccine_protected(Node* node) {
      return false;
    }

    void ascertain(Node* node, Node* src) {
      // check node as ascertained - if already ascertained at equal or lower level that source level + 1,
      // don't need to re-notify other neighbors.  *Might* to need to notify source. I.e., if L2 node connects to an L0 node
      // otherwise, proceed
      // if currently infectious, new hospitalized event
      // test neighbors (or reuse previous identifications)
      // register contact tracing events for found neighbors (ignoring source neighbor)
    }

    bool expose(Node* node, Node* src, double eventtime) {
        //cerr << "node: " << node->get_id() << endl;
        auto state = static_cast<StateType>(node->get_state());
        if (state == NUM_OF_STATE_TYPES) {
          cerr << "ERROR: Exposure of node in unsupported state.  Node state is " << state << endl;
          exit(-1);
        }

        bool didTransmissionOccur = (src->get_state() == INFECTIOUS) and
          (static_cast<DiseaseState>(node->get_state()) == SUSCEPTIBLE) and
          (not protected(node));

        // add infectiousness event
        if (didTransmissionOccur) {
          add_event(eventtime + time_to_event(EtoI_EVENT), EtoI_EVENT, node);
          community.update_state(node, EXPOSED);
        }

        return didTransmissionOccur;
    }

    void onset(Node* node) {
      double hospTime = 0;
      double recoverTime = 0;
      double deathTime = 0;
      if(ascertain_list[node->get_id()]) {
        // draw different hosp time
      } else {
        // draw other hosp time
      }
      // draw death time
      // draw recover time
      double tLimit = min(hospTime, recoverTime, deathTime);
      // draw contact times for each neighbor
      // add exposures for each less than tLimit
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
      if (ascertain_list[node->get_id()]) { // if already ascertained
        // trigger additional tracing?
      } else {
        // trigger initial tracing? addtional tracing?
      }
      community.update_state(node, HOSPITALIZED);
    }

    double verbose(double Now) {
      if (Now > day) {
        community.log(cerr, Now);
        day = static_cast<int>(Now);
      }
      return Now;
    }

    vector< vector<double> > run_simulation(
      vector<Event> evts, // initial list of events - e.g., background vaccination, incipient infection
      const double timelimit = numeric_limits<double>::max()
    ) {
        for (auto e : evts) add_event(e);
        while (verbose(next_event(timelimit)) < timelimit) continue;
        return log_data;
    }

    void logupdate(Event event, DiseaseState to) {
      log_data[event.n->get_id()][to] = event.time;
    }

    void process(Event event) {
      switch(event.type) {
        case EtoI_EVENT:
          if(onset(event.node)) logupdate(et, INFECTIOUS);
          break;
        case ItoR_EVENT:
          if(recover(event.node)) logupdate(et, RECOVERED);
          break;
        case ItoH_EVENT:
          if(hospitalize(event.node)) logupdate(et, HOSPITALIZED);
          break;
        case HtoD_EVENT:
          if(hospitaldeath(event.node)) logupdate(et, DEAD);
          break;
        case ItoD_EVENT:
          if(communitydeath(event.node)) logupdate(et, DEAD);
          break;
        case StoE_EVENT:
          if(expose(event.node, event.source, event.time)) {
            logupdate(et, EXPOSED);
            // TODO infector log
          }
          break;
        case V1_EVENT:
          if(vaccine_campaign(event.type)) {

          }
          break;
        case TRACING_EVENT:
          if(ascertain(event.node, event.source)) {
            // TODO ascertainment level log?
          }
          break;
        default:
            cerr << "ERROR: Unsupported event type: " << event.type << endl;
            exit(-2);
            break;          // superfluous after exit(), but included in case refactoring makes it necessary
      }
    }


// TODO incorporate log handling here by having all event handlers return bool
    double next_event(const double maxtime) {
        if ( EventQ.empty() ) {
            return maxtime;
        } else {
            auto event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            switch(event.type) {
                case EtoI_EVENT: if(onset(event.node)) {
                  log_data[et.node->get_id()][INFECTIOUS] = event.time;
                }
                break;
                case ItoR_EVENT: recover(event.node);                        break;
                case ItoH_EVENT: hospitalize(event.node);                    break;
                case HtoD_EVENT: hospitaldeath(event.node);                  break;
                case ItoD_EVENT: communitydeath(event.node);                 break;
                case StoE_EVENT: expose(event.node, event.source, event.time); break;
                case V1_EVENT:   vaccine_campaign(event.type);               break;
                case TRACING_EVENT: ascertain(event.node, event.source);     break;
                default:
                    cerr << "ERROR: Unsupported event type: " << event.type << endl;
                    exit(-2);
                    break;          // superfluous after exit(), but included in case refactoring makes it necessary
            }
            log(event, event.time);
            return event.time;
        }
    }

    void add_event( double time, EventType et, Node* node = nullptr, Node* source = nullptr) {
        add_event(Event(time, et, node, source));
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