// roughly
// bool prevented(time now, time last_vac_event_time, state last_vac_event_type)

// probably also needs to dole out events
// vac_event apply(Node n, time now)

class VaccineAbstract {
  public:
    virtual bool isProtect(double nw, Event lastvacevent, const mt19937& rng) = false;
    virtual Event operator() (const Node* n, const double tm) = Event(tm, NUM_OF_EVENT_TYPES, n);
}

class VaccineNonLeaky : public VaccineAbstract {
public:
  bool isProtect(double nw, Event lastvacevent, const mt19937& rng) override {
    return((lastvacevent == EFF_VACCINE) && (nw - lastvacevent.time > onset));
  }
  Event operator() override (const Node* n, const double tm, const ... rng) override {
    EventType et = runif(rng) > coverage ?
      MISSED_VACCINE : (runif(rng) < efficacy ? EFF_VACCINE : NON_VACCINE);
    return(Event(tm, et, n))
  }
  VaccineNonLeaky(double cov, double eff, double ons) {
    assert(0 <= cov && cov <= 1);
    assert(0 <= eff && eff <= 1);
    assert(0 <= ons);
    coverage = cov;
    efficacy = eff;
    onset    = ons;
  }
protected:
  double coverage;
  double efficacy;
  double onset;
}

class VaccineLeakyLinear : public VaccineAbstract {
public:
  bool isProtect(double nw, Event lastvacevent, const mt19937& rng) override {
    return((lastvacevent == EFF_VACCINE) && (runif(rng) < max((nw - lastvacevent.time)/onset, 1)*efficacy));
  }
  Event operator() override (const Node* n, const double tm, const ... rng) override {
    EventType et = runif(rng) > coverage ? MISSED_VACCINE : EFF_VACCINE;
    return(Event(tm, et, n))
  }
  VaccineNonLeaky(double cov, double eff, double ons) {
    assert(0 <= cov && cov <= 1);
    assert(0 <= eff && eff <= 1);
    assert(0 <= ons);
    coverage = cov;
    efficacy = eff;
    onset    = ons;
  }
protected:
  double coverage;
  double efficacy;
  double onset;
}

class Vaccine {
  public:
    Vaccine () : efficacy({1.0,1.0}), coverage({1.0,1.0}), isLeaky(true), timeToProtection(numeric_limits<double>::max()) {}
    vector<double> efficacy; // doses 1 and 2
    vector<double> coverage; // same
    bool isLeaky;
    double timeToProtection;
};
