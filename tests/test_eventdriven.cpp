#include <iostream>
#include <queue>
#include <string>
#include "../NetworkSimplate.h"

using namespace std;

class TC {
public:
  string name;
  double time;
  TC(const string nm, const double t) {
    name = nm; time = t;
  }
  bool operator<(const TC& right) const {
    return (time < right.time);
  }
  bool operator>(const TC& right) const {
    return (time > right.time);
  }
};

class TestSim : public EventDrivenSim<TC> {
  void process(TC event) { 
    cout << event.name << endl;
  }
  double verbose(const double Now) {
    cout << Now << endl;
    return EventDrivenSim<TC>::verbose(Now);
  }
};

int main(int argc, char* argv[]) {

  TestSim sim;
  sim.run({ TC("bob",0.5), TC("alice",0.1), TC("carl",0.3) });
  
  return 0;
}
