#include <iostream>
#include <queue>
#include <string>
#include "../NetworkSimplate.h"

using namespace std;

using TC = Event<string>;

class TestSim : public EventDrivenSim<TC> {
  void process(TC event) {
    cout << event.which << endl;
  }
  void verbose(const double Now) {
    cout << Now << endl;
    return EventDrivenSim<TC>::verbose(Now);
  }
};

int main() {

  TestSim sim;
  sim.run({ TC(0.5,"bob",nullptr), TC(0.1,"alice",nullptr), TC(0.3,"carl",nullptr) });

  return 0;
}
