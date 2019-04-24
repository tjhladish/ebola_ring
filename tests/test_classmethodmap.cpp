#include <iostream>
#include <map>
#include <functional>
#include <iostream>

using namespace std;
using namespace std::placeholders;

enum EgState {
  STATE_0, STATE_1, N_STATES
};

class Example {
public:

  Example() {}

  bool foo(EgState& eg) { return eg == STATE_1; } // false
  bool bar(EgState& eg) { return eg == STATE_1; } // true
  bool err(EgState& eg) { exit(5); return eg == N_STATES; } // true, but exit

  map<EgState,function<bool(EgState&)>> test{
    {STATE_0, bind(&Example::foo, this, _1)},
    {STATE_1, bind(&Example::bar, this, _1)},
    {N_STATES, bind(&Example::err, this, _1)}
  };

  bool process(EgState eg) {
    return test[eg](eg);
  }

};

int main() {
  auto exemplar = Example();
  for (auto eg : {STATE_0, STATE_1, N_STATES}) cout << (exemplar.process(eg) ? "TRUE" : "FALSE") << "\n";
  return 0;
}
