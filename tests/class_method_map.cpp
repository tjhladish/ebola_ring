#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include "../NetworkSimplate.h"

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

  pair<EgState, function<bool(EgState&)>> fpair(EgState e, bool (Example::*f)(EgState&)) {
    return { e, bind(f, this, _1) };
  }

  map<EgState,function<bool(EgState&)>> test{
    {STATE_0, bind(&Example::foo, this, _1)},
    {STATE_1, bind(&Example::bar, this, _1)},
    {N_STATES, bind(&Example::err, this, _1)}
  };

  bool process(EgState eg) {
    return test[eg](eg);
  }

  map<EgState,function<bool(EgState&)>> test2{
    fpair(STATE_0, &Example::foo),
    fpair(STATE_1, &Example::bar),
    fpair(N_STATES, &Example::err)
  };

  map<EgState,function<bool(EgState&)>> test3 = mapper<EgState,EgState,Example>(
    {STATE_0,       STATE_1,       N_STATES},
    {&Example::foo, &Example::bar, &Example::err},
    this
  );

  // TODO find a way to do binding
  // map<EgState,function<bool(EgState&)>> test3;
  // transform(test2.begin(), test2.end(), inserter(test3, test3.end()), [&](auto& p))
  //
  bool process2(EgState eg) {
    return test2[eg](eg);
  }

  bool process3(EgState eg) {
    return test3[eg](eg);
  }

};

int main() {
  auto exemplar = Example();
  for (auto eg : {STATE_0, STATE_1}) cout << (exemplar.process(eg) ? "TRUE" : "FALSE") << "\n";
  for (auto eg : {STATE_1, STATE_0}) cout << (exemplar.process2(eg) ? "2TRUE" : "2FALSE") << "\n";
  for (auto eg : {STATE_0, STATE_1}) cout << (exemplar.process2(eg) ? "3TRUE" : "3FALSE") << "\n";
  return 0;
}
