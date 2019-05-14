#include <iostream>
#include "../NetworkSimplate.h"

using namespace std;

enum EgState1 {
  A, B, C, N_STATES_ONE
};

enum EgState2 {
  D, E, F, N_STATES_TWO
};

int main() {
  auto exemplar1 = Event<EgState1>{ 0.0, A, nullptr }, exemplar3 = Event<EgState1>{ 0.0, C, nullptr };
  auto exemplar2 = Event<EgState2>{ 0.1, D, nullptr };
  // auto exemplar4 = Event{0.1, E, nullptr }; // error: missing template arguments before ‘{’ token
  // supposedly impute-able in later standards
  Event<EgState1> exemplar4{0.0, B, nullptr };
  // intentional failing check:
  // Event<EgState2> exemplar5{0.0, B, nullptr }; // error: no matching function for call to
  cout << exemplar1 << " & " << exemplar2 << " & " <<  exemplar3 << " & " <<  exemplar4 << endl;
  return 0;
}
