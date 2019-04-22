#include <iostream>
#include <queue>
#include <string>

using namespace std;

class TC {
public:
  string name;
  double time;
  TC(string nm, double t) : name(nm), time(t) {}
  bool operator<(const TC& right) const {
    return (time < right.time);
  }
  bool operator>(const TC& right) const {
    return (time > right.time);
  }
};

int main(int argc, char* argv[]) {
  auto TQ = priority_queue<TC,vector<TC>,greater<TC>>();
  TQ.push(TC("alice",0.1)); TQ.push(TC("bob",0.5)); TQ.push(TC("carl",0.3));
  
  while(!TQ.empty()) { cout << TQ.top().name << " "; TQ.pop(); }
  cout << '\n';
  
  return 0;
}
