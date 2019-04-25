#include <iostream>
#include <map>
#include <functional>
#include <random>

using namespace std;
using namespace std::placeholders;

int main() {
  mt19937 rng1, rng2, rng3;
  rng1.seed(1); rng2.seed(1); rng3.seed(1);
  uniform_int_distribution<int> sixsided(1,6);
  uniform_int_distribution<int> coin(1,2);

  std::cout << "Approach 1: " << '\n';
  for (auto i : {1,2,3}) std::cout << sixsided(rng1) << '\n';
  std::cout << coin(rng1) << '\n';

  // to make bind version work, have to force rng2 to be passed by reference
  std::cout << "Approach 2: " << '\n';
  auto roll = bind(sixsided, std::ref(rng2));
  auto toss = bind(coin, std::ref(rng2));
  for (auto i : {1,2,3}) std::cout << roll() << '\n';
  std::cout << toss() << '\n';

  std::cout << "Approach 3: " << '\n';
  vector<function<double()>> dice = {
    bind(sixsided, std::ref(rng3)),
    bind(sixsided, std::ref(rng3)),
    bind(sixsided, std::ref(rng3))
  };
  for (auto d : dice) std::cout << d() << '\n';
  std::cout << bind(coin, std::ref(rng3))() << '\n';

  return 0;
}
