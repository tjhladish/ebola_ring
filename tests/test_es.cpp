#include "../EbolaSim.h"

int main(int argc, char* argv[]) {
  Network n("test", Network::Undirected);
  n.read_edgelist(argv[1], ',');
  SimPars ps = {
    &n, {
      { EXPOSE,   [](mt19937& rng) { rng.discard(1); return 5.0; }},
      { INCUBATE, [](mt19937& rng) { rng.discard(1); return 5.0; }},
      { DIE, [](mt19937& rng) { rng.discard(1); return 10.0; }},
      { RECOVER, [](mt19937& rng) { rng.discard(1); return 5.0; }},
      { HOSPITAL, [](mt19937& rng) { rng.discard(1); return 2.0; }}
    }, 0
  };

  EbolaSim es(ps);
  es.run(es.defaultEvents(), 10.0);
  return 0;
}
