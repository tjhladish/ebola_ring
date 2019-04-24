#include "../EbolaSim.h"

int main(int argc, char* argv[]) {
  Network n("test", Network::Undirected);
  n.read_edgelist(argv[1], ',');
  n.graphviz(argv[2]);
  EbolaSim es(&n);
  es.run(es.defaultEvents(), 10.0);
  return 0;
}
