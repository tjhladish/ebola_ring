#include <iostream>
#include "../IDCommunity.h"

int main(int argc, char* argv[]) {
  Network n("test", Network::Undirected);
  n.read_edgelist(argv[1], ',');
  n.graphviz(argv[2]);
  IDCommunity idc(&n);
  
  cout << idc.current_epidemic_size() << '\n';
  return 0;
}
