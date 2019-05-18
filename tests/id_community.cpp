#include <iostream>
#include "../IDCommunity.h"

int main(int /* argc */, char* argv[]) {
  Network n("test", Network::Undirected);
  n.read_edgelist(argv[1], ',');
  n.graphviz(argv[2]);
  IDCommunity idc(&n);

  cout << "net edge count: "<< n.get_edges().size() << endl;

  cout << idc.header << endl;
  cout << idc << endl;

  for (auto tos : { EXPOSED, INFECTIOUS, RECOVERED }) {
    idc.update_state(0, tos);
    cout << idc << endl;
  }

  idc.reset();
  cout << idc << '\n';

  return 0;
}
