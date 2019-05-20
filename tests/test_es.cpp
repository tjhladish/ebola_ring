#include "../EbolaSim.h"

inline double gamma_alpha(double mean, double sd) { return pow(mean/sd, 2); }
inline double gamma_beta(double mean, double sd)  { return pow(sd,2)/mean; }
inline function<double(mt19937&)> dgamma(double mn, double sd) { return gamma_distribution<double>(gamma_alpha(mn,sd), gamma_beta(mn,sd)); }

int main(int /* argc */, char* argv[]) {
  Network n("test", Network::Undirected);
  n.read_edgelist(argv[1], ',');
  mt19937 sharedrng;
  SimPars ps = {
    &n, {
      { EXPOSE,    dgamma(4,2) }, // time to exposure
      { INCUBATE,  dgamma(5,2) }, // time to infectiousness, given exposure
      { RECOVER,   dgamma(9,4) }, // time to removed, given exposure
      { HOSPITAL,  dgamma(8,6) }  // time to removed, given exposure
    },
    1.0, // trace probability
    0.7, // background efficacy
    0.7,  // background coverage
    sharedrng
  };

  EbolaSim es(ps);
  cout << EbolaSim::outhead << endl;
  es.run(es.defaultEvents());
  cout << es;
  return 0;
}
