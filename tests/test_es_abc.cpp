#include "../EbolaSim.h"
#include "AbcSmc.h"

const gsl_rng* GSL_RNG = gsl_rng_alloc (gsl_rng_taus2); // RNG for AbcSmc

inline double gamma_alpha(double mean, double sd) { return pow(mean/sd, 2); }
inline double gamma_beta(double mean, double sd)  { return pow(sd,2)/mean; }
inline function<double(mt19937&)> dgamma(double mn, double sd) { return gamma_distribution<double>(gamma_alpha(mn,sd), gamma_beta(mn,sd)); }
inline double rcoverage(double efficacy, mt19937& rng) {
  double u = uniform_real_distribution<double>(0,1)(rng), inveff = 1.0/efficacy;
  return inveff - sqrt(pow(inveff,2.0)-2.0*inveff*u + u);
}


vector<double> simulator(vector<double> args, const unsigned long int rng_seed, const unsigned long int /*serial*/, const ABC::MPI_par* /*mp*/) {
  const int net_replicate   = (int) args[0];
  //const int epi_replicate   = (int) args[1];
  const double back_vac_eff = args[2];
  const double trace_prob   = args[3];

  Network n(Network::Undirected);
  string network_dir      = "./networks/";
  string network_filename = ABC::get_nth_line("network_filenames.txt", net_replicate);
//    cerr << "Line " << net_replicate << " was " << network_filename << endl;
  n.read_edgelist(network_dir + network_filename, ',');
  assert(n.size() > 0);

  mt19937 globalrng(rng_seed);

// backkground coverage calculation:
// P(coverage=X|infection) = P(infection|coverage=X)P(coverage=X)/P(infection)
// P(infection|coverage=X) = 1-P(!infection|coverage=X) = 1 - coverage*efficacy
// CDF(coverage=X|infection) = int_0^X (1-coverage*efficacy)P(coverage=X)/P(infection)
// assuming uniform probability of coverage (for initial pass)
// P(infection)/P(coverage=X) = C = int_0^1 1 - coverage*efficacy = 1-efficacy/2
// so:
// CDF(coverage=X) = (1/(1-efficacy/2))*x(1-(efficacy/2)*x)
// which means if y ~ U, then 0 = x^2 - (2/eff)*x + (2/eff - 1)y
// which has a quadratic rule soln
// x = (-b +/- sqrt(b^2-4ac))/2a = 1/eff +/- sqrt((1/eff)^2-(2/eff - 1)y)
// at the asymptotic values of y (0,1), it's clear that the - term is correct

  double background_coverage = rcoverage(back_vac_eff, globalrng);

  SimPars ps = {
    &n, {
      { EXPOSE,    dgamma(2,1) }, // time to exposure
      { INCUBATE,  dgamma(5,2) }, // time to infectiousness, given exposure
      { RECOVER,   dgamma(9,4) }, // time to removed, given exposure
      { HOSPITAL,  dgamma(8,6) }  // time to removed, given exposure
    },
    trace_prob,
    back_vac_eff,
    background_coverage,  // background coverage
    globalrng
  };

  EbolaSim es(ps);
  cout << EbolaSim::loghead << endl;
  es.run(es.defaultEvents(), 100.0);
  vector<double> metrics(1,0);
  return metrics;
}


void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --process --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            buffer_size = atoi(argv[++i]);
        } else {
            usage();
            exit(101);
        }
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(GSL_RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(GSL_RNG);
    }

    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
