#include "../EbolaSim.h"
#include "AbcSmc.h"

const gsl_rng* GSL_RNG = gsl_rng_alloc (gsl_rng_taus2); // RNG for AbcSmc

inline double gamma_alpha(double mean, double sd) { return pow(mean/sd, 2); }
inline double gamma_beta(double mean, double sd)  { return pow(sd,2)/mean; }
inline function<double(mt19937&)> dgamma(double mn, double sd) { return gamma_distribution<double>(gamma_alpha(mn,sd), gamma_beta(mn,sd)); }

enum VaccineMechanism {
  LEAKY,
  NONLEAKY,
  N_VACMECHS
};

vector<double> simulator(vector<double> args, const unsigned long int /* rng_seed */, const unsigned long int serial, const ABC::MPI_par* /*mp*/) {

  std::uniform_int_distribution<int> reseed(0,1000);
  std::uniform_real_distribution<double> runif(0,1);
  const int net_replicate   = (int) args[0];
  mt19937 globalrng(net_replicate);
  const int epi_replicate   = (int) args[1];

  globalrng.seed(reseed(globalrng)*(epi_replicate+1));

  //const double background_coverage = EbolaSim::rcoverage(args[2], globalrng);

  globalrng.seed(epi_replicate);

  const double trace_prob   = args[3];
  const double exp_sd   = args[4];
  const double vaccine_delay   = args[5];
  const VaccineMechanism back_vac_mech   = static_cast<VaccineMechanism>(args[6]);
  const double background_coverage = args[7];
  const int networksdir = static_cast<int>(args[8]);
  // if back_vac_mech == 0, leaky
  // if back_vac_mech == 1, non-leaky
  // use the efficacious efficacy for bias whatever the vaccine mechanism
//  double background_coverage = useBias ? EbolaSim::rcoverage(args[2], globalrng) : uniform_real_distribution<double>(0,1)(globalrng);

  // if the mechanism is nonleaky, deprecate the background coverage to only the efficaciously vaccinated individuals
  // if efficacy is zero, don't modify
  // if (back_vac_mech == NONLEAKY and args[2] != 0.0) background_coverage *= args[2];

  //const double back_vac_eff = (back_vac_mech == LEAKY or args[2] == 0.0) ? args[2] : 1.0;
  const double back_vac_eff = args[2];

  // at this point, all the first case to be potentially resisted;
  // if it is, quick return, otherwise continue
  bool first_vaccinated = runif(globalrng) < background_coverage;
  bool first_resisted = runif(globalrng) < back_vac_eff;
  vector<double> metrics;
  if (!(first_vaccinated and first_resisted)) {
    Network n(Network::Undirected);
    string network_dir      = networksdir == 0 ? "./networks/" : networksdir == 1 ? "./declustered/" : "./reclustered/";
    string network_filename = ABC::get_nth_line("network_filenames.txt", net_replicate);
  //    cerr << "Line " << net_replicate << " was " << network_filename << endl;
    n.read_edgelist(network_dir + network_filename, ',', false);
    assert(n.size() > 0);

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

    SimPars ps = {
      &n, {
        { EXPOSE,    dgamma(12, exp_sd) }, // time to exposure
        { INCUBATE,  dgamma(9.9, 5.5) }, // time to infectiousness, given exposure
        { RECOVER,   dgamma(8.9, 4) }, // time to removed, given exposure
        { HOSPITAL,  dgamma(8, 6) }  // time to removed, given exposure
      },
      trace_prob,
      back_vac_eff,
      background_coverage,  // background coverage
      globalrng,
      vaccine_delay,
      6.0, // offset look time limit
      back_vac_mech
    };

    EbolaSim es(ps, first_vaccinated);
    es.run(es.defaultEvents());
  //  EbolaSim::dump(cout, es, to_string(serial));
    metrics = {
      // background_coverage,
      es.countPreVax(), // INFECTIOUS AT VACCINE TIME
      es.count_at(true, true, 0), // vaccinated, EBOV pos, at vaccine-time
      es.count_at(true, false, 0), // vaccinated, EBOV neg, at vaccine-time
      es.count_at(false, true, 0), // not, EBOV pos, at vaccine-time
      es.count_at(false, false, 0)/*, // not, EBOV neg, at vaccine-time
      es.count_at(true, true, 6), // vaccinated, EBOV pos, at vaccine-time
      es.count_at(true, false, 6), // vaccinated, EBOV neg, at vaccine-time
      es.count_at(false, true, 6), // not, EBOV pos, at vaccine-time
      es.count_at(false, false, 6) // not, EBOV neg, at vaccine-time */
    };

    EbolaSim::dump(cout, es, to_string(serial));

  } else {
    metrics = {
      // background_coverage,
      0, // INFECTIOUS AT VACCINE TIME
      0, // vaccinated, EBOV pos, at vaccine-time
      0, // vaccinated, EBOV neg, at vaccine-time
      0, // not, EBOV pos, at vaccine-time
      0 /*, // not, EBOV neg, at vaccine-time
      es.count_at(true, true, 6), // vaccinated, EBOV pos, at vaccine-time
      es.count_at(true, false, 6), // vaccinated, EBOV neg, at vaccine-time
      es.count_at(false, true, 6), // not, EBOV pos, at vaccine-time
      es.count_at(false, false, 6) // not, EBOV neg, at vaccine-time */
    };
  }

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
