#include <fstream>
#include <iostream>
#include <standards/accumulator.hpp>
#include "percolation.hpp"

struct options {
  uint32_t seed, length;
  double prob;
  uint32_t sweeps;
  std::string file;
  bool valid;

  options(uint32_t argc, char *argv[], bool print = true) :
    // default parameters
    seed(29833), length(8), prob(0.5), sweeps(1 << 16), file(),
    valid(true) {
    for (uint32_t i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 's' :
          if (++i == argc) { usage(print); return; }
          seed = std::atoi(argv[i]); break;
        case 'l' :
          if (++i == argc) { usage(print); return; }
          length = std::atoi(argv[i]); break;
        case 'p' :
          if (++i == argc) { usage(print); return; }
          prob = std::atof(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(print); return; }
          sweeps = std::atoi(argv[i]); break;
        case 'c' :
          if (++i == argc) { usage(print); return; }
          file = argv[i]; break;
        case 'h' :
          usage(print, std::cout); return;
        default :
          usage(print); return;
        }
        break;
      default :
        usage(print); return;
      }
    }
    if (length == 0 || prob < 0) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "System Linear Size     = " << length << std::endl
                << "Bond Probability       = " << prob << std::endl
                << "MCS                    = " << sweeps << std::endl
                << "Output filename        = " << file << std::endl;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -l int    System Linear Size\n"
         << "  -p double Bond Probability\n"
         << "  -m int    MCS\n"
         << "  -c name   Output final configuration to file\n"
         << "  -h        this help\n";
    valid = false;
  }
};

int main(int argc, char* argv[]) {
  std::cout << "Bond Percolation on Square Lattice\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  percolation simulator(p.seed, p.length, p.prob);
  standards::accumulator nc("Number of Clusters"), wmax("Size of Largest Cluster");

  for (uint32_t mcs = 0; mcs < p.sweeps; ++mcs) {
    simulator.update();
    nc << simulator.num_clusters();
    wmax << 1.0 * simulator.max_size() / (p.length * p.length);
  }
  std::cout << nc << std::endl
            << wmax << std::endl;
  if (p.file.length()) {
    std::ofstream os(p.file);
    for (uint32_t s = 0; s < p.length * p.length; ++s) {
      os << simulator.spin(s) << ' ';
      if (s % p.length == p.length - 1) os << std::endl;
    }
    std::cout << 1.0 * simulator.max_size() / (p.length * p.length) << std::endl;
  }
}
