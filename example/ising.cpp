#include <fstream>
#include <iostream>
#include <standards/accumulator.hpp>
#include "ising.hpp"

struct options {
  uint32_t seed, length;
  double temperature;
  uint32_t sweeps, therm;
  std::string file;
  bool valid;

  options(uint32_t argc, char *argv[], bool print = true) :
    // default parameters
    seed(29833), length(8), temperature(2 / std::log(std::sqrt(2.0) + 1)), sweeps(1 << 16), therm(sweeps >> 3), file(),
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
        case 't' :
          if (++i == argc) { usage(print); return; }
          temperature = std::atof(argv[i]); break;
        case 'm' :
          if (++i == argc) { usage(print); return; }
          sweeps = std::atoi(argv[i]);
          therm = sweeps >> 3; break;
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
    if (length == 0 || temperature <= 0. || sweeps == 0) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "System Linear Size     = " << length << std::endl
                << "Temperature            = " << temperature << std::endl
                << "MCS for Thermalization = " << therm << std::endl
                << "MCS for Measurement    = " << sweeps << std::endl
                << "Output filename        = " << file << std::endl;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -l int    System Linear Size\n"
         << "  -t double Temperature\n"
         << "  -m int    MCS for Measurement\n"
         << "  -c name   Output final configuration to file"
         << "  -h        this help\n";
    valid = false;
  }
};

int main(int argc, char* argv[]) {
  std::cout << "Swendsen-Wang Cluster Algorithm for Square Lattice Ising Model\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  ising simulator(p.seed, p.length, 1.0 / p.temperature);
  standards::accumulator energy("Energy Density"),
    magnetization("Magnetization"), magnetization2("Magnetization^2");

  for (uint32_t mcs = 0; mcs < (p.therm + p.sweeps); ++mcs) {
    simulator.update();
    if (mcs >= p.therm) {
      energy << simulator.energy();
      double mag = simulator.magnetization();
      magnetization << mag;
      magnetization2 << mag * mag;
    }
  }
  std::cout << energy << std::endl
            << magnetization << std::endl
            << magnetization2 << std::endl;
  if (p.file.length()) {
    std::ofstream os(p.file);
    for (uint32_t s = 0; s < p.length * p.length; ++s) {
      os << simulator.spin(s) << ' ';
      if (s % p.length == p.length - 1) os << std::endl;
    }
    std::cout << simulator.magnetization() << std::endl;
  }
}
