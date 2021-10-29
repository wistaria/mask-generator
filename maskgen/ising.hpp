#pragma once

#include <cstdint>
#include <random>
#include <cluster/union_find.hpp>
#include "lattice.hpp"

using std::uint32_t;

class ising {
public:
  ising(uint32_t seed, uint32_t L, double beta) :
    eng_(seed), lattice_(L), spins_(lattice_.num_sites(), 1), fragments_(lattice_.num_sites()), flip_(lattice_.num_sites()) {
      if (beta < 0) beta = std::log(std::sqrt(2.0) + 1) / 2; /* Tc */
      prob_ = 1 - std::exp(-2 * beta);
  }

  void update() const {
    // initialize cluster information
    std::fill(fragments_.begin(), fragments_.end(), fragment_t());
      
    // cluster generation
    for (uint32_t b = 0; b < lattice_.num_bonds(); ++b) {
      if (spins_[lattice_.source(b)] == spins_[lattice_.target(b)] && r_uniform01_(eng_) < prob_)
        unify(fragments_, lattice_.source(b), lattice_.target(b));
    }
      
    // assign cluster id & accumulate cluster properties
    int nc = 0;
    for (auto& f : fragments_) {
      if (f.is_root()) f.set_id(nc++);
    }
    for (auto& f : fragments_) f.set_id(cluster_id(fragments_, f));
      
    // flip spins
    for (uint32_t c = 0; c < nc; ++c) flip_[c] = (r_uniform01_(eng_) < 0.5);
    for (uint32_t s = 0; s < lattice_.num_sites(); ++s)
      if (flip_[fragments_[s].id()]) spins_[s] ^= 1;
  }

  auto config() const { return spins_; }
  uint32_t spin(uint32_t s) const { return spins_[s]; } // return 0 or 1

  double energy() const {
    double ene = 0;
    for (uint32_t b = 0; b < lattice_.num_bonds(); ++b) {
      ene -= (spins_[lattice_.source(b)] == spins_[lattice_.target(b)] ? 1.0 : -1.0);
    }
    return ene / lattice_.num_sites();
  }
  double magnetization() const {
    double mag = 0;
    for (uint32_t s = 0; s < lattice_.num_sites(); ++s) mag += 2 * spins_[s] - 1;
    return mag / lattice_.num_sites();
  }        

private:
    mutable std::mt19937 eng_;
    mutable std::uniform_real_distribution<> r_uniform01_;
    square_lattice lattice_;
    double prob_;
    mutable std::vector<int> spins_;
    typedef cluster::union_find::node fragment_t;
    mutable std::vector<fragment_t> fragments_;
    mutable std::vector<int> flip_;
};
