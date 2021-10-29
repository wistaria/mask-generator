#pragma once

#include <cstdint>
#include <random>
#include <cluster/union_find.hpp>
#include "lattice.hpp"

using std::uint32_t;

class percolation {
public:
  percolation(uint32_t seed, uint32_t L, double prob) :
    eng_(seed), lattice_(L), fragments_(lattice_.num_sites()), prob_(prob) {
      if (prob_ < 0) prob_ = 0.5;
  }

  void update() {
    // initialize cluster information
    std::fill(fragments_.begin(), fragments_.end(), fragment_t());
      
    // cluster generation
    for (uint32_t b = 0; b < lattice_.num_bonds(); ++b) {
      // std::cout << b << ' ' << (b & 1) << ' ' << lattice_.source(b) << ' ' << lattice_.target(b) << std::endl;
      if (r_uniform01_(eng_) < prob_)
        unify(fragments_, lattice_.source(b), lattice_.target(b));
    }
      
    // assign cluster id & accumulate cluster properties
    int nc = 0;
    for (auto& f : fragments_) {
      if (f.is_root()) f.set_id(nc++);
    }
    for (auto& f : fragments_) f.set_id(cluster_id(fragments_, f));
      
    // accumulate cluster properties
    nc_ = 0;
    wmax_ = 0;
    imax_ = 0;
    for (auto& f : fragments_) {
      if (f.is_root()) {
        ++nc;
        auto w = f.weight();
        if (w >= wmax_) {
          wmax_ = w;
          imax_ = f.id();
        }
      }
    }
  }

  auto config() const {
    std::vector<uint32_t> conf;
    for (uint32_t s = 0; s < lattice_.num_sites(); ++s) conf.push_back(spin(s));
    return conf;
  }
  uint32_t spin(uint32_t s) const { return (fragments_[s].id() == imax_) ? 1 : 0; }

  double num_clusters() const { return nc_; }

  double max_size() const { return wmax_; }

private:
    std::mt19937 eng_;
    std::uniform_real_distribution<> r_uniform01_;
    square_lattice lattice_;
    double prob_;
    typedef cluster::union_find::node fragment_t;
    std::vector<fragment_t> fragments_;
    uint32_t nc_, wmax_, imax_;
};
