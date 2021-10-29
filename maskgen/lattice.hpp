#pragma once

#include <cstdint>

using std::uint32_t;

class square_lattice {
public:
  square_lattice(uint32_t L) : L_(L) {}
  uint32_t num_sites() const { return L_ * L_; }
  uint32_t num_bonds() const { return 2 * num_sites(); }
  uint32_t source(uint32_t b) const { return (b >> 1); }
  uint32_t target(uint32_t b) const {
    uint32_t s = source(b);
    uint32_t x = s % L_;
    uint32_t y = s / L_;
    if ((b & 1) == 0) {
      x = (x + 1) % L_;
    } else {
      y = (y + 1) % L_;
    }
    return L_ * y + x;
  }
private:
  uint32_t L_;
};
