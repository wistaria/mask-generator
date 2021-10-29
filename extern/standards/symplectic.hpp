/*
   Copyright (C) 2015-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#include <cmath>
#include <string>
#include <vector>

namespace standards {

class symplectic_euler {
public:
  static std::string name() { return "Symplectic Euler method"; }
  symplectic_euler(std::size_t dim) : dim_(dim), k_(dim) {}
  template<class VEC, class F>
  void step(double h, const FORCE& force, VEC* q, VEC* p) const {
    for (std::size_t i = 0; i < dim_; ++i) q[i] += h * p[i];
    force(q, k_);
    for (std::size_t i = 0; i < dim_; ++i) y[i] += h * k_[i];
  }
  template<class VEC, class F>
  void step_n(double h, std::size_t n, const FORCE& force, VEC* q, VEC* p) const {
    for (std::size_t s = 0; s < n; ++s) step(h, force, q, p);
  }
private:
  std::size_t dim_;
  mutable std::vector<double> k_;
};
  
class leap_frog {
public:
  static std::string name() { return "Leap frog method"; }
  leap_frog(std::size_t dim) : dim_(dim), k_(dim) {}
  template<class VEC, class F>
  void step(double h, const FORCE& force, VEC* q, VEC* p) const {
    const double h2 = h / 2;
    force(q, k_);
    for (std::size_t i = 0; i < dim_; ++i) p[i] += h2 * k_[i];
    for (std::size_t i = 0; i < dim_; ++i) q[i] += h * p[i];
    force(q, k_);
    for (std::size_t i = 0; i < dim_; ++i) p[i] += h2 * k_[i];
  }
  template<class VEC, class F>
  void step_n(double h, std::size_t n, const FORCE& force, VEC* q, VEC* p) const {
    if (n > 0) {
      const double h2 = h / 2;
      force(q, k_);
      for (std::size_t i = 0; i < dim_; ++i) p[i] += h2 * k_[i];
      for (std::size_t s = 0; s < n-1; ++s) {
        for (std::size_t i = 0; i < dim_; ++i) q[i] += h * p[i];
        force(q, k_);
        for (std::size_t i = 0; i < dim_; ++i) p[i] += h * k_[i];
      }
      for (std::size_t i = 0; i < dim_; ++i) q[i] += h * p[i];
      force(q, k_);
      for (std::size_t i = 0; i < dim_; ++i) p[i] += h2 * k_[i];
    }
  }
private:
  std::size_t dim_;
  mutable std::vector<double> k_;
};
  
class velocity_velret {
public:
  static std::string name() { return "Velocity Verlet method"; }
  velocity_velret(std::size_t dim) : dim_(dim), k1_(dim), k2_(dim) {}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const double h22 = h * h / 2;
    const double h2 = h / 2;
    const std::size_t dim2 = dim_ / 2;
    f(t, y, k1_);
    for (std::size_t i = 0; i < dim2; ++i) y[i] = y[i] +  h * k1_[i] + h22 * k1_[i + dim2];
    f(t + h, y, k2_);
    for (std::size_t i = dim2; i < dim_; ++i) y[i] = y[i] +  h2 * (k2_[i]+k1_[i]);
  }
private:
  std::size_t dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
};

// ref .American Journal of Physics 73, 938 (2005) 
//      or original paper : Phys. Lett. A 150, 262–268 (1990)

class yoshida_symplectic_4 {
public:
  static std::string name() { return "4th-order Yoshida symplectic method"; }
  yoshida_symplectic_4(std::size_t dim) : dim_(dim), k1_(dim), k2_(dim), k3_(dim){}
  template<class F>
  void step(double t, double h, std::vector<double>& y, F const& f) const {
    const int  dim2 = dim_ / 2;
    double c1 = 1.0/(4.0 - std::pow(2.0,4.0/3.0));  
    double c2 = (1.0 - std::pow(2.0,1.0/3.0))/(4.0 - std::pow(2.0,4.0/3.0));
    double b1 = 1.0 / (2.0 - std::pow(2.0,1.0/3.0));
    double b2 = 1.0 / (1.0 - std::pow(2.0,2.0/3.0));
    for (int i = 0; i < dim2; ++i) y[i] = y[i] +  h * c1 * y[i + dim2];
    f(t, y, k1_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b1 * k1_[i + dim2];
      y[i] = y[i] +  h * c2 * y[i + dim2];
    }
    f(t, y, k2_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b2 * k2_[i + dim2];
      y[i] = y[i] +  h * c2 * y[i + dim2];
    }
    f(t, y, k3_);
    for (int i = 0; i < dim2; ++i){
      y[i + dim2] = y[i + dim2] + h * b1 * k3_[i + dim2];
      y[i] = y[i] +  h * c1 * y[i + dim2];
    }
  }
private:
  std::size_t dim_;
  mutable std::vector<double> k1_;
  mutable std::vector<double> k2_;
  mutable std::vector<double> k3_;
};


} // namespace standards