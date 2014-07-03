// Copyright 2014 by Douwe Gelling
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef _DA_H_
#define _DA_H_

#include <cmath>
#include <utility>
#include <algorithm>

#include "model.h"

// m = trg len
// n = src len
// i = trg index
// j = src index
namespace DiagonalAlignment {

  // find diagonal for a given position. Clamp return values between 0 and n, as these are
  // the extreme values that could have been obtained without an offset
  unsigned Diagonal(const unsigned i, const unsigned m, const unsigned n, const double o) {
    return std::max(std::min(static_cast<unsigned>((double(i) * (n+1)) / (m+1) + o * (n+1)), n), 0u);
  }

};

// likely doesn't work atm
struct Squared {
  static constexpr double init_learningrate_lk = 40000;
  static constexpr double init_learningrate_o = .1;

  static constexpr double min_lk = -20;
  static constexpr double max_lk = 150;
  static constexpr double min_o = -1;
  static constexpr double max_o = 1;

  static constexpr double init_kappa = 70;
  static constexpr double init_lambda = 70;

  static inline double Transform(const double x) {
    return x * x;
  }

  static inline double dTransform(const double x) {
    return x;
  }

  static inline double Feature(const unsigned i, const unsigned j,
                               const unsigned m, const unsigned n, const double o) {
    return double(i) / (m+1) - double(j) / (n+1) + o;
  }

  static inline double UnnormalizedProb(const unsigned i, const unsigned j,
                                 const unsigned m, const unsigned n,
                                 const double kappa, const double offset) {
    return exp(Transform(Feature(i, j, m, n, offset)) * kappa * -1);
  }

  static std::pair<double, double> ComputeZs(const unsigned i, const unsigned m, const unsigned n,
                         const double kappa, const double lambda, const double offset) {
    const unsigned floor = DiagonalAlignment::Diagonal(i,m,n,offset);
    const unsigned ceil = floor + 1;
    double ezt = 0;
    double ezb = 0;

    for (size_t j = 1; j <= floor; ++j) {
      ezb += UnnormalizedProb(i,j,m,n,kappa, offset);
    }
    for (size_t j = ceil; j <= n; ++j) {
      ezt += UnnormalizedProb(i,j,m,n,lambda, offset);
    }
    return std::make_pair(ezb,ezt);
  }

  static double ComputeZ(const unsigned i, const unsigned m, const unsigned n,
                         const double kappa, const double lambda, const double offset) {
    auto zs = ComputeZs(i, m, n, kappa, lambda, offset);
    return zs.first + zs.second;
  }

  static Params ComputeDLogZs(const unsigned i, const unsigned m,
                             const unsigned n, const double kappa, const double lambda, const double offset) {
    const double z = ComputeZ(i, m, n, kappa, lambda, offset);
    const unsigned floor = DiagonalAlignment::Diagonal(i,m,n,offset);
    const unsigned ceil = floor + 1;
    double pct = 0;
    double pcb = 0;
    double pot = 0;
    double pob = 0;

    for (size_t j = 1; j <= floor; ++j) {
      double x = Feature(i,j,m,n,offset);
      double single = UnnormalizedProb(i,j,m,n,kappa,offset) * x;
      pcb += single * x ;
      pob += single;
    }
    pob *= kappa;
    for (size_t j = ceil; j <= n; ++j) {
      double x = Feature(i,j,m,n,offset);
      double single = UnnormalizedProb(i,j,m,n,lambda,offset) * x;
      pct += single * x;
      pot += single;
    }
    pot *= lambda;

    //std::cerr << "prior offset: " << (pob + pot)/z << std::endl;

    return {pcb/z, pct/z, (pob + pot)/z};
  }
};


struct Absolute {
  static constexpr double init_learningrate_lk = 1000;
  static constexpr double init_learningrate_o = 0.025;

  static constexpr double min_lk = -10;
  static constexpr double max_lk = 30;
  static constexpr double min_o = -1;
  static constexpr double max_o = 1;

  static constexpr double init_kappa = 8;
  static constexpr double init_lambda = 8;

  static inline double Transform(double x) {
    if (x < 0) return -x;
    else return x;
  }

  static inline double dTransform(double x) {
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
  }


  static inline double Feature(const unsigned i, const unsigned j,
                               const unsigned m, const unsigned n, const double o) {
    return double(i) / (m+1) - double(j) / (n+1) + o;
  }

  static inline double UnnormalizedProb(const unsigned i, const unsigned j,
                                 const unsigned m, const unsigned n,
                                 const double kappa, const double offset) {
    return exp(Transform(Feature(i, j, m, n, offset)) * kappa * -1);
  }

  static std::pair<double, double> ComputeZs(const unsigned i, const unsigned m, const unsigned n,
                         const double kappa, const double lambda, const double offset) {
    const unsigned floor = DiagonalAlignment::Diagonal(i,m,n,offset);
    const unsigned ceil = floor + 1;
    double ezt = 0;
    double ezb = 0;

    for (size_t j = 1; j <= floor; ++j) {
      ezb += UnnormalizedProb(i,j,m,n,kappa, offset);
    }
    for (size_t j = ceil; j <= n; ++j) {
      ezt += UnnormalizedProb(i,j,m,n,lambda, offset);
    }
    return std::make_pair(ezb,ezt);
  }

  static double ComputeZ(const unsigned i, const unsigned m, const unsigned n,
                         const double kappa, const double lambda, const double offset) {
    auto zs = ComputeZs(i, m, n, kappa, lambda, offset);
    return zs.first + zs.second;
  }

  static Params ComputeDLogZs(const unsigned i, const unsigned m,
                             const unsigned n, const double kappa, const double lambda, const double offset) {
    const double z = ComputeZ(i, m, n, kappa, lambda, offset);
    const unsigned floor = DiagonalAlignment::Diagonal(i,m,n,offset);
    const unsigned ceil = floor + 1;
    double pct = 0;
    double pcb = 0;
    double pot = 0;
    double pob = 0;

    //std::cerr << "position: " << i << std::endl;

    for (size_t j = 1; j <= floor; ++j) {
      double x = Feature(i,j,m,n,offset);
      double single = UnnormalizedProb(i,j,m,n,kappa,offset);
      pcb += single * Transform(x);
      pob += single * dTransform(x);
    }
    pob *= kappa;
    for (size_t j = ceil; j <= n; ++j) {
      double x = Feature(i,j,m,n,offset);
      double single = UnnormalizedProb(i,j,m,n,lambda,offset);
      pct += single * Transform(x);
      pot += single * dTransform(x);
    }
    pot *= lambda;

    //std::cerr << "contribution: " << (pob + pot) / z << std::endl;

    return {pcb/z, pct/z, (pob + pot) /z};
  }
};

#endif
