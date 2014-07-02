// Copyright 2013 by Chris Dyer
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

#ifndef _TTABLES_H_
#define _TTABLES_H_

#include <cmath>

#include <vector>
#include "corpus.h"

namespace Md {
inline double digamma(double x) {
  double result = 0, xx, xx2, xx4;
  for (; x < 7; ++x) result -= 1 / x;
  x -= 1.0 / 2.0;
  xx = 1.0 / x;
  xx2 = xx * xx;
  xx4 = xx2 * xx2;
  result += log(x) + (1. / 24.) * xx2 - (7.0 / 960.0) * xx4 +
            (31.0 / 8064.0) * xx4 * xx2 - (127.0 / 30720.0) * xx4 * xx4;
  return result;
}
}

class TTable {
  typedef std::unordered_map<Token, double> Word2Double;
  typedef std::vector<Word2Double> Word2Word2Double;
  Word2Word2Double ttable;

 public:

  TTable(): ttable() {}

  inline double prob(const Token& e, const Token& f) const {
    if (e < ttable.size()) {
      const Word2Double& cpd = ttable[e];
      const Word2Double::const_iterator it = cpd.find(f);
      if (it == cpd.end()) return 1e-9;
      return it->second;
    } else {
      return 1e-9;
    }
  }

  inline void Increment(const Token& e, const Token& f) { 
    if (e >= ttable.size()) {
      ttable.resize(e + 1);
    }
    ttable[e][f] += 1.0; 
  }
  
  inline void Increment(const Token& e, const Token& f, double x) {
    if (e >= ttable.size()) {
      ttable.resize(e + 1);
    }
    ttable[e][f] += x;
  }

  inline void clear() {
    ttable.clear();
  }

  void UpdateFrom(TTable& other) {
    ttable.swap(other.ttable);
    other.ttable.clear();
    Normalize();
  }

  void UpdateFromVB(TTable& other, const double alpha) {
    ttable.swap(other.ttable);
    other.ttable.clear();
    NormalizeVB(alpha);
  }

  void NormalizeVB(const double alpha) {
    for (auto &cpd: ttable) {
      double tot = 0;
      for (auto &it: cpd)
        tot += it.second + alpha;
      for (auto &it: cpd)
        it.second = exp(Md::digamma(it.second + alpha) - Md::digamma(tot));
    }
  }

  void Normalize() {
    for (auto &cpd: ttable) {
      double tot = 0;
      for (auto &it: cpd)
        tot += it.second;
      for (auto &it: cpd)
        it.second /= tot;
    }
  }

  // adds counts from another TTable - probabilities change!
  TTable& operator+=(const TTable& rhs) {
    if (ttable.size() < rhs.ttable.size()) {
      ttable.resize(rhs.ttable.size());
    }
    for (size_t ind = 0; ind < rhs.ttable.size(); ++ind) {
      Word2Double const &cpd = rhs.ttable[ind];
      Word2Double &tgt = ttable[ind];
      for (auto const &j: cpd) {
        tgt[j.first] += j.second;
      }
    }
    return *this;
  }

  void ExportToFile(std::string &filename, Dict& d) {
    std::ofstream file(filename);
    for (size_t ind = 0; ind < ttable.size(); ++ind) {
      const std::string& a = d.Convert(ind);
      Word2Double& cpd = ttable[ind];
      for (auto const &it: cpd) {
        const std::string& b = d.Convert(it.first);
        double c = log(it.second);
        file << a << '\t' << b << '\t' << c << std::endl;
      }
    }
    file.close();
  }
};

#endif
