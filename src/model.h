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

#ifndef _MODEL_H_
#define _MODEL_H_

#include <algorithm>                    // for fill
#include <string>                       // for string
#include <tuple>                        // for tuple
#include <unordered_map>                // for unordered_map
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "ttables.h"                    // for TTable
class Dict;


struct Params {
  double kappa;
  double lambda;
  double offset;

  Params(double k, double l, double o) : kappa(k), lambda(l), offset(o) {}
  Params() : kappa(0), lambda(0), offset(0) {}

  Params &operator+=(Params const &rhs) {
    kappa += rhs.kappa;
    lambda += rhs.lambda;
    offset += rhs.offset;
    return *this;
  }
};

// assumes that int and size_t are 64 bits. Depends on platform and implementation
struct TupHash {
  size_t operator()(const std::tuple<unsigned short, unsigned short, unsigned short> &x) const {
    return (size_t)std::get<0>(x) << 32ul | (size_t)std::get<1>(x) << 16ul | std::get<2>(x);
  }
};

typedef std::unordered_map<std::tuple<unsigned short, unsigned short, unsigned short>, unsigned, TupHash> tup_map;

inline tup_map &operator+=(tup_map &lhs, tup_map const &rhs) {
  for (auto &pair : rhs) {
    lhs[pair.first] += pair.second;
  }
  return lhs;
}

typedef std::pair<unsigned, unsigned> uintpair;
inline uintpair &operator+=(uintpair &lhs, uintpair const &rhs) {
  lhs.first += rhs.first;
  lhs.second += rhs.second;
  return lhs;
}

template <typename T>
inline std::vector<T> &addToVec(std::vector<T> &lhs, std::vector<T> const &rhs) {
  if (lhs.size() < rhs.size()) {
    lhs.resize(rhs.size());
  }
  for (size_t i = 0; i < rhs.size(); ++i) {
    lhs[i] += rhs[i];
  }
  return lhs;
}

// stats collected during E-M, for reporting and optimizing
struct Stats {
  TTable tcounts;
  std::vector<Params> emp_counts;
  std::vector<Params> mod_counts;
  std::vector<uintpair> toks;

  std::vector<tup_map> size_index_counts;
  double tot_len_ratio = 0;

  double likelihood = 0;
  int lc = 0;
  double c0 = 0;
  unsigned denom = 0;

  void reset() {
    fill(emp_counts.begin(), emp_counts.end(), Params{0,0,0});
    likelihood = 0;
    lc = 0;
    c0 = 0;
  }

  void add_count(Stats &other, bool first_iter = false) {
    tcounts += other.tcounts;
    other.tcounts.clear();
    addToVec(emp_counts, other.emp_counts);
    fill(other.emp_counts.begin(), other.emp_counts.end(), Params{0,0,0});
    likelihood += other.likelihood;
    other.likelihood = 0;
    if (first_iter) {
      tot_len_ratio += other.tot_len_ratio;
      other.tot_len_ratio = 0;

      addToVec(toks, other.toks);
      fill(other.toks.begin(), other.toks.end(), uintpair{0,0});
      addToVec(size_index_counts, other.size_index_counts);
      fill(other.size_index_counts.begin(), other.size_index_counts.end(), tup_map{});
    }
  }

  void print(Dict &posdict);
  void count_toks();

 private:
  bool first_iteration = true;
};

struct Model {

  TTable ttable;
  std::vector<Params> params;
};

// get param or make new one if not present
template <typename T>
T &getSafe(std::vector<T> &vec, unsigned id, T fillval = T()) {
  // make new and return ref
  if (vec.size() <= id)
  {
    vec.resize(id+1, fillval);
    return vec[id];
  }
  //return ref
  return vec[id];
}

void normalize_counts(std::vector<Params> &vec, std::vector<std::pair<unsigned,unsigned>> &toks);
void normalize_sum_counts(std::vector<Params> &vec, std::vector<std::pair<unsigned,unsigned>> &toks);
void print_stats(std::vector<Params> &vec, std::string prefix, unsigned denom = 0);
void print_table(std::vector<Params> &vec, Dict &posdict, std::string prefix);

#endif