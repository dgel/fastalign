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

#include "model.h"
#include <ext/alloc_traits.h>
#include <cmath>                        // for log, pow
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for operator<<, basic_ostream, etc
#include "corpus.h"                     // for Dict

using namespace std;

void Stats::count_toks() {
  denom = 0;
  for (auto &x : toks) {
    denom += x.first;
    denom += x.second;
  }
}

void Stats::print(Dict &posdict) {
  // log(e) = 1.0
  double base2_likelihood = likelihood / log(2);

  if (first_iteration) {
    double mean_srclen_multiplier = tot_len_ratio / lc;
    cerr << "\nexpected target length = source length * "
         << mean_srclen_multiplier << endl;
    first_iteration = false;
  }

  cerr << "  log_e likelihood: " << likelihood << endl;
  cerr << "  log_2 likelihood: " << base2_likelihood << endl;
  cerr << "     cross entropy: " << (-base2_likelihood / denom) << endl;
  cerr << "        perplexity: " << pow(2.0, -base2_likelihood / denom) << endl;
  cerr << "      posterior p0: " << c0 / denom << endl;
  cerr << endl;
}

void normalize_counts(vector<Params> &vec, vector<pair<unsigned,unsigned>> &toks) {
  for (unsigned i = 3; i < vec.size(); ++i) {
    auto &param = vec[i];
    auto &toknum = toks[i];
    param.kappa /= toknum.first;
    param.lambda /= toknum.second;
    param.offset /= toknum.first + toknum.second;
  }
}

void normalize_sum_counts(vector<Params> &vec, vector<pair<unsigned,unsigned>> &toks) {
  for (unsigned i = 3; i < vec.size(); ++i) {
    auto &param = vec[i];
    auto &toknum = toks[i];
    unsigned toks = toknum.first + toknum.second;
    param.kappa = (param.kappa + param.lambda) / toks;
    param.lambda = param.kappa;
    param.offset /= toks;
  }
}


void print_stats(vector<Params> &vec, string prefix, unsigned denom)
{
  cerr << prefix << "  kappa (sum): ";
  double sum = 0;
  unsigned nPos = vec.size();

  for (unsigned i = 0; i < nPos; ++i)
  {
    sum += vec[i].kappa;
  }
  if (denom != 0)
    sum /= denom;
  cerr << sum << endl;

  cerr << prefix << " lambda (sum): ";
  sum = 0;
  for (unsigned i = 0; i < nPos; ++i)
  {
    sum += vec[i].lambda;
  }
  if (denom != 0)
    sum /= denom;
  cerr << sum << endl;

  cerr << prefix << " offset (sum): ";
  sum = 0;
  for (unsigned i = 0; i < nPos; ++i)
  {
    sum += vec[i].offset;
  }
  if (denom != 0)
    sum /= denom;
  cerr << sum << endl;
}

void print_table(vector<Params> &vec, Dict &posdict, string prefix)
{
    unsigned full_w = 10;
    unsigned nPos = vec.size();

    bool haveTags = posdict.max() != 3;

    cerr << "         === " << prefix << " ===" << endl;
    if (haveTags)
    {
      cerr << setw(full_w) << "POS";
    }
    cerr << setw(full_w) << "kappa" << " " << setw(full_w) << "lambda" << " " << setw(full_w) << "offset" << endl;
    for (unsigned i = 3; i < nPos; ++i)
    {
      if (haveTags)
      {
        cerr << setw(full_w) << posdict.Convert(i);
      }
      cerr << setw(full_w)  << vec[i].kappa << " " << setw(full_w) << vec[i].lambda << " " << setw(full_w) << vec[i].offset << endl;
    }
}