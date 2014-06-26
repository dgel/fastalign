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

#include "fast_align_ng.h"
#include <ext/alloc_traits.h>
#include <system_error>                 // for system_error
#include <iomanip>                      // for operator<<, setw

using namespace std;

unsigned max_align(Options const &opts, vector<double> &probs) {
  double max_p = -1;
  unsigned max_index = 0;
  unsigned i = opts.use_null ? 0 : 1;
  for (; i < probs.size(); ++i) {
    if (probs[i] > max_p) {
      max_index = i;
      max_p = probs[i];
    }
  }
  return max_index;
}

void print_align(Options const &opts, unsigned i, unsigned j, OutputNode *out) {
  if (i > 0) {
    string &outputbuf = out->out.back();
    if (opts.is_reverse)
       outputbuf += to_string(j) + "-" + to_string(i - 1) + " ";
    else
      outputbuf += to_string(i - 1) + "-" + to_string(j) + " ";
  }
}

int main(int argc, char **argv) 
{
  ios_base::sync_with_stdio(false);

  Options opts;
  parseArgs(argc, argv, opts);

  if (opts.variational_bayes && opts.alpha <= 0.0) {
    cerr << "--alpha must be > 0\n";
    return 1;
  }

  if (opts.squared_model) {
    // most likely doesn't work atm
    return run<Squared>(opts);
  } else {
    return run<Absolute>(opts);
  }
}
