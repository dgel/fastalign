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

#include <iterator>
#include <fstream>
#include <memory>

#include "fast_align_ng.h"

 double Squared::learningrate_lk;
 double Squared::learningrate_o;
 double Absolute::learningrate_lk;
 double Absolute::learningrate_o;

using namespace std;

bool printProgress(size_t lc, bool flag) {
  if (lc % 1000 == 0) {
    cerr << '.';
    flag = true;
  }
  if (lc % 50000 == 0) {
    cerr << " [" << lc << "]" << endl;
    flag = false;
  }
  return flag;
}

void print_align(Options const &opts, vector<double> &probs, vector<unsigned> &src, unsigned j, OutputNode *out) {
  double max_p = -1;
  int max_index = -1;
  unsigned i = opts.use_null ? 0 : 1;
  for (; i <= src.size(); ++i) {
    if (probs[i] > max_p) {
      max_index = i;
      max_p = probs[i];
    }
  }
  if (max_index > 0) {
    string &outputbuf = out->out.back();
    if (opts.is_reverse)
       outputbuf += to_string(j) + "-" + to_string(max_index - 1) + " ";
    else
      outputbuf += to_string(max_index - 1) + "-" + to_string(j) + " ";
  }
}

void print_mat(vector<vector<double>> &probs){
  for (size_t j = 1; j < probs[0].size(); ++j)
  {
    for (size_t i = 0; i < probs.size(); ++i){
      cerr << setiosflags(ios::fixed) <<  setw(5) << setprecision(4) << log(probs[i][j]) << " ";
    }
    cerr << endl;
  }
  cerr << endl;
}

int main(int argc, char **argv) 
try
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
} catch (std::system_error &err) {
  cerr << "System Error!\n" << err.what() << endl;
}
