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

#ifndef _OPTS_H_
#define _OPTS_H_

#include <string>

struct Options {
  std::string input;           // corpus filename
  unsigned ITERATIONS;         // number of iterations
  unsigned n_threads;          // number of threads
  double alpha;                // parameter for variational bayes
  double prob_align_null;      // null alignment prob
  double prob_align_not_null;  // and complement
  std::string conditional_probability_filename;  // filename to output
                                                 // probability
                                                 // table
  std::string pos_filename;  // filename for pos-tags. same format as corpus,
                             // can
                             // leave line empty.
  double startkappa;         // lambda parameter for delta function from article
  double startlambda;

  bool is_reverse;        // swap source and target
  bool favor_diagonal;    // use parameterized version of Model 2 (otherwise use
                          // Model 1)
  bool optimize_tension;  // optimize lambda parameter from article
  bool variational_bayes;  // use variational bayes
  bool use_null;           // use null alignment
  bool nosplit;
  bool squared_model;

  size_t n_no_update_diag;
  size_t n_no_update_offset;

  unsigned kNULL;
};

void parseArgs(int argc, char **argv, Options &opts);

#endif