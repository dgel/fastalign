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

#include <tclap/CmdLine.h>
#include <stdlib.h>                     // for exit
#include <iostream>                     // for operator<<, cerr, ostream
#include <string>                       // for string

struct Options {
  std::string input;           // corpus filename
  unsigned ITERATIONS;         // number of iterations
  unsigned n_threads;          // number of threads
  unsigned batch_size;         // number of lines for a thread to process at a time
  
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
  bool absolute_model;

  size_t n_no_update_diag;
  size_t n_no_update_offset;

  unsigned kNULL;
};


template <typename ModelType>
void parseArgs(int argc, char **argv, Options &opts) {
  using namespace TCLAP;
  using namespace std;
  try {
    CmdLine cmd("fast_align_ng", ' ', "0.3");

    // flags with argument
    ValueArg<string> oInput("i", "input", "Input parallel corpus", true, "",
                            "string", cmd);
    ValueArg<string> oPos("P", "POS", "Input POS data, for predicted side (trg by default)", false, "",
                            "string", cmd);
    ValueArg<size_t> oIter("I", "iterations",
                           "number of iterations in EM training (default = 5)",
                           false, 5, "unsigned", cmd);
    ValueArg<size_t> oThreads("", "parallel",
                           "number of threads to use (default = 1)",
                           false, 1, "unsigned", cmd);
    ValueArg<size_t> oBatchSize("", "batch_size",
                           "number of lines for a thread to process at a time",
                           false, 20000, "unsigned", cmd);
    ValueArg<double> oP0("p", "p0", "p_null parameter (default = 0.08)", false,
                         0.08, "double", cmd);
    ValueArg<double> oKappa(
        "K", "kappa",
        "starting Kappa for diagonal distance parameter (default = 6)", false,
        6, "double", cmd);
    ValueArg<double> oLambda(
        "L", "lambda",
        "starting lambda for diagonal distance parameter (default = 6)", false,
        6, "double", cmd);
    ValueArg<double> oAlpha(
        "a", "alpha",
        "alpha parameter for optional Dirichlet prior (default = 0.01)", false,
        0.01, "double", cmd);
    ValueArg<string> oCond("c", "conditional_probabilities",
                           "Output conditional probability table", false, "",
                           "string", cmd);
    ValueArg<size_t> oNUpDiag("", "n_no_update_diag",
                           "number of iterations during which diagonal is not updated",
                           false, 1, "unsigned", cmd);
    ValueArg<size_t> oNUpOffset("", "n_no_update_offset",
                           "number of iterations during which offset is not updated",
                           false, 1, "unsigned", cmd);
    // switches
    SwitchArg sRev(
        "r", "reverse",
        "Run alignment in reverse (condition on target and predict source)",
        cmd, false);
    SwitchArg sDiag("", "no_favor_diagonal",
        "Don't favor alignment points close to the monotonic diagonal",
        cmd, true);
    SwitchArg sTens("", "no_optimize_tension",
                    "Don't optimize how close to the diagonal alignment points "
                    "should be",
                    cmd, true);
    SwitchArg sVar("", "no_variational_bayes",
                   "Don't use Dirichlet prior on lexical translation distributions ",
                   cmd, true);
    SwitchArg sNosplit("", "no_split_params",
                   "Don't split the parameters about the diagonal",
                   cmd, false);
    SwitchArg sNull("N", "no_null_word", "No null word", cmd, false);
    SwitchArg sAbsolute("", "absolute_model", "Use the absolute model", cmd, false);

    // parse arguments
    cmd.parse(argc, argv);


    // extract values
    opts.input = oInput.getValue();
    opts.pos_filename = oPos.getValue();
    opts.ITERATIONS = oIter.getValue();
    opts.n_threads = oThreads.getValue();
    opts.batch_size = oBatchSize.getValue();

    opts.prob_align_null = oP0.getValue();
    opts.prob_align_not_null = 1. - opts.prob_align_null;

    opts.startkappa = oKappa.isSet() ? oKappa.getValue() : ModelType::init_kappa;
    opts.startlambda = oLambda.isSet() ? oLambda.getValue() : ModelType::init_lambda;
    opts.alpha = oAlpha.getValue();
    opts.conditional_probability_filename = oCond.getValue();

    opts.n_no_update_diag = oNUpDiag.getValue();
    opts.n_no_update_offset = oNUpOffset.getValue();

    opts.is_reverse = sRev.getValue();
    opts.favor_diagonal = sDiag.getValue();
    opts.optimize_tension = sTens.getValue();
    opts.variational_bayes = sVar.getValue();
    opts.use_null = !sNull.getValue();
    opts.nosplit = sNosplit.getValue();
    opts.absolute_model = sAbsolute.getValue();

    if (opts.variational_bayes && opts.alpha <= 0.0) {
      cerr << "--alpha must be > 0\n";
      exit(1);
    }

  }

  // handle errors in arguments
  catch (ArgException &e) {
    cerr << "HIER" << endl;

    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    exit(1);
  }
}

#endif