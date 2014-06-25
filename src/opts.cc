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

#include <tclap/CmdLine.h>

#include "opts.h"

using namespace TCLAP;
using namespace std;

void parseArgs(int argc, char **argv, Options &opts) {
  try {
    TCLAP::CmdLine cmd("fast-align", ' ', "0.2");

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
                           false, 7500, "unsigned", cmd);
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
    SwitchArg sDiag(
        "d", "no_favor_diagonal",
        "Don't favor alignment points close to the monotonic diagonal",
        cmd, true);
    SwitchArg sTens("o", "no_optimize_tension",
                    "Don't optimize how close to the diagonal alignment points "
                    "should be",
                    cmd, true);
    SwitchArg sVar("v", "no_variational_bayes",
                   "Don't use Dirichlet prior on lexical translation distributions ",
                   cmd, true);
    SwitchArg sNosplit("n", "no_split_params",
                   "Don't split the parameters about the diagonal",
                   cmd, false);
    SwitchArg sNull("N", "no_null_word", "No null word", cmd, false);
    SwitchArg sSquared("s", "squared_model", "Use the squared model", cmd, false);

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

    opts.startkappa = oKappa.getValue();
    opts.startlambda = oLambda.getValue();
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
    opts.squared_model = sSquared.getValue();
  }

  // handle errors in arguments
  catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
}
