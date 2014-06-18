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


#ifndef _FAST_ALIGN_H_
#define _FAST_ALIGN_H_

#include <utility>
#include <memory>
#include <fstream>
#include <unordered_map>

#include "da.h"
#include "opts.h"

// forward declare 
bool printProgress(size_t lc, bool flag);
bool print_align(Options const &opts, std::vector<double> &probs, std::vector<unsigned> &src, unsigned j,
                 bool first_al);
void print_mat(std::vector<std::vector<double>> &probs);

template <typename ModelType>
double update_lk(double emp, double mod, double orig)
{
  orig += (mod - emp) * ModelType::learningrate_lk;
  if (orig < -10) orig = -10;
  if (orig > 30) orig = 30;
  return orig;
}

template <typename ModelType>
double update_o(double emp, double mod, double orig)
{
  orig += (mod - emp) * ModelType::learningrate_o;
  if (orig < -1) orig = -1;
  if (orig > 1) orig = 1;
  return orig;
}



template <typename ModelType>
double processWord(Model &model, Stats &stats, const Options &opts, std::vector<double> &probs,
                   std::vector<unsigned> &src, std::vector<unsigned> &trg, unsigned j, unsigned posNum,
                   bool final_iteration) {
  auto const &params = getSafe(model.params, posNum, Params{opts.startkappa, opts.startlambda, 0});
  auto &emp_counts = getSafe(stats.emp_counts, posNum);

  const unsigned f_j = trg[j];
  double sum = 0;
  double prob_a_i =
      1.0 / (src.size() + (opts.use_null ? 1 : 0));  // uniform (model 1)

  // compute probability of alignment positions (including null if
  // applicable)
  if (opts.use_null) {
    if (opts.favor_diagonal) prob_a_i = opts.prob_align_null;
    probs[0] = model.ttable.prob(opts.kNULL, f_j) * prob_a_i;
    sum += probs[0];
  }


  // find diagonal -- i and j have switched meaning compared to paper
    //                  same for m and n
  unsigned diag = 0;
  if (opts.favor_diagonal) {
    diag = DiagonalAlignment::Diagonal(j+1, trg.size(), src.size(), params.offset);
  }
  // need 2 loops for favor_diagonal case
  if (opts.favor_diagonal) {
    double const az = ModelType::ComputeZ(j + 1, trg.size(), src.size(),
                                     params.kappa, params.lambda, params.offset) /
         opts.prob_align_not_null;


    // up to diag uses kappa, beyond uses lambda;
    for (unsigned i = 1; i <= diag; ++i) {
      prob_a_i = ModelType::UnnormalizedProb(
                       j + 1, i, trg.size(), src.size(), params.kappa, params.offset) /
                   az;
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      sum += probs[i];
    }
    for (unsigned i = diag + 1; i <= src.size(); ++i) {
      prob_a_i = ModelType::UnnormalizedProb(
                       j + 1, i, trg.size(), src.size(), params.lambda, params.offset) /
                   az;
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      sum += probs[i];
    }
  // can keep single loop if not
  } else {
    for (unsigned i = 1; i <= src.size(); ++i) {
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      sum += probs[i];
    }
  }
 
  // update counts if necessary
  if (!final_iteration) {
    double p = 0;
    // add to counts
    if (opts.use_null) {
      p = probs[0] / sum;
      stats.c0 += p;
      stats.tcounts.Increment(opts.kNULL, f_j, p);
    }
    double tmpoffset = 0;
    // empirical counts
    for (unsigned i = 1; i <= diag; ++i) {
      p = probs[i] / sum;
      stats.tcounts.Increment(src[i - 1], f_j, p);
      double x = ModelType::Feature(j+1, i, trg.size(), src.size(), params.offset);
      emp_counts.kappa += ModelType::Transform(x) * p;
      tmpoffset += ModelType::dTransform(x) * p;
      //std::cerr << "position " << j << ", " << i << " prob: " << p << " (tprob: " << model.ttable.prob(src[i - 1], f_j) << ") feature: " << x << " dTransform: " << ModelType::dTransform(x) << std::endl;
    }
    emp_counts.offset += tmpoffset * params.kappa;
    tmpoffset = 0;
    for (unsigned i = diag + 1; i <= src.size(); ++i) {
      p = probs[i] / sum;
      stats.tcounts.Increment(src[i - 1], f_j, p);
      double x = ModelType::Feature(j+1, i, trg.size(), src.size(), params.offset);
      emp_counts.lambda += ModelType::Transform(x) * p;
      tmpoffset += ModelType::dTransform(x) * p;
      //std::cerr << "position " << j << ", " << i << " prob: " << p << " (tprob: " << model.ttable.prob(src[i - 1], f_j) << ") feature: "<< x << " dTransform: " << ModelType::dTransform(x) << std::endl;
    }
    emp_counts.offset += tmpoffset * params.lambda;
    //std::cerr << "posterior offset: " << of << std::endl;
  }
  return sum;
}


template <typename ModelType>
void count(std::istream &in, Dict &d, std::unique_ptr<std::ifstream>  &posdata, Dict &posdict,
           Model &model, Stats &stats, const Options &opts, unsigned iter) {
  std::string line;
  std::vector<double> probs;

  std::string posline;
  std::vector<unsigned> src, trg, posint;
  bool flag = false;

  bool first_iter = iter == 0;
  const bool final_iteration = (iter == (opts.ITERATIONS - 1));

  while (std::getline(in, line)) {
    flag = printProgress(++stats.lc, flag);
    src.clear();
    trg.clear();
    d.ParseLine(line, src, trg);
    posint.clear();

    if (opts.is_reverse) swap(src, trg);
    if (src.empty() || trg.empty()) {
      std::cerr << "Error in line " << stats.lc << "\n" << line << std::endl;
      exit(1);
    }
    if (posdata != nullptr) {
      if (!std::getline(*posdata, posline)) {
        std::cerr << "POS file has too few lines" << std::endl;
        std::exit(1);
      }
      size_t delimpos = posline.find("|||");

      if (opts.is_reverse) {
        posdict.ConvertWhitespaceDelimitedLine(posline, posint, 0, delimpos - 1);
      } else {
        posdict.ConvertWhitespaceDelimitedLine(posline, posint, delimpos + 4);
      }
      if (trg.size() != posint.size()) {
        std::cerr << "POS line has fewer tags than trg line has words.\n"
                "Perhaps src and trg should be swapped?"
             << std::endl;
        exit(1);
      }
    }

    if (first_iter)
      stats.tot_len_ratio +=
          static_cast<double>(trg.size()) / static_cast<double>(src.size());

    probs.resize(src.size() + 1);

    bool first_al = true;  // used when printing alignments

    unsigned posNum = 3;
    // for every target location
    for (unsigned j = 0; j < trg.size(); ++j) {
      if (posdata != nullptr)
      {
        posNum = posint[j];
      }

      if (first_iter)
      {
        // std::cerr << posNum << " " << posdict.Convert(posNum) << std::endl;
        auto &curtoks = getSafe(stats.toks, posNum);
        auto const &params = getSafe(model.params, posNum, Params{opts.startkappa, opts.startlambda, 0});
        auto diag = DiagonalAlignment::Diagonal(j+1, trg.size(), src.size(), params.offset);
        curtoks.first += diag;
        curtoks.second += src.size() - diag;

        getSafe(stats.size_index_counts, posNum)[std::tuple<unsigned short, unsigned short, unsigned short>(trg.size(), src.size(), j+1)] += 1;
      }

      // find probability for every position
      double sum =
          processWord<ModelType>(model, stats, opts, probs, src, trg, j, posNum, final_iteration);
      if (final_iteration) {
        // print alignments
        // /cerr << "printing an alignment" << std::endl;
        first_al = print_align(opts, probs, src, j, first_al);
      }
      stats.likelihood += log(sum);
    }
    if (final_iteration) std::cout << std::endl;
  }
  if (flag) {
    std::cerr << std::endl;
  }
}

template <typename ModelType>
void update(Model &model, Stats &stats, Dict &posdict, Options opts, std::istream * const posdata, size_t iter) {
  bool update_tension = opts.favor_diagonal && opts.optimize_tension && (iter >= opts.n_no_update_diag || iter >= opts.n_no_update_offset);
  if (update_tension) {


    // keep kappa and lambda same for all values
    if (opts.nosplit)
    {
      normalize_sum_counts(stats.emp_counts, stats.toks);
    } else {
      normalize_counts(stats.emp_counts, stats.toks);
    }

    ModelType::learningrate_lk = ModelType::init_learningrate_lk;
    ModelType::learningrate_o = ModelType::init_learningrate_o;

    // 8 update steps
    for (int ii = 0; ii < 8; ++ii) {
      std::fill(stats.mod_counts.begin(), stats.mod_counts.end(), Params{0,0,0});

      for (size_t ind = 3; ind < stats.size_index_counts.size(); ++ind)
      {
        auto const &params = getSafe(model.params, ind);
        auto &mod_counts = getSafe(stats.mod_counts, ind);
        for (auto &pair : stats.size_index_counts[ind])
        {
          auto &tup = pair.first;
          unsigned count = pair.second;
          //std::cerr << "sent pos: " << std::get<2>(tup) << std::endl;
          auto dzs = ModelType::ComputeDLogZs(std::get<2>(tup)
                                      , std::get<0>(tup), std::get<1>(tup), params.kappa, params.lambda, params.offset);
          mod_counts.kappa += dzs.kappa * count;
          mod_counts.lambda += dzs.lambda * count;
          mod_counts.offset += dzs.offset * count;
        }
      }

      print_stats(stats.mod_counts, boost::lexical_cast<std::string>(ii) + " model al-", stats.denom);
      

      if (opts.nosplit)
      {
        normalize_sum_counts(stats.mod_counts, stats.toks);
      } else {
        normalize_counts(stats.mod_counts, stats.toks);
      }

      unsigned nPos = stats.mod_counts.size();
      for (unsigned i = 3; i < nPos; ++i)
      {
        auto &params = model.params[i];
        auto const &emp = stats.emp_counts[i];
        auto const &mod = stats.mod_counts[i];

        // std::cerr << "updating parameter for POS: " << posdict.Convert(i) << std::endl;
        // std::cerr << "original: " << params.kappa << " " << params.lambda << " " << params.offset << std::endl;
        // std::cerr << "empirical: " << emp.kappa << " " << emp.lambda << " " << emp.offset << std::endl;
        // std::cerr << "model: " << mod.kappa << " " << mod.lambda << " " << mod.offset << std::endl;
        if (iter >= opts.n_no_update_diag) {
          params.kappa = update_lk<ModelType>(emp.kappa, mod.kappa, params.kappa);
          params.lambda = update_lk<ModelType>(emp.lambda, mod.lambda, params.lambda);
        }
        if (iter >= opts.n_no_update_offset)
          params.offset = update_o<ModelType>(emp.offset, mod.offset, params.offset);
      }
      ModelType::learningrate_lk *= 0.9;
      ModelType::learningrate_o *= 0.9;
      // std::cerr << std::endl;
      // print_table(model.params, posdict, "  final ");
      // std::cerr << std::endl;
    }
    
    std::cerr << std::endl;
    print_table(model.params, posdict, "  final ");
    std::cerr << std::endl;
    
  }
  if (opts.variational_bayes)
    model.ttable.UpdateFromVB(stats.tcounts, opts.alpha);
  else
    model.ttable.UpdateFrom(stats.tcounts);
}

template <typename ModelType>
int run(Options &opts) {
  Dict d;  // integerization map
  Dict posdict;

  opts.kNULL = d.kNULL_;
  Model model;
  Stats stats;

  for (unsigned iter = 0; iter < opts.ITERATIONS; ++iter) {
    const bool final_iteration = (iter == (opts.ITERATIONS - 1));
    std::cerr << "ITERATION " << (iter + 1) << (final_iteration ? " (FINAL)" : "")
         << std::endl;
    std::ifstream in(opts.input);
    if (!in) {
      std::cerr << "Can't read " << opts.input << std::endl;
      return 1;
    }

    // gets deleted automatically
    std::unique_ptr<std::ifstream> posdata;
    if (opts.pos_filename != "")
    {
      posdata = std::unique_ptr<std::ifstream>(new std::ifstream(opts.pos_filename));
      if (!(*posdata)) {
        std::cerr << "Can't read " << opts.pos_filename << std::endl;
        return 1;
      }
    } else {
      posdict.Convert("Default");
    }

    stats.reset();

    // process current iteration
    // collect counts
    count<ModelType>(in, d, posdata, posdict, model, stats, opts, iter);
    // make update if necessary
    stats.count_toks();
    stats.print(posdict);
    print_stats(stats.emp_counts, "     posterior ", stats.denom);
    if (!final_iteration) {
      update<ModelType>(model, stats, posdict, opts, posdata.get(), iter);
    }
  }
  // after last iteration, write probabilities
  if (!opts.conditional_probability_filename.empty()) {
    std::cerr << "conditional probabilities: " << opts.conditional_probability_filename
         << std::endl;
    model.ttable.ExportToFile(opts.conditional_probability_filename.c_str(), d);
  }
  return 0;
}

#endif
