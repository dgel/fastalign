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

#include <math.h>                       // for log
#include <algorithm>                    // for fill
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <mutex>                        // for mutex, lock_guard
#include <string>                       // for allocator, operator+, etc
#include <thread>                       // for thread
#include <utility>                      // for get, pair
#include <vector>                       // for vector, vector<>::iterator
#include "da.h"                         // for Diagonal
#include "model.h"                      // for Stats, Params, Model, etc
#include "opts.h"                       // for Options
#include "reader.h"                     // for PosLines, PairLines, Reader
#include "threaded_output.h"            // for OutputNode, ThreadedOutput
#include "ttables.h"                    // for TTable
#include "utils.h"                      // for par_map, n_ranges, etc
class Dict;

unsigned max_align(Options const &opts, std::vector<double> &probs);
void print_align(Options const &opts, unsigned i, unsigned j, OutputNode *out);

template <typename ModelType>
double update(double emp, double mod, double orig,
              double learnrate, double min, double max)
{
  orig += (mod - emp) * learnrate;
  if (orig < min) orig = min;
  if (orig > max) orig = max;
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
    }
    emp_counts.offset += tmpoffset * params.kappa;
    tmpoffset = 0;
    for (unsigned i = diag + 1; i <= src.size(); ++i) {
      p = probs[i] / sum;
      stats.tcounts.Increment(src[i - 1], f_j, p);
      double x = ModelType::Feature(j+1, i, trg.size(), src.size(), params.offset);
      emp_counts.lambda += ModelType::Transform(x) * p;
      tmpoffset += ModelType::dTransform(x) * p;
    }
    emp_counts.offset += tmpoffset * params.lambda;
  }
  return sum;
}

template <typename ModelType>
void count(PairLines &in, PosLines &posdata, Model &model, Stats &stats, const Options &opts, unsigned iter, OutputNode *out) {
  std::vector<double> probs;

  bool first_iter = iter == 0;
  bool final_iteration = iter == (opts.ITERATIONS - 1);

  size_t lineno = 0;
  PosLines::iterator positer;
  if (!posdata.empty()) {
    positer = posdata.begin();
  }

  for (std::pair<std::vector<unsigned>, std::vector<unsigned>> &lines : in) {
    std::vector<unsigned> &src = std::get<0>(lines);
    std::vector<unsigned> &trg = std::get<1>(lines);

    if (first_iter)
      stats.tot_len_ratio +=
          static_cast<double>(trg.size()) / static_cast<double>(src.size());

    probs.resize(src.size() + 1);
    out->out.emplace_back();

    unsigned posNum = 3;
    // for every target location
    for (unsigned j = 0; j < trg.size(); ++j) {
      if (!posdata.empty())
      {
        posNum = (*positer)[j];
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
        print_align(opts, max_align(opts, probs),  j, out);
      }
      stats.likelihood += log(sum);
    }
    ++lineno;
    if (!posdata.empty()) {
      ++positer;
    }
  }
}

template <typename ModelType>
void update(Model &model, Stats &stats, Dict &posdict, Options opts, size_t iter) {
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
        auto const& map = stats.size_index_counts[ind];
        // assign each thread a range of buckets
        std::vector<std::pair<size_t, size_t>> ranges = n_ranges(map.bucket_count(), opts.n_threads);
        // sum DlogZs in each range in parallel
        auto vals = par_map(ranges, [&](std::pair<size_t, size_t> &range) {
          Params accum = {0,0,0};
          // loop over buckets
          for (size_t i = range.first; i < range.second; ++i) {
            // loop over key,value pairs in bucket
            for (auto iter = map.cbegin(i), end = map.cend(i); iter != end; ++iter) {
              auto &tup = iter->first;
              unsigned count = iter->second;
              
              auto dzs = ModelType::ComputeDLogZs(std::get<2>(tup)
                                          , std::get<0>(tup), std::get<1>(tup), params.kappa, params.lambda, params.offset);
              accum.kappa += dzs.kappa * count;
              accum.lambda += dzs.lambda * count;
              accum.offset += dzs.offset * count;
            }
          }
          return accum;
        });

        // add to relevant mod_counts
        for (Params &p: vals) {
          mod_counts += p;
        }
      }

      print_stats(stats.mod_counts, std::to_string(ii) + " model al-", stats.denom);
      
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

        if (iter >= opts.n_no_update_diag) {
          params.kappa = update<ModelType>(emp.kappa, mod.kappa, params.kappa, 
                                           ModelType::learningrate_lk, ModelType::min_lk, ModelType::max_lk);
          params.lambda = update<ModelType>(emp.lambda, mod.lambda, params.lambda,
                                            ModelType::learningrate_lk, ModelType::min_lk, ModelType::max_lk);
        }
        if (iter >= opts.n_no_update_offset)
          params.offset = update<ModelType>(emp.offset, mod.offset, params.offset,
                                            ModelType::learningrate_o, ModelType::min_o, ModelType::max_o);
      }
      ModelType::learningrate_lk *= 0.9;
      ModelType::learningrate_o *= 0.9;
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
void countWorker(Reader &reader, ThreadedOutput &outp, Model &model, 
                 Stats &stats, std::mutex &statsmut, Options &opts, unsigned iter, unsigned N) {
  PairLines lines;
  PosLines pos;
  Stats localstats;
  bool done = false;
  // collect counts N lines at a time
  while (!done) {
    localstats.reset();
    OutputNode *out;
    done = reader.read_n_lines_threaded(lines, pos, N, outp, out);

    if (lines.size() == 0) {
      out->done = true;
      return;
    }
    count<ModelType>(lines, pos, model, localstats, opts, iter, out);
    out->done = true;
    {
      std::lock_guard<std::mutex> statslock(statsmut);
      stats.add_count(localstats, iter == 0);
    }
  }
}

template <typename ModelType>
int run(Options &opts) {
  Model model;

  Stats stats;
  std::mutex statsmut;
  
  Reader reader(opts.input, opts.pos_filename, stats, opts.is_reverse);
  opts.kNULL = reader.dict.kNULL_;

  ThreadedOutput outp(reader.iomut);
  
  std::thread outputthread(handleOutput, std::ref(outp), 5);

  for (unsigned iter = 0; iter < opts.ITERATIONS; ++iter) {
    const bool final_iteration = (iter == (opts.ITERATIONS - 1));
    std::cerr << "ITERATION " << (iter + 1) << (final_iteration ? " (FINAL)" : "")
         << std::endl;
    
    stats.reset();
    reader.rewind();

    run_n(opts.n_threads, countWorker<ModelType>, reader, outp, model, 
          stats, statsmut, opts, iter, opts.batch_size);
    
    // make update if necessary
    stats.count_toks();
    stats.print(reader.posdict);
    print_stats(stats.emp_counts, "     posterior ", stats.denom);
    if (!final_iteration) {
      update<ModelType>(model, stats, reader.posdict, opts, iter);
    }
  }
  outp.finished();

  // after last iteration, write probabilities
  if (!opts.conditional_probability_filename.empty()) {
    std::cerr << "conditional probabilities: " << opts.conditional_probability_filename
         << std::endl;
    model.ttable.ExportToFile(opts.conditional_probability_filename.c_str(), reader.dict);
  }
  outputthread.join();
  return 0;
}

#endif
