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

size_t max_align(Options const &opts, std::vector<double> &probs);
void print_align(Options const &opts, Token i, Token j, OutputNode *out);

template <typename ModelType>
double update(double emp, double mod, double orig,
              double learnrate, double min, double max, double reg_param, double reg_target)
{
  orig += (mod - emp) * learnrate + reg_param * (reg_target - orig);
  if (orig < min) orig = min;
  if (orig > max) orig = max;
  return orig;
}

template <typename ModelType>
double processWord(Model &model, Stats &stats, const Options &opts, std::vector<double> &probs,
                   std::vector<unsigned> &src, std::vector<unsigned> &trg, unsigned j, Token posNum,
                   bool final_iteration) {
  auto const &params = getSafe(model.params, posNum, Params{opts.startkappa, opts.startlambda, 0});
  auto &emp_counts = getSafe(stats.emp_counts, posNum);

  const Token f_j = trg[j];
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
  size_t diag = 0;
  // need 2 loops for favor_diagonal case
  if (opts.favor_diagonal) {
    diag = DiagonalAlignment::Diagonal(j+1, trg.size(), src.size(), params.offset);
    double az = 0;
    for (size_t i = 1; i <= diag; ++i) {
      prob_a_i = ModelType::UnnormalizedProb(
                       j + 1, i, trg.size(), src.size(), params.kappa, params.offset);
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      az += prob_a_i;
    }
    for (size_t i = diag + 1; i <= src.size(); ++i) {
      prob_a_i = ModelType::UnnormalizedProb(
                       j + 1, i, trg.size(), src.size(), params.lambda, params.offset);
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      az += prob_a_i;
    }
    az = opts.prob_align_not_null / az;
    for (size_t i = 1; i <= src.size(); ++i) {
      probs[i] *= az;
      sum += probs[i];
    }
  // can keep single loop if not
  } else {
    for (size_t i = 1; i <= src.size(); ++i) {
      probs[i] = model.ttable.prob(src[i - 1], f_j) * prob_a_i;
      sum += probs[i];
    }
  }

  // update counts if necessary
  if (!final_iteration) {
    double p;
    // add to counts
    if (opts.use_null) {
      p = probs[0] / sum;
      stats.c0 += p;
      stats.tcounts.Increment(opts.kNULL, f_j, p);
    }
    double tmpoffset = 0;
    tmpoffset = 0;
    double tmp;
    // empirical counts
    for (size_t i = 1; i <= diag; ++i) {
      p = probs[i] / sum;
      stats.tcounts.Increment(src[i - 1], f_j, p);
      // optimization that happens to work for both modeltypes.
      double x = ModelType::Feature(j+1, i, trg.size(), src.size(), params.offset);
      tmp = ModelType::dTransform(x) * p;
      emp_counts.lambda += x * tmp ;
      tmpoffset += tmp;
    }
    emp_counts.offset += tmpoffset * params.kappa;

    tmpoffset = 0;
    for (size_t i = diag + 1; i <= src.size(); ++i) {
      p = probs[i] / sum;
      stats.tcounts.Increment(src[i - 1], f_j, p);
      double x = ModelType::Feature(j+1, i, trg.size(), src.size(), params.offset);
      // optimization that happens to work for both modeltypes.
      tmp = ModelType::dTransform(x) * p;
      emp_counts.lambda += x * tmp ;
      tmpoffset += tmp;
    }

    emp_counts.offset += tmpoffset * params.lambda;

  }
  return sum;
}

template <typename ModelType>
void count(PairLines &in, PosLines &posdata, Model &model, Stats &stats, const Options &opts, size_t iter, OutputNode *out) {
  std::vector<double> probs;

  bool first_iter = iter == 0;
  bool final_iteration = iter == (opts.ITERATIONS - 1);

  size_t lineno = 0;
  PosLines::iterator positer;
  if (!posdata.empty()) {
    positer = posdata.begin();
  }

  for (std::pair<std::vector<Token>, std::vector<Token>> &lines : in) {
    std::vector<Token> &src = std::get<0>(lines);
    std::vector<Token> &trg = std::get<1>(lines);

    if (first_iter)
      stats.tot_len_ratio +=
          static_cast<double>(trg.size()) / static_cast<double>(src.size());

    probs.resize(src.size() + 1);

    if (final_iteration)
      out->out.emplace_back();

    Token posNum = 3;
    // for every target location
    for (size_t j = 0; j < trg.size(); ++j) {
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

        stats.get_size_index_count(posNum)[std::tuple<unsigned short, unsigned short, unsigned short>(trg.size(), src.size(), j+1)] += 1;
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
void count_single_thread(Reader &reader, ThreadedOutput &outp, Model &model, Stats &stats, const Options &opts, size_t iter) {
  std::vector<double> probs;

  bool first_iter = iter == 0;
  bool final_iteration = iter == (opts.ITERATIONS - 1);

  PairLine lines;
  Line posdata;

  OutputNode *out = nullptr;

  size_t lineno = 0;

  while (reader.read_1_line(lines, posdata)) {
    if (final_iteration && lineno % opts.batch_size == 0) {
      if (out) {
        out->done = true;
      }
      out = outp.getOutput();
    }

    std::vector<Token> &src = std::get<0>(lines);
    std::vector<Token> &trg = std::get<1>(lines);

    if (first_iter)
      stats.tot_len_ratio +=
          static_cast<double>(trg.size()) / static_cast<double>(src.size());

    probs.resize(src.size() + 1);

    if (final_iteration)
      out->out.emplace_back();

    Token posNum = 3;
    // for every target location
    for (size_t j = 0; j < trg.size(); ++j) {
      if (!posdata.empty())
      {
        posNum = posdata[j];
      }

      if (first_iter)
      {
        // std::cerr << posNum << " " << posdict.Convert(posNum) << std::endl;
        auto &curtoks = getSafe(stats.toks, posNum);
        auto const &params = getSafe(model.params, posNum, Params{opts.startkappa, opts.startlambda, 0});
        auto diag = DiagonalAlignment::Diagonal(j+1, trg.size(), src.size(), params.offset);
        curtoks.first += diag;
        curtoks.second += src.size() - diag;

        stats.get_size_index_count(posNum)[std::tuple<unsigned short, unsigned short, unsigned short>(trg.size(), src.size(), j+1)] += 1;
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
  }
  if (out) {
    out->done = true;
  }
}




template <typename ModelType>
void update(Model &model, Stats &stats, Dict &posdict, Options opts, size_t iter) {

  bool update_tension = opts.favor_diagonal
                      && opts.optimize_tension
                      && (iter >= opts.n_no_update_diag || iter >= opts.n_no_update_offset);
  if (update_tension) {


    // keep kappa and lambda same for all values
    if (opts.nosplit)
    {
      normalize_sum_counts(stats.emp_counts, stats.toks);
    } else {
      normalize_counts(stats.emp_counts, stats.toks);
    }

    double learningrate_lk = ModelType::init_learningrate_lk;
    double learningrate_o = ModelType::init_learningrate_o;

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
            //double tmp = 0;
            // loop over key,value pairs in bucket
            for (auto iter = map.begin(i), end = map.end(i); iter != end; ++iter) {
              auto &tup = iter->first;
              size_t count = iter->second;

              auto dzs = ModelType::ComputeDLogZs(std::get<2>(tup)
                                          , std::get<0>(tup), std::get<1>(tup), params.kappa, params.lambda, params.offset);
              accum.kappa += dzs.kappa * count;
              accum.lambda += dzs.lambda * count;
              accum.offset += dzs.offset * count;
              //tmp += dzs.offset * count;
            }
            //std::cerr << "intermediate: " << tmp << std::endl;

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
                                            learningrate_lk, ModelType::min_lk, ModelType::max_lk, 0.05, 80);
           params.lambda = update<ModelType>(emp.lambda, mod.lambda, params.lambda,
                                             learningrate_lk, ModelType::min_lk, ModelType::max_lk, 0.05, 80);
        }
        if (iter >= opts.n_no_update_offset)
          params.offset = update<ModelType>(emp.offset, mod.offset, params.offset,
                                            learningrate_o, ModelType::min_o, ModelType::max_o, 0.03, 0);
      }
      learningrate_lk *= 0.9;
      learningrate_o *= 0.9;
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
    auto status = reader.read_n_lines_threaded(lines, pos, N, outp, out);
    if (status == ReadStatus::NOTHING_READ) {
      return;
    } else if (status == ReadStatus::FINISHED) {
      done = true;
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

    if (opts.n_threads > 0) {
      run_n(opts.n_threads, countWorker<ModelType>, reader, outp, model,
            stats, statsmut, opts, iter, opts.batch_size);
    } else {
      count_single_thread<ModelType>(reader, outp, model, stats, opts, iter);
    }

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
    model.ttable.ExportToFile(opts.conditional_probability_filename, reader.dict);
  }
  outputthread.join();
  return 0;
}

#endif
