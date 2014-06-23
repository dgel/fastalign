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
#include <queue>
#include <thread>
#include <mutex>
#include <chrono>

#include "da.h"
#include "opts.h"

struct OutputNode {
  bool done = false;
  std::vector<std::string> out;
};

struct ThreadedOutput {
  bool done = false;
  std::queue<OutputNode> que;
};

typedef std::vector<std::pair<std::vector<unsigned>,std::vector<unsigned>>> PairLines;
typedef std::vector<std::vector<unsigned>> PosLines;

class Reader {
  std::ifstream data;
  std::unique_ptr<std::ifstream> posdata;
  Stats &stats;

  bool is_reverse;

public:
  Dict dict;
  Dict posdict;

  Reader(std::string &datafile, std::string &posfile, Stats &stats, bool is_reverse) 
  : data(datafile), stats(stats), is_reverse(is_reverse)
  {
    if (!data) {
      std::cerr << "Can't read " << datafile << std::endl;
      exit(1);
    }

    // gets deleted automatically
    if (posfile != "")
    {
      posdata = std::unique_ptr<std::ifstream>(new std::ifstream(posfile));
      if (!(*posdata)) {
        std::cerr << "Can't read " << posfile << std::endl;
        exit(1);
      }
    } else {
      posdict.Convert("Default");
    }
  }



  bool readNlines(PairLines &lines, PosLines &poslines, size_t N) {
    std::cerr.flush();
    lines.clear();
    poslines.clear();

    std::string line;
    size_t nread = 0;
    while (std::getline(data, line)) {
      ++stats.lc;
      lines.emplace_back();
      if (is_reverse) {
        dict.ParseLine(line, lines.back().second, lines.back().first);
      } else {
        dict.ParseLine(line, lines.back().first, lines.back().second);
      }

      if(lines.back().first.empty() || lines.back().second.empty()) {
        std::cerr << "Error in line " << stats.lc << "\n" << line << std::endl;
        exit(1);
      }
      if (posdata != nullptr) {
        if (!std::getline(*posdata, line)) {
          std::cerr << "POS file has too few lines" << std::endl;
          std::exit(1);
        }

        size_t delimpos = line.find("|||");
        poslines.emplace_back();

        if (is_reverse) {
          posdict.ConvertWhitespaceDelimitedLine(line, poslines.back(), 0, delimpos - 1);
        } else {
          posdict.ConvertWhitespaceDelimitedLine(line, poslines.back(), delimpos + 4);
        }
        if (lines.back().second.size() != poslines.back().size()) {
          std::cerr << "POS line has fewer tags than trg line has words.\n"
                       "Perhaps src and trg should be swapped?"
               << std::endl;
          exit(1);
        }    
      }
      if (++nread == N)
        break;
    }

    // signal if done processing corpus
    if (nread < N) {
      return true;
    }
    return false;
  }

  void rewind() {
    data.clear();
    data.seekg(0);
    if (posdata) {
      posdata->clear();
      posdata->seekg(0);
    }
  }
};

bool printProgress(size_t lc, bool flag);
void print_align(Options const &opts, std::vector<double> &probs, std::vector<unsigned> &src, unsigned j, OutputNode *out);
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
  std::cerr << "counting " << in.size() << " sentences" << std::endl;
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
        // print alignments
        print_align(opts, probs, src, j, out);
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

void handleOutput(ThreadedOutput &out, std::mutex &mut, size_t sleeptime) {
  while (true) {
    if (mut.try_lock()) {
      while (!out.que.empty() && out.que.front().done) {
        auto &head = out.que.front().out;
        for (std::string &line : head) {
          if (!line.empty())
            line.pop_back();
          std::cout << line << std::endl;
        }
        std::cerr << std::endl;
        out.que.pop();
      }
      mut.unlock();
      if (out.done) {
        return;
      }
    } 
    std::this_thread::sleep_for (std::chrono::seconds(sleeptime));
  }
}



template <typename F, typename... Args>
void runN(size_t N, F func, Args &... args) {
  std::vector<std::thread> threads;
  for (size_t i = 0; i < N; ++i) {
    threads.emplace_back(func, std::ref(args)...);
  }
  for (std::thread &t : threads) {
    t.join();
  }
}

template <typename ModelType>
void countWorker(Reader &reader, ThreadedOutput &outp, std::mutex &iomut, Model &model, Stats &stats, std::mutex &statsmut, Options &opts, unsigned iter) {
  PairLines lines;
  PosLines pos;
  Stats localstats;
  bool done = false;
  // collect counts 1000 lines at a time
  while (!done) {
    localstats.reset();
    OutputNode *out;
    // create scope for lock
    {
      std::lock_guard<std::mutex> io_lock(iomut);
      done = reader.readNlines(lines, pos, 3500);
      if (lines.size() == 0)
        return;
      outp.que.emplace();
      out = &outp.que.back();
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
  ThreadedOutput outp;
  std::mutex iomut;
  
  std::thread outputthread(handleOutput, std::ref(outp), std::ref(iomut), 30);

  for (unsigned iter = 0; iter < opts.ITERATIONS; ++iter) {
    const bool final_iteration = (iter == (opts.ITERATIONS - 1));
    std::cerr << "ITERATION " << (iter + 1) << (final_iteration ? " (FINAL)" : "")
         << std::endl;
    
    stats.reset();
    reader.rewind();

    runN(opts.n_threads, countWorker<ModelType>, reader, outp, iomut, model, stats, statsmut, opts, iter);
    
    // make update if necessary
    stats.count_toks();
    stats.print(reader.posdict);
    print_stats(stats.emp_counts, "     posterior ", stats.denom);
    if (!final_iteration) {
      update<ModelType>(model, stats, reader.posdict, opts, iter);
    }
  }
  {
    std::lock_guard<std::mutex> lock(iomut);
    outp.done = true;
  }
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
