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

#ifndef _FA_UTILS_H_
#define _FA_UTILS_H_

#include <thread>


// map a function over a vector, spawning a thread for each function call
template <typename T, typename Func>
auto par_map(std::vector<T> &in, Func f) -> std::vector<decltype(f(in[0]))> {
  typedef decltype(f(in[0])) Res;
  std::vector<Res> result(in.size());
  std::vector<std::thread> threads;

  for (size_t i = 0; i < in.size(); ++i) {
    threads.emplace_back([&,i] {result[i] = f(in[i]);});
  }
  for (auto &t : threads) {
    t.join();
  }
  return result;
}

// run a function func in N threads, with the supplied arguments
template <typename F, typename... Args>
void run_n(size_t N, F func, Args &... args) {
  std::vector<std::thread> threads;
  
  for (size_t i = 0; i < N; ++i) {
    threads.emplace_back(func, std::ref(args)...);
  }
  for (std::thread &t : threads) {
    t.join();
  }
}

// divide the range [0..size) into N ranges, as equally distributed as possible.
std::vector<std::pair<size_t, size_t>> n_ranges(size_t size, size_t N) {
  size_t incr = size / N;
  size_t rem = size % N;
  std::vector<std::pair<size_t, size_t>> result;

  for (size_t i = 0, low = 0; i < N; ++i) {
    result.push_back({low, low + incr + (i < rem ? 1 : 0)});
    low += incr + (i < rem ? 1 : 0);
  }

  return result;
}


#endif
