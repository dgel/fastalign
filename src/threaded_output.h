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

#ifndef _THREADED_OUTPUT_H_
#define _THREADED_OUTPUT_H_

#include <mutex>                        // for mutex
#include <queue>                        // for queue
#include <string>                       // for string
#include <vector>                       // for vector


// used for coordinating output, each thread writes into the out vector of an OutputNode
// OutputNodes are flushed in correct order by a separate thread when writing thread is
// done with it
struct OutputNode {
  bool done = false;
  std::vector<std::string> out;
};

struct ThreadedOutput {
  bool done = false;
  std::mutex &iomut;
  std::queue<OutputNode> que;

  ThreadedOutput(std::mutex &mut);
  void finished();
};


void handleOutput(ThreadedOutput &out, size_t sleeptime);


#endif