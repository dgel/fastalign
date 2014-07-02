#include "threaded_output.h"
#include <iostream>                     // for endl, ostream, etc
#include <thread>                       // for this_thread::sleep_for
#include <chrono>                       // for chrono::seconds


ThreadedOutput::ThreadedOutput(std::mutex &mut)
    : iomut(mut)
{};

void ThreadedOutput::finished()
{
  std::lock_guard<std::mutex> lock(iomut);
  done = true;
}

OutputNode *ThreadedOutput::getOutput() {
  std::lock_guard<std::mutex> lock(iomut);
  que.emplace();
  return &que.back();
}

void handleOutput(ThreadedOutput &out, size_t sleeptime) {
  while (true) {
    if (out.iomut.try_lock()) {
      while (!out.que.empty() && out.que.front().done) {
        auto &head = out.que.front().out;
        for (std::string &line : head) {
          if (!line.empty())
            line.pop_back();
          std::cout << line << std::endl;
        }
        out.que.pop();
      }
      out.iomut.unlock();
      if (out.done) {
        return;
      }
    } 
    std::this_thread::sleep_for (std::chrono::seconds(sleeptime));
  }
}

