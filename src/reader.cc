#include "reader.h"
#include <ext/alloc_traits.h>
#include <stdlib.h>                     // for exit
#include <iostream>                     // for cerr
#include <queue>                        // for queue
#include "model.h"                      // for Stats
#include "threaded_output.h"            // for OutputNode, ThreadedOutput

Reader::Reader(std::string &datafile, std::string &posfile, Stats &stats,
               bool is_reverse)
    : data(datafile), stats(stats), is_reverse(is_reverse), print_newline(false), iomut() {

  if (!data) {
    std::cerr << "Can't read " << datafile << std::endl;
    exit(1);
  }

  // gets deleted automatically
  if (posfile != "") {
    posdata = std::unique_ptr<std::ifstream>(new std::ifstream(posfile));
    if (!(*posdata)) {
      std::cerr << "Can't read " << posfile << std::endl;
      exit(1);
    }
  } else {
    posdict.Convert("Default");
  }
}

bool print_progress(size_t lc) {
  if (lc % 1000 == 0) {
    std::cerr << '.';
  }
  if (lc % 50000 == 0) {
    std::cerr << " [" << lc << "]" << std::endl;
    return false;
  }
  return true;
}

bool Reader::read_n_lines(PairLines &lines, PosLines &poslines, size_t N) {
  lines.clear();
  poslines.clear();

  size_t nread = 0;
  while (std::getline(data, readbuf)) {
    // give user some idea of progress
    print_newline = print_progress(++stats.lc);
    lines.emplace_back();
    read_line(lines.back());
    if (posdata != nullptr) {
      poslines.emplace_back();
      read_pos_line(poslines.back(), is_reverse ? lines.back().first.size() : lines.back().second.size());
    }
    if (++nread == N) break;
  }
  // signal if done processing corpus
  if (nread < N) {
    if (print_newline) {
      std::cerr << std::endl;
      print_newline = false;
    }
    return true;
  }
  return false;
}

void Reader::read_line(PairLine &pair) {
  if (is_reverse) {
    dict.ParseLine(readbuf, &pair.second, &pair.first);
  } else {
    dict.ParseLine(readbuf, &pair.first, &pair.second);
  }

  if (pair.first.empty() || pair.second.empty()) {
    std::cerr << "Error in line " << stats.lc << "\n" << readbuf << std::endl;
    exit(1);
  }
}

void Reader::read_pos_line(Line &posline, size_t n_tokens) {
  if (!std::getline(*posdata, readbuf)) {
    std::cerr << "POS file has too few lines" << std::endl;
    std::exit(1);
  }
  if (is_reverse) {
    posdict.ParseLine(readbuf, &posline, nullptr);
  } else {
    posdict.ParseLine(readbuf, nullptr, &posline);
  }
  if (n_tokens != posline.size()) {
    std::cerr << "Error: POS line has fewer tags than trg line has words.\n"
                 "    Perhaps src and trg should be swapped?" << std::endl;
    std::cerr << "In line: " << stats.lc << std::endl;
    exit(1);
  }
}

ReadStatus Reader::read_n_lines_threaded(PairLines &lines, PosLines &pos, size_t N, ThreadedOutput &outp, OutputNode *&out) {
  std::lock_guard<std::mutex> io_lock(iomut);
  bool done = read_n_lines(lines, pos, N);
  if (lines.size() == 0) {
    return ReadStatus::NOTHING_READ;
  }
  outp.que.emplace();
  out = &outp.que.back();
  return done? ReadStatus::FINISHED : ReadStatus::CONTINUE;
}

void Reader::rewind() {
  data.clear();
  data.seekg(0);
  if (posdata) {
    posdata->clear();
    posdata->seekg(0);
  }
  print_newline = true;
}

