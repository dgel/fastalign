#include "reader.h"
#include <ext/alloc_traits.h>
#include <stdlib.h>                     // for exit
#include <iostream>                     // for cerr
#include <queue>                        // for queue
#include "model.h"                      // for Stats
#include "threaded_output.h"            // for OutputNode, ThreadedOutput

Reader::Reader(std::string &datafile, std::string &posfile, Stats &stats,
               bool is_reverse)
    : data(datafile), stats(stats), is_reverse(is_reverse), iomut() {

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
  bool print_newline = false;

  std::string line;
  size_t nread = 0;
  while (std::getline(data, line)) {
    ++stats.lc;
    // give user some idea of progress
    print_newline = print_progress(++stats.lc);
    
    read_line(lines, line);
    read_pos_line(lines, poslines, line);
    
    if (++nread == N) break;
  }

  // signal if done processing corpus
  if (nread < N) {
    if (print_newline)
      std::cerr << std::endl;
    return true;
  }
  return false;
}

void Reader::read_line(PairLines &lines, std::string &line) {
  lines.emplace_back();
  if (is_reverse) {
    dict.ParseLine(line, lines.back().second, lines.back().first);
  } else {
    dict.ParseLine(line, lines.back().first, lines.back().second);
  }

  if (lines.back().first.empty() || lines.back().second.empty()) {
    std::cerr << "Error in line " << stats.lc << "\n" << line << std::endl;
    exit(1);
  }
}

void Reader::read_pos_line(PairLines const &lines, PosLines &poslines, std::string &line) {
  if (posdata != nullptr) {
    if (!std::getline(*posdata, line)) {
      std::cerr << "POS file has too few lines" << std::endl;
      std::exit(1);
    }

    size_t delimpos = line.find("|||");
    poslines.emplace_back();

    if (is_reverse) {
      posdict.ConvertWhitespaceDelimitedLine(line, poslines.back(), 0,
                                             delimpos - 1);
    } else {
      posdict.ConvertWhitespaceDelimitedLine(line, poslines.back(),
                                             delimpos + 4);
    }
    if (lines.back().second.size() != poslines.back().size()) {
      std::cerr << "Error: POS line has fewer tags than trg line has words.\n"
                   "    Perhaps src and trg should be swapped?" << std::endl;
      std::cerr << "In line: " << stats.lc << std::endl;
      exit(1);
    }
  }
}

void Reader::rewind() {
  data.clear();
  data.seekg(0);
  if (posdata) {
    posdata->clear();
    posdata->seekg(0);
  }
}

bool Reader::read_n_lines_threaded(PairLines &lines, PosLines &pos, size_t N, ThreadedOutput &outp, OutputNode *&out) {
  std::lock_guard<std::mutex> io_lock(iomut);
  bool done = read_n_lines(lines, pos, N);
  outp.que.emplace();
  out = &outp.que.back();
  return done;
}