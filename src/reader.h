#ifndef _READER_H_
#define _READER_H_

#include <fstream>                      // for ifstream
#include <memory>                       // for unique_ptr
#include <mutex>                        // for mutex
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "corpus.h"                     // for Dict
struct OutputNode;
struct Stats;
struct ThreadedOutput;

typedef std::vector<std::pair<std::vector<unsigned>,std::vector<unsigned>>> PairLines;
typedef std::vector<std::vector<unsigned>> PosLines;


// Handles reading the parallel data, pos data if present and conversion to unsigned
class Reader {
  std::ifstream data;
  std::unique_ptr<std::ifstream> posdata;
  Stats &stats;

  bool is_reverse;

  void read_line(PairLines &lines, std::string &line);
  void read_pos_line(PairLines const &lines, PosLines &poslines, std::string &line);

  // read N lines into the supplied lines vector and poslines vector
  bool read_n_lines(PairLines &lines, PosLines &poslines, size_t N);

public:
  Dict dict;
  Dict posdict;
  std::mutex iomut;

  Reader(std::string &datafile, std::string &posfile, Stats &stats, bool is_reverse);

  bool read_n_lines_threaded(PairLines &lines, PosLines &pos, size_t N, ThreadedOutput &outp, OutputNode *&out);
  
  void rewind();
};

bool print_progress(size_t lc);

#endif