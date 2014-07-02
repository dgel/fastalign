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

typedef std::vector<Token> Line;
typedef std::pair<Line, Line> PairLine;
typedef std::vector<PairLine> PairLines;
typedef std::vector<Line> PosLines;


enum class ReadStatus {
  CONTINUE,
  FINISHED,
  NOTHING_READ
};

// Handles reading the parallel data, pos data if present and conversion to Token
class Reader {
  std::ifstream data;
  std::unique_ptr<std::ifstream> posdata;
  Stats &stats;

  std::string readbuf;

  bool is_reverse;
  bool print_newline;

  void read_line(PairLine &lines);
  void read_pos_line(Line &posline, size_t n_tokens);

  // read N lines into the supplied lines vector and poslines vector
  bool read_n_lines(PairLines &lines, PosLines &poslines, size_t N);

public:
  Dict dict;
  Dict posdict;
  std::mutex iomut;

  Reader(std::string &datafile, std::string &posfile, Stats &stats, bool is_reverse);

  ReadStatus read_n_lines_threaded(PairLines &lines, PosLines &pos, size_t N, ThreadedOutput &outp, OutputNode *&out);
  bool read_1_line(PairLine &lines, Line &pos);
  
  void rewind();
};

bool print_progress(size_t lc);

#endif