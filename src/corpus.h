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


#ifndef _CPYPDICT_H_
#define _CPYPDICT_H_

#include <string>
#include <locale>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <functional>
#include <cassert>
    
class Dict {
  typedef std::unordered_map<std::string, unsigned,
                                  std::hash<std::string> > Map;

 public:
  Dict() : b0_("<bad0>"), eps_("<eps>") {
    words_.reserve(1000);
    words_.push_back("<bad0>");
    kNULL_ = Convert(eps_);
    kDIV_ = Convert("|||");
  }

  unsigned max() const { return words_.size(); }

  void ConvertWhitespaceDelimitedLine(const std::string &line,
                                      std::vector<unsigned> &out, size_t cur = 0, size_t endpos = 0) {
    size_t last = 0;
    if (endpos == 0) {
      endpos = line.size();
    }
    bool inToken = 0;
    out.clear();
    while (cur < endpos) {
      if (isspace(line[cur++])) {
        if (!inToken) continue;
        out.push_back(Convert(line.substr(last, cur - last - 1)));
        inToken = false;
      } else {
        if (inToken) continue;
        last = cur - 1;
        inToken = true;
      }
    }
    if (inToken) out.push_back(Convert(line.substr(last, cur - last)));
  }

  inline void ParseLine(const std::string &line, std::vector<unsigned> &src, std::vector<unsigned> &trg) {
    static std::vector<unsigned> tmp;
    src.clear();
    trg.clear();
    ConvertWhitespaceDelimitedLine(line, tmp);

    unsigned i = 0;
    while (i < tmp.size() && tmp[i] != kDIV_) {
      src.push_back(tmp[i]);
      ++i;
    }
    if (i < tmp.size() && tmp[i] == kDIV_) {
      ++i;
      for (; i < tmp.size(); ++i) trg.push_back(tmp[i]);
    }
  }

  unsigned Convert(const std::string &word, bool frozen = false) {
    Map::iterator i = d_.find(word);
    if (i == d_.end()) {
      if (frozen) return 0;
      unsigned word_num = words_.size();
      d_[word] = word_num;
      words_.push_back(word);
      return word_num;
    } else {
      return i->second;
    }
  }

  const std::string &Convert(const unsigned id) const { return words_[id]; }

  unsigned kNULL_;
  unsigned kDIV_;

 private:
  std::string b0_;
  std::string eps_;
  std::vector<std::string> words_;
  Map d_;
};

inline void ReadFromFile(const std::string &filename, Dict *d,
                  std::vector<std::vector<unsigned> > *src,
                  std::set<unsigned> *src_vocab) {
  src->clear();
  std::cerr << "Reading from " << filename << std::endl;
  std::ifstream in(filename.c_str());
  assert(in);
  std::string line;
  // size_t lc = 0;
  while (getline(in, line)) {
    //++lc;
    src->push_back(std::vector<unsigned>());
    d->ConvertWhitespaceDelimitedLine(line, src->back());
    for (unsigned i = 0; i < src->back().size(); ++i)
      src_vocab->insert(src->back()[i]);
  }
}

#endif
