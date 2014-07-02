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
#include <functional>
#include <cassert>

#include "map_type.h"

typedef std::uint32_t Token;
    
class Dict {
  typedef MapType<std::string, Token,
                                  std::hash<std::string> > Map;

 public:
  Dict() : b0_("<bad0>"), eps_("<eps>") {
    #ifdef __NEED_SET_EMPTY_KEY__
      d_.set_empty_key("<<e>>");
      d_.set_deleted_key("<<d>>");
    #endif

    words_.reserve(1000);
    words_.push_back("<bad0>");
    kNULL_ = Convert(eps_);
    kDIV_ = Convert("|||");
  }

  unsigned max() const { return words_.size(); }

  size_t Skip(const std::string &line, size_t cur) {
    size_t last = 0;
    bool inToken = false;
    Token tmp;

    while (cur < line.size()) {
      if (isspace(line[cur++])) {
        if (!inToken) continue;
        tmp = Convert(line.substr(last, cur - last - 1));
        if (tmp == kDIV_) {
          return cur;
        } 
      } else {
        if (inToken) continue;
        last = cur - 1;
        inToken = true;
      }
    }
    return line.size();
  }


  size_t ConvertWhitespaceDelimitedLinePart(const std::string &line, std::vector<Token> *out, size_t cur) {
    if (!out) {
      return Skip(line, cur);
    }

    size_t last = cur;
    bool inToken = false;
    Token tmp;

    out->clear();
    while (cur < line.size()) {
      if (isspace(line[cur++])) {
        if (!inToken) continue;
        tmp = Convert(line.substr(last, cur - last - 1));
        
        if (tmp == kDIV_) {
          return cur;
        } else {
          out->push_back(tmp);
        }
        inToken = false;
      } else {
        if (inToken) continue;
        last = cur - 1;
        inToken = true;
      }
    }
    if (inToken) {
      tmp = Convert(line.substr(last, cur - last));
      if (tmp != kDIV_) {
        out->push_back(tmp);
      }
    }
    return line.size();
  }

  inline void ParseLine(const std::string &line, std::vector<Token> *src, std::vector<Token> *trg) {
    size_t midpoint = ConvertWhitespaceDelimitedLinePart(line, src, 0);
    ConvertWhitespaceDelimitedLinePart(line, trg, midpoint);
  }

  Token Convert(const std::string &word, bool frozen = false) {
    Map::iterator i = d_.find(word);
    if (i == d_.end()) {
      if (frozen) return 0;
      Token word_num = words_.size();
      d_[word] = word_num;
      words_.push_back(word);
      return word_num;
    } else {
      return i->second;
    }
  }

  const std::string &Convert(const unsigned id) const { return words_[id]; }

  Token kNULL_;
  Token kDIV_;

 private:
  std::string b0_;
  std::string eps_;
  std::vector<std::string> words_;
  Map d_;
};


#endif
