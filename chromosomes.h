//
// chromosomes.h
//
// Copyright 2015 Peter Andrews @ CSHL
//

#include <istream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "./error.h"

bool read_this(const std::string & needed, std::istream & in) {
  for (unsigned int i = 0; i != needed.size(); ++i) {
    const char read = in.get();
    if (read != needed[i]) {
      return false;
    }
  }
  return true;
}

class ChromosomeInfo {
 public:
  ChromosomeInfo(const std::string & sam_file_name,
                 const bool simple_only = false) {
    std::ifstream in(sam_file_name.c_str());
    if (!in) throw paa::Error("Problem opening sam file") << sam_file_name;
    uint64_t offset_ = 0;
    if (in.peek() != '@') {
      throw paa::Error("Is the sam header missing?") << sam_file_name;
    }
    while (in && in.peek() == '@') {
      if (read_this("@SQ\tSN:", in)) {
        std::string chromosome;
        unsigned int length_ = 0;
        in >> chromosome;
        if (read_this("\tLN:", in)) {
          in >> length_;
          if (!in) {
            throw paa::Error("Parse error in");
          }
        } else {
          throw paa::Error("Parse Problem read_this LN:");
        }
        in.ignore(1000, '\n');
        if (simple_only && chromosome.find_first_of("_M") !=
            std::string::npos) {
          continue;
        }
        lookup[chromosome] = names.size();
        names.push_back(chromosome);
        lengths.push_back(length_);
        offsets.push_back(offset_);
        offset_ += length_;
        ends.push_back(offset_);
      } else {
        in.ignore(1000, '\n');
      }
    }
    lookup["*"] = names.size();
  }

  ChromosomeInfo(const std::vector<std::string> & snames,
                 const std::vector<uint64_t> & sizes,
                 const bool simple_only = false,
                 const bool one_strand = true) {
    uint64_t offset_ = 0;
    for (unsigned int c2 = 0; c2 < snames.size();
         (one_strand ? ++c2 : c2 += 2)) {
      const std::string & chromosome = snames[c2];
      const uint64_t length_ = sizes[c2];
      if (simple_only &&
          chromosome.find_first_of("_M") != std::string::npos) {
        continue;
      }
      lookup[chromosome] = names.size();
      names.push_back(chromosome);
      lengths.push_back(length_);
      offsets.push_back(offset_);
      offset_ += length_;
      ends.push_back(offset_);
    }
    lookup["*"] = names.size();
  }


  const std::string & name(const unsigned int i) const {
    return names[i];
  }
  unsigned int offset(const unsigned int i) const {
    return offsets[i];
  }
  unsigned int end(const unsigned int i) const {
    return ends[i];
  }
  unsigned int length(const unsigned int i) const {
    return lengths[i];
  }
  unsigned int index(const std::string & chr) const {
    auto found = lookup.find(chr);
    if (found == lookup.end()) throw paa::Error("Unknown chromosome") << chr;
    return found->second;
  }
  unsigned int abspos(const std::string & chr, const unsigned int pos) const {
    return offsets[index(chr)] + pos;
  }
  unsigned int abspos(const std::pair<std::string,
                      unsigned int> & chr_pos) const {
    return offsets[index(chr_pos.first)] + chr_pos.second;
  }
  std::pair<std::string, unsigned int> chrpos(const unsigned int pos) const {
    unsigned int chromosome = 0;
    while (chromosome != names.size()) {
      if (pos < ends[chromosome]) {
        return std::make_pair(names[chromosome], pos - offsets[chromosome]);
      }
      ++chromosome;
    }
    throw paa::Error("chrpos not located for") << pos;
  }
  unsigned int size() const {
    return names.size();
  }
  bool contains(const std::string & chr) const {
    auto found = lookup.find(chr);
    if (found == lookup.end()) {
      return false;
    } else {
      return true;
    }
  }
  void output(std::ostream & out) const {
    for (unsigned int i = 0; i != size(); ++i) {
      out << i << '\t'
          << names[i] << '\t'
          << lengths[i] << '\t'
          << offsets[i] << '\t'
          << ends[i]
          << std::endl;
    }
  }

 private:
  std::map<std::string, unsigned int> lookup;
  std::vector<std::string> names;
  std::vector<unsigned int> offsets;
  std::vector<unsigned int> ends;
  std::vector<unsigned int> lengths;
};



