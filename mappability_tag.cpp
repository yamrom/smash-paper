//
// mappability_tag
//
// Copyright 2015 Peter Andrews @ CSHL
//

#include <algorithm>
using std::find_first_of;
using std::min;
using std::max_element;

#include <numeric>
using std::accumulate;

#include <exception>
using std::exception;

#include <functional>
using std::ref;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <map>
using std::map;

#include <sstream>
using std::istringstream;
using std::ostringstream;

#include <string>
using std::string;

#include <utility>
using std::move;
using std::pair;

#include <vector>
using std::vector;

#include "./chromosomes.h"
#include "./error.h"
#include "./util.h"

int main(int argc, char ** argv) try {
  --argc;
  if (argc != 2) throw paa::Error("usage: mappability_tag fasta_file in.sam");

  const string ref_name = argv[1];
  const Mappability mappability(ref_name + ".bin/map.bin");
  const ChromosomeInfo all_chr(ref_name + ".bin/sam_header.txt", false);

  const string input_name = argv[2];
  ifstream input(input_name.c_str());
  if (!input) throw paa::Error("Could not open smash sam file") << input_name;

  string line;
  const bool output_sam = true;
  unsigned int n = 0;
  while (getline(input, line)) {
    ++n;

    if (line[0] == '@') {
      if (output_sam) cout << line << endl;
    } else {
      istringstream linestr(line);
      string name;
      int flag;
      string chr;
      uint32_t pos;
      int qual;
      string cigar;
      linestr >> name >> flag >> chr >> pos >> qual >> cigar;
      const bool small_chr = chr.find("_gl000") != string::npos ||
          chr.find("chrM") != string::npos;
      if (output_sam) cout << line;

      if (!linestr) cerr << "parse error";
      istringstream cigarstr(cigar);
      uint32_t count;
      char code;
      int offset = 0;
      int uindex = 0;
      string optional;
      const auto abspos = all_chr.abspos(chr, pos);
      if (cigar != "*") {
        while (cigarstr >> count >> code) {
          if (code == '=') {
            ostringstream info;
            const auto left_m = mappability.left(abspos + offset + count - 1);
            const auto left = left_m ? left_m - 1 : 255;
            const auto right_m = mappability.right(abspos + offset - 1);
            const auto right = right_m ? right_m : 255;
            if (uindex < 10) {
              info << '\t' << 'L' << uindex << ":i:" << left
                   << '\t' << 'R' << uindex << ":i:" << right;
              optional += info.str();
            }
            if (left > count && !small_chr) {
              cerr << left_m << " " << right_m << " "
                   << line << '\t' << optional << endl;
              throw paa::Error("left mappability too big") << left;
            }
            if (right > count && !small_chr)
              throw paa::Error("right mappability too big") << right;

            ++uindex;
            // cout << n << " " << chr << " " << pos << " " << count
            // << " " << left << " " << right << endl;
          } else if (code != 'S' && code != 'M') {
            throw paa::Error("unexpected cigar") << code;
          }
          offset += count;
        }
      }
      if (output_sam) cout << optional << endl;
    }
  }

  return 0;
} catch(exception & e) {
  cerr << e.what() << endl;
  return 1;
} catch(...) {
  cerr << "Some exception was caught." << endl;
  return 1;
}

