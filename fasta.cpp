/* Modifications of sparseMEM Copyright Peter Andrews 2013 CSHL */

#include "./fasta.h"

#include <sys/mman.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "./util.h"
#include "./error.h"
using paa::Error;

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ostringstream;

// Return the reverse complement of sequence. This allows searching
// the plus strand of instances on the minus strand.
void reverse_complement(string * const seq_rc) {
  // Reverse in-place.
  reverse(seq_rc->begin(), seq_rc->end());
  for (uint64_t i = 0; i != seq_rc->size(); ++i) {
    // Adapted from Kurtz code in MUMmer v3.
    char & ch = (*seq_rc)[i];
    switch (ch) {
    case 'a': ch = 't'; break;
    case 'c': ch = 'g'; break;
    case 'g': ch = 'c'; break;
    case 't': ch = 'a'; break;
    case 'r': ch = 'y'; break; /* a or g */
    case 'y': ch = 'r'; break; /* c or t */
    case 'm': ch = 'k'; break; /* a or c */
    case 'k': ch = 'm'; break; /* g or t */
    case 'b': ch = 'v'; break; /* c, g or t */
    case 'd': ch = 'h'; break; /* a, g or t */
    case 'h': ch = 'd'; break; /* a, c or t */
    case 'v': ch = 'b'; break; /* a, c or g */
    case 'A': ch = 'T'; break;
    case 'C': ch = 'G'; break;
    case 'G': ch = 'C'; break;
    case 'T': ch = 'A'; break;
    case 'R': ch = 'Y'; break; /* a or g */
    case 'Y': ch = 'R'; break; /* c or t */
    case 'M': ch = 'K'; break; /* a or c */
    case 'K': ch = 'M'; break; /* g or t */
    case 'B': ch = 'V'; break; /* c, g or t */
    case 'D': ch = 'H'; break; /* a, g or t */
    case 'H': ch = 'D'; break; /* a, c or t */
    case 'V': ch = 'B'; break; /* a, c or g */
      default:
        break;
    }
  }
}

// Trim a string, giving start and end in trimmed version.
// NOTE: Assumes line.length() > 0!!!!
void trim(const string &line, uint64_t &start, uint64_t &end) {
  // Trim leading spaces.
  for (uint64_t i = start; i < line.size(); ++i) {
    if (line[i] != ' ') {
      start = i;
      break;
    }
  }
  // Trim trailing spaces.
  for (uint64_t i = line.size() ; i != 1; --i) {
    if (line[i-1] != ' ') {
      end = i;
      break;
    }
  }
}
Sequence::~Sequence() {
  if (memory_mapped && using_mapping) {
    if (munmap(seq, N * sizeof(*seq)))
      throw Error("sequence memory unmap failure");
  } else {
    if (using_mapping) free(seq);
  }
}
Sequence::Sequence(const RefArgs & arguments)
    : RefArgs(arguments), using_mapping(false) {
  const time_t start_time = time(nullptr);

  // Reference cache filename
  ostringstream saved_reference_stream;
  saved_reference_stream << ref_fasta << ".bin";
  bin_base = saved_reference_stream.str();
  mkdir(bin_base);
  saved_reference_stream << "/rc" << rcref << ".ref";
  bin_base = saved_reference_stream.str();
  saved_reference_stream << ".bin";
  const string saved_reference = saved_reference_stream.str();

  const uint64_t fasta_size = file_size(ref_fasta);

  // Load or create reference
  if (readable(saved_reference)) {
    if (verbose) cerr << "# loading reference binary" << endl;

    FILE * reference = fopen(saved_reference.c_str(), "rb");
    if (!reference) throw Error("could not open reference bin file")
                        << saved_reference << "for reading";

    uint64_t fasta_saved_size;
    bread(reference, fasta_saved_size, "fasta_size");
    if (fasta_size != fasta_saved_size)
      throw Error("")
          << "reference fasta size has changed\n"
          << "maybe the reference has changed?\n"
          << "If so, you will need to delete the current reference to proceed";
    bread(reference, N, "N");
    using_mapping = true;
    bread(bin_base + ".seq.bin", seq, "seq", N);
    uint64_t descr_size;
    bread(reference, descr_size, "descr_size");
    startpos.resize(descr_size);
    sizes.resize(descr_size);
    descr.resize(descr_size);
    for (unsigned int i = 0; i != descr_size; ++i) {
      bread(reference, startpos[i], "startpos");
      bread(reference, sizes[i], "sizes");
      uint64_t string_size;
      bread(reference, string_size, "string_size");
      descr[i].resize(string_size);
      bread(reference, descr[i][0], "description string", string_size);
    }
    bread(reference, maxdescrlen, "maxdescrlen");
    if (fclose(reference) != 0) throw Error("problem closing reference file");
  } else {
    if (verbose) cerr << "# loading reference from fasta" << endl;

    string meta, line;
    uint64_t length = 0;

    // Everything starts at zero.
    startpos.push_back(0);

    ifstream data(ref_fasta);

    if (!data.is_open()) throw Error("unable to open") << ref_fasta;

    while (!data.eof()) {
      getline(data, line);  // Load one line at a time.
      if (!data.eof() && line.size() == 0) continue;

      uint64_t start = 0, end = line.size();

      // Meta tag line and start of a new sequence.
      if (data.eof() || line[0] == '>') {
        // Save previous sequence and meta data.
        if (length > 0) {
          const uint64_t this_start = startpos.back();
          if (verbose) cerr << "# " << meta << " " << length
                            << " " << this_start << endl;
          descr.push_back(meta);
          if (rcref || !data.eof()) {
            seq_vec.push_back('`');  // ` character used to separate strings
            startpos.push_back(seq_vec.size());
          }
          sizes.push_back(length);
          if (rcref) {
            descr.push_back(meta);
            sizes.push_back(length);
            string reversed(seq_vec.begin() + this_start,
                            seq_vec.begin() + this_start + length);
            reverse_complement(&reversed);
            seq_vec.insert(seq_vec.end(), reversed.begin(), reversed.end());
            if (!data.eof()) {
              seq_vec.push_back('`');  // ` character used to separate strings
              startpos.push_back(seq_vec.size());
            }
          }
          if (data.eof()) break;
        }
        // Reset parser state.
        start = 1;
        meta = "";
        length = 0;
      }
      trim(line, start, end);
      // Collect meta data.
      if (line[0] == '>') {
        for (uint64_t i = start; i != end; ++i) {
          if (line[i] == ' ') break;
          meta += line[i];
        }
      } else {  // Collect sequence data.
        length += end - start;
        for (uint64_t i = start; i != end; ++i) {
          seq_vec.push_back(tolower(line[i]));
        }
      }
    }
    seq_vec.push_back('$');
    if (verbose) cerr << "# seq_vec.length=" << seq_vec.size() << endl;

    N = seq_vec.size();
    seq = &seq_vec[0];

    // Get maximum reference sequence description length.
    maxdescrlen = 0;
    for (uint64_t i = 0; i != descr.size(); ++i) {
      if (maxdescrlen < descr[i].size())  maxdescrlen = descr[i].size();
    }

    if (verbose) cerr << "# saving reference binary" << endl;

    FILE * reference = fopen(saved_reference.c_str(), "wb");
    if (reference == nullptr)
      throw Error("could not open reference")
          << saved_reference << "for writing";
    bwrite(reference, fasta_size, "fasta_size");
    bwrite(reference, N, "N");
    bwrite(bin_base + ".seq.bin", seq[0], "seq", N);
    const uint64_t descr_size = descr.size();
    bwrite(reference, descr_size, "descr_size");
    for (unsigned int i = 0; i != descr_size; ++i) {
      bwrite(reference, startpos[i], "startpos");
      bwrite(reference, sizes[i], "sizes");
      const uint64_t string_size = descr[i].size();
      bwrite(reference, string_size, "string_size");
      bwrite(reference, descr[i][0], "description string", string_size);
    }
    bwrite(reference, maxdescrlen, "maxdescrlen");
    if (fclose(reference) != 0)
      throw Error("problem closing reference file");
  }

  const time_t end_time = time(nullptr);
  if (verbose) cerr << "# constructed reference in "
                    << end_time - start_time << " seconds" << endl;
}

string Sequence::sam_header() const {
  ostringstream out;
  out << "@HD\tVN:1.0\tSO:unsorted" << endl;
  for (unsigned int chr = 0; chr < sizes.size();
       chr += rcref ? 2 : 1) {
    out << "@SQ\tSN:" << descr[chr] << "\tLN:" << sizes[chr] << endl;
  }
  out << "@PG\tID:longMEM\tPN:longMEM\tVN:0.5" << endl;
  return out.str();
}
