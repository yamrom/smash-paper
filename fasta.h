/* Modifications of sparseMEM Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_FASTA_H_
#define LONGMEM_FASTA_H_

#include <string>
#include <vector>

#include "./size.h"

void reverse_complement(std::string * const seq_rc);
void trim(const std::string & line, uint64_t & start, uint64_t & end);

class RefArgs {
 public:
  RefArgs() : ref_fasta(nullptr), rcref(false), verbose(false) {}
  const char * ref_fasta;
  bool rcref;
  bool verbose;
 private:
  RefArgs & operator=(const RefArgs & disabled_assignment_operator);
};

class Sequence : public RefArgs {
 public:
  explicit Sequence(const RefArgs & arguments);
  ~Sequence();
  char operator[] (const uint64_t n) const { return seq[n]; }
  std::string sam_header() const;
  uint64_t N;  // !< Length of the sequence.
  std::vector<char> seq_vec;
  char * seq;
  std::vector<std::string> descr;
  std::vector<uint64_t> startpos;
  std::vector<uint64_t> sizes;
  uint64_t maxdescrlen;
  std::string bin_base;
 private:
  bool using_mapping;
  Sequence(const Sequence & disabled_copy_constructor);
  Sequence & operator=(const Sequence & disabled_assignment_operator);
};

#endif  // LONGMEM_FASTA_H_
