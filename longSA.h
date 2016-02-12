/* Modifications of sparseMEM Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_LONGSA_H_
#define LONGMEM_LONGSA_H_

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <iosfwd>

#include "./size.h"
#include "./fasta.h"

// Stores the LCP array in an unsigned char (0-255).  Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP;
struct vec_uchar {
  struct item_t {
    item_t() {}
    item_t(size_t i, ANINT v) {
      idx = i;
      val = v;
    }
    size_t idx;
    ANINT val;
    bool operator < (const item_t & t) const { return idx < t.idx;  }
  };
  vec_uchar() : N_vec(0), vec(nullptr), cap_M(0), N_M(0), M(nullptr),
                using_mapping(false) {}
  ~vec_uchar();

  // Vector X[i] notation to get LCP values.
  ANINT operator[] (const size_t idx) const {
    if (vec[idx] == std::numeric_limits<unsigned char>::max())
      return std::lower_bound(M, M + N_M, item_t(idx, 0))->val;
    else
      return vec[idx];
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(const size_t idx, const ANINT v);
  // Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init();
  void load(const std::string & base, FILE * index);
  void save(const std::string & base, FILE * index) const;
  void resize(const size_t N);

 private:
  uint64_t N_vec;
  unsigned char * vec;
  uint64_t cap_M;
  uint64_t N_M;
  item_t * M;
  bool using_mapping;
  // std::vector<unsigned char> vec;  // LCP values from 0-65534
  // std::vector<item_t> M;
  vec_uchar(const vec_uchar & disabled_copy_constructor);
  vec_uchar & operator=(const vec_uchar & disabled_assignment_operator);
};

// depth : [start...end]
struct interval_t {
  interval_t() : depth(-1), start(1), end(0) { }
  interval_t(const uint64_t s, const uint64_t e,
             const uint64_t d) : depth(d), start(s), end(e) {}
  void reset(const uint64_t e) {
    depth = 0;
    start = 0;
    end = e;
  }
  uint64_t size() const { return end - start + 1; }
  uint64_t depth, start, end;
};

// Match find by findMEM.
struct match_t {
  match_t() {
    ref = 0;
    query = 0;
    len = 0;
  }
  match_t(const uint64_t r, const uint64_t q, const uint64_t l) {
    ref = r;
    query = q;
    len = l;
  }
  uint64_t ref;  // position in reference sequence
  uint64_t query;  // position in query
  uint64_t len;  // length of match
};

class Args;
class SAArgs {
 public:
  SAArgs() : verbose(false), mappability(false), ref_args() {}
  operator const RefArgs & () const { return ref_args; }
  bool verbose;
  bool mappability;
 private:
  RefArgs ref_args;
  SAArgs & operator=(const SAArgs & disabled_assignment_operator);
  friend class Args;
};

class Aligner;
struct longSA : public SAArgs {
  uint64_t size() const { return N; }

  bool using_mapping;

  const Sequence ref;
  const uint64_t N;  // !< Length of the sequence.

  const uint64_t logN;  // ceil(log(N))
  const uint64_t Nm1;  // N - 1

  //  std::vector<ANINT> SA;  // Suffix array.
  //  std::vector<ANINT> ISA;  // Inverse suffix array.
  ANINT * SA;
  ANINT * ISA;
  vec_uchar LCP;  // Simulates a vector<int> LCP.

  // Constructor builds suffix array.
  explicit longSA(const SAArgs & arguments);
  ~longSA();

  // Modified Kasai et all for LCP computation.
  void computeLCP();

  // Binary search for left boundry of interval.
  inline uint64_t bsearch_left(const char c, const uint64_t i,
                                    uint64_t l, uint64_t r) const;
  // Binary search for right boundry of interval.
  inline uint64_t bsearch_right(const char c, const uint64_t i,
                                     uint64_t l, uint64_t r) const;

  // Simple suffix array search.
  inline bool search(const std::string &P, uint64_t &start,
                     uint64_t &end) const;

  // Simple top down traversal of a suffix array.
  inline bool top_down(const char c, const uint64_t i,
                       uint64_t &start, uint64_t &end) const;
  inline bool top_down_faster(const char c, const uint64_t i,
                              uint64_t &start, uint64_t &end) const;

  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(const std::string &P, const uint64_t prefix,
                       interval_t &cur, const ANINT min_len) const;

  // Simulate a suffix link.
  inline bool suffixlink(interval_t * m) const;

  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(interval_t * link) const {
    const uint64_t thresh = 2 * link->depth * logN;
    uint64_t exp = 0;  // Threshold link expansion.
    uint64_t start = link->start;
    uint64_t end = link->end;
    while (LCP[start] >= link->depth) {
      if (++exp >= thresh) return false;
      --start;
    }
    while (end < Nm1 && LCP[end+1] >= link->depth) {
      if (++exp >= thresh) return false;
      ++end;
    }
    link->start = start;
    link->end = end;
    return true;
  }

  // Given a position i in S, finds a left maximal match of minimum length
  inline void find_Lmaximal(Aligner & query, const uint64_t prefix,
                            const uint64_t i,
                            const uint64_t len) const;

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  void collectMEMs(Aligner & query, const uint64_t prefix,
                   const interval_t mli, interval_t xmi) const;

  // Find all MEMs
  void findMEM(Aligner & query) const;

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  // NOTE: min_len must be > 1
  void MAM(Aligner & query) const;
  inline bool is_leftmaximal(const std::string &P, const uint64_t p1,
                             const uint64_t p2) const;

  // Find Maximal Exact Matches (MEMs)
  void MEM(Aligner & query) const;

  // Maximal Unique Match (MUM)
  void MUM(Aligner & query) const;

  void show_mappability(const std::string & filename) const;

  template <class Out>
  void show(Out & output, const bool bin) const;

 private:
  longSA(const longSA & disabled_copy_constructor);
  longSA & operator=(const longSA & disabled_assignment_operator);
};

#endif  // LONGMEM_LONGSA_H_

