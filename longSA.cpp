/* Modifications of sparseMEM Copyright Peter Andrews 2013 CSHL */

#include "./longSA.h"

#include <math.h>
#include <limits.h>
#include <sys/mman.h>

#include <algorithm>
using std::max;
using std::sort;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <limits>
using std::numeric_limits;

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ostream;
using std::ofstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "./query.h"
#include "./util.h"
#include "./error.h"
using paa::Error;

// LS suffix sorter (integer alphabet).
void suffixsort(ANINT *x, ANINT *p,
                const ANINT n, const ANINT k, const ANINT l);

void vec_uchar::set(const size_t idx, const ANINT v) {
  if (v >= numeric_limits<unsigned char>::max()) {
    vec[idx] = numeric_limits<unsigned char>::max();
    if (N_M == cap_M) {
      if (cap_M == 0) cap_M = 1024;
      cap_M *= 2;
      if ((M = reinterpret_cast<item_t *>(realloc(M, sizeof(item_t) * cap_M)))
          == nullptr) throw Error("M realloc failed");
    }
    M[N_M++] = item_t(idx, v);
  } else {
    vec[idx] = (unsigned char)v;
  }
}

void vec_uchar::init() {
  sort(M, M + N_M);
}

void vec_uchar::load(const string & base, FILE * index) {
  using_mapping = true;
  bread(index, N_vec, "N_vec");
  bread(base + ".lcp.vec.bin", vec, "vec", N_vec);
  bread(index, N_M, "N_M");
  bread(base + ".lcp.m.bin", M, "M", N_M);
}

void vec_uchar::save(const string & base, FILE * index) const {
  bwrite(index, N_vec, "N_vec");
  bwrite(base + ".lcp.vec.bin", vec[0], "vec", N_vec);
  bwrite(index, N_M, "N_M");
  bwrite(base + ".lcp.m.bin", M[0], "M", N_M);
}
vec_uchar::~vec_uchar() {
  if (memory_mapped && using_mapping) {
    if (munmap(vec, N_vec * sizeof(unsigned char)))
      throw Error("vec memory unmap failure");
    if (N_M)
      if (munmap(M, N_M * sizeof(item_t)))
        throw Error("M memory unmap failure");
  } else {
    free(vec);
    free(M);
  }
}
void vec_uchar::resize(const size_t N) {
  N_vec = N;
  if ((vec = (unsigned char *)malloc(sizeof(unsigned char) * N)) == nullptr)
    throw Error("malloc error for lcp");
}

longSA::longSA(const SAArgs & arguments)
    : SAArgs(arguments), using_mapping(false),
      ref(arguments), N(ref.N),  // S(ref.seq),
      logN((uint64_t)(ceil(log(N) / log(2.0)))), Nm1(N - 1) {
  const time_t start_time = time(nullptr);

  // Index cache filename
  ostringstream saved_index_stream;
  saved_index_stream << ref.ref_fasta << ".bin";
  saved_index_stream << "/rc" << ref.rcref;
  saved_index_stream << ".i" << sizeof(ANINT) << ".index";
  const string bin_base = saved_index_stream.str();
  saved_index_stream << ".bin";
  const string saved_index = saved_index_stream.str();

  const uint64_t fasta_size = file_size(ref.ref_fasta);

  // Load or create index
  if (readable(saved_index)) {
    if (verbose) cerr << "# loading index binary" << endl;

    FILE * index = fopen(saved_index.c_str(), "rb");
    if (index == nullptr)
      throw Error("could not open index") << saved_index << "for reading";

    uint64_t fasta_saved_size;
    bread(index, fasta_saved_size, "fasta_size");
    if (fasta_size != fasta_saved_size)
      throw Error("saved fasta size used for index does not"
                  "match current fasta size\n"
                  "maybe the reference has changed?\n"
                  "you may need to delete the current index to proceed");
    uint64_t dummy;
    bread(index, dummy, "logN");
    bread(index, dummy, "logN");
    uint64_t SA_size;
    bread(index, SA_size, "SA_size");

    using_mapping = true;
    bread(bin_base + ".sa.bin", SA, "SA", SA_size);
    bread(bin_base + ".isa.bin", ISA, "ISA", SA_size);
    LCP.load(bin_base, index);
    if (fclose(index) != 0) throw Error("problem closing index file");
  } else {
    if (verbose) cerr << "# creating index from reference" << endl;

    if ((SA = reinterpret_cast<ANINT *>(malloc(sizeof(ANINT) * N))) == nullptr)
      throw Error("SA malloc error");
    if ((ISA = reinterpret_cast<ANINT *>(malloc(sizeof(ANINT) * N))) == nullptr)
      throw Error("ISA malloc error");

    int64_t char2int[UCHAR_MAX+1];  // Map from char to integer alphabet.

    // Zero char2int mapping.
    for (ANINT i = 0; i <= UCHAR_MAX; ++i) char2int[i] = 0;

    // Determine which characters are used in the string S.
    for (uint64_t i = 0; i != N; ++i) char2int[(ANINT)ref[i]]=1;

    // Count the size of the alphabet.
    ANINT alphasz = 0;
    for (ANINT i = 0; i <= UCHAR_MAX; ++i) {
      if (char2int[i])
        char2int[i] = alphasz++;
      else
        char2int[i] = -1;
    }

    // Remap the alphabet.
    for (uint64_t i = 0; i != N; ++i) ISA[i] = (ANINT)ref[i];
    for (uint64_t i = 0; i != N; ++i) ISA[i]=char2int[ISA[i]] + 1;
    // First "character" equals 1 because of above plus one,
    // l=1 in suffixsort().
    const ANINT alphalast = alphasz + 1;

    // Use LS algorithm to construct the suffix array.
    suffixsort(&ISA[0], &SA[0], N - 1, alphalast, 1);

    // Use algorithm by Kasai et al to construct LCP array.
    LCP.resize(N);
    computeLCP();   // SA + ISA -> LCP
    LCP.init();

    if (verbose) cerr << "# saving index" << endl;

    FILE * index = fopen(saved_index.c_str(), "wb");
    if (index == nullptr)
      throw Error("could not open index") << saved_index << "for writing";
    bwrite(index, fasta_size, "fasta_size");
    bwrite(index, logN, "logN");
    bwrite(index, Nm1, "Nm1");
    bwrite(index, N, "SA_size");
    bwrite(bin_base + ".sa.bin", SA[0], "SA", N);
    bwrite(bin_base + ".isa.bin", ISA[0], "ISA", N);
    LCP.save(bin_base, index);
    if (fclose(index) != 0)
      throw Error("problem closing index file");
  }
  const time_t end_time = time(nullptr);
  if (verbose) cerr << "# constructed index in "
                    << end_time - start_time << " seconds" << endl;

  const bool exercise = false;
  if (exercise) {
    if (verbose) cerr << "# exercising memory mapped file" << endl;
    uint64_t temp = 0;
    for (uint64_t i = 0; i < N; i += 4096)
      temp += ref[i];
    for (uint64_t i = 0; i < N; i += 4096)
      temp += SA[i];
    for (uint64_t i = 0; i < N; i += 4096)
      temp += ISA[i];
    for (uint64_t i = 0; i < N; i += 4096)
      temp += LCP[i];
    if (verbose) cerr << "# done exercising memory mapped file" << endl;
  }
}

longSA::~longSA() {
  if (memory_mapped && using_mapping) {
    if (munmap(SA, N * sizeof(ANINT))) throw Error("SA Memory unmap failure");
    if (munmap(ISA, N * sizeof(ANINT))) throw Error("ISA Memory unmap failure");
  } else {
    free(SA);
    free(ISA);
  }
}

// Uses the algorithm of Kasai et al 2001 which was described in
// Manzini 2004 to compute the LCP array.
void longSA::computeLCP() {
  uint64_t h = 0;
  for (uint64_t i = 0; i < N; ++i) {
    uint64_t m = ISA[i];
    if (m == 0) {
      LCP.set(m, 0);  // LCP[m]=0;
    } else {
      const uint64_t j = SA[m-1];
      while (i+h < N && j+h < N && ref[i+h] == ref[j+h]) ++h;
      LCP.set(m, h);  // LCP[m] = h;
    }
    h = max(0L, (int64_t)h - 1L);
  }
}

// Binary search for left boundry of interval.
uint64_t longSA::bsearch_left(const char c, const uint64_t i,
                                   uint64_t l, uint64_t r) const {
  if (c == ref[SA[l]+i]) return l;
  while (r > l + 1) {
    const uint64_t m = (l+r) / 2;
    if (c <= ref[SA[m] + i])
      r = m;
    else
      l = m;
  }
  return r;
}

// Binary search for right boundry of interval.
uint64_t longSA::bsearch_right(const char c, const uint64_t i,
                                    uint64_t l, uint64_t r) const {
  if (c == ref[SA[r]+i]) return r;
  while (r - l > l + 1) {
    const uint64_t m = (l+r) / 2;
    if (c < ref[SA[m] + i])
      r = m;
    else
      l = m;
  }
  return l;
}

// Simple top down traversal of a suffix array.
bool longSA::top_down(const char c, const uint64_t i,
                      uint64_t &start, uint64_t &end) const {
  if (c < ref[SA[start]+i]) return false;
  if (c > ref[SA[end]+i]) return false;
  const uint64_t l = bsearch_left(c, i, start, end);
  const uint64_t l2 = bsearch_right(c, i, start, end);
  start = l;
  end = l2;
  return l <= l2;
}

// Top down traversal of the suffix array to match a pattern.  NOTE:
// NO childtab as in the enhanced suffix array (ESA).
bool longSA::search(const string &P, uint64_t &start,
                    uint64_t &end) const {
  start = 0;
  end = N - 1;
  uint64_t i = 0;
  while (i < P.length()) {
    if (top_down(P[i], i, start, end) == false) {
      return false;
    }
    ++i;
  }
  return true;
}

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void longSA::traverse(const string &P, const uint64_t prefix,
                      interval_t &cur, const ANINT min_len) const {
  if (cur.depth >= min_len) return;

  while (prefix+cur.depth < P.length()) {
    uint64_t start = cur.start;
    uint64_t end = cur.end;
    // If we reach a mismatch, stop.
    if (top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false)
      return;

    // Advance to next interval.
    cur.depth += 1;
    cur.start = start;
    cur.end = end;

    // If we reach min_len, stop.
    if (cur.depth == min_len) return;
  }
}

// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool longSA::top_down_faster(const char c, const uint64_t i,
                             uint64_t &start, uint64_t &end) const {
  uint64_t l, r, m, r2 = end, l2 = start;
  int64_t vgl;
  bool found = false;
  const int64_t cmp_with_first = (int64_t)c - (int64_t)ref[SA[start]+i];
  const int64_t cmp_with_last = (int64_t)c - (int64_t)ref[SA[end]+i];
  if (cmp_with_first < 0) {
    l = start + 1;
    l2 = start;  // pattern doesn't occur!
  } else if (cmp_with_last > 0) {
    l = end + 1;
    l2 = end;
    // pattern doesn't occur!
  } else {
    // search for left border:
    l = start;
    r = end;
    if (cmp_with_first == 0) {
      found = true;
      r2 = r;
    } else {
      while (r > l + 1) {
        m = (l+r) / 2;
        vgl = (int64_t)c - (int64_t)ref[SA[m] + i];
        if (vgl <= 0) {
          if (!found && vgl == 0) {
            found = true;
            l2 = m;
            r2 = r;  // search interval for right border
          }
          r = m;
        } else {
          l = m;
        }
      }
      l = r;
    }
    // search for right border (in the range [l2:r2])
    if (!found) {
      l2 = l - 1;  // pattern not found => right border to the left of 'l'
    }
    if (cmp_with_last == 0) {
      l2 = end;  // right border is the end of the array
    } else {
      while (r2 > l2 + 1) {
        m = (l2 + r2) / 2;
        vgl = (int64_t)c - (int64_t)ref[SA[m] + i];
        if (vgl < 0)
          r2 = m;
        else
          l2 = m;
      }
    }
  }
  start = l;
  end = l2;
  return l <= l2;
}

// Suffix link simulation using ISA/LCP heuristic.
bool longSA::suffixlink(interval_t * m) const {
  if (m->depth <= 1) {
    m->depth = 0;
    return false;
  }
  --m->depth;
  m->start = ISA[SA[m->start] + 1];
  m->end = ISA[SA[m->end] + 1];
  return expand_link(m);
}

// For a given offset in the prefix k, find all MEMs.
void longSA::findMEM(Aligner & query) const {
  const string & P = query();
  // Offset all intervals at different start points.
  uint64_t prefix = 1;
  interval_t mli(0, N - 1, 0);  // min length interval
  interval_t xmi(0, N - 1, 0);  // max match interval

  // Right-most match used to terminate search.
  while (prefix <= P.length()) {
    // Traverse until minimum length matched.
    traverse(P, prefix, mli, query.min_len);
    if (mli.depth > xmi.depth) xmi = mli;
    if (mli.depth <= 1) {
      mli.reset(N - 1);
      xmi.reset(N - 1);
      ++prefix;
      continue;
    }

    if (mli.depth >= query.min_len) {
      traverse(P, prefix, xmi, P.length());  // Traverse until mismatch.
      collectMEMs(query, prefix, mli, xmi);  // Using LCP to find MEM length.
      // When using ISA/LCP trick, depth = depth - 1. prefix += 1.
      ++prefix;
      if ( suffixlink(&mli) == false ) {
        mli.reset(N - 1);
        xmi.reset(N - 1);
        continue;
      }
      suffixlink(&xmi);
    } else {
      // When using ISA/LCP trick, depth = depth - 1. prefix += 1.
      ++prefix;
      if ( suffixlink(&mli) == false ) {
        mli.reset(N - 1);
        xmi.reset(N - 1);
        continue; }
      xmi = mli;
    }
  }
}

// Finds left maximal matches given a right maximal match at position i.
inline void longSA::find_Lmaximal(
    Aligner & query, const uint64_t prefix,
    const uint64_t i, const uint64_t len) const {
  const string & P = query();

  // Is this function broken??? Probably

  // Advance to the left 1 step.
  // If we reach the end and the match is long enough, print.
  if (prefix == 0 || i == 0) {
    if (len >= query.min_len)
      query.process_match(match_t(i, prefix, len));
    return;  // Reached mismatch, done.
  } else if (P[prefix-1] != ref[i - 1]) {
    // If we reached a mismatch, print the match if it is long enough.
    if (len >= query.min_len)
      query.process_match(match_t(i, prefix, len));
    return;  // Reached mismatch, done.
  }
}

// Use LCP information to locate right maximal matches.
// Test each for left maximality.
void longSA::collectMEMs(Aligner & query, const uint64_t prefix,
                         const interval_t mli, interval_t xmi) const {
  // All of the suffixes in xmi's interval are right maximal.
  for (uint64_t i = xmi.start; i <= xmi.end; ++i)
    find_Lmaximal(query, prefix, SA[i], xmi.depth);

  if (mli.start == xmi.start && mli.end == xmi.end) return;

  while (xmi.depth >= mli.depth) {
    // Attempt to "unmatch" xmi using LCP information.
    if (xmi.end+1 < N)
      xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
    else
      xmi.depth = LCP[xmi.start];

    // If unmatched XMI is > matched depth from mli, then examine rmems.
    if (xmi.depth >= mli.depth) {
      // Scan RMEMs to the left, check their left maximality..
      while (LCP[xmi.start] >= xmi.depth) {
        --xmi.start;
        find_Lmaximal(query, prefix, SA[xmi.start], xmi.depth);
      }
      // Find RMEMs to the right, check their left maximality.
      while (xmi.end+1 < N && LCP[xmi.end+1] >= xmi.depth) {
        ++xmi.end;
        find_Lmaximal(query, prefix, SA[xmi.end], xmi.depth);
      }
    }
  }
}

struct by_ref {
  bool operator() (const match_t &a, const match_t &b) const {
  if (a.ref == b.ref)
    return a.len > b.len;
  else
    return a.ref < b.ref;
  }
};

// Finds maximal almost-unique matches (MAMs) These can repeat in the
// given query pattern P, but occur uniquely in the indexed reference S.
void longSA::MAM(Aligner & query) const {
  const string &P = query();
  interval_t cur(0, N - 1, 0);
  uint64_t prefix = 0;
  while (prefix < P.length()) {
    // Traverse SA top down until mismatch or full string is matched.
    traverse(P, prefix, cur, P.length());
    if (cur.depth <= 1) {
      cur.depth = 0;
      cur.start = 0;
      cur.end = N - 1;
      ++prefix;
      continue;
    }
    if (cur.size() == 1 && cur.depth >= query.min_len) {
      if (is_leftmaximal(P, prefix, SA[cur.start])) {
        // Yes, it's a MAM.
        query.process_match(match_t(SA[cur.start], prefix, cur.depth));
      }
    }
    do {
      cur.depth = cur.depth-1;
      cur.start = ISA[SA[cur.start] + 1];
      cur.end = ISA[SA[cur.end] + 1];
      ++prefix;
      if ( cur.depth == 0 || expand_link(&cur) == false ) {
        cur.depth = 0;
        cur.start = 0;
        cur.end = N - 1;
        break;
      }
    } while (cur.depth > 0 && cur.size() == 1);
  }
}

// Returns true if the position p1 in the query pattern and p2 in the
// reference is left maximal.
bool longSA::is_leftmaximal(const string &P, const uint64_t p1,
                            const uint64_t p2) const {
  if (p1 == 0 || p2 == 0)
    return true;
  else
    return P[p1-1] != ref[p2-1];
}

// Maximal Unique Match (MUM)
void longSA::MUM(Aligner & query) const {
  // Find unique MEMs.
  query.set_print(false);
  MAM(query);
  vector<match_t> matches;
  query.forget(matches);
  query.set_print(true);

  // Adapted from Stephan Kurtz's code in cleanMUMcand.c in MUMMer v3.20.
  uint64_t currentright, dbright = 0;
  bool ignorecurrent, ignoreprevious = false;
  sort(matches.begin(), matches.end(), by_ref());
  for (uint64_t i = 0; i != matches.size(); ++i) {
    ignorecurrent = false;
    currentright = matches[i].ref + matches[i].len - 1;
    if (dbright > currentright) {
      ignorecurrent = true;
    } else {
      if (dbright == currentright) {
        ignorecurrent = true;
        if (!ignoreprevious && matches[i-1].ref == matches[i].ref)
          ignoreprevious = true;
      } else {
        dbright = currentright;
      }
    }
    if (i > 0 && !ignoreprevious) {
      query.process_match(matches[i-1]);
    }
    ignoreprevious = ignorecurrent;
  }
  if (!ignoreprevious) {
    if (matches.size() > 0) {
      query.process_match(matches[matches.size()-1]);
    }
  }
}

void longSA::MEM(Aligner & query) const {
  if (query.min_len < 1) return;
  findMEM(query);
}

class BinWriter {
 public:
  explicit BinWriter(const string & filename) {
    outfile = fopen(filename.c_str(), "wb");
    if (outfile == nullptr)
      throw Error("could not open outfile") << filename << "for writing";
  }
  ~BinWriter() {
    if (fclose(outfile) != 0)
      Error("problem closing BinWriter outfile");
  }
  FILE * outfile;
};
template <class Out>
BinWriter & operator<<(BinWriter & writer, const Out out) {
  bwritec(writer.outfile, &out, "some binary value", 1);
  return writer;
}


template <class Out>
void longSA::show(Out & output, const bool bin = false) const {
  if (N > 10000000000) cerr << "trying to allocate "
                            << N * sizeof(uint64_t) / 1000000000
                            << " GB - make take a while!" << endl;
  output << "test" << "\n";
  vector<uint64_t> min_lengths(N, 0);

  const uint64_t line = 72;
  const uint64_t notify = std::max(ref.sizes[0] / line, 1UL);

  if (verbose) cerr << "# computing mappability (each dot represents "
                    << notify << " bases and line length is first chromosome)"
                    << endl;
  uint64_t ref_seq = 0;
  if (verbose) cerr << "# " << ref.descr[ref_seq] << " ";
  for (uint64_t i = 0; i != N; ++i) {
    min_lengths[i] = LCP[i] + 1;
    if (i) {
      min_lengths[i - 1] = max(min_lengths[i - 1], min_lengths[i]);
      if (i % notify == 0) {
        while ((ref_seq + 1) != ref.descr.size() &&
               i > ref.startpos[ref_seq + 1]) {
          ++ref_seq;
          if (verbose) cerr << endl << "# " << ref.descr[ref_seq] << " ";
        }
        if (verbose) cerr << "." << std::flush;
      }
    }
  }
  if (verbose) cerr << endl;

  const bool show_all = false;
  if (!bin) output << "chrom\tpos\tlmin\trmin\n";
  if (show_all) cerr << "# pos\tchrom\ti\tref\tsa\tisa\tlcp\tlcp+\tlmin\trmin"
                     << endl;
  if (verbose) cerr << "# outputting mappability" << endl;

  for (uint64_t chrom = 0; chrom != ref.startpos.size(); chrom += 2) {
    const string & name = ref.descr[chrom];
    const uint64_t startpos = ref.startpos[chrom];
    const uint64_t size_ = ref.sizes[chrom];
    if (verbose) cerr << "# " << name << " ";
    for (uint64_t i = 0; i != size_; ++i) {
      if (verbose && i && i % notify == 0) {
        cerr << "." << std::flush;
      }
      const uint64_t pos = i + startpos;
      const uint64_t sapos = ISA[pos];
      const uint64_t rcsapos = ISA[startpos + 2 * size_ - i];
      if (1) {
        if (sapos >= ref.N) throw Error("range error sa");
        if (rcsapos >= ref.N) throw Error("range error rcsa");
      }
      if (min_lengths[sapos] + i >= size_) min_lengths[sapos] = 0;
      if (min_lengths[rcsapos] >= i) min_lengths[rcsapos] = 0;
      if (show_all)
        cerr << "# " << pos << "\t" << name << "\t" << i + 1 << "\t"
             << ref[pos] << "\t" << SA[pos] << "\t" << ISA[pos] << "\t"
             << LCP[sapos] << "\t"
             << (sapos + 1 < ref.N ? LCP[sapos+1] : 0UL) << "\t"
             << min_lengths[rcsapos] << "\t"
             << min_lengths[sapos] << endl;
      if (bin) {
        output << static_cast<char>(
            min_lengths[rcsapos] < 255 ? min_lengths[rcsapos] : 255);
        output << static_cast<char>(
            min_lengths[sapos] < 255 ? min_lengths[sapos] : 255);
      } else {
        output << name << "\t" << i + 1 << "\t" << min_lengths[rcsapos] << "\t"
               << min_lengths[sapos] << "\n";
        cerr << name << "\t" << i + 1 << "\t" << min_lengths[rcsapos] << "\t"
               << min_lengths[sapos] << "\n";
      }
    }
    if (verbose) cerr << "\n";
  }
  exit(0);
}

void longSA::show_mappability(const string & filename) const {
  if (filename == "-") {
    show(cout);
  } else if (filename.find(".bin") != string::npos) {
    BinWriter binout(filename);
    show(binout, true);
  } else {
    cerr << "output to output " << filename << endl;
    ofstream output(filename.c_str());
    ostream & out = output;
    out << "output to output" << endl;
    show(out);
    output << "done output to output" << endl;
    output << endl;
  }
}
