/* Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_QUERY_H_
#define LONGMEM_QUERY_H_

#include <algorithm>
#include <string>
#include <vector>

#include "./locked.h"
#include "./longSA.h"
#include "./util.h"
#include "./memsam.h"

// #define use_ctr 1

class Aligner;
struct Alignment {
  int64_t rcpos;  // position in ref + rcref of unflipped query start
  int64_t pos;  // position in chromosome of query start
  int64_t qpos;  // position in query of first hit start
  uint64_t seq_index;  // index of chromosome sequence
  uint64_t prefix;  // length before match
  uint64_t length;  // length of match
  uint64_t suffix;  // length after match
  uint64_t n_matches;  // how many MAMs in read
  uint64_t n_unique_bases;  // n bases in MAMs that match the ref
  uint64_t n_matched_bases;  // how many bases match reference
  uint64_t alignment_index;  // index of alignment
  Alignment * prev_alignment;  // linked list of alignments
  Alignment * next_alignment;  // linked list of alignments
  Alignment * best_mate;  // not the best, just preferred
  std::vector<char> cigar;  // cigar string
  bool rc;  // was query mapped to reversed reference?
#ifdef use_ctr
  Alignment();
#else
  void reset();
  void partial_reset();
#endif
  void resolve(const match_t & input_match,
               const uint64_t query_length,
               const Sequence & ref);
  void print(const Aligner &);
};

struct Query {
  void clear();
  void swap(Query & other);
  std::string name;
  std::string query;
  std::string original;
  std::string errors;
  std::string optional;
};

class NewQueryArgs {
 public:
  NewQueryArgs() : nucleotides_only(false) {}
  bool nucleotides_only;
 private:
  NewQueryArgs & operator=(const NewQueryArgs & disabled_assignment_operator);
};

class NewQuery : public Query, public NewQueryArgs {
 public:
  explicit NewQuery(const NewQueryArgs & args);
  void reset(std::string & name_);
  void extend(const std::string & line);
  void set_errors(std::string & line);
  void add_optional(std::string & opt);
  bool complete() const;

 private:
  NewQuery & operator=(const NewQuery & disabled_assignment_operator);
};


class OutputSorter {
 public:
  OutputSorter(const std::string header_,
               const uint64_t buffer_size_ = 500000000,
               const uint64_t max_line_size_ = 10000) :
      header(header_), buffer_size(buffer_size_), max_line_size(max_line_size_),
      buffer(buffer_size) {}
  void end_line() {
    buffer[end++] = '\0';
    // cerr << "end line " << begin << " " << end << endl;
    output_reads.emplace_back(&buffer[begin]);
    // Prepare for next line
    begin = end;
    const uint64_t needed_size = end + max_line_size;
    if (needed_size > buffer_size) {
      flush();
    }
  }
  void flush();
  void printf(const char * format, ...) {
    va_list args;
    va_start(args, format);
    end += vsprintf(&buffer[end], format, args);
    va_end(args);
  }

 private:
  std::string header;
  uint64_t buffer_size;
  uint64_t max_line_size;
  uint64_t begin = 0;
  uint64_t end = 0;
  uint32_t file_sequence = 0;
  std::vector<char> buffer;
  std::vector<MemSam> output_reads;
};

class OutputArgs {
 public:
  OutputArgs() : sam_out(false), nomap(false) {}
  bool sam_out;
  bool nomap;

 private:
  OutputArgs & operator=(const OutputArgs & disabled_assignment_operator);
};

enum mum_t { MUM, MAM, MEM };
class AlignerArgs : public OutputArgs {
 public:
  AlignerArgs() : OutputArgs(), type(MAM), min_len(20), min_block(20) {}
  mum_t type;
  unsigned int min_len;
  unsigned int min_block;
 private:
  AlignerArgs & operator=(const AlignerArgs & disabled_assignment_operator);
};

class Aligner : public Query, public AlignerArgs {
 public:
  // Used by PairRunner class
  Aligner(const AlignerArgs & args, const longSA & sa_);
  Aligner(const Aligner & other);
  void clear();
  void reset(NewQuery & new_query);
  void set_nomap();
  void run();
  void print_matches(OutputSorter & output);
  bool has_mate(const Aligner & read2) const;
  void set_mate(const Aligner & other);
  // Used by longSA class
  const std::string & operator()() const { return query; }
  void process_match(const match_t & match);
  void forget(std::vector<match_t> & matches_);
  void set_print(const bool print_);

 private:
  void prepare_matches();
  const longSA & sa;
  std::string rcquery;
  bool print;
  unsigned int read_flag;
  Alignment * best_alignment;
  std::vector<match_t> matches;
  std::vector<Alignment> alignments;
  std::vector<Alignment *> sorted_alignments;
  uint64_t n_alignments;
  Aligner & operator=(const Aligner & disabled_assignment_operator);
};

class PairArgs : public AlignerArgs, public NewQueryArgs {
 public:
  PairArgs() : AlignerArgs(), NewQueryArgs() {}
 private:
  PairArgs & operator=(const PairArgs & disabled_assignment_operator);
};

class Pair : public PairArgs {
 public:
  Pair(const PairArgs & args, const longSA & sa_);
  Pair(const Pair & other);
  static void * runner_thread(void * obj);
  RingBuffer<NewQuery> queue;
  uint64_t n_queries;
 private:
  void run();
  Aligner read1;
  Aligner read2;
  OutputSorter output;
  Pair & operator=(const Pair & disabled_assignment_operator);
};

class PairsArgs : public PairArgs {
 public:
  PairsArgs() : PairArgs(), max_n_threads(2), n_threads(1), verbose(false) {}
  unsigned int max_n_threads;
  unsigned int n_threads;
  bool verbose;
 private:
  PairArgs & operator=(const PairArgs & disabled_assignment_operator);
};

class Pairs : public PairsArgs {
 public:
  Pairs(const PairsArgs & args, const longSA & sa);
  ~Pairs();
  Pair * get_pair();
  void switch_pair(Pair * & pair);
  void release_pair(Pair * pair);
  void done();
 private:
  const time_t start_time;
  std::vector<pthread_t> thread_ids;
  std::vector<Pair> pairs;
  RingBuffer<Pair *> available;
  Pairs(const Pairs & disabled_copy_constructor);
  Pairs & operator=(const Pairs & disabled_assignment_operator);
};

class ReaderArgs {
 public:
  ReaderArgs() : query_input(nullptr), fastq(false),
                 sam_in(false), verbose(false) {}
  const char * query_input;
  bool fastq;
  bool sam_in;
  bool verbose;
 private:
  ReaderArgs & operator=(const ReaderArgs & disabled_assignment_operator);
};

class QueryReader : public ReaderArgs {
 public:
  QueryReader(const ReaderArgs & args, Pairs & pairs_,
              const char * const query_input_);
  ~QueryReader();
  void run();
  void run_in_thread();
  void join();
  void print_sam_header();
  uint64_t n_sequences;
 private:
  static void * fasta_thread(void * obj);
  Pairs & pairs;
  uint64_t n_full;
  pthread_attr_t attr;
  pthread_t thread_id;
  QueryReader & operator=(const QueryReader & disabled_assignment_operator);
};

class ReadersArgs : public ReaderArgs {
 public:
  ReadersArgs() : n_input(0), input(nullptr) {}
  unsigned int n_input;
  char * * input;
 private:
  ReadersArgs & operator=(const ReadersArgs & disabled_assignment_operator);
};

class Readers : public ReadersArgs {
 public:
  Readers(const ReadersArgs & args, Pairs & pairs_);
  ~Readers();
 private:
  const time_t start_time;
  Pairs & pairs;
  std::vector<QueryReader> readers;
};

#endif  // LONGMEM_QUERY_H_
