/* Copyright Peter Andrews 2013 CSHL */

#include "./query.h"

#include <pthread.h>

#include <exception>
using std::exception;

#include <vector>
using std::vector;

#include <algorithm>
#include <string>
using std::string;
using std::swap;

#include <sstream>
using std::istringstream;
using std::ostringstream;

#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "./util.h"
#include "./error.h"
using paa::Error;

#ifdef use_ctr

// Alignment
Alignment::Alignment()
    : rcpos(0), pos(0), qpos(0), seq_index(0), prefix(0), length(0), suffix(0),
      n_matches(0), n_unique_bases(0), n_matched_bases(0), alignment_index(0),
      prev_alignment(0), next_alignment(0), best_mate(0), cigar("*"),
      rc(false) { }
#else
void Alignment::partial_reset() {
  n_matches = 0;
  n_unique_bases = 0;
  n_matched_bases = 0;
  alignment_index = 0;
  prev_alignment = nullptr;
  next_alignment = nullptr;
  best_mate = nullptr;
}
void Alignment::reset() {
  partial_reset();
  rcpos = 0;
  pos = 0;
  qpos = 0;
  seq_index = 0;
  prefix = 0;
  length = 0;
  suffix = 0;
  cigar.resize(2);
  cigar[0] = '*';
  cigar[1] = 0;
  rc = false;
}
#endif

inline void Alignment::resolve(const match_t & input_match,
                               const uint64_t query_length,
                               const Sequence & ref) {
#ifndef use_ctr
  partial_reset();
#endif
  vector<uint64_t>::const_iterator it =
      upper_bound(ref.startpos.begin(), ref.startpos.end(), input_match.ref);
  seq_index = it - ref.startpos.begin() - 1;
  rcpos = input_match.ref - input_match.query;
  pos = rcpos - *--it;
  const unsigned int extra = query_length - input_match.len - input_match.query;
  if (ref.rcref && ((seq_index % 2) == 1)) {
    seq_index -= 1;
    pos = ref.sizes[seq_index] - pos;  // Check this!
    pos -= query_length;
    prefix = extra;
    suffix = input_match.query;
    rc = true;
  } else {
    prefix = input_match.query;
    suffix = extra;
    rc = false;
  }
  qpos = input_match.query;
  length = input_match.len;
  cigar.resize(query_length * 5);
  cigar[0] = '*';
  cigar[1] = 0;
}

inline void Alignment::print(const Aligner &) { }

inline void Query::clear() {
  name.clear();
  query.clear();
  original.clear();
  errors.clear();
  optional.clear();
}
inline void Query::swap(Query & other) {
  name.swap(other.name);
  original.swap(other.original);
  query.swap(other.query);
  errors.swap(other.errors);
  optional.swap(other.optional);
}

// NewQuery
NewQuery::NewQuery(const NewQueryArgs & args)
    : NewQueryArgs(args) {}

inline void NewQuery::reset(string & name_) {
  clear();
  name.swap(name_);
}

inline void NewQuery::extend(const string & line) {
  uint64_t end = line.size();
  while (end && line[end - 1] == ' ') --end;
  for (uint64_t i = 0; i != end; ++i) {
    if (line[i] == ' ') continue;
    const char c = tolower(line[i]);
    if (nucleotides_only) {
      switch (c) {
        case 'a': case 't': case 'g': case 'c':
          query.push_back(c);
          break;
        default:
          query.push_back('~');
      }
    } else {
      query.push_back(c);
    }
    original.push_back(line[i]);
  }
}

inline void NewQuery::set_errors(string & line) {
  if (line.empty()) throw Error("empty errors");
  errors.swap(line);
}

inline void NewQuery::add_optional(string & opt) {
  optional += "\t";
  optional += opt;
}

inline bool NewQuery::complete() const {
  return query.size();
}

// Aligner
Aligner::Aligner(const AlignerArgs & args, const longSA & sa_)
    : Query(), AlignerArgs(args), sa(sa_), rcquery(""),
      print(true), read_flag(0), best_alignment(nullptr), matches(0),
      alignments(0), sorted_alignments(0), n_alignments(0) {}
Aligner::Aligner(const Aligner & other)
    : Query(other), AlignerArgs(other),
      sa(other.sa), rcquery(other.rcquery),
      print(other.print), read_flag(other.read_flag),
      best_alignment(other.best_alignment), matches(other.matches),
      alignments(other.alignments),
      sorted_alignments(other.sorted_alignments),
      n_alignments(other.n_alignments) {}

inline void Aligner::clear() {
  Query::clear();
  rcquery.clear();
  print = true;
  read_flag = 0;
  best_alignment = nullptr;
  matches.clear();
  alignments.clear();
  sorted_alignments.clear();
  n_alignments = 0;
}
inline void Aligner::reset(NewQuery & new_query) {
  clear();
  Query::swap(new_query);
  rcquery.swap(new_query.errors);
  if (name.size() >= 2) {  // Look for :0 or :1 for mate info
    const uint64_t pos = name.size() - 2;
    if (name[pos] == ':') {
      if (name[pos + 1] == '0') {
        name.resize(pos);
        read_flag = is_paired | is_first;
      } else if (name[pos + 1] == '1') {
        name.resize(pos);
        read_flag = is_paired | is_second;
      }
    }
  }
}

struct to_merge {
  bool operator() (const Alignment * a, const Alignment * b) const {
    if (a->rc == b->rc) {
      if (a->seq_index == b->seq_index) {
        if (a->pos == b->pos) {
          return a->prefix < b->prefix;
        } else {
          return a->pos < b->pos;
        }
      } else {
        return a->seq_index < b->seq_index;
      }
    } else {
      return a->rc < b->rc;
    }
  }
};

struct to_print {
  bool operator() (const Alignment * a, const Alignment * b) const {
    if (a->qpos == b->qpos) {
      return a->rc < b->rc;
    } else {
      return a->qpos < b->qpos;
    }
  }
};

inline void Aligner::prepare_matches() {
  n_alignments = 0;
  best_alignment = nullptr;
  if (matches.size()) {
    alignments.resize(matches.size());
    for (uint64_t i = 0; i != alignments.size(); ++i) {
      alignments[i].resolve(matches[i], query.size(), sa.ref);
    }
    for (uint64_t i = alignments.size() - 1; ; --i) {
      // Remove off-chromosome mappings with negative positions
      if (alignments[i].pos < 0) {
        matches.erase(matches.begin() + i);
        alignments.erase(alignments.begin() + i);
      }
      if (i == 0) break;
    }
    for (uint64_t i = 0; i != alignments.size(); ++i) {
      sorted_alignments.push_back(&alignments[i]);
    }
    if (matches.size() && sam_out) {
      sort(sorted_alignments.begin(), sorted_alignments.end(), to_merge());
      uint64_t cigar_end = 0;
      uint64_t last_end = 0;
      for (uint64_t i = 0; i != alignments.size(); ++i) {
        Alignment * const a = sorted_alignments[i];
        Alignment * const na = (i + 1 == alignments.size()) ?
            nullptr : sorted_alignments[i + 1];
        ++a->n_matches;
        a->n_unique_bases += a->length;
        if (a->prefix)
          cigar_end += sprintf(&a->cigar[cigar_end], "%lu%c",
                               a->prefix - last_end,
                               last_end ? 'M' : 'S');
        cigar_end += sprintf(&a->cigar[cigar_end], "%lu=", a->length);
        if (!na || na->pos != a->pos || na->seq_index != a->seq_index ||
            na->rc != a->rc) {
          if (a->suffix)
            cigar_end += sprintf(&a->cigar[cigar_end], "%luS", a->suffix);

          for (uint64_t j = 0; j != query.size(); ++j) {
            const int64_t ref_pos = a->rcpos + j;
            if (ref_pos >= 0L && ref_pos < (int64_t)sa.size() &&
                sa.ref[ref_pos] == query[j]) ++a->n_matched_bases;
          }

          a->cigar.resize(cigar_end);
          cigar_end = 0;
          last_end = 0;
        } else {
          last_end = a->prefix + a->length;
          a->cigar.swap(na->cigar);
          // ::swap(a->qpos, na->qpos);
          na->qpos = std::min(a->qpos, na->qpos);
          ::swap(a->n_matches, na->n_matches);
          a->n_matches = 0;
          ::swap(a->n_unique_bases, na->n_unique_bases);
          a->n_matched_bases = 0;
        }
      }
      sort(sorted_alignments.begin(), sorted_alignments.end(), to_print());
      best_alignment = sorted_alignments.front();
      Alignment * previous_alignment = nullptr;
      for (uint64_t i = 0; i != alignments.size(); ++i) {
        Alignment * const a = sorted_alignments[i];
        if (a->n_matches) {
          a->alignment_index = n_alignments++;
          if (previous_alignment) {
            a->prev_alignment = previous_alignment;
            previous_alignment->next_alignment = a;
          }
          previous_alignment = a;
        }
      }
    }
  }
}

inline void Aligner::set_nomap() {
  if (n_alignments == 0 && sam_out && nomap) {
    ++n_alignments;
    read_flag = read_flag | is_unmapped;
    alignments.resize(1);
    Alignment & a = alignments.back();
#ifndef use_ctr
    a.reset();
#endif
    sorted_alignments.push_back(&a);
    // best_alignment = sorted_alignments.front();
  }
}

inline void Aligner::run() {
  if (errors.empty()) errors.assign(query.size(), '!');
  if (type == MAM) sa.MAM(*this);
  else if (type == MUM) sa.MUM(*this);
  else if (type == MEM) sa.MEM(*this);
  prepare_matches();
  set_nomap();
}

inline void Aligner::print_matches(OutputSorter & output) {
  if (alignments.size()) {
    for (uint64_t i = 0; i != alignments.size(); ++i) {
      const Alignment * a = sorted_alignments[i];
      if (sorted_alignments.size() == 0 ||
          alignments.size() != sorted_alignments.size()) {
        cerr << "Size problem" << endl;
        exit(1);
      }
      if (a == nullptr) {
        cerr << "Alignment problem" << endl;
        exit(1);
      }
      if (sam_out) {
        if (a->n_matches || (read_flag & is_unmapped)) {
          if (read_flag & is_unmapped)
            output.printf("%s\t%u\t%s\t%ld\t0\t%s",
                          name.c_str(), read_flag,
                          (a->best_mate ?
                           sa.ref.descr[a->best_mate->seq_index].c_str() :
                           "*"),
                          (a->best_mate ? a->best_mate->pos + 1 : 0), "*");
          else
            output.printf("%s\t%u\t%s\t%ld\t50\t%s",
                          name.c_str(), read_flag | (a->rc ? is_reversed : 0) |
                          (a->alignment_index ? is_not_primary : 0),
                          sa.ref.descr[a->seq_index].c_str(),
                          a->pos + 1, &a->cigar[0]);
          if (a->best_mate)
            output.printf("\t%s\t%ld\t0",
                          sa.ref.descr[a->best_mate->seq_index].c_str(),
                          a->best_mate->pos + 1);
          else
            output.printf("\t*\t0\t0");
          if (a->rc) {
            if (rcquery.empty()) {
              rcquery = original;
              reverse_complement(&rcquery);
            }
            output.printf("\t%s\t", rcquery.c_str());
            for (uint64_t j = errors.size(); j != 0 ; )
              output.printf("%c", errors[--j]);
          } else {
            output.printf("\t%s\t%s",
                           original.c_str(), errors.c_str());
          }
          if (a->n_matches) {
            output.printf("\tXM:i:%lu\tXU:i:%lu\tXE:i:%lu"
                          "\tXS:A:%c\tNH:i:%lu\tHI:i:%lu",
                          a->n_matches, a->n_unique_bases,
                          a->n_matched_bases,
                          a->rc ? '-' : '+', n_alignments, a->alignment_index);
          } else {
            output.printf("\tXM:i:0\tNH:i:0");
          }
          if (a->prev_alignment && a->prev_alignment != a) {
            const Alignment * const prev = a->prev_alignment;
            output.printf("\tcc:Z:%s\tcp:i:%ld\txo:A:%c\txc:Z:%s",
                          sa.ref.descr[prev->seq_index].c_str(),
                          prev->pos + 1, prev->rc == a->rc ? '=' : '!',
                          &prev->cigar[0]);
          }
          if (a->next_alignment && a->next_alignment != a) {
            const Alignment * const next = a->next_alignment;
            output.printf("\tCC:Z:%s\tCP:i:%ld\tXO:A:%c\tXC:Z:%s",
                          sa.ref.descr[next->seq_index].c_str(),
                          next->pos + 1, next->rc == a->rc ? '=' : '!',
                          &next->cigar[0]);
          }
          if (optional.size()) output.printf("%s", optional.c_str());
          output.printf("\n");
          output.end_line();
        }
      } else {
        output.printf("> %s\n  %s", name.c_str(),
                       sa.ref.descr[a->seq_index].c_str());
        for (uint64_t j = 0; j != sa.ref.maxdescrlen -
                sa.ref.descr[a->seq_index].length() + 1; ++j)
          output.printf(" ");
        output.printf(" %8ld  %8ld  %8ld %d\n", a->pos + 1,
                       a->prefix + 1, a->length, a->rc ? 1 : 0);
      }
    }
  }
}

inline bool Aligner::has_mate(const Aligner & read2) const {
  return (read_flag & is_first) && (read2.read_flag & is_second);
}

inline void Aligner::set_mate(const Aligner & other) {
  if (n_alignments && other.n_alignments) {
    if (other.best_alignment) {
      for (uint64_t i = 0; i != alignments.size(); ++i) {
        alignments[i].best_mate = other.best_alignment;
      }
    } else {
      read_flag = read_flag | is_mate_unmapped;
      for (uint64_t i = 0; i != alignments.size(); ++i) {
        alignments[i].best_mate = best_alignment;
      }
    }
  }
}

void Aligner::process_match(const match_t & match) {
  matches.push_back(match);
}

void Aligner::forget(vector<match_t> & matches_) {
  matches.swap(matches_);
}

void Aligner::set_print(const bool print_) { print = print_; }


// OutputSorter
void OutputSorter::flush() {
  if (end) {
    mkdir("mapout");
    sort(output_reads.begin(), output_reads.end());
    ostringstream file_name;
    file_name << "mapout/mapout" << (uint64_t)this << "."
              << ++file_sequence << ".txt";
    FILE * out = fopen(file_name.str().c_str(), "wb");
    if (out == nullptr) {
      throw Error("Problem opening out file");
    }
    fputs(header.c_str(), out);
    for (const auto read : output_reads) {
      fputs(read.data, out);
    }
    fclose(out);
    begin = 0;
    end = 0;
    output_reads.clear();
  }
}

// Pair
Pair::Pair(const PairArgs & args, const longSA & sa_)
    : PairArgs(args), queue(1000UL, NewQuery(args), true, false, true),
      n_queries(0), read1(args, sa_), read2(args, sa_),
      output(sa_.ref.sam_header()) {
}
Pair::Pair(const Pair & other)
    : PairArgs(other), queue(other.queue), n_queries(other.n_queries),
      read1(other.read1), read2(other.read2), output(other.output) {
}

inline void Pair::run() {
  NewQuery * new_query = nullptr;
  uint64_t n_print = 0;
  uint64_t not_ready = 0;
  while (1) {
    while ((new_query = queue.first()) == nullptr) ++not_ready;
    const bool second = n_queries++ % 2;
    Aligner & read = second ? read2 : read1;
    if (new_query->complete()) {
      if (0) cerr << "got query " << new_query->name
                  << " with size " << new_query->query.size() << endl;
      read.reset(*new_query);
      queue.yield();
      read.run();
    } else {
      // cerr << "# incomplete query - finished" << endl;
      --n_queries;
      queue.yield();
      break;
    }
    if (second) {
      if (read1.has_mate(read2)) {
        read1.set_mate(read2);
        read2.set_mate(read1);
      }
      read1.print_matches(output);
      read2.print_matches(output);
      read1.clear();
      read2.clear();
      n_print += 2;
    }
  }
  if ((n_queries % 2) == 1) {
    read1.print_matches(output);
    read1.clear();
    ++n_print;
  }
  output.flush();
  if (not_ready) cerr << "No query available " << not_ready << " times" << endl;
}

void * Pair::runner_thread(void * obj) {
  try {
    reinterpret_cast<Pair *>(obj)->run();
  }
  catch(exception & e) {
    cerr << e.what() << endl;
    exit(1);
  }
  catch(...) {
    cerr << "Some exception was caught." << endl;
    exit(1);
  }
  pthread_exit(nullptr);
}

Pairs::Pairs(const PairsArgs & args, const longSA & sa)
    : PairsArgs(args), start_time(time(nullptr)), thread_ids(n_threads),
      pairs(n_threads, Pair(args, sa)),
      available(n_threads, nullptr, true, false, true) {
  if (verbose)
    cerr << "# running " << n_threads << " thread"
         << (n_threads > 1 ? "s" : "") << " to answer queries" << endl;
  // if (sam_out) sa.ref.print_sam_header();
  // Set up chromosome map for absolute position determination
  uint64_t offset = 0;
  for (unsigned int i = 0; i < sa.ref.descr.size(); i += 2) {
    // cerr << sa.ref.descr[i] << " " << offset << endl;
    MemSam::chromosomes[sa.ref.descr[i]] = offset;
    offset += sa.ref.sizes[i];
  }
  MemSam::chromosomes["*"] = offset;
  // Initialize pairs and start threads
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  for (unsigned int thread = 0; thread != n_threads; ++thread) {
    if (pthread_create(&thread_ids[thread], &attr, &Pair::runner_thread,
                       &pairs[thread]))
      throw Error("Problem creating runner_thread") << thread;
    release_pair(&pairs[thread]);
  }
}

Pairs::~Pairs() {
  uint64_t n_processed = 0;
  for (unsigned int t = 0; t != n_threads; ++t) {
    pthread_join(thread_ids[t], nullptr);
    n_processed += pairs[t].n_queries;
  }
  if (verbose) cerr << "# ran " << n_processed << " queries in "
                    << time(nullptr) - start_time << " seconds" << endl;
}

inline Pair * Pairs::get_pair() {
  return available.pop();
}

inline void Pairs::switch_pair(Pair * & pair) {
  Pair * temp = available.pop();
  available.push(pair);
  pair = temp;
}

inline void Pairs::release_pair(Pair * pair) {
  available.push(pair);
}

inline void Pairs::done() {
  NewQuery * query = nullptr;
  for (vector<Pair>::iterator pair = pairs.begin(); pair != pairs.end();
       ++pair) {
    while ((query = pair->queue.next(0)) == nullptr) continue;
    query->clear();
    pair->queue.commit();
  }
}
QueryReader::QueryReader(const ReaderArgs & args, Pairs & pairs_,
                         const char * const query_input_ = nullptr) :
    ReaderArgs(args), n_sequences(0), pairs(pairs_), n_full(0) {
  query_input = query_input_;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
}

QueryReader::~QueryReader() {
  pthread_attr_destroy(&attr);
  if (query_input && verbose) {
    cerr << "# query reader for " << query_input << " processed "
         << n_sequences << " sequences" << endl;
  }
}

inline void QueryReader::run() {
  ifstream data(query_input);
  if (!data) throw Error("unable to open") << query_input;

  Pair * pair = pairs.get_pair();
  NewQuery * query = nullptr;
  const char start_char = fastq ? '@' : '>';
  string meta;
  string line;
  string name, ref, pos, mapq, cigar, mref, mpos, tlen, seq, errors, optional;
  unsigned int flag;
  while (data) {
    getline(data, line);
    if (line.size()) {
      // Get next query from buffer
      if ((n_sequences++ % 2) == 0) {
        while ((query = pair->queue.next()) == nullptr) {
          // Get available pair to own if pair query buffer is full
          pairs.switch_pair(pair);
          ++n_full;
        }
      } else {
        ++query;  // Assumes NewQuery ring buffer has even size
      }
      // Start of new sequence
      if (sam_in) {
        istringstream input(line);
        input >> name >> flag >> ref >> pos >> mapq >> cigar
              >> mref >> mpos >> tlen >> seq >> errors;
        if (flag & is_first) name += ":0";
        else if (flag & is_second) name += ":1";
        query->reset(name);
        query->extend(seq);
        query->set_errors(errors);
        while (input >> optional) query->add_optional(optional);
      } else {
        uint64_t start = 0, end = line.length();
        if (line[0] != start_char)
          throw Error("missing query start character") << start_char <<
              "in input line" << line;

        // Process query name
        trim(line, start = 1, end);
        for (uint64_t i = start; start < end && i < end; ++i) {
          if (line[i] == ' ') {  // done with name
            if (i + 1 != end) {  // look for illumina mate info
              if (line[i+1] == '1') {
                meta += ":0";
              } else if (line[i+1] == '2') {
                meta += ":1";
              }
            }
            break;
          }
          meta += line[i];
        }
        query->reset(meta);
        // Collect sequence data.
        getline(data, line);
        if (line.empty()) throw Error("empty sequence");
        query->extend(line);
        if (fastq) {  // Assumes one line for errors, no spaces
          getline(data, line);
          getline(data, line);
          query->set_errors(line);
        }
      }
      pair->queue.commit();
    }
  }
  pairs.release_pair(pair);

  if (!n_sequences) cerr << "# no reads processed" << endl;
}

void * QueryReader::fasta_thread(void * obj) {
  try {
    reinterpret_cast<QueryReader *>(obj)->run();
  }
  catch(exception & e) {
    cerr << e.what() << endl;
    exit(1);
  }
  catch(...) {
    cerr << "Some exception was caught." << endl;
    exit(1);
  }
  pthread_exit(nullptr);
}

void QueryReader::run_in_thread() {
  if (pthread_create(&thread_id, &attr, &fasta_thread, this))
    throw Error("Problem creating fasta_thread");
}

void QueryReader::join() {
  pthread_join(thread_id, nullptr);
}

Readers::Readers(const ReadersArgs & args, Pairs & pairs_)
    : ReadersArgs(args), start_time(time(nullptr)), pairs(pairs_),
      readers(args.n_input, QueryReader(args, pairs_)) {
  if (verbose) cerr << "# running " << readers.size() << " query reader"
                    << (readers.size() > 1 ? "s" : "") << endl;
  for (unsigned int r = 0; r != readers.size(); ++r) {
    readers[r].query_input = args.input[r];
    readers[r].run_in_thread();
  }
}
Readers::~Readers() {
  // Wait for reading threads to finish
  uint64_t n_sequences = 0;
  for (unsigned int r = 0; r != readers.size(); ++r) {
    readers[r].join();
    n_sequences += readers[r].n_sequences;
  }
  pairs.done();
  readers.clear();
  const time_t elapsed_time = time(nullptr) - start_time;
  if (verbose && n_sequences)
    cerr << "# read " << n_sequences << " queries in "
         << elapsed_time << " seconds";
  if (verbose && elapsed_time) cerr << " for a rate of "
                                    << 60E-6 * n_sequences / elapsed_time
                                    << " million queries per minute";
  cerr << endl;
}
