/* Modifications of sparseMEM Copyright Peter Andrews 2013 CSHL */

#include <unistd.h>
#include <limits.h>
#include <getopt.h>

#include <exception>
using std::exception;

#include <iostream>
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include "./longSA.h"
#include "./query.h"
#include "./util.h"
#include "./error.h"
using paa::Error;

extern char verbose;

// Class to handle options parsing and to deliver options to classes
class Args : public SAArgs, public PairsArgs, public ReadersArgs {
 public:
  Args(int argc_, char * argv_[], char * envp_[]);
  bool verbose;
 private:
  void check_integer_sizes() const;
  void usage(const string & prog) const;
  int argc;
  char * * argv;
  char * * envp;
  Args(const Args & disabled_copy_constructor);
  Args & operator=(const Args & disabled_assignment_operator);
};

// mummer
int main(int argc, char* argv[], char * envp[]) {
  try {
    // Process options and args.
    const Args args(argc, argv, envp);

    // Create suffix array
    const longSA sa(args);

    // Optionally compute mappability
    if (args.mappability) {
      sa.show_mappability(args.input[0]);
      return 0;
    }

    // Pairs manages Pair worker threads
    Pairs pairs(args, sa);

    // Readers read queries and pass them off to a Pair in Pairs
    Readers readers(args, pairs);
  }
  catch(exception & e) {
    cerr << "Error" << endl;
    cerr << e.what() << endl;
    return 1;
  }
  catch(...) {
    cerr << "Some exception was caught." << endl;
    return 1;
  }
  return 0;
}

Args::Args(int argc_, char * argv_[], char * envp_[])
    : SAArgs(), PairsArgs(), ReadersArgs(),
      verbose(false), argc(argc_), argv(argv_), envp(envp_) {
  // Collect arguments from the command line. These options are allowed.
  struct option long_options[] = {
    {"l", 1, nullptr, 0},  // 0
    {"mumreference", 0, nullptr, 0},  // 1
    {"maxmatch", 0, nullptr, 0},  // 2
    {"mum", 0, nullptr, 0},  // 3
    {"mumcand", 0, nullptr, 0},   // 4
    {"n", 0, nullptr, 0},  // 5
    {"qthreads", 1, nullptr, 0},  // 6
    {"samout", 0, nullptr, 0},  // 7
    {"verbose", 0, nullptr, 0},  // 8
    {"nomap", 0, nullptr, 0},  // 9
    {"rcref", 0, nullptr, 0},  // 10
    {"fastq", 0, nullptr, 0},  // 11
    {"samin", 0, nullptr, 0},  // 12
    {"mappability", 0, nullptr, 0},  // 13
    {"cached", 0, nullptr, 0},  // 14
    {"normalmem", 0, nullptr, 0},  // 15
    {"minblock", 1, nullptr, 0},  // 16
    {nullptr, 0, nullptr, 0}
  };
  while (1) {
    int longindex = -1;
    int c = getopt_long_only(argc, argv, "", long_options, &longindex);
    if (c == -1) {
      break;  // Done parsing flags.
    } else if (c == '?') {  // If the user entered junk, let him know.
      cerr << "Invalid arguments." << endl;
      usage(argv[0]);
    } else {
      // Branch on long options.
      switch (longindex) {
        case 0: min_len = atol(optarg); break;
        case 1: type = MAM; break;
        case 2: type = MEM; break;
        case 3: type = MUM; break;
        case 4: type = MAM; break;
        case 5: nucleotides_only = true; break;
        case 6: max_n_threads = atoi(optarg) ; break;
        case 7: sam_out = true; break;
        case 8:  // This is messy - a design flaw
          ::verbose = true;
          ref_args.verbose = true;
          SAArgs::verbose = true;
          PairsArgs::verbose = true;
          ReadersArgs::verbose = true;
          verbose = true;
          break;
        case 9: nomap = true; break;
        case 10: ref_args.rcref = true; break;
        case 11: fastq = true; break;
        case 12: sam_in = true; break;
        case 13: mappability = true; break;
        case 14: read_ahead = false; break;
        case 15: memory_mapped = false; break;
        case 16: min_block = atoi(optarg); break;
        default: break;
      }
    }
  }

  // Validate arguments
  argc -= optind;
  if (argc < 2) {
    cerr << "There are too few arguments" << endl;
    usage(argv[0]);
  }
  if (fastq && sam_in) throw Error("-fastq cannot be used with -samin");
  if (nomap && !sam_out) throw Error("-nomap can only be used with -sam_out");
  if (mappability && !ref_args.rcref)
    throw Error("-mappability requires -rcref");
  char * * args = argv + optind;
  ref_args.ref_fasta = *args;
  n_input = argc - 1;
  input = args + 1;
  check_integer_sizes();
  n_threads = max_n_threads;  // do something better later
}

// Make sure the compiled size is good for reference length
void Args::check_integer_sizes() const {
  if (sizeof(ANINT) < 8 || sizeof(SINT) < 8) {
    string program(argv[0]);
    const uint64_t fasta_size = file_size(ref_args.ref_fasta);
    const uint64_t reference_size = fasta_size * (ref_args.rcref ? 2 : 1);
    if (reference_size <= INT_MAX) return;
    remove(program, "-long");
    remove(program, "-medium");
    if (reference_size > UINT_MAX - 100000) {
      program += "-long";
    } else if (reference_size > INT_MAX - 100000 && sizeof(SINT) != 8) {
      program += "-medium";
    }
    if (program != argv[0]) {
      if (verbose)
        cerr << "# switching program to " << program
             << " because reference is too long for " << argv[0] << endl
             << "# you can run " << program
             << " directly instead if desired" << endl;
      // See if program name is absolute - and if not fix that - todo
      execve(program.c_str(), argv, envp);
      throw Error("Problem running program") << program;
    }
  } else if (sizeof(ANINT) < 4) {
    cerr << "# int size other than 4 and 8 are not planned for "
         << "- use program at own risk." << endl;
  }
}

// Display proper command line usage
void Args::usage(const string & prog) const {
  cerr << "Usage: " << prog <<
      " [options] <reference-file> <query-file> ...\n"
      "Implemented MUMmer v3 options:\n"
      "-mum         "
      "  compute maximal matches that are unique in both sequences\n"
      "-mumreference"
      "  compute maximal matches that are unique in the reference-\n"
      "             "
      "  sequence but not necessarily in the query-sequence (default)\n"
      "-mumcand       same as -mumreference\n"
      "-maxmatch    "
      "  compute all maximal matches regardless of their uniqueness\n"
      "-l             set the minimum length of a match\n"
      "               if not set, the default value is 20\n"
      "-n             match only the characters a, c, g, or t\n"
      "\n"
      "Additional options:\n"
      "-verbose       output diagnostics and progress to stderr\n"
      "-samin         input in SAM format\n"
      "-samout        output in basic SAM format\n"
      "-qthreads      number of threads to use for queries\n"
      "-nomap         output unmapped reads too (only when -samout)\n"
      "-rcref         reverse complement reference\n"
      "-fastq         fastq input\n"
      "-mappability   output mappability measures only\n"
      "-minblock      with -samout, after merge of mapped segments\n"
      "               by read-start position, ensures that a mapped block\n"
      "               is of a minimum length\n"
      "-cached        shorter real time for subsequent runs only\n"
      "-normalmem    turn off memory mapping" << endl;
  exit(1);
}


