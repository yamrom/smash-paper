/* Copyright Peter Andrews 2013 CSHL */

#include "./util.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
using std::cerr;
using std::endl;

#include <sstream>
using std::stringstream;

#include <string>
using std::string;

#include "./error.h"
using paa::Error;

bool read_ahead = true;
bool memory_mapped = true;

void remove(string & input, const string & search) {
  const size_t pos = input.find(search);
  if (pos != string::npos) input.erase(pos, search.size());
}
void replace(string & input, const char a, const char b) {
  size_t pos = 0;
  while ((pos = input.find(a, pos)) != string::npos) {
    input[pos] = b;
    ++pos;
  }
}
void replace_substring(string & str,
                       const string & oldstr, const string & newstr) {
  const size_t pos = str.find(oldstr);
  if (pos != string::npos) {
    str.replace(pos, oldstr.size(), newstr);
    // replace(str, oldstr, newstr);
  }
}
void replace_substring(std::string & str,
                       const char * const oldstr, const char * const newstr) {
  return replace_substring(str, string(oldstr), string(newstr));
}


void mkdir(const std::string & dir_name) {
  ::mkdir(dir_name.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
}

uint64_t file_size(const string & file_name) {
    struct stat st;
    stat(file_name.c_str(), &st);
    return static_cast<uint64_t>(st.st_size);
}

bool readable(const string & file) {
  return !access(file.c_str(), R_OK);
}

warn::warn(const string & message) {
  cerr << message << endl;
}

void bwritec(FILE * output, const void * data, const std::string & name,
             const uint64_t count) {
  const uint64_t written = fwrite(data, 1, count, output);
  if (written != count) {
    perror(nullptr);
    throw Error("problem writing") << count << "bytes at" << name
                                   << "only wrote" << written;
  }
}

void bwritec(const std::string & filename, const void * data,
             const std::string & name,
             const uint64_t count) {
  FILE * output = fopen(filename.c_str(), "wb");
  if (output == nullptr)
    throw Error("could not open output") << filename << "for writing";
  bwritec(output, data, name, count);
  if (fclose(output) != 0)
    throw Error("problem closing output file") << filename;
}

void breadc(FILE * input, void * data, const std::string & name,
            const uint64_t count) {
  const uint64_t read = fread(data, 1, count, input);
  if (read != count)
    throw Error("problem reading") << count << "elements at" << name;
}

void breadc(const std::string & filename, void * & data,
            const std::string & name,
            const uint64_t count) {
  if (memory_mapped) {
    int input = open(filename.c_str(), O_RDONLY);
    if (input == -1)
      throw Error("could not open input") << filename << "for reading";
    if (count) {
      // cerr << "read ahead " << read_ahead << " for " << name << endl;
      if ((data = mmap64(nullptr, count, PROT_READ, MAP_SHARED |
                         (read_ahead ? MAP_POPULATE : 0),
                         input, 0)) == MAP_FAILED) {
        throw Error("Memory mapping error for") << name;
      }
    }
    if (close(input) == -1)
      throw Error("problem closing input file") << filename;
  } else {
    FILE * input = fopen(filename.c_str(), "rb");
    if ((data = malloc(count)) == nullptr)
      throw Error("malloc error for") << name;
    breadc(input, data, name, count);
    if (fclose(input) != 0)
      throw Error("problem closing input file") << filename;
  }
}

MappedFile::MappedFile() {}

void MappedFile::load(const std::string & file_name_) {
  file_name = file_name_;
  const size_t input_size = file_size(file_name);
  page = static_cast<size_t>(getpagesize());
  // cerr << "Page size is " << page << endl;
  if (!input_size) {
    throw Error("Zero input file size for") << file_name;
  }
  int input = open(file_name.c_str(), O_RDONLY);
  if (input == -1) {
    throw Error("could not open input") << file_name << "for reading";
  }
  file_begin = static_cast<char * const>(
      mmap64(nullptr, input_size, PROT_READ, MAP_PRIVATE, input, 0));
  if (file_begin == MAP_FAILED) {
    throw Error("Memory mapping error for") << file_name;
  }
  close(input);
  file_end = file_begin + input_size;
  sequential();
  // needed();
}

MappedFile::MappedFile(const std::string & file_name_) {
  load(file_name_);
}
void MappedFile::unmap() {
  if (munmap(file_begin, size())) {
    perror("munmap error");
    throw Error("unmap error for") << file_name;
  }
}
int MappedFile::advise(const char * start, size_t length,
                       const int advice) const {
  if (start == nullptr && length == 0) {
    start = begin();
    length = size();
  } else {
    start = (const char *)(((uint64_t)start / page_size()) * page_size());
  }
  return madvise(const_cast<char *>(start), length, advice);
}
void MappedFile::sequential(const char * start, const size_t length) const {
  if (advise(start, length, MADV_SEQUENTIAL)) {
    throw Error("sequential madvise error");
  }
}
void MappedFile::random(const char * start, const size_t length) const {
  if (advise(start, length, MADV_RANDOM)) {
    throw Error("random madvise error");
  }
}
void MappedFile::unneeded(const char * start, const size_t length) const {
  if (advise(start, length, MADV_DONTNEED)) {
    throw Error("unneeded madvise error");
  }
}
void MappedFile::needed(const char * start, const size_t length) const {
  if (advise(start, length, MADV_WILLNEED)) {
    throw Error("needed madvise error");
  }
}

