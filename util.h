/* Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_UTIL_H_
#define LONGMEM_UTIL_H_

#include <stdint.h>

#include <exception>
#include <sstream>
#include <string>
#include <fstream>

extern bool read_ahead;
extern bool memory_mapped;

template <class Val>
Val sqr(const Val val) {
  return val * val;
}

// Remove a portion of string, if found
void remove(std::string & input, const std::string & search);
void replace(std::string & input, const char a, const char b);
void replace_substring(std::string & str,
                       const std::string & oldstr, const std::string & newstr);
void replace_substring(std::string & str,
                       const char * const oldstr, const char * const newstr);


// Files and directories
void mkdir(const std::string & dir_name);
uint64_t file_size(const std::string & file_name);
bool readable(const std::string & file);

// Binary data read and write
void bwritec(FILE * output, const void * data, const std::string & name,
             const uint64_t count);

template<class T>
void bwrite(FILE * output, const T & data, const std::string & name,
            const uint64_t count = 1) {
  bwritec(output, &data, name, sizeof(T) * count);
}

void bwritec(const std::string & filename, const void * data,
             const std::string & name,
             const uint64_t count);

template<class T>
void bwrite(const std::string & filename, const T & data,
            const std::string & name,
            const uint64_t count = 1) {
  bwritec(filename, &data, name, sizeof(T) * count);
}

void breadc(FILE * input, void * data, const std::string & name,
            const uint64_t count);

template<class T>
void bread(FILE * input, T & data, const std::string & name,
          const uint64_t count = 1) {
  breadc(input, &data, name, sizeof(T) * count);
}

void breadc(const std::string & filename, void * & data,
            const std::string & name,
            const uint64_t count);

template<class T>
void bread(const std::string & filename, T * & data,
           const std::string & name,
           const uint64_t count = 1) {
  breadc(filename, (void*&)data, name, count * sizeof(T));
}

class warn {
 public:
  explicit warn(const std::string & message);
};

class MappedFile {
 public:
  MappedFile();
  explicit MappedFile(const std::string & file_name_);
  // MappedFile(const MappedFile &) = delete;
  // MappedFile operator=(const MappedFile &) = delete;
  void load(const std::string & file_name_);
  const std::string & name() const {
    return file_name;
  }
  void unmap();
  int advise(const char * start, size_t length, const int advice) const;
  void sequential(const char * start = nullptr,
                  const size_t length = 0) const;
  void random(const char * start = nullptr,
              const size_t length = 0) const;
  void unneeded(const char * start = nullptr,
                const size_t length = 0) const;
  void needed(const char * start = nullptr,
              const size_t length = 0) const;
  char * begin() const {
    return file_begin;
  }
  char * end() const {
    return file_end;
  }
  const char * cbegin() const {
    return file_begin;
  }
  const char * cend() const {
    return file_end;
  }
  size_t size() const {
    return static_cast<size_t>(file_end - file_begin);
  }
  size_t page_size() const {
    return page;
  }

 private:
  std::string file_name;
  char * file_begin;
  char * file_end;
  size_t page;
};

enum class Dir : unsigned int {
  left = 0, right = 1
};

class Mappability {
 public:
  explicit Mappability(const std::string & file_name) : data(file_name) { }
  unsigned int operator()(const unsigned int n, const Dir dir) const {
    return *(reinterpret_cast<unsigned char *>(2 + data.begin() + n * 2UL +
                                               static_cast<unsigned int>(dir)));
  }
  unsigned int left(const unsigned int n) const {
    return *(reinterpret_cast<unsigned char *>(2 + data.begin() + n * 2UL));
  }
  unsigned int right(const unsigned int n) const {
    return *(reinterpret_cast<unsigned char *>(2 + data.begin() + n * 2UL + 1));
  }
  unsigned int size() const {
    return data.size() / 2;
  }
 private:
  MappedFile data;
};

class Reference {
 public:
  explicit Reference(const std::string & file_name) :
    reference(cache_name(file_name)) {}
  inline const char & operator[] (const unsigned int n) const {
    return *(reference.begin() + n);
  }
  void random() const {
    reference.random();
  }
  unsigned int size() const {
    return reference.size();
  }

 private:
  static std::string cache_name(const std::string & file_name) {
    // cerr << "Loading ref" << endl;
    std::ifstream input(file_name.c_str());
    std::string mod_name = file_name;
    replace(mod_name, '/', '+');
    std::string out_name = "/data/temp/paa/ref-cache-" + mod_name + ".txt";
    std::ifstream test(out_name);
    // Should lock file
    if (!test) {
      std::ofstream output(out_name);
      output << "X";
      std::string line;
      while (input >> line) {
        if (line[0] != '>')
          output << line;
      }
      // cerr << "Done loading and caching ref" << endl;
    }
    return out_name;
  }
  MappedFile reference;
};

#endif  // LONGMEM_UTIL_H_
