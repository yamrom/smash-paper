/* Copyright Peter Andrews 2015 CSHL */

#ifndef LONGMEM_FILES_H_
#define LONGMEM_FILES_H_

#include <cstdint>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "./error.h"

namespace paa {

extern bool read_ahead;
extern bool memory_mapped;

// Files and directories
uint64_t file_size(const std::string & file_name);
void mkdir(const std::string & dir_name);
bool readable(const std::string & file);
void unlink(const std::string & file);

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
            const std::string & name, const uint64_t count);

template<class T>
void bread(const std::string & filename, T * & data,
           const std::string & name,
           const uint64_t count = 1) {
  breadc(filename, (void*&)data, name, count * sizeof(T));
}

class MappedFile {
 public:
  // deleted
  MappedFile(const MappedFile &) = delete;
  MappedFile operator=(const MappedFile &) = delete;
  MappedFile & operator=(MappedFile &&) = delete;

  MappedFile();
  explicit MappedFile(const std::string & file_name_,
                      const bool warn_empty = true);
  MappedFile(MappedFile && other);
  ~MappedFile();
  void load(const std::string & file_name_,
            const bool warn_empty = true);
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
  char * begin() {
    return file_begin;
  }
  char * end() {
    return file_end;
  }
  const char * begin() const {
    return file_begin;
  }
  const char * end() const {
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

//
// a potentially dangerous class - exists in a build state or a read state
// and some member functions can only be used in one of the states
//
template<class Type>
class MappedVector {
 public:
  typedef Type * iterator;
  typedef const Type * const_iterator;

  // deleted
  MappedVector(const MappedVector &) = delete;
  MappedVector & operator=(const MappedVector &) = delete;
  MappedVector & operator=(MappedVector &&) = delete;

  // use these member functions during building only
  MappedVector() : MappedVector{initial_size} {}
  explicit MappedVector(const uint64_t start_size) : capacity{start_size},
    data{static_cast<Type *>(::operator new(sizeof(Type) * capacity))} { }
  MappedVector(MappedVector && old) :
      file{std::move(old.file)}, mapped{old.mapped},
    n_elem{old.n_elem}, capacity{old.capacity}, data{old.data} {
      old.mapped = false;
      old.capacity = 0;
    }
  template<typename... VArgs>
  void emplace_back(VArgs && ... args) {
    expand_if_needed();
    new (data + n_elem++) Type(std::forward<VArgs>(args) ...);
  }
  void push_back(const Type val) {
    expand_if_needed();
    data[n_elem++] = val;
  }
  template <class P>
  void insert_back(P b, P e) {
    while (b != e) {
      push_back(*(b++));
    }
  }
  Type * begin() { return data; }
  Type * end() { return data + n_elem; }
  void clear() { n_elem = 0; }
  void reduce_size(const uint64_t new_size) { n_elem = new_size; }
  void save(const std::string & file_name) const {
    bwritec(file_name, data, "MappedVector", n_elem * sizeof(Type));
  }
  void write(std::FILE * out_file) const {
    bwritec(out_file, data, "MappedVector", n_elem * sizeof(Type));
  }
  void write_n(std::FILE * out_file, const uint64_t n) const {
    bwritec(out_file, data, "MappedVector", n * sizeof(Type));
  }

  // use these member function during reading only
  explicit MappedVector(const std::string & file_name,
                        const bool warn_empty = true) :
      file{file_name, warn_empty},
    mapped{file.begin() == nullptr ? false : true},
    n_elem{file.size() / sizeof(Type)},
    capacity{n_elem},
    data{reinterpret_cast<Type *>(file.begin())} { }
  const std::string & name() const { return file.name(); }
  void load(const std::string & file_name_, const bool warn_empty = true) {
    file.load(file_name_, warn_empty);
    mapped = file.begin() == nullptr ? false : true;
    n_elem = file.size() / sizeof(Type);
    data = reinterpret_cast<Type *>(file.begin());
  }

  // these member functions can be used anytime
  uint64_t bytes() const { return capacity * sizeof(Type); }
  uint64_t size() const { return n_elem; }
  const Type * begin() const { return data; }
  const Type * end() const { return data + n_elem; }
  Type operator[](const uint64_t index) const { return data[index]; }
  Type & operator[](const uint64_t index) { return data[index]; }
  const Type back() const { return data[n_elem - 1]; }
  Type & back() { return data[n_elem - 1]; }
  ~MappedVector() { if (!mapped && capacity) delete data; }

 private:
  void expand_if_needed() {
    if (n_elem + 1 > capacity) {
      capacity *= 2;
      Type * new_data_location =
          static_cast<Type *>(::operator new(sizeof(Type) * capacity));
      memcpy(new_data_location, data, n_elem * sizeof(Type));
      delete data;
      data = new_data_location;
    }
  }
  static constexpr uint64_t initial_size{1000};
  MappedFile file{};
  bool mapped{false};
  uint64_t n_elem{0};
  uint64_t capacity{0};
  Type * data{nullptr};
};

template<class BigInt, class SmallInt>
class CompressedInts {
 public:
  explicit CompressedInts(const uint64_t capacity = 0,
                          const uint64_t n_lookup = 0) :
      lookup{n_lookup}, small{capacity} { }
  explicit CompressedInts(const std::string & dir,
                          const unsigned int small_start = 0,
                          const unsigned int lookup_start = 0) :
      lookup{dir + "/lookup.bin"},
    small{dir + "/counts.bin"},
    big{dir + "/over.bin"},
    small_position{small_start},
    big_position{lookup[lookup_start]} {}
  BigInt next_int() const {
    if (small_position == small.size()) {
      throw paa::Error("Tried to read past end of CompressedInts");
      return std::numeric_limits<BigInt>::max();
    } else {
      const SmallInt s = small[small_position++];
      if (s == std::numeric_limits<SmallInt>::max()) {
        return big[big_position++];
      } else {
        return s;
      }
    }
  }
  void add_int(unsigned int b) {
    if (b > std::numeric_limits<BigInt>::max()) {
      std::cerr << "add_int encountered big int " << b << std::endl;
      b = std::numeric_limits<BigInt>::max();
    }
    if (b >= std::numeric_limits<SmallInt>::max()) {
      small.push_back(std::numeric_limits<SmallInt>::max());
      big.push_back(static_cast<BigInt>(b));
    } else {
      small.push_back(static_cast<SmallInt>(b));
    }
  }
  void print_savings() const {
    std::cerr << "Saved " << small.size() * (sizeof(BigInt)-sizeof(SmallInt)) -
        big.size() * sizeof(BigInt) << " bytes" << std::endl;
  }
  void add_lookup_entry() {
    lookup.push_back(big.size());
  }

  void save(const std::string & dir_name) const {
    mkdir(dir_name);
    std::ostringstream counts_filename;
    counts_filename << dir_name << "/counts.bin";
    std::FILE * counts_file = fopen(counts_filename.str().c_str(), "wb");
    if (counts_file == nullptr) throw paa::Error("Problem opening counts file")
                                    << counts_filename.str();
    if (fwrite(small.begin(), sizeof(SmallInt), small.size(), counts_file) !=
        small.size())
      throw paa::Error("Problem writing in counts file")
          << counts_filename.str();
    if (fclose(counts_file)) throw paa::Error("Problem closing counts file")
                                 << counts_filename.str();
    std::ostringstream over_filename;
    over_filename << dir_name << "/over.bin";
    std::FILE * over_file = fopen(over_filename.str().c_str(), "wb");
    if (over_file == nullptr) throw paa::Error("Problem opening over file")
                                  << over_filename.str();
    if (fwrite(big.begin(), sizeof(BigInt), big.size(), over_file) !=
        big.size())
      throw paa::Error("Problem writing in over file") << over_filename.str();
    if (fclose(over_file)) throw paa::Error("Problem closing over file")
                               << over_filename.str();
    std::ostringstream lookup_filename;
    lookup_filename << dir_name << "/lookup.bin";
    std::FILE * lookup_file = fopen(lookup_filename.str().c_str(), "wb");
    if (lookup_file == nullptr) throw paa::Error("Problem opening lookup file")
                                    << lookup_filename.str();
    if (fwrite(lookup.begin(), sizeof(unsigned int),
               lookup.size(), lookup_file) != lookup.size())
      throw paa::Error("Problem writing in lookup file")
          << lookup_filename.str();
    if (fclose(lookup_file)) throw paa::Error("Problem closing lookup file")
                                 << lookup_filename.str();
    print_savings();
  }
  void relocate(const uint64_t new_small, const uint64_t new_lookup) const {
    small_position = new_small;
    big_position = lookup[new_lookup];
  }
  uint64_t size() const {
    return small.size();
  }

 private:
  MappedVector<unsigned int> lookup;
  MappedVector<SmallInt> small;
  MappedVector<BigInt> big;
  mutable uint64_t small_position{0};
  mutable uint64_t big_position{0};
};

typedef CompressedInts<uint16_t, uint8_t> Compressed;

}  // namespace paa

#endif  // LONGMEM_FILES_H_
