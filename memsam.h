/* Copyright Peter Andrews 2013 CSHL */

#ifndef LONGMEM_MEMSAM_H_
#define LONGMEM_MEMSAM_H_

#include <unistd.h>

#include <map>
#include <string>
#include "./util.h"
#include "./error.h"

enum MapFlag {
  is_paired = 1 << 0,
  is_proper = 1 << 1,
  is_unmapped = 1 << 2,
  is_mate_unmapped = 1 << 3,
  is_reversed = 1 << 4,
  is_mate_reversed = 1 << 5,
  is_first = 1 << 6,
  is_second = 1 << 7,
  is_not_primary = 1 << 8,
  is_bad_vendor_quality = 1 << 9,
  is_a_duplicate = 1 << 10
};

enum field_index {
  name_field = 0,
  flag_field,
  chromosome_field,
  position_field,
  mapq_field,
  cigar_field,
  mate_chromosome_field,
  mate_position_field,
  template_length_field,
  bases_field,
  errors_field,
  extra_field
};

template <class CHAR>  // template just to be in header
const CHAR * next_field(const CHAR * field) {
  while (*field++ != '\t') {}
  return field;
}
template <class CHAR>
const CHAR * next_read(const CHAR * field) {
  while (*field++ != '\n') {}
  return field;
}

class MemSam {
 public:
  static std::map <std::string, uint64_t> chromosomes;
  explicit MemSam(const char * in)
      : data{in}, absolute_position_{read_absolute_position()} { }
  const char * jump_to(int field) const {
    auto byte = data;
    while (field--) {
      byte = next_field(byte);
    }
    return byte;
  }
  const std::string to_string(int field) const {
    auto string_start = jump_to(field);
    auto string_end = next_field(string_start) - 1;
    return std::string(string_start, string_end);
  }
  uint64_t position() const {
    return atol(jump_to(position_field));
  }
  uint64_t mate_position() const {
    return atol(jump_to(mate_position_field));
  }
  uint64_t absolute_position() const {
    return absolute_position_;
  }
  uint64_t read_absolute_position() const {
    auto chromosome_start = jump_to(chromosome_field);
    auto position_start = next_field(chromosome_start);
    auto chr_ = std::string(chromosome_start, position_start - 1);
    // std::cerr << "returning " << chr_ << std::endl;
    return atol(position_start) +
        chromosomes.at(chr_);
  }
  uint64_t mate_absolute_position() const {
    auto chromosome_start = jump_to(mate_chromosome_field);
    auto position_start = next_field(chromosome_start);
    return atol(position_start) +
        chromosomes.at(std::string(chromosome_start, position_start - 1));
  }
  const std::string chromosome() const {
    return to_string(chromosome_field);
  }
  const std::string mate_chromosome() const {
    return to_string(mate_chromosome_field);
  }
  const std::string name() const {
    return to_string(name_field);
  }
  const std::string bases() const {
    return to_string(bases_field);
  }
  const std::string cigar() const {
    return to_string(cigar_field);
  }
  const std::string rest() const {
    auto start = jump_to(errors_field);
    return std::string(start, next_read(start) - 1);
  }
  uint16_t flag() const {
    return atoi(jump_to(flag_field));
  }
  uint16_t mapq() const {
    return atoi(jump_to(mapq_field));
  }
  bool is_secondary() const {
    return rest().find("HI:i:0") == std::string::npos;
  }
  bool is_better_than(const MemSam & other) const {
    auto this_flag = flag();
    auto other_flag = other.flag();
    if ((this_flag & is_mate_unmapped) > (other_flag & is_mate_unmapped))
      return false;
    if ((this_flag & is_mate_unmapped) == (other_flag & is_mate_unmapped)) {
      if (this_flag & is_mate_unmapped) {
        // Compare errors - best is smallest
      } else {
        // Compare matches - best is ?
      }
    }
    return true;
  }

  bool operator<(const MemSam & right) const {
    auto lpos = absolute_position();
    auto rpos = right.absolute_position();
    if (lpos == rpos) {
      auto lname = name();
      auto rname = right.name();
      if (lname == rname) {
        const unsigned int mate_info = flag() &
            (is_first | is_second | is_reversed);
        const unsigned int r_mate_info = right.flag() &
            (is_first | is_second | is_reversed);
        if (mate_info == r_mate_info) {
          // return false;
          throw paa::Error("flags equal") << "\n" << (*this) << "\n" << right;
        }
        return mate_info < r_mate_info;
      } else {
        return lname < rname;
      }
    } else {
      return lpos < rpos;
    }
  }

  template <class Out>
  void out(Out & o) const {
    o << name() << '\t' << flag() << '\t'
      << chromosome() << '\t' << position() << '\t'
      << mapq() << '\t' << cigar() << '\t'
      << mate_chromosome() << '\t' << mate_position() << "\t0\t"
      << bases() << '\t' << rest();
  }

  template <class Out>
  void dupe_out(Out & o) const {
    o << name() << '\t' << (flag()|is_a_duplicate) << '\t'
      << chromosome() << '\t' << position() << '\t'
      << mapq() << '\t' << cigar() << '\t'
      << mate_chromosome() << '\t' << mate_position() << "\t0\t"
      << bases() << '\t' << rest();
  }

  const char * data;

 private:
  uint64_t absolute_position_;
};

template <class Out>
Out & operator<<(Out & out, const MemSam & sam) {
  sam.out(out);
  return out;
}

#endif  // LONGMEM_MEMSAM_H_
