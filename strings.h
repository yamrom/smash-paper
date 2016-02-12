/* Copyright Peter Andrews 2015 CSHL */

#ifndef LONGMEM_STRINGS_H_
#define LONGMEM_STRINGS_H_

#include <string>

namespace paa {

std::string to_lower(const std::string & in);
std::string to_upper(const std::string & in);
void remove(std::string & input, const std::string & search);
void replace(std::string & input, const char a, const char b);
void replace_substring(std::string & str,
                       const std::string & oldstr, const std::string & newstr);
void replace_substring(std::string & str,
                       const char * const oldstr, const char * const newstr);

}  // namespace paa

#endif  // LONGMEM_STRINGS_H_
