/* Copyright Peter Andrews 2015 CSHL */

#include "./strings.h"

#include <algorithm>
#include <string>

namespace paa {

std::string to_lower(const std::string & up) {
  std::string lower(up.size(), 'X');
  std::transform(up.begin(), up.end(), lower.begin(), ::tolower);
  return lower;
}
std::string to_upper(const std::string & low) {
  std::string upper(low.size(), 'X');
  std::transform(low.begin(), low.end(), upper.begin(), ::toupper);
  return upper;
}

void remove(std::string & input, const std::string & search) {
  const size_t pos = input.find(search);
  if (pos != std::string::npos) input.erase(pos, search.size());
}
void replace(std::string & input, const char a, const char b) {
  size_t pos = 0;
  while ((pos = input.find(a, pos)) != std::string::npos) {
    input[pos] = b;
    ++pos;
  }
}
void replace_substring(std::string & str,
                       const std::string & oldstr, const std::string & newstr) {
  const size_t pos = str.find(oldstr);
  if (pos != std::string::npos) {
    str.replace(pos, oldstr.size(), newstr);
    // replace(str, oldstr, newstr);
  }
}
void replace_substring(std::string & str,
                       const char * const oldstr, const char * const newstr) {
  return replace_substring(str, std::string(oldstr), std::string(newstr));
}

}  // namespace paa
