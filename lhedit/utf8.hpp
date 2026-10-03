// utf8.hpp -- the little UTF-8 awareness the editor needs.
//
// The l-hyphen input files are mostly ASCII, but their comments are often
// written in French ("# contact adhésif"). The text stays stored as bytes and
// the cursor as a byte offset; these helpers keep the cursor on character
// boundaries and turn byte offsets into screen columns (one per code point,
// which is right for accented Latin letters).

#ifndef UTF8_HPP
#define UTF8_HPP

#include <string>

namespace utf8 {

// A continuation byte (10xxxxxx) never starts a character.
inline bool isContinuation(char c) { return (static_cast<unsigned char>(c) & 0xC0) == 0x80; }

// Byte length of the character starting at byte b (1 for a stray byte).
inline size_t length(const std::string& s, size_t b) {
  size_t e = b + 1;
  while (e < s.size() && isContinuation(s[e])) e++;
  return e - b;
}

// Screen column of byte offset b.
inline int column(const std::string& s, size_t b) {
  if (b > s.size()) b = s.size();
  int col = 0;
  for (size_t i = 0; i < b; ++i) {
    if (!isContinuation(s[i])) col++;
  }
  return col;
}

inline int columns(const std::string& s) { return column(s, s.size()); }

// Byte offset of screen column col, clamped to the end of the string.
inline size_t byteAt(const std::string& s, int col) {
  size_t b = 0;
  while (b < s.size() && col > 0) {
    b += length(s, b);
    col--;
  }
  return b;
}

// Moves a byte offset back onto the start of the character it falls in.
inline size_t snap(const std::string& s, size_t b) {
  if (b > s.size()) return s.size();
  while (b > 0 && b < s.size() && isContinuation(s[b])) b--;
  return b;
}

// The first n columns of s, never cutting a character in half.
inline std::string clip(const std::string& s, int n) { return s.substr(0, byteAt(s, n)); }

}  // namespace utf8

#endif /* end of include guard: UTF8_HPP */
