#ifndef NULL_SIZE_T_HPP
#define NULL_SIZE_T_HPP

#include <cstddef>
#include <iostream>
#include <limits>
#include <string>

// pour éviter d'utiliser des pointeurs (avec les problèmes de copies)
// on utilise plutôt un indice (size_t), et pour traiter le cas "pas d'indice"
// on défini le null_size_t (doit être interprété comme l'équivalent d'un nullptr)
const size_t null_size_t = std::numeric_limits<size_t>::max();

template <char nullChar> void put_Ptr_size_t(std::ostream &os, const size_t value) {
  if (value == null_size_t)
    os << nullChar;
  else
    os << value;
}

template <char nullChar> size_t get_Ptr_size_t(std::istream &is) {
  std::string str;
  is >> str;
  if (str[0] == nullChar)
    return null_size_t;
  else
    return std::stoul(str.c_str());
}

#endif /* end of include guard: NULL_SIZE_T_HPP */
