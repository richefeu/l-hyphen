#ifndef NULL_SIZE_T_HPP
#define NULL_SIZE_T_HPP

#include <cstddef>
#include <limits>

// pour éviter d'utiliser des pointeurs (avec les problèmes de copies)
// on utilise plutôt un indice (size_t), et pour traiter le cas "pas d'indice"
// on défini le null_size_t (doit être interprété comme l'équivalent d'un nullptr)
const size_t null_size_t = std::numeric_limits<size_t>::max();

#endif /* end of include guard: NULL_SIZE_T_HPP */
