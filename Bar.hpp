#ifndef BAR_HPP
#define BAR_HPP

#include "vec2.hpp"

/**
 *  Une barre extensible entre deux noeuds
 *
 */
class Bar {
public:
  size_t i, j; ///< le numéro des noeuds dans Cell::nodes
  double l0;   ///< longueur initiale de la barre [m]
  double kn;   ///< raideur suivant l'axe [N/m]
  double fn;   ///< force axiale [N]

  Bar(size_t t_i, size_t t_j);

  void init(double t_kn, vec2r &posi, vec2r &posj);
};

#endif /* end of include guard: BAR_HPP */
