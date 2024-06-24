#ifndef NEIGHBOR_HPP
#define NEIGHBOR_HPP

#include "vec2.hpp"

class Neighbor;

struct Connect {
  size_t ic; // indice de la première cellule dans Sample::cells
  size_t jc; // indice de la seconde cellule dans Sample::cells (ic <= jc)
  size_t in; // indice du noeud de la première cellule dans Sample::cells[ic]
  size_t jn; // indice du noeud de la seconde cellule dans Sample::cells[jc]
  Neighbor *woami{nullptr};
  Connect(size_t t_ic, size_t t_jc, size_t t_in, size_t t_jn, Neighbor *t_woami) : ic(t_ic), jc(t_jc), in(t_in), jn(t_jn) {
    woami = t_woami;
  }
};

/**
 * Contact potentiel (voisin ou interaction)
 *
 */
class Neighbor {
public:
  size_t ic; // indice de la première cellule dans Sample::cells
  size_t jc; // indice de la seconde cellule dans Sample::cells (ic <= jc)

  size_t in; // indice du noeud de la première cellule dans Sample::cells[ic]
  size_t jn; // indice du noeud de la seconde cellule dans Sample::cells[jc]
  // jn code pour la barre qui commence par jn

  vec2r n; // vecteur normal (de j vers i)

  Neighbor *brother{nullptr};

  // contact
  int contactState;
  double fn; // force normale de contact
  double ft; // force de frottement

  // cohesion
  int glueState; // 0 = pas de colle, 1 = colle (avec 3eme corps)
  double fn_coh; // force normale de cohésion solide
  double ft_coh; // force tangentielle de cohésion solide
  double kn_coh;
  double kt_coh;
  double fn_coh_max;
  double ft_coh_max;
  double yieldPower;

  Neighbor(size_t t_ic, size_t t_jc, size_t t_in, size_t t_jn);
};

namespace std {
// par defaut, la comparaison std::less est utilisée pour assurer l'ordre dans un std::set
template <> struct less<Neighbor> {
  bool operator()(const Neighbor &lhs, const Neighbor &rhs) const {
    if (lhs.jc < rhs.jc)
      return true;
    if (lhs.jc == rhs.jc && lhs.in < rhs.in)
      return true;
    if (lhs.jc == rhs.jc && lhs.in == rhs.in && lhs.jn < rhs.jn)
      return true;
    return false;
  }
};
} // namespace std

#endif /* end of include guard: NEIGHBOR_HPP */
