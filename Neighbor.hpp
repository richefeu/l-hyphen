//  Copyright or © or Copr. l-hyphen
//
//  This software is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#pragma once

#include "vec2.hpp"

#define NOT_TOUCHING 0
#define TOUCHING 1

#define GLUE_NONE 0
#define GLUE_FORCE_THRESHOLD 1
#define GLUE_GC 2

/// Contact potentiel (voisin ou interaction)
///
class Neighbor {
public:
  size_t ic{0}; // indice de la première cellule dans Sample::cells
  size_t jc{0}; // indice de la seconde cellule dans Sample::cells (ic <= jc)

  size_t in{0}; // indice du noeud de la première cellule dans Sample::cells[ic]
  size_t jn{0}; // indice du noeud de la seconde cellule dans Sample::cells[jc]
                // jn code pour la barre qui commence par jn

  vec2r n; // vecteur normal (de j vers i)
  vec2r T;

  Neighbor *brother{nullptr};

  // contact
  int contactState{-1};
  double fn{0.0}; // force normale de contact
  double ft{0.0}; // force de frottement

  // cohesion
  int glueState{0};   // 0 = pas de colle, 1 = colle, 2, colle Gc
  double fn_coh{0.0}; // force normale de cohésion solide
  double ft_coh{0.0}; // force tangentielle de cohésion solide
  double kn_coh{0.0};
  double kt_coh{0.0};
  double fn_coh_max{0.0};
  double ft_coh_max{0.0};

  // rupture avec critère en energie
  double length{0.0}; // longueur d'interface collée
  double yieldPower{0.0};
  double Gc{0.0};

  Neighbor(size_t t_ic, size_t t_jc, size_t t_in, size_t t_jn);
};

namespace std {
// par defaut, la comparaison std::less est utilisée pour assurer l'ordre dans un std::set
template <> struct less<Neighbor> {
  bool operator()(const Neighbor &lhs, const Neighbor &rhs) const {
    if (lhs.jc < rhs.jc) {
      return true;
    }
    if (lhs.jc == rhs.jc && lhs.in < rhs.in) {
      return true;
    }
    if (lhs.jc == rhs.jc && lhs.in == rhs.in && lhs.jn < rhs.jn) {
      return true;
    }
    return false;
  }
};
} // namespace std
