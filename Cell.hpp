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

#include "Bar.hpp"
#include "Neighbor.hpp"
#include "Node.hpp"
#include "null_size_t.hpp"

#include <set>
#include <vector>

/// Une cellule contenant une boucle (fermée ou non) de barres extensibles
/// 
class Cell {
public:
  std::vector<Node> nodes;               ///< les noeuds
  std::vector<Bar> bars;                 ///< les barres
  std::set<Neighbor> neighbors;          ///< les voisins
  std::vector<Neighbor *> vec_neighbors; ///< vecteur de pointeurs de voisins pour un accès rapide

  double radius{0.0};   ///< un seul rayon pour toute la cellule
  double surface{0.0};  ///< surface (volume) intérieure
  double surface0{0.0}; ///< surface (volume) intérieure initial
  bool close{false};    ///< flag pour dire si la cellule est fermée ou pas
  vec2r center;         ///< "centre" de la cellule
  double p_int{0.0};    ///< pression interieure (éventuelle) pour cellule fermée

  Cell();

  void reorderNodes(); // réordonne la position des noeuds de la cellule en fonction de leur angle
  void insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR); // ajoute ou supprime un voisin

  void
  connectOrderedNodes(double width, double t_kn, double t_kr, double t_mz_max, double t_p_int,
                      bool closed = true);     // ajoute des barres entre les noeuds et défini leurs longueurs initiales
  void CellSurface();                          // mesure la surface exposée à la pression
  void CellCenter();                           // mesure le centre de la cellule
  void CellForce(vec2r &force);                // mesure la force axiale
  double getElasticNRJ(double compressFactor); // mesure la force NRJ
};
