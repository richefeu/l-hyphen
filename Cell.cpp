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

#include "Cell.hpp"

Cell::Cell() : radius(0.0), surface(0.0), surface0(0.0), close(false) {}

/**
 *  Reorders the position of the nodes in the cell based on their angle.
 */
void Cell::reorderNodes() {
  if (nodes.empty()) {
    std::cout << "@Cell::reorderNodes, no node in the cell\n";
    return;
  }

  // Structure to sort nodes in clockwise order
  struct sortedNode {
    double x0, y0;
    double x, y;
    sortedNode(double X0, double Y0, double X, double Y) : x0(X0), y0(Y0), x(X), y(Y) {}
    bool operator<(const sortedNode &rhs) const {
      double angle = atan2(y - y0, x - x0);
      double rhsAngle = atan2(rhs.y - rhs.y0, rhs.x - rhs.x0);
      return angle < rhsAngle;
    }
  };

  // Calculate the centroid of the node positions
  double xc = 0.0;
  double yc = 0.0;
  for (size_t i = 0; i < nodes.size(); i++) {
    xc += nodes[i].pos.x;
    yc += nodes[i].pos.y;
  }
  xc /= static_cast<double>(nodes.size());
  yc /= static_cast<double>(nodes.size());

  // Create a multiset of sorted nodes
  std::multiset<sortedNode> sortedNodes;
  for (size_t i = 0; i < nodes.size(); i++) {
    sortedNodes.insert(sortedNode(xc, yc, nodes[i].pos.x, nodes[i].pos.y));
  }

  // Update the nodes with the sorted positions
  size_t currentId = 0;
  for (const auto &s : sortedNodes) {
    nodes[currentId].pos.x = s.x;
    nodes[currentId].pos.y = s.y;
    currentId++;
  }
}

/**
 *   Ajoute ou supprime un voisin. isNear doit être fourni et la methode déterminera si il faut ajouter ou
 *   bien supprimer (c'est basé sur le fait que la 'liste' de voisin est un std::set, et donc ordonnée et le
 *   'find' est rapide)
 *
 *   @param ci      cellule i
 *   @param cj      cellule j
 *   @param in      noeud n dans la cellule i
 *   @param jn      noeud n dans la cellule j
 *   @param isNEAR  à tester avant l'appel
 */
void Cell::insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR) {
  Neighbor toFind(ci, cj, in, jn);

  auto exist_it = neighbors.find(toFind);
  bool NEW = (exist_it == neighbors.end());

  if (isNEAR && NEW) {
    neighbors.insert(toFind);
  } else if (!isNEAR && !NEW) {
    neighbors.erase(exist_it);
  }
}

/**
 *  Ajoute des barres entre les noeuds et défini leurs longueurs initiales
 *  ainsi que les angles initaux entre les barres adjacentes.
 *  on suppose que les noeuds sont numérotés suivant un sens de rotation quelconque
 *  mais qu'ils se suivent.
 *
 *  @param width     épaisseur des barres
 *  @param t_kn      raideur axiale des barres
 *  @param t_kr      raideur angulaire entre les barres adjacentes
 *  @param t_mz_max  moment seuil plastique
 *  @param t_closed  spécifie si les barres forment une boucle (true par defaut)
 */
void Cell::connectOrderedNodes(double width, double t_kn, double t_kr, double t_mz_max, double t_p_int, bool t_closed) {
  if (nodes.empty()) {
    std::cout << "@Cell::connectOrderedNodes, Cannot connect the nodes because there's no node!\n";
    return;
  }

  close = t_closed;
  radius = 0.5 * width;
  bars.clear();

  for (size_t i = 1; i < nodes.size(); ++i) {
    bars.push_back(Bar(i - 1, i));
  }

  if (t_closed == true) {
    bars.push_back(Bar(nodes.size() - 1, 0));
    p_int = t_p_int;
  }

  // init bars
  for (size_t b = 0; b < bars.size(); ++b) {
    bars[b].init(t_kn, nodes[bars[b].i].pos, nodes[bars[b].j].pos);
  }

  // init nodes
  if (nodes.size() > 2) {
    if (t_closed == true) {
      nodes[0].init(t_kr, t_mz_max, nodes.size() - 1, 1);
    } else {
      nodes[0].init(t_kr, t_mz_max, null_size_t, 1);
    }

    for (size_t n = 1; n < nodes.size() - 1; n++) {
      nodes[n].init(t_kr, t_mz_max, n - 1, n + 1);
    }
    if (t_closed == true) {
      nodes.back().init(t_kr, t_mz_max, nodes.size() - 2, 0);
    } else {
      nodes.back().init(t_kr, t_mz_max, nodes.size() - 2, null_size_t);
    }
  }
}

/**
 *  Calcul la surface exposée à la pression par la méthode des produits vectoriels
 *
 */
void Cell::CellSurface() {
  surface = 0.0;
  if (close == false) {
    return;
  }

  for (size_t i = 1; i < nodes.size() - 1; ++i) {
    surface += 0.5 * (cross((nodes[i].pos - nodes[0].pos), (nodes[i + 1].pos - nodes[0].pos)));
  }
  surface = fabs(surface);
}

/**
 *  Calculates the center of a closed cell.
 *
 */
void Cell::CellCenter() {
  center.reset();
  for (size_t i = 0; i < nodes.size(); i++) {
    center += nodes[i].pos;
  }
  center /= double(nodes.size());
}

/**
 *  Calculates the total force acting on the Cell.
 *
 *  @param force  reference to the vec2r object to store the calculated force
 *
 */
void Cell::CellForce(vec2r &force) {
  force.reset();
  for (size_t i = 0; i < nodes.size(); i++) {
    force += nodes[i].force;
  }
}

/**
 *  Calculates the elastic energy of the cell.
 *
 *  @param compressFactor  the compression factor
 *
 *  @return the elastic energy of the cell
 */
double Cell::getElasticNRJ(double compressFactor) {
  if (close == false) {
    return 0.0;
  }
  double NRJ = 0.0;
  for (size_t b = 0; b < bars.size(); b++) {
    NRJ += 0.5 * bars[b].fn * bars[b].fn / bars[b].kn;
  }
  for (size_t n = 0; n < nodes.size(); n++) {
    NRJ += 0.5 * nodes[n].mz * nodes[n].mz / nodes[n].kr;
  }
  NRJ += 0.5 * p_int * p_int / compressFactor;
  return NRJ;
}
