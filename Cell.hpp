#ifndef CELL_HPP
#define CELL_HPP

#include "Node.hpp"
#include "Bar.hpp"
#include "Neighbor.hpp"

#include <vector>
#include <set>

/**
 * @brief une cellule contenant une boucle (fermée ou non) de barres extensibles
 *
 */
class Cell {
public:
  std::vector<Node> nodes;
  std::vector<Bar> bars;
  std::set<Neighbor> neighbors;

  double radius; // un seul rayon pour toute la cellule [m]

  Cell();
 
  void reorderNodes();
  void insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR);
  void connectOrderedNodes(double width, double kn_, double kr_, double mz_max_, bool closed = true);
};


#endif /* end of include guard: CELL_HPP */
