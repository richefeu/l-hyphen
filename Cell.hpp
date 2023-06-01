#ifndef CELL_HPP
#define CELL_HPP

#include "Bar.hpp"
#include "Neighbor.hpp"
#include "Node.hpp"
#include "null_size_t.hpp"

#include <set>
#include <vector>

/**
 * @brief une cellule contenant une boucle (fermée ou non) de barres extensibles
 *
 */
class Cell {
public:
  std::vector<Node> nodes;
  std::vector<Bar> bars;
  std::set<Neighbor> neighbors;

  double radius;   // un seul rayon pour toute la cellule
  double surface;  // surface (volume) intérieure
  double surface0; // surface (volume) intérieure initial
  bool close;      //  flag pour dire si la cellule est fermée ou pas
  vec2r center;    // "centre" de la cellule
  double p_int;    // pression interieure (éventuelle) pour cellule fermée

  Cell();

  void reorderNodes();
  void insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR);
  void connectOrderedNodes(double width, double kn_, double kr_, double mz_max_, double pint, bool closed = true);
  void CellSurface(); // mesure la surface exposée à la pression
  void CellCenter();  //
  void CellForce(vec2r & force);
	double getElasticNRJ(double compressFactor_);
};

#endif /* end of include guard: CELL_HPP */
