#ifndef NODE_HPP
#define NODE_HPP

#include "vec2.hpp"

/**
 * @brief un point ayant une masse (et d'autre propriétés)
 *
 */
class Node {
public:
  vec2r pos; // position
  vec2r vel; // vitesse
  vec2r acc; // acceleration

  vec2r force; // force résultante

  double mass; // masse

  size_t ictrl; // null_size_t si pas de control

  size_t prevNode; // numéro du noeud précédant
  size_t nextNode; // numéro du noeud suivant (null_size_t si pas de barre)

  double kr;     // raideur angulaire [N.m(/radian)]
  double mz;     // moment entre les barres connectées
  double mz_max; // moment-seuil de plastification parfaite

  // Ctor
  Node(double x, double y);

  void init(double t_kr, double t_mz_max, size_t prev, size_t next);
};

#endif /* end of include guard: NODE_HPP */
