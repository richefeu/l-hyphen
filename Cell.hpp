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

  /**
   * @brief ré-ordonne la position des noeuds de la cellule en fonction de leur angle
   *
   */
  void reorderNodes();

  /**
   * @brief ajoute ou bien supprime un voisin. isNear doit être fourni et la methode déterminera si il faut ajouter ou
   *        bien supprimer (c'est basé sur le fait que la 'liste' de voisin est un std::set, et donc ordonnée et le
   *        'find' est rapide)
   *
   * @param ci cellule i
   * @param cj cellule j
   * @param in noeud n dans la cellule i
   * @param jn noeud n dans la cellule j
   * @param isNEAR à tester avant l'appel
   */
  void insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR);

  /**
   * @brief ajoute des barres entre les noeuds et défini leurs longueurs initiales
   *        ainsi que les angles initaux entre les barres adjacentes.
   *        on suppose que les noeuds sont numérotés suivant un sens de rotation quelconque
   *        mais qu'ils se suivent.
   *
   * @param width épaisseur des barres
   * @param kn_ raideur axiale des barres
   * @param kr_ raideur angulaire entre les barres adjacentes
   * @param mz_max_ moment seuil plastique
   * @param closed spécifie si les barres forment une boucle (true par defaut)
   */
  void connectOrderedNodes(double width, double kn_, double kr_, double mz_max_, bool closed = true);
};


#endif /* end of include guard: CELL_HPP */
