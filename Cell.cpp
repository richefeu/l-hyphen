#include "Cell.hpp"

Cell::Cell() : radius(0.0) {}

/**
 * @brief ré-ordonne la position des noeuds de la cellule en fonction de leur angle
 *
 */
void Cell::reorderNodes() {
  if (nodes.empty()) {
    std::cout << "@Cell::reorderNodes, no node in the cell\n";
    return;
  }

  // cette structure permet de trier les noeuds dans le sens trigo 
  struct sortedNode {
    double x0, y0;
    double x, y;
    sortedNode(double X0, double Y0, double X, double Y) : x0(X0), y0(Y0), x(X), y(Y) {}
    bool operator<(const sortedNode &rhs) const {
      double Angle = atan2(y - y0, x - x0);
      double rhsAngle = atan2(rhs.y - rhs.y0, rhs.x - rhs.x0);
      if (Angle < rhsAngle)
        return true;
      return false;
    }
  };

  double xc = 0.0;
  double yc = 0.0;
  for (size_t i = 0; i < nodes.size(); i++) {
    xc += nodes[i].pos.x;
    yc += nodes[i].pos.y;
  }
  xc /= (double)(nodes.size());
  yc /= (double)(nodes.size());

  std::multiset<sortedNode> sortedNodes;

  for (size_t i = 0; i < nodes.size(); i++) {
    sortedNodes.insert(sortedNode(xc, yc, nodes[i].pos.x, nodes[i].pos.y));
  }

  size_t currentId = 0;
  for (auto s : sortedNodes) {
    nodes[currentId].pos.x = s.x;
    nodes[currentId].pos.y = s.y;
    currentId++;
  }
}

/**
 * @brief ajoute ou supprime un voisin. isNear doit être fourni et la methode déterminera si il faut ajouter ou
 *        bien supprimer (c'est basé sur le fait que la 'liste' de voisin est un std::set, et donc ordonnée et le
 *        'find' est rapide)
 *
 * @param ci cellule i
 * @param cj cellule j
 * @param in noeud n dans la cellule i
 * @param jn noeud n dans la cellule j
 * @param isNEAR à tester avant l'appel
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
void Cell::connectOrderedNodes(double width, double kn_, double kr_, double mz_max_, bool closed) {
  if (nodes.empty()) {
    std::cout << "@Cell::connectOrderedNodes, Cannot connect the nodes because there's no node!\n";
    return;
  }

  radius = 0.5 * width;
  bars.clear();

  for (size_t i = 1; i < nodes.size(); i++) {
    bars.push_back(Bar(i - 1, i));
  }
  if (closed == true)
    bars.push_back(Bar(nodes.size() - 1, 0));

  // init bars
  for (size_t b = 0; b < bars.size(); b++) {
    bars[b].init(kn_, nodes[bars[b].i].pos, nodes[bars[b].j].pos);
  }

  // init nodes
  if (nodes.size() > 2) {
    nodes[0].init(kr_, mz_max_, nodes.size() - 1, 1);
    for (size_t n = 1; n < nodes.size() - 1; n++) {
      nodes[n].init(kr_, mz_max_, n - 1, n + 1);
    }
    if (closed == true)
      nodes.back().init(kr_, mz_max_, nodes.size() - 2, 0);
  }
}