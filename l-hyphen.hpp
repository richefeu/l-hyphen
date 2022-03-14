// l-hyphen (lodge hyphen) is a multi-purpose system made of cells that are themeself composed of bars
//
// Zoom cell -- Z-cells -- H-cells
// 12-cells
// cells made of bars (CMB)
// bar-cells
// P-cells
// cells of bars
// CIC (cells in cells)
// bars for cells (B4C) *****
// P-cells (P = poly, polygonal, plain, ... )
// Soft cells
// grenade - granada (espagnol)
// grain na damaged
// Joinable, l-hyphen
//
// ============================================
// Ce fichier contient tout le code du modèle !
// Pour l'utiliser il suffit de l'inclure

#ifndef MP_CBN_HPP
#define MP_CBN_HPP

#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

#include "svgtools.hpp"
#include "vec2.hpp"

// pour éviter d'utiliser des pointeurs (avec les problèmes de copies)
// on utilise plutôt un indice (size_t), et pour traiter le cas "pas d'indice"
// on défini le null_size_t (doit être interprété comme l'équivalent d'un nullptr)
const size_t null_size_t = std::numeric_limits<size_t>::max();

/**
 * @brief contact potentiel (voisin ou interaction)
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

  Neighbor(size_t ic_, size_t jc_, size_t in_, size_t jn_)
      : ic(ic_), jc(jc_), in(in_), jn(jn_), n(), contactState(0), fn(0.0), ft(0.0), glueState(0), fn_coh(0.0),
        ft_coh(0.0), kn_coh(0.0), kt_coh(0.0), fn_coh_max(0.0), ft_coh_max(0.0), yieldPower(2.0) {}
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

#define VELOCITY_CONTROL 0
#define FORCE_CONTROL 1
/**
 * @brief Un control contient des flags pour dire si une force ou une vitesse est
 *        imposée selon les direction x et y. Il concerne uniquement les noeuds
 *
 */
class Control {
public:
  int xmode;
  double xvalue;
  int ymode;
  double yvalue;

  /**
   * @brief Construct a new Control object
   *
   */
  Control(int xmode_, double xvalue_, int ymode_, double yvalue_) {
    xmode = xmode_;
    xvalue = xvalue_;
    ymode = ymode_;
    yvalue = yvalue_;
  }
};

/**
 * @brief un point ayant une masse (et d'autre propriétés)
 *
 */
class Node {
public:
  vec2r pos;
  vec2r vel;
  vec2r acc;

  vec2r force; // force résultante

  double mass;

  size_t ictrl; // null_size_t si pas de control

  size_t prevNode; // numéro du noeud précédant
  size_t nextNode; // numéro du noeud suivant (null_size_t si pas de barre)

  double kr;     // raideur angulaire [N.m(/radian)]
  double mz;     // moment entre les barres connectées
  double mz_max; // moment-seuil de plastification parfaite

  // Ctor
  Node(double x, double y)
      : pos(x, y), vel(), acc(), force(), mass(1.0), ictrl(null_size_t), prevNode(null_size_t), nextNode(null_size_t),
        kr(0.0), mz(0.0), mz_max(0.0) {}

  void init(double kr_, double mz_max_, size_t prev, size_t next) {
    prevNode = prev;
    nextNode = next;
    kr = kr_;
    mz_max = mz_max_;
  }
};

/**
 * @brief une barre extensible entre deux noeuds
 *
 */
class Bar {
public:
  size_t i, j; // le numéro des noeuds dans Cell::nodes
  double l0;   // longueur initiale de la barre [m]
  double kn;   // raideur suivant l'axe [N/m]
  double fn;   // force axiale [N]

  Bar(size_t i_, size_t j_) : i(i_), j(j_), l0(0.0), kn(0.0) {}

  void init(double kn_, vec2r &posi, vec2r &posj) {
    vec2r branch = posj - posi;
    l0 = branch.normalize();
    kn = kn_;
  }
};

/**
 * @brief une cellule contenant un système de barres extensibles
 *
 */
class Cell {
public:
  std::vector<Node> nodes;
  std::vector<Bar> bars;
  std::set<Neighbor> neighbors;

  double radius; // un seul rayon pour toute la cellule [m]
  
  Cell() : radius(0.0) {}

  /**
   * @brief ré-ordonne la position des noeuds de la cellule en fonction de leur angle
   *
   */
  void reorderNodes() {
    if (nodes.empty()) {
      std::cout << "@Cell::reorderNodes, no node in the cell\n";
      return;
    }

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
  void insertOrRemove(size_t ci, size_t cj, size_t in, size_t jn, bool isNEAR) {
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
  void connectOrderedNodes(double width, double kn_, double kr_, double mz_max_, bool closed = true) {
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
};

/**
 * @brief cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *        lors de la définition de certaine cellules
 *
 */
struct RegularCellDataset {
  int nbFaces;
  double x;
  double y;
  double rot;
  double Rext;
  double barWidth;
};

/**
 * @brief cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *        lors de la définition de certaine cellules
 *
 */
struct TwoPointsDataset {
  double xo;
  double yo;
  double xe;
  double ye;
  double barWidth;
};

/**
 * @brief cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *        lors de la définition de certaine cellules
 *
 */
struct CellProperties {
  double kn;
  double kr;
  double mz_max;
};

/**
 * @brief Le système entier
 *
 */
class Sample {
public:
  std::vector<Cell> cells;
  std::vector<Control> controls;

  // pour les sorties SVG
  double xmin, xmax;
  double ymin, ymax;

  double dt;
  double dt_2;
  double dt2_2;

  double globalViscosity;
  double numericalDissipation;

  vec2r gravity;

  double distVerlet; // entre noeuds et barres

  // parametres mécaniques d'interactions entre les cellules
  double kn; // raideur normale de contact
  double kt; // raideur tangentielle de contact
  // double kn_coh; // raideur normale de cohésion
  // double kt_coh; // raideur tengentielle de cohésion
  double mu;   // coefficient de frottement (entre les cellules)
  double fadh; // force normale d'adhésion au contact

  int nstep;
  int nstepPeriodVerlet;
  int nstepPeriodSVG;
  int nstepPeriodRecord;
  int nstepPeriodConf;
  int isvg;
  int iconf;

  /**
   * @brief Construct a new Sample object
   *
   */
  Sample()
      : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), dt(0.0), globalViscosity(0.0), numericalDissipation(0.0), gravity(),
        distVerlet(0.0), kn(1000.0), kt(1000.0), mu(0.0), fadh(0.0), nstep(1000), nstepPeriodVerlet(1),
        nstepPeriodSVG(10), nstepPeriodRecord(1), nstepPeriodConf(10), isvg(0), iconf(0) {}

  /**
   * @brief affiche un petit entete sympatique
   * 
   */
  void head() {
    std::cout << "    L-HYPHEN\n";
    std::cout << "    _    _\n";
    std::cout << "   / \\__/ \\\n";
    std::cout << "  |  o  <  |\n";
    std::cout << "   \\_\\ _/_/\n";
    std::cout << "     /_/\n";
    std::cout << "    /_/\n\n";
  }

  /**
   * @brief trouver les limites de la zone dessinée dans les fichiers svg
   *
   * @param factor facteur multiplicateur
   */
  void findDisplayArea(double factor = 1.0) {
    if (cells.empty()) {
      std::cout << "@Sample::findDisplayArea, cannot find the scene limits because no cell has been set!\n";
      return;
    }

    xmin = xmax = cells[0].nodes[0].pos.x;
    ymin = ymax = cells[0].nodes[0].pos.y;
    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {

        if (cells[c].nodes[n].pos.x < xmin)
          xmin = cells[c].nodes[n].pos.x;
        if (cells[c].nodes[n].pos.x > xmax)
          xmax = cells[c].nodes[n].pos.x;
        if (cells[c].nodes[n].pos.y < ymin)
          ymin = cells[c].nodes[n].pos.y;
        if (cells[c].nodes[n].pos.y > ymax)
          ymax = cells[c].nodes[n].pos.y;
      }
    }

    double halfL = 0.5 * (xmax - xmin);
    double halfH = 0.5 * (ymax - ymin);
    double x0 = 0.5 * (xmin + xmax);
    double y0 = 0.5 * (ymin + ymax);
    xmin = x0 - factor * halfL;
    xmax = x0 + factor * halfL;
    ymin = y0 - factor * halfH;
    ymax = y0 + factor * halfH;
  }

  /**
   * @brief ajoute une cellule de forme polygonale régulière
   *
   * @param h les données pour définir un polyèdre régulier
   * @param p les propriétés mécaniques
   */
  void addRegularPolygonalCell(RegularCellDataset h, CellProperties p) {
    Cell C;
    double a0 = 2.0 * M_PI / (double)(h.nbFaces);
    double Rcell = h.Rext - 0.5 * h.barWidth;
    for (double a = 0.0; a < 2.0 * M_PI - 0.99 * a0; a += a0) {
      double theta = a0 + h.rot + a;
      C.nodes.emplace_back(h.x + Rcell * cos(theta), h.y + Rcell * sin(theta));
    }
    C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, true);
    cells.push_back(C);
  }

  /**
   * @brief ajoute une boite à partir de deux points
   *
   * @param h les deux points qui definissent les dimensions extérieures de la boite
   * @param p les propriétés mécaniques
   */
  void addBoxCell(TwoPointsDataset h, CellProperties p) {
    Cell C;
    double r = 0.5 * h.barWidth;
    C.nodes.emplace_back(h.xo + r, h.yo + r);
    C.nodes.emplace_back(h.xe - r, h.yo + r);
    C.nodes.emplace_back(h.xe - r, h.ye - r);
    C.nodes.emplace_back(h.xo + r, h.ye - r);
    C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, true);
    cells.push_back(C);
  }

  /**
   * @brief ajoute une ligne (cellule non fermée) entre deux points
   *
   * @param h les deux points
   * @param p les propriétés mécaniques
   */
  void addLine(TwoPointsDataset h, CellProperties p) {
    Cell C;
    C.nodes.emplace_back(h.xo, h.yo);
    C.nodes.emplace_back(h.xe, h.ye);
    C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, false);
    cells.push_back(C);
  }

  /**
   * @brief
   *
   * @param nx
   * @param ny
   * @param horizontalDistance
   * @param verticalDistance
   * @param h
   * @param p
   */
  void addRegularPolygonalCellsOnTriangularGrid(int nx, int ny, double horizontalDistance, double verticalDistance,
                                                RegularCellDataset h, CellProperties p) {

    double x0 = h.x;
    double y0 = h.y;

    int odd = 0;
    double xshift = 0.0;
    for (int iy = 0; iy < ny; iy++) {
      double y = y0 + iy * verticalDistance;
      for (double ix = 0; ix < nx - odd; ix++) {
        double x = x0 + xshift + ix * horizontalDistance;
        h.x = x;
        h.y = y;
        addRegularPolygonalCell(h, p);
      }
      odd = 1 - odd;
      if (odd == 1)
        xshift = 0.5 * horizontalDistance;
      else
        xshift = 0.0;
    }
  }

  /**
   * @brief
   *
   * @param nx
   * @param ny
   * @param horizontalDistance
   * @param xleft
   * @param ybottom
   * @param barWidth
   * @param p
   */
  void addSquareBrickWallCells(int nx, int ny, double horizontalDistance, double xleft, double ybottom, double barWidth,
                               CellProperties p) {
    addRegularPolygonalCellsOnTriangularGrid(
        nx, ny, horizontalDistance, horizontalDistance,
        {.nbFaces = 4,
         .x = xleft - horizontalDistance,
         .y = ybottom - horizontalDistance,
         .rot = M_PI / 4.,
         .Rext = 0.5 * (horizontalDistance - barWidth) * sqrt(2.0) + 0.5 * barWidth,
         .barWidth = barWidth},
        p);
  }

  /**
   * @brief crée une structure en nid d'abeille
   *
   * @param nx nombre de cellule suivant x
   * @param ny nombre de lignes
   * @param CellExternWidth largeur externe d'une cellule hexagonale
   * @param xleft  position x la plus à gauche
   * @param ybottom position y la plus en bas
   * @param barWidth épaisseur des barres
   * @param p paramètres mécaniques de toutes les cellules
   */
  void addHoneycombCells(int nx, int ny, double CellExternWidth, double xleft, double ybottom, double barWidth,
                         CellProperties p) {
    double verticalDistance = 1.5 * CellExternWidth / sqrt(3.0);
    addRegularPolygonalCellsOnTriangularGrid(nx, ny, CellExternWidth, verticalDistance,
                                             {.nbFaces = 6,
                                              .x = xleft + 0.5 * CellExternWidth + 0.5 * barWidth,
                                              .y = ybottom + 0.5 * verticalDistance + 0.5 * barWidth,
                                              .rot = M_PI / 6.0,
                                              .Rext = (CellExternWidth - barWidth) / sqrt(3.0) + 0.5 * barWidth,
                                              .barWidth = barWidth},
                                             p);
  }

  /**
   * @brief Utiliser cette méthode pour definir le pas de temps (comme ça dt/2 et dt^2/2 seront pré-calculées)
   *
   * @param dt_ pas de temps
   */
  void setTimeStep(double dt_) {
    dt = dt_;
    dt_2 = 0.5 * dt;
    dt2_2 = dt_2 * dt;
  }

  void setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max, double yieldPower) {

    for (size_t ci = 0; ci < cells.size(); ci++) {
      for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
           ++InterIt) {

        Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));
        if (Inter->glueState == 1) {
          Inter->kn_coh = kn_coh;
          Inter->kt_coh = kt_coh;
          Inter->fn_coh_max = fn_coh_max;
          Inter->ft_coh_max = ft_coh_max;
          Inter->yieldPower = yieldPower;
        }
      }
    }
  }

  /**
   * @brief ajoute un control au noeud n de la cellule c
   *
   * @param c numéro de la cellule
   * @param n numéro du noeud au sein de la cellule c
   * @param xmode VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
   * @param xvalue valeur suivant x (une force ou une vitesse selon xmode)
   * @param ymode VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
   * @param yvalue valeur suivant y (une force ou une vitesse selon ymode)
   */
  void setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue) {
    Control C(xmode, xvalue, ymode, yvalue);
    controls.push_back(C);
    cells[c].nodes[n].ictrl = controls.size() - 1;
  }

  /**
   * @brief ajoute un control à tous les noeuds de la cellule c
   *
   * @param c numéro de la cellule
   * @param xmode VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
   * @param xvalue valeur suivant x (une force ou une vitesse selon xmode)
   * @param ymode VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
   * @param yvalue valeur suivant y (une force ou une vitesse selon ymode)
   */
  void setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue) {
    Control C(xmode, xvalue, ymode, yvalue);
    controls.push_back(C);
    size_t ictrl = controls.size() - 1;
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].ictrl = ictrl;
    }
  }

  /**
   * @brief définir un même control à tous les noeuds qui se trouvent dans une zone rectangulaire
   *
   * @param xmin left limit of the box
   * @param xmax right limit of the box
   * @param ymin bottom limit of the box
   * @param ymax top limit of the box
   * @param xmode imposed mode in the x-direction (VELOCITY_CONTROL or FORCE_CONTROL)
   * @param xvalue imposed value in the x-direction
   * @param ymode imposed mode in the y-direction (VELOCITY_CONTROL or FORCE_CONTROL)
   * @param yvalue imposed value in the y-direction
   */
  void setNodeControlInBox(double xmin, double xmax, double ymin, double ymax, int xmode, double xvalue, int ymode,
                           double yvalue) {
    Control C(xmode, xvalue, ymode, yvalue);

    controls.push_back(C);
    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        vec2r pos = cells[c].nodes[n].pos;
        if (pos.x >= xmin && pos.x <= xmax && pos.y >= ymin && pos.y <= ymax) {
          cells[c].nodes[n].ictrl = controls.size() - 1;
        }
      }
    }
  }

  /**
   * @brief vérifie si un noeud (cell ci, node in) et une barre (cell cj, barre commancant par le noeud jn) sont
   *        proches. Selon le cas, la paire dans la liste de voisins est soit ajoutée soit retirée
   *
   * @param ci cellule i
   * @param cj cellule j
   * @param in numero du noeud dans la cellule ci
   * @param jn numero du noeud du début d'une barre dans la cellule cj
   * @param epsilonEnds longueur à ignoré aux extremités de la barre
   */
  void addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds = 0.0) {
    // jn est vu ici comme l'indice d'une barre (noeud du debut en fait)

    if (cells[cj].nodes[jn].nextNode == null_size_t) { // cas sans barre (une extrémité)

      vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
      double sqrDist = norm2(branch);
      double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
      bool isNear = sqrDist < sumR * sumR;
      cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

    } else { // la barre reliant jn à jnext existe

      size_t jnext = cells[cj].nodes[jn].nextNode;
      vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
      vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
      double u_length = u.normalize();
      double proj = b * u;

      if (proj <= epsilonEnds) { // début de la barre

        double sqrDist = norm2(b);
        double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
        bool isNear = sqrDist < sumR * sumR;
        cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

      } else if (proj >= u_length - epsilonEnds) { // fin de la barre

        b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
        double sqrDist = norm2(b);
        double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
        bool isNear = sqrDist < sumR * sumR;
        cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

      } else { // sur la barre

        vec2r t(-u.y, u.x);
        double dist = fabs(b * t) - (cells[ci].radius + cells[cj].radius);
        bool isNear = dist < distVerlet;
        cells[ci].insertOrRemove(ci, cj, in, jn, isNear);
      }
    }
  }

  /**
   * @brief cette fonction permet de lire un fichier contenant une liste de positions x,y avec numéro
   *        de cellule. Peut importe les numéros du moment qu'ils sont différents pour chaque cellule.
   *        Il faut faire un pré-nettoyage pour qu'il y ait au moins 3 noeuds par cellule.
   *        Toutes les cellules sont fermées.
   *
   * @param name nom du fichier
   * @param barWidth épaisseur de toutes les barres
   */
  void readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max) {
    std::ifstream file(name);

    struct data {
      double x, y;
      size_t id;
      bool operator<(const data &rhs) const {
        if (id < rhs.id)
          return true;
        return false;
      }
    };

    std::multiset<data> dataset;

    data D;
    while (file.good()) {
      file >> D.x >> D.y >> D.id;
      if (file.eof())
        break;
      dataset.insert(D);
    }

    int cellId = -1;
    Cell C;
    for (auto i : dataset) {
      if (cellId != (int)(i.id)) {
        if (cellId != -1) {
          cells.push_back(C);
          C.nodes.clear();
        }
        cellId = (int)(i.id);
      }
      C.nodes.emplace_back(i.x, i.y);
    }
    if (!C.nodes.empty())
      cells.push_back(C);
    /*
    for (size_t i = 0; i < cells.size(); i++) {
      std::cout << i << "  " << cells[i].nodes.size() << '\n';
    }
    */

    for (size_t i = 0; i < cells.size(); i++) {
      cells[i].reorderNodes();
      cells[i].connectOrderedNodes(barWidth, Kn, Kr, Mz_max, true);
    }

    /*
    for (size_t i = 0; i < cells.size(); i++) {
      std::cout << i << " ----- \n";
      for (size_t n = 0; n < cells[i].nodes.size();n++) {
        std::cout << cells[i].nodes[n].pos << '\n';
      }
    }
    */
  }

  /**
   * @brief met à jour la liste des voisins
   *
   */
  void updateNeighbors() {
    // TODO : regarder d'abort la proximité des cellules (?)

    for (size_t ci = 0; ci < cells.size(); ci++) {
      // TODO: cas même cellule ici (?)
      for (size_t cj = ci + 1; cj < cells.size(); cj++) {

        for (size_t in = 0; in < cells[ci].nodes.size(); in++) {
          for (size_t jn = 0; jn < cells[cj].nodes.size(); jn++) {
            addNodeToBarNeighbor(ci, cj, in, jn);
          }
        }

        for (size_t jn = 0; jn < cells[cj].nodes.size(); jn++) {
          for (size_t in = 0; in < cells[ci].nodes.size(); in++) {
            addNodeToBarNeighbor(cj, ci, jn, in, cells[cj].radius);
          }
        }
      }
    }
  }

  /**
   * @brief met un point de colle aux endroits où la distance est suffisement proche.
   *
   * @param epsilonDist distane maximale
   */
  void glue(double epsilonDist) {
    updateNeighbors();

    for (size_t ci = 0; ci < cells.size(); ci++) {
      for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
           ++InterIt) {

        Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));
        size_t cj = Inter->jc;
        size_t in = Inter->in;
        size_t jn = Inter->jn;

        if (cells[cj].nodes[jn].nextNode == null_size_t) { // pas de barre -> ctc disque-disque

          vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
          double sqrDist = norm2(branch);
          double sumR = cells[ci].radius + cells[cj].radius + epsilonDist;
          if (sqrDist - sumR * sumR < 0.0) {
            Inter->glueState = 1;
            Inter->fn_coh = 0.0;
            Inter->ft_coh = 0.0;
          }

        } else {

          size_t jnext = cells[cj].nodes[jn].nextNode; // jnext est la fin de la barre (jn c'est le début)
          vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
          vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
          double u_length = u.normalize();
          double proj = b * u;

          if (proj <= 0.0) {

            double sqrDist = norm2(b);
            double sumR = cells[ci].radius + cells[cj].radius + epsilonDist;
            if (sqrDist - sumR * sumR < 0.0) {
              Inter->glueState = 1;
              Inter->fn_coh = 0.0;
              Inter->ft_coh = 0.0;
            }

          } else if (proj >= u_length) {

            b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
            double sqrDist = norm2(b);
            double sumR = cells[ci].radius + cells[cj].radius + epsilonDist;
            if (sqrDist - sumR * sumR < 0.0) {
              Inter->glueState = 1;
              Inter->fn_coh = 0.0;
              Inter->ft_coh = 0.0;
            }

          } else {

            vec2r t(-u.y, u.x);
            double sumR = cells[ci].radius + cells[cj].radius;
            double dist = b * t; // dot product
            double dn = fabs(dist) - sumR;
            if (dn < epsilonDist) {
              Inter->glueState = 1;
              Inter->fn_coh = 0.0;
              Inter->ft_coh = 0.0;
            }
          }
        }
      }
    }
  }

  /**
   * @brief calcul des forces d'interaction (entre cellules)
   *
   */
  void computeInteractionForces() {

    vec2r vrel, T;

    for (size_t ci = 0; ci < cells.size(); ci++) {
      for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
           ++InterIt) {

        // il faut convertir l'iterateur vers un set, car sinon on ne pourra pas modifier les valeurs.
        Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));

        size_t cj = Inter->jc;
        size_t in = Inter->in;
        size_t jn = Inter->jn;

        if (cells[cj].nodes[jn].nextNode == null_size_t) {
          // =============================================================================
          // On est dans le cas où la barre dans cj n'existe pas.
          // En d'autre terme, il s'agit d'une interaction entre un disque de ci
          // et un disque de cj (qui correspond à une extrémité de barre)

          vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
          double sqrDist = norm2(branch);
          double sumR = cells[ci].radius + cells[cj].radius;

          if (sqrDist - sumR * sumR < 0.0) { // overlap
            double b_length = sqrt(sqrDist);
            double dn = b_length - sumR;
            Inter->n = branch / b_length;
            Inter->fn = -kn * dn - fadh;
          } else {
            if (Inter->glueState == 1) {
              double b_length = sqrt(sqrDist);
              Inter->n = branch / b_length;
            }
            Inter->fn = 0.0;
            Inter->ft = 0.0;
          }

          if (Inter->glueState == 1) {
            // ...
            vec2r vrel = cells[cj].nodes[jn].vel - cells[ci].nodes[in].vel;
            // vec2r T(-n.y, n.x);
            Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
            // Inter->ft_coh += -kt_coh * vrel * T * dt;

            // gerer la rupture
          }

          // transfert des force vers les noeuds concernés
          vec2r finc = (Inter->fn + Inter->fn_coh) * Inter->n;
          cells[ci].nodes[in].force += finc;
          cells[cj].nodes[jn].force -= finc;

        } else {
          // =============================================================================
          // un disque (ci,in) avec une barre (cj, jn -> jnext)

          size_t jnext = cells[cj].nodes[jn].nextNode;
          // jnext est le numéro du noeud à la fin de la barre dans la cellule cj (jn c'est le début)
          vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
          vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
          double u_length = u.normalize();
          double proj = b * u;

          // on doit gérer le calcul des forces (contact, frottement, cohesion...)
          // pour chacun des 3 cas suivants :
          //   disque (ci, in)   ---   disque (cj, jn) début de la barre
          //   disque (ci, in)   ---   disque (cj, jnext) fin de la barre
          //   disque (ci, in)   ---   barre (cj, jn--jnext)
          if (proj <= 0.0) { // ====================== disque j du début

            double sqrDist = norm2(b);
            double sumR = cells[ci].radius + cells[cj].radius;
            if (sqrDist - sumR * sumR < 0.0) {
              Inter->contactState = 1;
              double b_lenght = sqrt(sqrDist);
              double dn = b_lenght - sumR;
              Inter->n = b / b_lenght;
              Inter->fn = -kn * dn - fadh;

            } else {
              Inter->contactState = 0;
            }

            // transfert des force vers les noeuds concernés
            vec2r finc = Inter->fn * Inter->n;
            cells[ci].nodes[in].force += finc;
            cells[cj].nodes[jn].force -= finc;

          } else if (proj >= u_length) { // ====================== disque jnext (de fin)

            b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
            double sqrDist = norm2(b);
            double sumR = cells[ci].radius + cells[cj].radius;
            if (sqrDist - sumR * sumR < 0.0) {
              Inter->contactState = 1;
              double b_lenght = sqrt(sqrDist);
              double dn = b_lenght - sumR;
              Inter->n = b / b_lenght;
              Inter->fn = -kn * dn - fadh;

              vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jnext].vel;
              T.set(-Inter->n.y, Inter->n.x);
              Inter->ft -= kt * (vrel * T) * dt;
              double threshold = mu * Inter->fn;
              if (Inter->ft > threshold) {
                Inter->ft = threshold;
              } else if (Inter->ft < -threshold) {
                Inter->ft = -threshold;
              }

              if (Inter->glueState == 0)
                Inter->fn -= fadh;

              // transfert des force vers les noeuds concernés
              vec2r finc = Inter->fn * Inter->n;
              cells[ci].nodes[in].force += finc;
              cells[cj].nodes[jnext].force -= finc;

            } else {
              Inter->contactState = 0;
              Inter->fn = 0.0;
              Inter->ft = 0.0;
            }

            if (Inter->glueState == 1) {
              if (Inter->contactState == 0) {
                double b_lenght = sqrt(sqrDist);
                Inter->n = b / b_lenght;
                vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jnext].vel;
                T.set(-Inter->n.y, Inter->n.x);
              }

              Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
              Inter->ft_coh += -Inter->kt_coh * vrel * T * dt;

              if (Inter->fn_coh > 0.0)
                Inter->fn_coh = 0.0;

              // rupture
              double zeta = -Inter->fn_coh / Inter->fn_coh_max +
                            pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
              if (zeta > 0.0) {
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
              } else {
                // transfert des force vers les noeuds concernés
                vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * T;
                cells[ci].nodes[in].force += finc;
                cells[cj].nodes[jnext].force -= finc;
              }
            }

          } else { // ====================== sur la barre

            vec2r urot(-u.y, u.x); // on tourne u de 90°
            double wend = 0.0;
            double wbeg = 0.0;
            double sumR = cells[ci].radius + cells[cj].radius;
            double dist = b * urot; // dot product
            double dn = fabs(dist) - sumR;

            if (dn < 0.0) { // overlap
              Inter->contactState = 1;
              Inter->n = urot;
              if (dist < 0.0) {
                Inter->n *= -1.0;
              }
              // ici le vecteur normal est orienté de la barre vers le disque

              Inter->fn = -kn * dn;

              // vitesse relative noeud par rapport à barre
              wend = proj / u_length;
              wbeg = 1.0 - wend;
              // on neglige l'épaisseur des barres
              // à voir si c'est ok car on fait une approximation
              vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
              T.set(-Inter->n.y, Inter->n.x);
              Inter->ft -= kt * (vrel * T) * dt;
              double threshold = mu * Inter->fn;
              if (Inter->ft > threshold) {
                Inter->ft = threshold;
              } else if (Inter->ft < -threshold) {
                Inter->ft = -threshold;
              }

              if (Inter->glueState == 0)
                Inter->fn -= fadh;

              // transfert des force vers les noeuds concernés
              // n est orienté de la barre vers le noeud
              // donc fn positif correspond à une force de la barre qui pousse le noeud

              vec2r finc = Inter->fn * Inter->n + Inter->ft * T;
              cells[ci].nodes[in].force += finc;

              cells[cj].nodes[jn].force -= wbeg * finc;
              cells[cj].nodes[jnext].force -= wend * finc;

            } else {
              Inter->contactState = 0;
              Inter->fn = 0.0;
              Inter->ft = 0.0;
            }

            if (Inter->glueState == 1) {

              // si on se trouve dans un cas sans contact, certaines variables n'ont pas été calculés
              if (Inter->contactState == 0) {
                Inter->n = urot;
                if (dist < 0.0) {
                  Inter->n *= -1.0;
                }

                wend = proj / u_length;
                wbeg = 1.0 - wend;
                vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
                T.set(-Inter->n.y, Inter->n.x);
              }

              Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
              Inter->ft_coh += -Inter->kt_coh * vrel * T * dt;

              if (Inter->fn_coh > 0.0)
                Inter->fn_coh = 0.0;

              // rupture
              double zeta = -Inter->fn_coh / Inter->fn_coh_max +
                            pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
              if (zeta > 0.0) {
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
              } else {
                // transfert des force vers les noeuds concernés
                vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * T;
                cells[ci].nodes[in].force += finc;

                cells[cj].nodes[jn].force -= wbeg * finc;
                cells[cj].nodes[jnext].force -= wend * finc;
              }
            } // fin glueState
          }   // fin "sur la barre"
        }

      } // boucle sur les voisins dans ci
    }   // boucle sur les cellules ci
  }

  /**
   * @brief Calcul des forces internes aux cellules
   *
   */
  void computeNodeForces() {

    // Forces normales dans les barres
    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t b = 0; b < cells[c].bars.size(); b++) {
        size_t i = cells[c].bars[b].i;
        size_t j = cells[c].bars[b].j;
        vec2r posi = cells[c].nodes[i].pos;
        vec2r posj = cells[c].nodes[j].pos;

        vec2r branch = posj - posi;
        vec2r n = branch;
        double l = n.normalize();            // n est orienté de i vers j
        double dn = l - cells[c].bars[b].l0; // dn positif = allongement
        double fn = -cells[c].bars[b].kn * dn;

        // TODO: ajouter plasticité dans les barres

        vec2r finc = fn * n; // en cas d'allongement finc est orienté de j vers i
        cells[c].nodes[i].force -= finc;
        cells[c].nodes[j].force += finc;
      }
    }

    // Moment au niveau des noeuds (entre les barres d'une même cellule).
    // Puisque les noeuds n'ont pas de ddl en rotation, imposer un moment est fait ici
    // en imposant une force à chaque extrémité de barre
    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        if (cells[c].nodes[n].prevNode == null_size_t || cells[c].nodes[n].nextNode == null_size_t) {
          continue;
        }
        size_t prev = cells[c].nodes[n].prevNode;
        size_t next = cells[c].nodes[n].nextNode;

        double omegaNext = getRotationVelocityBar(cells[c].nodes[n].pos, cells[c].nodes[next].pos,
                                                  cells[c].nodes[n].vel, cells[c].nodes[next].vel);
        double omegaPrev = getRotationVelocityBar(cells[c].nodes[n].pos, cells[c].nodes[prev].pos,
                                                  cells[c].nodes[n].vel, cells[c].nodes[prev].vel);
        cells[c].nodes[n].mz += -cells[c].nodes[n].kr * (omegaNext - omegaPrev) * dt;

        // TODO: ajouter la possibilité de désactiver la plasticité
        if (cells[c].nodes[n].mz > cells[c].nodes[n].mz_max) {
          cells[c].nodes[n].mz = cells[c].nodes[n].mz_max;
        } else if (cells[c].nodes[n].mz < -cells[c].nodes[n].mz_max) {
          cells[c].nodes[n].mz = -cells[c].nodes[n].mz_max;
        }

        vec2r nextVector = cells[c].nodes[next].pos - cells[c].nodes[n].pos;
        vec2r prevVector = cells[c].nodes[n].pos - cells[c].nodes[prev].pos;

        // TODO: ici c'est optimisable....
        vec2r nextDir = nextVector;
        vec2r prevDir = prevVector;
        double nextL = nextDir.normalize();
        double prevL = prevDir.normalize();
        vec2r prevDirRot(-prevDir.y, prevDir.x);

        // mz sur next
        vec2r nextDirRot(-nextDir.y, nextDir.x);
        double F = cells[c].nodes[n].mz / nextL;
        vec2r finc = F * nextDirRot;
        cells[c].nodes[next].force += finc;
        cells[c].nodes[n].force -= finc;

        // -mz sur prev
        F = cells[c].nodes[n].mz / prevL;
        finc = F * prevDirRot;
        cells[c].nodes[prev].force += finc;
        cells[c].nodes[n].force -= finc;
      }
    }

    // visco globale
    if (globalViscosity > 0.0) {
      for (size_t c = 0; c < cells.size(); c++) {
        for (size_t n = 0; n < cells[c].nodes.size(); n++) {
          cells[c].nodes[n].force -= globalViscosity * cells[c].nodes[n].vel;
        }
      }
    }
  }

  /**
   * @brief compute the node accelerations
   *
   */
  void NodeAccelerations() {

    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        cells[c].nodes[n].acc = gravity;
        cells[c].nodes[n].force.reset();
        if (cells[c].nodes[n].ictrl != null_size_t) {
          Control *ctrl = &(controls[cells[c].nodes[n].ictrl]);
          if (ctrl->xmode == FORCE_CONTROL)
            cells[c].nodes[n].force.x += ctrl->xvalue;
          if (ctrl->ymode == FORCE_CONTROL)
            cells[c].nodes[n].force.y += ctrl->yvalue;
        }
      }
    }

    computeNodeForces();
    computeInteractionForces();

    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        cells[c].nodes[n].acc += cells[c].nodes[n].force / cells[c].nodes[n].mass;
      }
    }
  }

  /**
   * @brief un pas en velocity verlet
   *
   */
  void SingleStep() {

    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        if (cells[c].nodes[n].ictrl == null_size_t) {
          cells[c].nodes[n].pos += cells[c].nodes[n].vel * dt + cells[c].nodes[n].acc * dt2_2;
          cells[c].nodes[n].vel += cells[c].nodes[n].acc * dt_2;
        } else {
          Control *ctrl = &(controls[cells[c].nodes[n].ictrl]);
          if (ctrl->xmode == FORCE_CONTROL) {
            cells[c].nodes[n].pos.x += cells[c].nodes[n].vel.x * dt + cells[c].nodes[n].acc.x * dt2_2;
            cells[c].nodes[n].vel.x += cells[c].nodes[n].acc.x * dt_2;
          } else if (ctrl->xmode == VELOCITY_CONTROL) {
            cells[c].nodes[n].vel.x = ctrl->xvalue;
            cells[c].nodes[n].force.x = 0.0;
            cells[c].nodes[n].pos.x += cells[c].nodes[n].vel.x;
          }

          if (ctrl->ymode == FORCE_CONTROL) {
            cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y * dt + cells[c].nodes[n].acc.y * dt2_2;
            cells[c].nodes[n].vel.y += cells[c].nodes[n].acc.y * dt_2;
          } else if (ctrl->ymode == VELOCITY_CONTROL) {
            cells[c].nodes[n].vel.y = ctrl->yvalue;
            cells[c].nodes[n].force.y = 0.0;
            cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y;
          }
        }
      }
    }

    NodeAccelerations();

    for (size_t c = 0; c < cells.size(); c++) {
      for (size_t n = 0; n < cells[c].nodes.size(); n++) {
        if (cells[c].nodes[n].ictrl == null_size_t) {
          cells[c].nodes[n].vel += cells[c].nodes[n].acc * dt_2;
        } else {
          Control *ctrl = &(controls[cells[c].nodes[n].ictrl]);
          if (ctrl->xmode == FORCE_CONTROL) {
            cells[c].nodes[n].vel.x += cells[c].nodes[n].acc.x * dt_2;
          }
          if (ctrl->ymode == FORCE_CONTROL) {
            cells[c].nodes[n].vel.y += cells[c].nodes[n].acc.y * dt_2;
          }
        }
      }
    }

    // une dissipation purement numérique
    if (numericalDissipation > 0.0) {
      for (size_t c = 0; c < cells.size(); c++) {
        for (size_t n = 0; n < cells[c].nodes.size(); n++) {
          cells[c].nodes[n].vel *= (1.0 - numericalDissipation);
        }
      }
    }
  }

  /**
   * @brief run the simulation!
   *
   */
  void run() {
    for (int step = 0; step < nstep; step++) {
      if (step % nstepPeriodVerlet == 0) {
        updateNeighbors();
      }

      if (step % nstepPeriodSVG == 0) {
        saveSVG(isvg);
        isvg++;
      }

      SingleStep();
    }
  }

  /**
   * @brief Save the current configuration in a file
   *
   * @param fname Name of the file
   */
  void saveCONF(const char *fname) {
    std::ofstream file(fname);

    file << "kn " << kn << '\n';
    file << "kt " << kt << '\n';
    file << "mu " << mu << '\n';
    file << "gravity " << gravity << '\n';
    file << "numericalDissipation " << numericalDissipation << '\n';
    file << "globalViscosity " << globalViscosity << '\n';
    file << "distVerlet " << distVerlet << '\n';
    file << "dt " << dt << '\n';

    file << "nstep " << nstep << '\n';
    file << "nstepPeriodVerlet " << nstepPeriodVerlet << '\n';
    file << "nstepPeriodSVG " << nstepPeriodSVG << '\n';
    file << "nstepPeriodRecord " << nstepPeriodRecord << '\n';
    file << "nstepPeriodConf " << nstepPeriodConf << '\n';
    file << "isvg " << isvg << '\n';
    file << "iconf " << iconf << '\n';
    file << "SVG_area " << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << '\n';

    file << "cells " << cells.size() << '\n';
    for (size_t i = 0; i < cells.size(); i++) {
      /*
      file << grains[i].radius << ' ' << grains[i].mass << ' ' << grains[i].inertia << ' ' << grains[i].young << ' '
           << grains[i].waterContent << ' ' << grains[i].initialVolume << ' ' << grains[i].pos << ' ' <<
      grains[i].vel
           << ' ' << grains[i].rot << ' ' << grains[i].vrot << '\n';
      */
    }

    /*
        size_t nbNeighbors = 0;
        for (size_t i = 0; i < grains.size(); i++) {
          for (std::set<Neighbor>::iterator N = grains[i].neighbors.begin(); N != grains[i].neighbors.end(); ++N) {
            if (N->dn >= 0.0)
              continue;
            nbNeighbors++;
          }
        }
        file << "neighbors " << nbNeighbors << '\n';
        for (size_t i = 0; i < grains.size(); i++) {
          for (std::set<Neighbor>::iterator N = grains[i].neighbors.begin(); N != grains[i].neighbors.end(); ++N) {
            if (N->dn >= 0.0)
              continue;
            file << N->i << ' ' << N->j << ' ' << N->glueState << ' ' << N->pos << ' ' << N->nji << ' ' << N->dn <<
       ' '
       << N->dt
                 << ' ' << N->fnji << ' ' << N->ftji << '\n';
          }
        }
    */
  }

  void saveCONF(int ifile) {
    char fname[256];
    sprintf(fname, "conf%d", ifile);
    std::cout << "save CONF file: " << fname << '\n';
    saveCONF(fname);
  }

  /**
   * @brief Load a configuration that has been saved in a file
   *
   * @param fname The name of the file
   */
  void loadCONF(const char *fname) {
    std::ifstream file(fname);

    std::string token;
    file >> token;
    while (file.good()) {
      if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
        getline(file, token); // ignore the rest of the current line
        file >> token;        // next token
        continue;
      } else if (token == "gravity") {
        file >> gravity;
      } else if (token == "numericalDissipation") {
        file >> numericalDissipation;
      } else if (token == "globalViscosity") {
        file >> globalViscosity;
      } else if (token == "distVerlet") {
        file >> distVerlet;
      } else if (token == "dt") {
        double timestep;
        file >> timestep;
        setTimeStep(timestep);
      } else if (token == "nstep") {
        file >> nstep;
      } else if (token == "nstepPeriodVerlet") {
        file >> nstepPeriodVerlet;
      } else if (token == "nstepPeriodSVG") {
        file >> nstepPeriodSVG;
      } else if (token == "nstepPeriodRecord") {
        file >> nstepPeriodRecord;
      } else if (token == "nstepPeriodConf") {
        file >> nstepPeriodConf;
      } else if (token == "isvg") {
        file >> isvg;
      } else if (token == "iconf") {
        file >> iconf;
      } else if (token == "kn") {
        file >> kn;
      } else if (token == "kt") {
        file >> kt;
      } else if (token == "setGlueSameProperties") {
        double kn_coh;
        double kt_coh;
        double fn_coh_max;
        double ft_coh_max;
        double yieldPower;
        file >> kn_coh >> kt_coh >> fn_coh_max >> ft_coh_max >> yieldPower;
        setGlueSameProperties(kn_coh, kt_coh, fn_coh_max, ft_coh_max, yieldPower);
      } else if (token == "mu") {
        file >> mu;
      } else if (token == "fadh") {
        file >> fadh;
      } else if (token == "cells") {
        /*
        grains.clear();
        size_t nbGrains;
        file >> nbGrains;
        Grain G;
        for (size_t i = 0; i < nbGrains; i++) {
          file >> G.radius >> G.mass >> G.inertia >> G.young >> G.waterContent >> G.initialVolume >> G.pos >> G.vel
        >> G.rot >> G.vrot; grains.push_back(G);
        }
        */
      } else if (token == "neighbors") {
        /*
        for (size_t i = 0; i < grains.size(); i++) {
          grains[i].neighbors.clear();
        }

        size_t nbNeighbors;
        file >> nbNeighbors;
        Neighbor N;
        for (size_t k = 0; k < nbNeighbors; k++) {
          file >> N.i >> N.j >> N.glueState >> N.pos >> N.nji >> N.dn >> N.dt >> N.fnji >> N.ftji;
          grains[N.i].neighbors.insert(N);
        }
        */
      } else if (token == "addRegularPolygonalCell") {
        RegularCellDataset h;
        CellProperties p;
        file >> h.nbFaces >> h.x >> h.y >> h.rot >> h.Rext >> h.barWidth;
        file >> p.kn >> p.kr >> p.mz_max;
        addRegularPolygonalCell(h, p);
      } else if (token == "addSquareBrickWallCells") {
        CellProperties p;
        int nx, ny;
        double hdist;
        double xleft, ybottom;
        double barWidth;
        file >> nx >> ny >> hdist >> xleft >> ybottom >> barWidth;
        file >> p.kn >> p.kr >> p.mz_max;
        addSquareBrickWallCells(nx, ny, hdist, xleft, ybottom, barWidth, p);
      } else if (token == "findDisplayArea") {
        double d;
        file >> d;
        findDisplayArea(d);
      } else if (token == "glue") {
        double d;
        file >> d;
        glue(d);
      } else if (token == "setCellControl") {
        size_t icell;
        int xmode, ymode;
        double xvalue, yvalue;
        file >> icell >> xmode >> xvalue >> ymode >> yvalue;
        setCellControl(icell, xmode, xvalue, ymode, yvalue);
      } else if (token == "readNodeFile") {
        std::string fileName;
        double barWidth, Kn, Kr, Mz_max;
        file >> fileName >> barWidth >> Kn >> Kr >> Mz_max;
        readNodeFile(fileName.c_str(), barWidth, Kn, Kr, Mz_max);
      } else {
        std::cout << "@Sample::loadCONF, this token is not known: " << token << '\n';
      }

      file >> token;
    }
  }

  void loadCONF(int ifile) {
    char fname[256];
    sprintf(fname, "conf%d", ifile);
    loadCONF(fname);
  }

  /**
   * @brief fonction très basique pour dessiner les cellules avec des lignes bleues
   *        TODO: pour d'autre sorties graphiques, on verra après :)
   *
   * @param num numero du fichier
   * @param nameBase nom de base (dans lequel sera inséré le numéro)
   * @param Canvaswidth largeur de canvas (par affichae rapide)
   *
   * @see https://bestofcpp.com/repo/tomkwok-svgasm
   */
  void saveSVG(int num, const char *nameBase = "sample%04d.svg", int Canvaswidth = 400) {
    char name[256];
    sprintf(name, nameBase, num);
    std::cout << "save file: " << name << '\n';
    std::ofstream ofs(name);
    SVGfile svg(ofs);

    int CanvasHeight = Canvaswidth * (ymax - ymin) / (xmax - xmin);
    svg.begin(Canvaswidth, CanvasHeight);
    viewZone vz(0, 0, Canvaswidth, CanvasHeight);
    double delta = (xmax - xmin) * 0.05; // c'est une bordure
    vz.adjustRange(xmin - delta, xmax + delta, ymin - delta, ymax + delta);

    char opt[256];

    for (size_t c = 0; c < cells.size(); c++) {

      for (size_t b = 0; b < cells[c].bars.size(); b++) {
        sprintf(opt, "stroke:blue;fill:none;stroke-linecap:round;stroke-width:%g", 2 * cells[c].radius * vz.scalex);
        size_t i = cells[c].bars[b].i;
        size_t j = cells[c].bars[b].j;
        svg.line(vz, cells[c].nodes[i].pos.x, cells[c].nodes[i].pos.y, cells[c].nodes[j].pos.x, cells[c].nodes[j].pos.y,
                 opt);
      }
    }

    svg.end();
  }

private:
  /**
   * @brief Get the Rotation Velocity Bar object
   *
   * @param a Position point A
   * @param b Position point B
   * @param va Vitesse point A
   * @param vb Vitesse point B
   * @return double vitesse de rotation de solide rigide
   */
  double getRotationVelocityBar(vec2r &a, vec2r &b, vec2r &va, vec2r &vb) {
    vec2r u = b - a;
    double l = u.normalize();
    vec2r t(-u.y, u.x);
    return ((vb - va) * t / l);
  }
};

#endif /* MP_CBN_HPP */
