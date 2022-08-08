// l-hyphen (lodge hyphen) is a multi-purpose system made of cells that are themeself composed of bars
//
// ============================================
// Ce fichier contient tout le code du modèle !
// Pour l'utiliser il suffit de l'inclure

#ifndef L_HYPHEN_HPP
#define L_HYPHEN_HPP

#include <fstream>
#include <iostream>

#include <set>
#include <vector>

#include "svgtools.hpp"
#include "vec2.hpp"

#include "Bar.hpp"
#include "Cell.hpp"
#include "Control.hpp"
#include "Neighbor.hpp"
#include "Node.hpp"

/**
 * @brief cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *        lors de la définition de certaine cellules
 *
 */
struct named_arg_RegularCellDataset {
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
struct named_arg_TwoPointsDataset {
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
struct named_arg_CellProperties {
  double kn;
  double kr;
  double mz_max;
};

struct CellNodeID {
  size_t c, n;
  CellNodeID(size_t c_, size_t n_) : c(c_), n(n_) {}
};

struct CapturedNodes {
  std::vector<CellNodeID> cellNodeIDs;
  std::string filename;
  CapturedNodes(const char *name) : filename(name) {}
};

/**
 * @brief Le système entier
 *
 */
class Lhyphen {
public:
  std::vector<Cell> cells;
  std::vector<Control> controls;

  std::vector<CapturedNodes> capturedNodes;

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
   * @brief Construct a new Lhyphen object
   *
   */
  Lhyphen();

  /**
   * @brief affiche un petit entete sympatique
   *
   */
  void head();

  /**
   * @brief trouver les limites de la zone dessinée dans les fichiers svg
   *
   * @param factor facteur multiplicateur
   */
  void findDisplayArea(double factor = 1.0);

  /**
   * @brief ajoute une cellule de forme polygonale régulière
   *
   * @param h les données pour définir un polyèdre régulier
   * @param p les propriétés mécaniques
   */
  void addRegularPolygonalCell(named_arg_RegularCellDataset h, named_arg_CellProperties p);

  /**
   * @brief ajoute une boite à partir de deux points
   *
   * @param h les deux points qui definissent les dimensions extérieures de la boite
   * @param p les propriétés mécaniques
   */
  void addBoxCell(named_arg_TwoPointsDataset h, named_arg_CellProperties p);

  /**
   * @brief ajoute une ligne (cellule non fermée) entre deux points
   *
   * @param h les deux points
   * @param p les propriétés mécaniques
   */
  void addLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p);
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
                                                named_arg_RegularCellDataset h, named_arg_CellProperties p);

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
                               named_arg_CellProperties p);
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
                         named_arg_CellProperties p);

  /**
   * @brief Utiliser cette méthode pour definir le pas de temps (comme ça dt/2 et dt^2/2 seront pré-calculées)
   *
   * @param dt_ pas de temps
   */
  void setTimeStep(double dt_);

  void setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max, double yieldPower);

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
  void setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue);

  /**
   * @brief ajoute un control à tous les noeuds de la cellule c
   *
   * @param c numéro de la cellule
   * @param xmode VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
   * @param xvalue valeur suivant x (une force ou une vitesse selon xmode)
   * @param ymode VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
   * @param yvalue valeur suivant y (une force ou une vitesse selon ymode)
   */
  void setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue);

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
                           double yvalue);

  /**
   * @brief vérifie si un noeud (cell ci, node in) et une barre (cell cj, barre commancant par le noeud jn) sont
   *        proches. Selon le cas, la paire dans la liste de voisins est soit ajoutée soit retirée
   *
   * @param ci cellule i
   * @param cj cellule j
   * @param in numero du noeud dans la cellule ci
   * @param jn numero du noeud du début d'une barre dans la cellule cj
   * @param epsilonEnds longueur à ignorer aux extremités de la barre
   */
  void addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds = 0.0);

  /**
   * @brief cette fonction permet de lire un fichier contenant une liste de positions x,y avec numéro
   *        de cellule. Peut importe les numéros du moment qu'ils sont différents pour chaque cellule.
   *        Il faut faire un pré-nettoyage pour qu'il y ait au moins 3 noeuds par cellule.
   *        Toutes les cellules sont fermées.
   *
   * @param name nom du fichier
   * @param barWidth épaisseur de toutes les barres
   */
  void readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max);

  /**
   * @brief met à jour la liste des voisins
   *
   */
  void updateNeighbors();

  /**
   * @brief met un point de colle aux endroits où la distance est suffisement proche.
   *
   * @param epsilonDist distance maximale
   */
  void glue(double epsilonDist);

  /**
   * @brief calcul des forces d'interaction (entre cellules)
   *
   */
  void computeInteractionForces();

  /**
   * @brief Calcul des forces internes aux cellules
   *
   */
  void computeNodeForces();

  /**
   * @brief compute the node accelerations
   *
   */
  void NodeAccelerations();

  /**
   * @brief un pas en velocity verlet
   *
   */
  void SingleStep();

  /**
   * @brief run the simulation!
   *
   */
  void run();

  /**
   * @brief Save the current configuration in a file
   *
   * @param fname Name of the file
   */
  void saveCONF(const char *fname);

  void saveCONF(int ifile);

  /**
   * @brief Load a configuration that has been saved in a file
   *
   * @param fname The name of the file
   */
  void loadCONF(const char *fname);

  void loadCONF(int ifile);

  /**
   * @brief fonction très basique pour dessiner les cellules avec des lignes bleues,
   *        TODO: pour d'autre sorties graphiques, on verra après :)
   *
   * @param num numero du fichier
   * @param nameBase nom de base (dans lequel sera inséré le numéro)
   * @param Canvaswidth largeur de canvas (par affichae rapide)
   *
   * @see https://bestofcpp.com/repo/tomkwok-svgasm
   */
  void saveSVG(int num, const char *nameBase = "sample%04d.svg", int Canvaswidth = 400);

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
  double getRotationVelocityBar(vec2r &a, vec2r &b, vec2r &va, vec2r &vb);
};

#endif /* L_HYPHEN_HPP */
