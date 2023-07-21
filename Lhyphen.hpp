// l-hyphen (lodge hyphen) is a multi-purpose system made of cells that are themeself composed of bars
//
// ============================================
// Ce fichier contient tout le code du modèle !
// Pour l'utiliser il suffit de l'inclure

#ifndef L_HYPHEN_HPP
#define L_HYPHEN_HPP

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <utility>
#include <vector>

#include "AABB_2D.hpp"
#include "ColorTable.hpp"
#include "profiler.hpp"
#include "svgtools.hpp"
#include "vec2.hpp"

#include "Bar.hpp"
#include "Cell.hpp"
#include "Control.hpp"
#include "Neighbor.hpp"
#include "Node.hpp"

/**
 *  Cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *  lors de la définition de certaine cellules
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
 *  Cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *  lors de la définition de certaine cellules
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
 *  Cette structure sert simplement à obtenir des arguments només pour une interface plus compréhensible
 *  lors de la définition de certaine cellules
 *
 */
struct named_arg_CellProperties {
  double kn;
  double kr;
  double mz_max;
  double p_int;
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

#define CELL_CONTAIN_NOTHING 0
#define CELL_CONTAIN_GAS 1
#define CELL_CONTAIN_LIQUID 2

/**
 *  Le système entier
 *
 */
class Lhyphen {
public:
  std::vector<Cell> cells;       ///< cellules
  std::vector<Control> controls; ///< contrôles

  std::vector<CapturedNodes> capturedNodes; ///< noeuds capturés

  std::vector<size_t> followedCells; ///< cellules suivies

  // Limites de dessin pour les sorties SVG
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  double t;               ///< temps courant
  double cyclicVelPeriod; ///< durée d'un cycle (inversion de la vitesse de chargement)

  double dt;    ///< pas de temps
  double dt_2;  ///< pas de temps divisé par deux
  double dt2_2; ///< pas de temps au carré divisé par deux

  double globalViscosity;      ///< viscosité globale
  double numericalDissipation; ///< coefficizent de dissipation numérique

  vec2r gravity; ///< gravité

  double distVerlet; ///< distance de Verlet entre noeuds et barres

  // int enablePressures;
  int cellContent;
  double compressFactor;

  // parametres mécaniques d'interactions (contact frottant avec ou sans adhésion) entre les cellules
  double kn;                  ///< raideur normale de contact
  double kt;                  ///< raideur tangentielle de contact
  double mu;                  ///< coefficient de frottement (entre les cellules)
  double fadh;                ///< force normale d'adhésion au contact
  int adaptativeStiffness{0}; ///< adaptativeStiffness

  int nstep;             ///< nombre de pas
  int nstepPeriodVerlet; ///< nombre de pas entre mise à jour des voisins
  int nstepPeriodSVG;    ///< nombre de pas entre sauvegardes SVG
  int nstepPeriodRecord; ///< nombre de pas entre chaque ligne d'enregistrement de données
  int nstepPeriodConf;   ///< nombre de pas entre sauvegardes de configuration
  int isvg;              ///< numéro actuel de sauvegarde SVG
  int iconf;             ///< numéro actuel de sauvegarde de configuration

  int SVG_colorCells; // 0=rien, 1=pressure, 2=NRJ elast
  int SVG_colorTableRescale;
  double SVG_colorTableMin;
  double SVG_colorTableMax;
  int SVG_cellForces; // 0=rien, 1=rouge/bleu

  ColorTable ctNeg, ctPos;

  int reorder{1}; ///< flag pour ré-ordonner les noeud dans la fonction ReadNodeFile

  Lhyphen();

  void head(); ///< affiche un petit entete sympatique

  void addRegularPolygonalCell(named_arg_RegularCellDataset h, named_arg_CellProperties p);
  void addBoxCell(named_arg_TwoPointsDataset h, named_arg_CellProperties p);
  void addLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p);
  void addMultiLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p, int n);
  void addRegularPolygonalCellsOnTriangularGrid(int nx, int ny, double horizontalDistance, double verticalDistance,
                                                named_arg_RegularCellDataset h, named_arg_CellProperties p);
  void addSquareBrickWallCells(int nx, int ny, double horizontalDistance, double xleft, double ybottom, double barWidth,
                               named_arg_CellProperties p);
  void addHoneycombCells(int nx, int ny, double CellExternWidth, double xleft, double ybottom, double barWidth,
                         named_arg_CellProperties p);
  void setTimeStep(double dt_);
  void setCellWallDensities(double rho, double thickness = 1.0);
  void setCellDensities(double rho, double thickness = 1.0);
  void setCellMasses(double m);
  void setNodeMasses(double m);
  void setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max, double yieldPower);
  void setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue);
  void setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue);
  void setNodeControlInBox(double xmin, double xmax, double ymin, double ymax, int xmode, double xvalue, int ymode,
                           double yvalue);
  void addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds = 0.0);
  void readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max, double p_int);

  void updateNeighbors();
  void glue(double epsilonDist);
  void computeInteractionForces();
  void computeNodeForces();
  void NodeAccelerations();
  void SingleStep();
  void integrate();

  void saveCONF(const char *fname);
  void saveCONF(int ifile);
  void loadCONF(const char *fname);
  void loadCONF(int ifile);
  void initCellInitialSurfaces();

  void findDisplayArea(double factor = 1.0);
  void saveSVG(int num, const char *nameBase = "sample%04d.svg", int Canvaswidth = 400);
  void InternalGasPressureForce();
  void InternalLiquidPressureForce();
  void setCellInternalPressure(size_t c, double p);
};

#endif /* L_HYPHEN_HPP */
