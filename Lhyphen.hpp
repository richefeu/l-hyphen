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

#include <fstream>
#include <iostream>

#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "AABB_2D.hpp"
#include "ColorTable.hpp"
#include "linkCells2D.hpp"
#include "profiler.hpp"
#include "svgtools.hpp"
#include "vec2.hpp"

#include "Bar.hpp"
#include "Cell.hpp"
#include "Control.hpp"
#include "Neighbor.hpp"
#include "Node.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define LH_ENABLED 1
#define LH_DISABLED 0

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
  CellNodeID(size_t t_c, size_t t_n) : c(t_c), n(t_n) {}
};

/**
 *  Une boite qui sert à capturer des sommets pour faire des sorties de post-traitement
 */
struct CapturedNodes {
  std::vector<CellNodeID> cellNodeIDs;
  std::string filename;
  double xmin, xmax, ymin, ymax;
  CapturedNodes(const char *name, double t_xmin, double t_xmax, double t_ymin, double t_ymax)
      : filename(name), xmin(t_xmin), xmax(t_xmax), ymin(t_ymin), ymax(t_ymax) {}
};

/**
 *  Une boite qui sert à capturer des sommets pour imposer une force ou une vitesse suivant x et y
 */
struct ControlBoxArea {
  double xmin, xmax, ymin, ymax;
  int xmode;
  double xvalue;
  int ymode;
  double yvalue;
  ControlBoxArea(double t_xmin, double t_xmax, double t_ymin, double t_ymax, int t_xmode, double t_xvalue, int t_ymode,
                 double t_yvalue)
      : xmin(t_xmin), xmax(t_xmax), ymin(t_ymin), ymax(t_ymax), xmode(t_xmode), xvalue(t_xvalue), ymode(t_ymode),
        yvalue(t_yvalue) {}
};

#define CELL_EMPTY 0
#define CELL_CONSTANT_PV 1
#define CELL_ELASTIC_PV 2

/**
 *  Le système entier
 *
 */
class Lhyphen {
public:
  std::vector<Cell> cells;       ///< cellules
  std::vector<Control> controls; ///< contrôles

  std::vector<CapturedNodes> capturedNodes; ///< noeuds capturés
  std::vector<ControlBoxArea> controlBoxAreas;

  std::vector<size_t> followedCells; ///< cellules suivies

  int nbThreads{1};

  // Limites de dessin pour les sorties SVG
  double xmin{0.0};
  double xmax{0.0};
  double ymin{0.0};
  double ymax{0.0};

  double t{0.0};               ///< temps courant
  double cyclicVelPeriod{0.0}; ///< durée d'un cycle (inversion de la vitesse de chargement)

  double dt{0.0};    ///< pas de temps
  double dt_2{0.0};  ///< pas de temps divisé par deux
  double dt2_2{0.0}; ///< pas de temps au carré divisé par deux

  double globalViscosity{0.0};      ///< viscosité globale
  double numericalDissipation{0.0}; ///< coefficizent de dissipation numérique

  vec2r gravity; ///< gravité

  double distVerlet{0.0}; ///< distance de Verlet entre noeuds et barres

  double linkCells_lx{0.0};
  double linkCells_ly{0.0};

  int cellContent;
  double compressFactor{0.0};

  // parametres mécaniques d'interactions (contact frottant avec ou sans adhésion) entre les cellules
  double kn{0.0};             ///< raideur normale de contact
  double kt{0.0};             ///< raideur tangentielle de contact
  double mu{0.0};             ///< coefficient de frottement (entre les cellules)
  double fadh{0.0};           ///< force normale d'adhésion au contact
  int adaptativeStiffness{0}; ///< adaptativeStiffness

  int nstep{0};             ///< nombre de pas
  int nstepPeriodVerlet{0}; ///< nombre de pas entre mise à jour des voisins
  int nstepPeriodSVG{0};    ///< nombre de pas entre sauvegardes SVG
  int nstepPeriodRecord{0}; ///< nombre de pas entre chaque ligne d'enregistrement de données
  int nstepPeriodConf{0};   ///< nombre de pas entre sauvegardes de configuration
  int isvg{0};              ///< numéro actuel de sauvegarde SVG
  int iconf{0};             ///< numéro actuel de sauvegarde de configuration

  int SVG_colorCells{0}; // 0=rien, 1=pressure, 2=NRJ elast
  int SVG_colorTableRescale{0};
  double SVG_colorTableMin{0.0};
  double SVG_colorTableMax{0.0};
  int SVG_cellForces{0}; // 0=rien, 1=rouge/bleu

  ColorTable ctNeg, ctPos; // ??????? c'est pas le bon endroit pour mettre ça

  int reorder{LH_ENABLED}; ///< flag pour ré-ordonner les noeud dans la fonction ReadNodeFile

  std::ofstream breakHistory;
  std::ofstream breakEvol;
  double cumulatedG{0.0};
  double cumulatedL{0.0};

  Lhyphen();

  void head(); ///< affiche un petit entete sympatique

  // pre-processing functions
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

  void setTimeStep(double t_dt);
  void setCellWallDampingRates(double alpha);
  void setCellWallDensities(double rho, double thickness = 1.0);
  void setCellDensities(double rho, double thickness = 1.0);
  void setCellMasses(double cellMass);
  void setNodeMasses(double nodeMass);
  void setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max, double yieldPower);
  void setGcGlueSameProperties(double kn_coh, double kt_coh, double Gc);
  void setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue);
  void setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue);
  void setNodeControlInBox(double t_xmin, double t_xmax, double t_ymin, double t_ymax, int xmode, double xvalue,
                           int ymode, double yvalue);
  void addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds = 0.0);
  double getMinimumNodeDistance();
  void readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max, double p_int);

  std::function<void()> updateNeighbors;
  void updateNeighbors_brute_force();
  void updateNeighbors_linkCells();
  void copy_neighbors_set_to_vec();

  int getPosition(size_t ci, size_t cj, size_t in, size_t jn, vec2r &pos);
  void glue(double epsilonDist, int modelGc = 1);
  void associateGlue(int modelGc = 1, double activationLength = 0.0);
  void glue_breakage(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jn, size_t jnext, double wbeg,
                     double wend);
  void computeInteractionForces();
  void computeNodeForces();
  void nodeAccelerations();
  void SingleStep();
  void integrate();

  void saveCONF(const char *fname);
  void saveCONF(int ifile);
  void loadCONF(const char *fname);
  void loadCONF(int ifile);
  void computeCellInitialSurfaces();

  void findDisplayArea(double factor = 1.0);
  void saveSVG(int num, const char *nameBase = "sample%04d.svg", int Canvaswidth = 400);
  void InternalGasPressureForce();
  void InternalLiquidPressureForce();
  void setCellInternalPressure(size_t c, double p);
};
