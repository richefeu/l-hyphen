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
  
  int enablePressures;

  // parametres mécaniques d'interactions (contact frottant avec ou sans adhésion) entre les cellules
  double kn;   // raideur normale de contact
  double kt;   // raideur tangentielle de contact
  double mu;   // coefficient de frottement (entre les cellules)
  double fadh; // force normale d'adhésion au contact

  int nstep;
  int nstepPeriodVerlet;
  int nstepPeriodSVG;
  int nstepPeriodRecord;
  int nstepPeriodConf;
  int isvg;
  int iconf;

  Lhyphen();

  void head();

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
  void setCellDensities(double rho);
  void setCellMasses(double m);
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

  void findDisplayArea(double factor = 1.0);
  void saveSVG(int num, const char *nameBase = "sample%04d.svg", int Canvaswidth = 400);
  void InternalPressureForce();
  void setCellInternalPressure(size_t c, double p);

private:
  double getRotationVelocityBar(vec2r &a, vec2r &b, vec2r &va, vec2r &vb);
};

#endif /* L_HYPHEN_HPP */
