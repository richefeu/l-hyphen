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

#include "Lhyphen.hpp"
#include "null_size_t.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

/// Construct a new Lhyphen object
///
Lhyphen::Lhyphen()
    : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), t(0.0), cyclicVelPeriod(0.0), dt(0.0), globalViscosity(0.0),
      numericalDissipation(0.0), gravity(), distVerlet(0.0), cellContent(CELL_EMPTY), compressFactor(0.0), kn(1000.0),
      kt(1000.0), mu(0.0), fadh(0.0), nstep(1000), nstepPeriodVerlet(1), nstepPeriodSVG(10), nstepPeriodRecord(1),
      nstepPeriodConf(10), isvg(0), iconf(0) {

  SVG_colorCells = 2;
  SVG_colorTableRescale = 0;
  SVG_colorTableMin = -10000.0;
  SVG_colorTableMax = 10000.0;

  SVG_cellForces = 1;

  updateNeighbors = [this]() { this->updateNeighbors_brute_force(); };

  linkCells_lx = 1.0;
  linkCells_ly = 1.0;

  ctNeg.setSize(128);
  ctPos.setSize(128);
  ctNeg.rebuild_interp_rgba({0, 127}, {{204, 204, 230, 255}, {0, 0, 255, 255}});
  ctPos.rebuild_interp_rgba({0, 127}, {{204, 204, 230, 255}, {255, 0, 0, 255}});
}

Lhyphen::~Lhyphen() {
  if (exprParser != nullptr) {
    delete exprParser;
  }
}

///  Affiche d'un petit entete sympatique
///
void Lhyphen::head() {
  std::cout << "    L-HYPHEN\n";
  std::cout << "    _    _\n";
  std::cout << "   / \\__/ \\\n";
  std::cout << "  |  o  <  |\n";
  std::cout << "   \\_\\ _/_/\n";
  std::cout << "     /_/\n";
  std::cout << "    /_/\n" << std::endl;
}

/// Affiche un diagnostic des paramètres de simulation avant l'intégration.
///
/// Calcule et affiche les rapports dt_critique/dt pour les contacts et les liens collés,
/// ainsi qu'une analyse des paramètres de la recherche de voisins de Verlet.
/// Sauvegarde le rapport dans un fichier "diagnostic.txt".
///
void Lhyphen::diagnostics() {
  const double INF = std::numeric_limits<double>::max();

  // ---- Comptages ----
  size_t totalNodes = 0, totalBars = 0, totalNeighbors = 0;
  for (auto &c : cells) {
    totalNodes += c.nodes.size();
    totalBars += c.bars.size();
    totalNeighbors += c.neighbors.size();
  }

  // ---- Masses ----
  double mMin = INF, mMax = 0.0;
  for (auto &c : cells)
    for (auto &n : c.nodes)
      if (n.mass > 0.0) {
        mMin = std::min(mMin, n.mass);
        mMax = std::max(mMax, n.mass);
      }

  // ---- Raideurs des barres ----
  double knBarMin = INF, knBarMax = 0.0;
  for (auto &c : cells)
    for (auto &b : c.bars)
      if (b.kn > 0.0) {
        knBarMin = std::min(knBarMin, b.kn);
        knBarMax = std::max(knBarMax, b.kn);
      }

  // ---- Rayons des cellules ----
  // ***** VR: c'est pas bon ce truc !!!!!!!!! <<<<<<<<<<<<<
  double rMin = INF, rMax = 0.0, rMean = 0.0;
  size_t nRad = 0;
  for (auto &c : cells)
    if (c.radius > 0.0) {
      rMin = std::min(rMin, c.radius);
      rMax = std::max(rMax, c.radius);
      rMean += c.radius;
      ++nRad;
    }
  if (nRad > 0)
    rMean /= (double)nRad;

  // ---- dt critique contact : dt_crit = sqrt(2 * m_min / kn_contact) ----
  // Cas le plus défavorable : deux noeuds de masse minimale en contact.
  double dtCritContact = INF;
  if (kn > 0.0 && mMin < INF)
    dtCritContact = std::sqrt(2.0 * mMin / kn);

  // ---- dt critique cohésion ----
  double knCohMax = 0.0;
  for (auto &c : cells)
    for (auto &inter : c.neighbors)
      knCohMax = std::max(knCohMax, inter.kn_coh);
  double dtCritCoh = INF;
  if (knCohMax > 0.0 && mMin < INF)
    dtCritCoh = std::sqrt(2.0 * mMin / knCohMax);

  // ---- Longueur moyenne des barres ----
  double lBarMean = 0.0;
  size_t nBarCount = 0;
  for (auto &c : cells)
    for (auto &b : c.bars)
      if (b.l0 > 0.0) {
        lBarMean += b.l0;
        ++nBarCount;
      }
  if (nBarCount > 0)
    lBarMean /= (double)nBarCount;

  // ====================================================================
  // Préparation du texte de diagnostic
  // ====================================================================

  std::ostringstream oss;

  // Affiche dt_crit/dt avec qualificatif :
  //   >= 50  : dt peut etre augmente
  //   >= 20  : OK
  //   >= 10  : proche de la limite
  //   >= 1   : tres proche de la limite (instabilite imminente)
  //   <  1   : INSTABLE
  auto formatDtLine = [&](const char *label, double dtCrit) -> std::string {
    if (dtCrit >= INF * 0.5)
      return "";
    double ratio = (dt > 0.0) ? dtCrit / dt : 0.0;
    std::ostringstream line;
    line << "    " << std::left << std::setw(10) << label << " : dt_crit = " << std::scientific << std::setprecision(3)
         << dtCrit << "   dt_crit/dt = " << std::fixed << std::setprecision(1) << std::right << std::setw(8) << ratio;
    if (ratio < 1.0)
      line << "  <<< INSTABLE !!";
    else if (ratio < 10.0)
      line << "  << tres proche de la limite";
    else if (ratio < 20.0)
      line << "  < proche de la limite";
    else if (ratio >= 50.0)
      line << "     (dt peut etre augmente)";
    return line.str();
  };

  oss << "\n";
  oss << "╔════════════════════════════════════════════════════════════════════════════╗\n";
  oss << "║                           DIAGNOSTIC SIMULATION                            ║\n";
  oss << "╚════════════════════════════════════════════════════════════════════════════╝\n";
  oss << "\n";
  oss << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  GEOMETRIE ET PROPRIETES\n";
  oss << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  Cellules           : " << cells.size() << "\n";
  oss << "  Noeuds             : " << totalNodes << "\n";
  oss << "  Barres             : " << totalBars << "\n";
  oss << "  Paires de voisins  : " << totalNeighbors / 2 << "\n";

  if (mMin < INF) {
    oss << "\n  Masses des noeuds  : min = " << std::scientific << std::setprecision(3) << mMin << "   max = " << mMax
        << "\n";
  }

  if (knBarMin < INF) {
    oss << "  Raideur des barres : min = " << std::scientific << std::setprecision(3) << knBarMin
        << "   max = " << knBarMax << "\n";
  }

  if (lBarMean > 0.0) {
    oss << "\n  Longueur de barre (moyenne) : " << std::fixed << std::setprecision(6) << lBarMean << "\n";
  }

  if (rMax > 0.0) {
    oss << "  Rayon des barres            : min = " << std::fixed << std::setprecision(6) << rMin << "   max = " << rMax
        << "\n";
  }

  // --- Pas de temps ---
  oss << "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  PAS DE TEMPS\n";
  oss << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  dt                   : " << std::scientific << std::setprecision(3) << dt << " s\n";
  oss << "  Nombre de pas        : " << nstep << "\n";
  oss << "  Durée totale         : " << std::fixed << std::setprecision(6) << (double)nstep * dt << " s\n";

  oss << "\n  Stabilité (dt_crit/dt >= 20 recommandé, < 10 critique) :\n";
  std::string contactLine = formatDtLine("contact", dtCritContact);
  std::string cohesionLine = formatDtLine("cohesion", dtCritCoh);
  if (!contactLine.empty())
    oss << contactLine << "\n";
  if (!cohesionLine.empty())
    oss << cohesionLine << "\n";

  // --- Sorties ---
  oss << "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  FICHIERS DE SORTIE\n";
  oss << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  if (nstepPeriodConf > 0) {
    oss << "  Configurations (conf) : tous les " << std::setw(6) << nstepPeriodConf << " pas  ->  " << std::setw(5)
        << (nstep / nstepPeriodConf) << " fichiers"
        << "  (dt_phys = " << std::scientific << std::setprecision(3) << nstepPeriodConf * dt << " s)\n";
  } else {
    oss << "  Configurations (conf) : désactivées\n";
  }

  if (nstepPeriodSVG > 0) {
    oss << "  Visualisations (SVG)  : tous les " << std::setw(6) << nstepPeriodSVG << " pas  ->  " << std::setw(5)
        << (nstep / nstepPeriodSVG) << " fichiers"
        << "  (dt_phys = " << std::scientific << std::setprecision(3) << nstepPeriodSVG * dt << " s)\n";
  } else {
    oss << "  Visualisations (SVG)  : désactivées\n";
  }

  // --- Voisins Verlet ---
  oss << "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  ALGORITHME DE RECHERCHE DE VOISINS\n";
  oss << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
  oss << "  Principe : Regroupe toutes les paires dans une zone élargie.\n";

  double Trecon = (dt > 0.0 && nstepPeriodVerlet > 0) ? (double)nstepPeriodVerlet * dt : 0.0;

  oss << "\n  distVerlet                  : " << std::fixed << std::setprecision(6) << distVerlet << "\n";
  oss << "  nstepPeriodVerlet           : " << nstepPeriodVerlet << "\n";

  if (Trecon > 0.0) {
    oss << "  temps entre mises à jour : " << std::scientific << std::setprecision(3) << Trecon << " s\n";
    oss << "\n  Budget de déplacement relatif entre mises à jour :\n";
    oss << "    v_rel_max = " << std::scientific << std::setprecision(3) << distVerlet / Trecon << " m/s\n";
  }

  // Checks pertinents uniquement en présence de contacts inter-cellules
  if (cells.size() > 1 && kn > 0.0) {
    oss << "\n  Algorithme de recherche : ";
    if (linkCells_lx > 0.0) {
      double minSz = 2.0 * rMax + distVerlet;
      oss << "link-cells\n";
      // oss << "    lx = " << std::fixed << std::setprecision(6) << linkCells_lx
      //     << "   ly = " << linkCells_ly << "\n";
      oss << "    Condition de validité : lx (et ly) >= 2 * r_max + distVerlet\n";
      oss << "      2 * r_max + distVerlet = 2 * " << std::fixed << std::setprecision(6) << rMax << " + " << distVerlet
          << " = " << minSz << "\n";
      bool lxOK = (linkCells_lx >= minSz);
      bool lyOK = (linkCells_ly >= minSz);
      oss << "      lx = " << linkCells_lx << (lxOK ? "  ✓ OK" : "  ✗ ATTENTION: paires possiblement manquées") << "\n";
      oss << "      ly = " << linkCells_ly << (lyOK ? "  ✓ OK" : "  ✗ ATTENTION: paires possiblement manquées") << "\n";
    } else {
      oss << "recherche brute O(N²)\n";
      if (totalNodes > 500) {
        oss << "    (" << totalNodes << " noeuds : envisager linkCells pour accélérer)\n";
      }
    }
  }

  std::string diagnosticText = oss.str();

  // Affichage sur la console
  std::cout << diagnosticText;

  // Sauvegarde dans un fichier
  std::ofstream diagFile("diagnostic.txt");
  if (diagFile.is_open()) {
    diagFile << diagnosticText;
    diagFile.close();
  }
}

// ======================================================================================================
// Adding cells (or packing of cells)
// ======================================================================================================

/// Ajoute une cellule de forme polygonale régulière
///
/// @param h les données pour définir un polyèdre régulier
/// @param p les propriétés mécaniques
///
void Lhyphen::addRegularPolygonalCell(named_arg_RegularCellDataset h, named_arg_CellProperties p) {
  Cell C;
  C.p_int = p.p_int;
  C.close = true;
  double a0 = 2.0 * M_PI / (double)(h.nbFaces);
  double Rcell = h.Rext - 0.5 * h.barWidth;
  for (double a = 0.0; a < 2.0 * M_PI - 0.99 * a0; a += a0) {
    double theta = a0 + h.rot + a;
    C.nodes.emplace_back(h.x + Rcell * cos(theta), h.y + Rcell * sin(theta));
  }
  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, p.p_int, true);
  cells.push_back(C);
}

/// Ajoute une cellule-boite à partir de deux points
///
/// @param h  les deux points qui definissent les dimensions extérieures de la boite
/// @param p  les propriétés mécaniques
///
void Lhyphen::addBoxCell(named_arg_TwoPointsDataset h, named_arg_CellProperties p) {
  Cell C;
  C.p_int = p.p_int;
  C.close = true;
  double r = 0.5 * h.barWidth;
  C.nodes.emplace_back(h.xo + r, h.yo + r);
  C.nodes.emplace_back(h.xe - r, h.yo + r);
  C.nodes.emplace_back(h.xe - r, h.ye - r);
  C.nodes.emplace_back(h.xo + r, h.ye - r);
  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, p.p_int, true);
  cells.push_back(C);
}

///  Ajoute une ligne (cellule non fermée) entre deux points
///
///  @param h  les deux points
///  @param p  les propriétés mécaniques
///
void Lhyphen::addLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p) {
  Cell C;
  C.close = false;
  C.nodes.emplace_back(h.xo, h.yo);
  C.nodes.emplace_back(h.xe, h.ye);
  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, 0.0, false);
  cells.push_back(C);
}

/// Ajoute une ligne formée de plusieurs barres (cellule non fermée) entre deux points
///
/// @param h  les deux points
/// @param p  les propriétés mécaniques
/// @param n  nombre de barres
///
void Lhyphen::addMultiLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p, int n) {
  double dx = (h.xe - h.xo) / (double)n;
  double dy = (h.ye - h.yo) / (double)n;

  Cell C;
  C.close = false;
  C.nodes.emplace_back(h.xo, h.yo);

  for (int i = 1; i <= n; i++) {
    C.nodes.emplace_back(h.xo + i * dx, h.yo + i * dy);
  }

  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, 0.0, false);

  cells.push_back(C);
}

///  Adds regular polygonal cells on a triangular grid.
///
///  @param nx                  number of cells in the x-direction
///  @param ny                  number of cells in the y-direction
///  @param horizontalDistance  horizontal distance between cells
///  @param verticalDistance    vertical distance between cells
///  @param h                   regular cell dataset
///  @param p                   cell properties
///
void Lhyphen::addRegularPolygonalCellsOnTriangularGrid(int nx, int ny, double horizontalDistance,
                                                       double verticalDistance, named_arg_RegularCellDataset h,
                                                       named_arg_CellProperties p) {

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
    if (odd == 1) {
      xshift = 0.5 * horizontalDistance;
    } else {
      xshift = 0.0;
    }
  }
}

/// Adds square brick wall cells to the given grid.
///
/// @param nx                  number of cells in the x direction
/// @param ny                  number of cells in the y direction
/// @param horizontalDistance  horizontal distance between cells
/// @param xleft               x-coordinate of the leftmost cell
/// @param ybottom             y-coordinate of the bottommost cell
/// @param barWidth            width of the bar between cells
/// @param p                   properties of the cells
///
void Lhyphen::addSquareBrickWallCells(int nx, int ny, double horizontalDistance, double xleft, double ybottom,
                                      double barWidth, named_arg_CellProperties p) {
  named_arg_RegularCellDataset h;
  h.nbFaces = 4;
  h.x = xleft - horizontalDistance;
  h.y = ybottom - horizontalDistance;
  h.rot = M_PI / 4.0;
  h.Rext = 0.5 * (horizontalDistance - barWidth) * sqrt(2.0) + 0.5 * barWidth;
  h.barWidth = barWidth;
  addRegularPolygonalCellsOnTriangularGrid(nx, ny, horizontalDistance, horizontalDistance, h, p);
}

///  Crée une structure en nid d'abeille
///
///  @param nx               nombre de cellule suivant x
///  @param ny               nombre de lignes
///  @param CellExternWidth  largeur externe d'une cellule hexagonale
///  @param xleft            position x la plus à gauche
///  @param ybottom          position y la plus en bas
///  @param barWidth         épaisseur des barres
///  @param p                paramètres mécaniques de toutes les cellules
///
void Lhyphen::addHoneycombCells(int nx, int ny, double CellExternWidth, double xleft, double ybottom, double barWidth,
                                named_arg_CellProperties p) {
  double verticalDistance = 1.5 * CellExternWidth / sqrt(3.0);
  named_arg_RegularCellDataset h;
  h.nbFaces = 6;
  h.x = xleft + 0.5 * CellExternWidth + 0.5 * barWidth;
  h.y = ybottom + 0.5 * verticalDistance + 0.5 * barWidth;
  h.rot = M_PI / 6.0;
  h.Rext = (CellExternWidth - barWidth) / sqrt(3.0) + 0.5 * barWidth;
  h.barWidth = barWidth;
  addRegularPolygonalCellsOnTriangularGrid(nx, ny, CellExternWidth, verticalDistance, h, p);
}

// ======================================================================================================
// Reading preprocessed file with nodes
// ======================================================================================================

double Lhyphen::getMinimumNodeDistance() {
  struct cellPair {
    size_t i;
    size_t j;
  };

  // Compute the expanded-AABB of each cell
  std::vector<AABB_2D> aabbs(cells.size());
  for (size_t ci = 0; ci < cells.size(); ci++) {
    aabbs[ci].set_single(cells[ci].nodes[0].pos);
    for (size_t in = 1; in < cells[ci].nodes.size(); in++) {
      aabbs[ci].add(cells[ci].nodes[in].pos);
    }
    aabbs[ci].enlarge(cells[ci].radius + 0.5 * distVerlet);
  }

  // Find intersecting expanded-AABB
  std::vector<cellPair> cellPairs;
  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (size_t cj = ci + 1; cj < cells.size(); cj++) {

      if (aabbs[ci].intersect(aabbs[cj])) {
        cellPair P;

        P.i = ci;
        P.j = cj;
        cellPairs.push_back(P);
      }
    }
  }

  double d2min = 1e20;
  for (size_t ip = 0; ip < cellPairs.size(); ++ip) {

    size_t ci = cellPairs[ip].i;
    size_t cj = cellPairs[ip].j;

    // test chaque sommet de j avec chaque barre de i
    for (size_t jn = 0; jn < cells[cj].nodes.size(); ++jn) {
      for (size_t in = 0; in < cells[ci].nodes.size(); ++in) {
        double d2 = norm2(cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos);
        if (d2 < d2min) {
          d2min = d2;
        }
      }
    }
  }

  return sqrt(d2min);
}

///  This function reads a file containing a list of x,y positions with cell numbers.
///  The numbers don't matter as long as they are different for each cell.
///  Pre-cleaning must be done so that there are at least 3 nodes per cell. All cells are closed.
///
///  @param name      name of the node file to read
///  @param barWidth  width of the bar (negative for auto-evaluated value)
///  @param Kn        value of Kn
///  @param Kr        value of Kr
///  @param Mz_max    maximum value of Mz
///  @param p_int     value of p_int
///
void Lhyphen::readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max, double p_int) {
  std::ifstream file(name);

  struct data {
    double x, y;
    size_t id;
    bool operator<(const data &rhs) const {
      if (id < rhs.id) {
        return true;
      }
      return false;
    }
  };

  std::multiset<data> dataset;

  data D;
  while (file.good()) {
    file >> D.x >> D.y >> D.id;
    if (file.eof()) {
      break;
    }
    dataset.insert(D);
  }

  size_t i0 = cells.size(); // in case some cells have already been added
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
  if (!C.nodes.empty()) {
    cells.push_back(C);
  }

  // If the barWidth is given as a negative value, it is defined as half the minimum distance between nodes of different
  // cells
  if (barWidth < 0.0) {
    barWidth = getMinimumNodeDistance();
    std::cout << "@Lhyphen::readNodeFile, barWidth has been set to " << barWidth << std::endl;
  }

  for (size_t i = i0; i < cells.size(); i++) {
    cells[i].reorderNodes();
    if (reorder == 1) {
      cells[i].connectOrderedNodes(barWidth, Kn, Kr, Mz_max, p_int, true);
    }
  }
}

// ======================================================================================================
// Setting parameters
// ======================================================================================================

///  Utiliser cette méthode pour definir le pas de temps (comme ça dt/2 et dt^2/2 seront pré-calculées)
///
///  @param t_dt  time step
///
void Lhyphen::setTimeStep(double t_dt) {
  dt = t_dt;
  dt_2 = 0.5 * dt;
  dt2_2 = dt_2 * dt;
}

// ne pas oublier de donner des masses et raideurs aux noeuds avant !!!!!
void Lhyphen::setCellWallDampingRates(double alpha_s, double alpha_b) {

  // Précalcul Ieff moyen
  double L = 0.0;
  double m = 0.0;
  size_t nbBars = 0;
  size_t nbNodes = 0;
  for (size_t c = 0; c < cells.size(); c++) {

    for (size_t b = 0; b < cells[c].bars.size(); b++) {
      L += cells[c].bars[b].l0;
      nbBars++;
    }

    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      m += cells[c].nodes[n].mass;
      nbNodes++;
    }
  }
  if (nbBars > 0) {
    L /= static_cast<double>(nbBars);
  }
  if (nbNodes > 0) {
    m /= static_cast<double>(nbNodes);
  }
  double Ieff = (m * L * L) / 3.0;

  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t b = 0; b < cells[c].bars.size(); b++) {

      size_t in = cells[c].bars[b].i;
      size_t jn = cells[c].bars[b].j;
      double imass = cells[c].nodes[in].mass;
      double jmass = cells[c].nodes[jn].mass;
      double meff = imass * jmass / (imass + jmass);

      cells[c].bars[b].visc = 2.0 * alpha_s * sqrt(meff * cells[c].bars[b].kn);
    }

    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].visc = 2.0 * alpha_b * sqrt(Ieff * cells[c].nodes[n].kr);
    }
  }
}

void Lhyphen::setCellWallDampings(double nu_s, double nu_b) {
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t b = 0; b < cells[c].bars.size(); b++) {
      cells[c].bars[b].visc = nu_s;
    }

    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].visc = nu_b;
    }
  }
}

/// Sets the cell wall densities for all cells in the Lhyphen.
///
/// @param rho The density of the cell walls.
/// @param thickness The thickness of the cell walls.
///
/// @throws None.
///
void Lhyphen::setCellWallDensities(double rho, double thickness) {
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass = 0.0;
    }
  }

  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t b = 0; b < cells[c].bars.size(); b++) {
      double barMass = rho * (cells[c].bars[b].l0 * 2.0 * cells[c].radius * thickness);
      cells[c].nodes[cells[c].bars[b].i].mass += 0.5 * barMass;
      cells[c].nodes[cells[c].bars[b].j].mass += 0.5 * barMass;
    }
  }
}

///  Sets the cell densities for the Lhyphen class.
///
///  @param rho        density of the cells
///  @param thickness  thickness of the cells
///
void Lhyphen::setCellDensities(double rho, double thickness) {
  setCellWallDensities(rho, thickness);
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass *=
          0.5; // parce que ce sera comptabilisé dans la surface (qui chevauche la barre de moitié)
    }
    cells[c].CellSurface();
    double massNodeInc = rho * (cells[c].surface * thickness);
    massNodeInc /= (double)cells[c].nodes.size();
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass += massNodeInc;
    }
  }
}

///  Sets the masses of all cells in the Lhyphen object. Mass of each node is cellMass / nbNodes
///
///  @param cellMass  mass value to be set for each cell
///
void Lhyphen::setCellMasses(double cellMass) {
  for (size_t c = 0; c < cells.size(); c++) {
    double nodeMass = cellMass / (double)cells[c].nodes.size();
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass = nodeMass;
    }
  }
}

///  Set the same mass to all nodes in the Lhyphen object.
///
///  @param nodeMass  mass to set for all nodes
///
void Lhyphen::setNodeMasses(double nodeMass) {
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass = nodeMass;
    }
  }
}

///  Set the same 'force' glue properties for all interactions.
///
///  @param kn_coh      cohesion normal stiffness
///  @param kt_coh      cohesion tangential stiffness
///  @param fn_coh_max  maximum cohesion normal force
///  @param ft_coh_max  maximum cohesion tangential force
///  @param yieldPower  yield power
///
void Lhyphen::setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max,
                                    double yieldPower) {

  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {

      Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));
      if (Inter->glueState == GLUE_FORCE_THRESHOLD) {
        Inter->kn_coh = kn_coh;
        Inter->kt_coh = kt_coh;
        Inter->fn_coh_max = fn_coh_max;
        Inter->ft_coh_max = ft_coh_max;
        Inter->yieldPower = yieldPower;
      }
    }
  }
}

///  Set the same 'Gc' glue properties for all interactions.
///
///  @param kn_coh      cohesion normal stiffness
///  @param kt_coh      cohesion tangential stiffness
///  @param Gc          fracture energy (by unit length)
///
void Lhyphen::setGcGlueSameProperties(double kn_coh, double kt_coh, double Gc) {
  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {

      Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));
      if (Inter->glueState == GLUE_GC) {
        Inter->kn_coh = kn_coh;
        Inter->kt_coh = kt_coh;
        Inter->Gc = Gc;
      }
    }
  }
}

///  Ajoute un control au noeud n de la cellule c
///
///  @param c       numéro de la cellule
///  @param n       numéro du noeud au sein de la cellule c
///  @param xmode   VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
///  @param xvalue  valeur suivant x (une force ou une vitesse selon xmode)
///  @param ymode   VELOCITY_CONTROL (1) ou FORCE_CONTROL (0)
///  @param yvalue  valeur suivant y (une force ou une vitesse selon ymode)
///
void Lhyphen::setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue) {
  Control C(xmode, xvalue, ymode, yvalue);
  controls.push_back(C);
  cells[c].nodes[n].ictrl = controls.size() - 1;
}

///  Ajoute un control à tous les noeuds de la cellule c
///
///  @param c       numéro de la cellule
///  @param xmode   VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
///  @param xvalue  valeur suivant x (une force ou une vitesse selon xmode)
///  @param ymode   VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
///  @param yvalue  valeur suivant y (une force ou une vitesse selon ymode)
///
void Lhyphen::setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue) {
  Control C(xmode, xvalue, ymode, yvalue);
  controls.push_back(C);
  size_t ictrl = controls.size() - 1;
  for (size_t n = 0; n < cells[c].nodes.size(); n++) {
    cells[c].nodes[n].ictrl = ictrl;
  }
}

///  Définir un même control à tous les noeuds qui se trouvent dans une zone rectangulaire
///
///  @param t_xmin  left limit of the box
///  @param t_xmax  right limit of the box
///  @param t_ymin  bottom limit of the box
///  @param t_ymax  top limit of the box
///  @param xmode   imposed mode in the x-direction (VELOCITY_CONTROL or FORCE_CONTROL)
///  @param xvalue  imposed value in the x-direction
///  @param ymode   imposed mode in the y-direction (VELOCITY_CONTROL or FORCE_CONTROL)
///  @param yvalue  imposed value in the y-direction
///
void Lhyphen::setNodeControlInBox(double t_xmin, double t_xmax, double t_ymin, double t_ymax, int xmode, double xvalue,
                                  int ymode, double yvalue) {
  Control C(xmode, xvalue, ymode, yvalue);

  controls.push_back(C);
  size_t ictrl = controls.size() - 1;
  int nbAdded = 0;
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      vec2r pos = cells[c].nodes[n].pos;
      if (pos.x >= t_xmin && pos.x <= t_xmax && pos.y >= t_ymin && pos.y <= t_ymax) {
        cells[c].nodes[n].ictrl = ictrl;
        nbAdded++;
        // std::cout << "@Lhyphen::setNodeControlInBox, node pos = " << pos << "\n";
      }
    }
  }
  std::cout << "@Lhyphen::setNodeControlInBox, nb Added nodes = " << nbAdded << std::endl;
}

///  Sets the internal pressure of a cell.
///
///  @param c  number of the cell
///  @param p  new internal pressure of the cell
///
void Lhyphen::setCellInternalPressure(size_t c, double p) { cells[c].p_int = p; }

// ======================================================================================================
// Updating the neighbor list
// ======================================================================================================

///  Function to copy the neighbors set of each cell into a vector of neighbor pointers
///
///
void Lhyphen::copy_neighbors_set_to_vec() {
  for (auto &cell : cells) {
    size_t nbv = cell.neighbors.size();
    cell.vec_neighbors.clear();
    cell.vec_neighbors.resize(nbv);

    size_t index = 0;
    for (auto &neighbor : cell.neighbors) {
      cell.vec_neighbors[index++] = (Neighbor *)&neighbor;
    }
  }
}

///  Checks whether a node (cell ci, node in) and a bar (cell cj, bar starting with node jn) are close.
///  Depending on the case, the pair in the list of neighbours is either added or removed.
///
///  @param ci           cell i
///  @param cj           cell j
///  @param in           node number in cell ci
///  @param jn           node number at the start of a bar in cell cj
///  @param epsilonEnds  length to be ignored at the ends of the bar
///
void Lhyphen::addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds) {
  // jn is seen here as the index of a bar (start node in fact)

  if (cells[cj].nodes[jn].nextNode == null_size_t) { // case without bar (bar-end)

    vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
    double sqrDist = norm2(branch);
    double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
    bool isNear = sqrDist < sumR * sumR;
    cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

  } else { // the bar linking jn to jnext exists

    size_t jnext = cells[cj].nodes[jn].nextNode;
    vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
    vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
    double u_length = u.normalize();
    double proj = b * u;

    if (proj <= epsilonEnds) { // start of bar

      double sqrDist = norm2(b);
      double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
      bool isNear = sqrDist < sumR * sumR;
      cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

    } else if (proj >= u_length - epsilonEnds) { // end of bar

      b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
      double sqrDist = norm2(b);
      double sumR = cells[ci].radius + cells[cj].radius + distVerlet;
      bool isNear = sqrDist < sumR * sumR;
      cells[ci].insertOrRemove(ci, cj, in, jn, isNear);

    } else { // on the bar

      vec2r T(-u.y, u.x);
      double dist = fabs(b * T) - (cells[ci].radius + cells[cj].radius);
      bool isNear = dist < distVerlet;
      cells[ci].insertOrRemove(ci, cj, in, jn, isNear);
    }
  }
}

/// Updates the list of neighbours
///
///
void Lhyphen::updateNeighbors_brute_force() {
  START_TIMER("updateNeighbors_brute_force");

  struct cellPair {
    size_t i;
    size_t j;
  };

  // Compute the expanded-AABB of each cell
  std::vector<AABB_2D> aabbs(cells.size());
  for (size_t ci = 0; ci < cells.size(); ci++) {
    aabbs[ci].set_single(cells[ci].nodes[0].pos);
    for (size_t in = 1; in < cells[ci].nodes.size(); in++) {
      aabbs[ci].add(cells[ci].nodes[in].pos);
    }
    aabbs[ci].enlarge(cells[ci].radius + 0.5 * distVerlet);
  }

  // Find intersecting expanded-AABB
  std::vector<cellPair> cellPairs;
  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (size_t cj = ci + 1; cj < cells.size(); cj++) {

      if (aabbs[ci].intersect(aabbs[cj])) {
        cellPair P;

        P.i = ci;
        P.j = cj;
        cellPairs.push_back(P);
      }
    }
  }

  // builds the list of neighbours
  for (size_t ip = 0; ip < cellPairs.size(); ++ip) {

    size_t ci = cellPairs[ip].i;
    size_t cj = cellPairs[ip].j;

    // test chaque sommet de i avec chaque barre de j
    for (size_t in = 0; in < cells[ci].nodes.size(); ++in) {
      for (size_t jn = 0; jn < cells[cj].nodes.size(); ++jn) {
        addNodeToBarNeighbor(ci, cj, in, jn);
      }
    }

    // test chaque sommet de j avec chaque barre de i
    for (size_t jn = 0; jn < cells[cj].nodes.size(); ++jn) {
      for (size_t in = 0; in < cells[ci].nodes.size(); ++in) {
        addNodeToBarNeighbor(cj, ci, jn, in, cells[ci].radius);
      }
    }
  }

  copy_neighbors_set_to_vec();
}

// En cours de devel. !!!!!!!!!! (ne pas utiliser)
// en plus il y a eu du changement, donc la procédure est devenue inadaptée
void Lhyphen::updateNeighbors_linkCells() {
  START_TIMER("updateNeighbors_linkCells");

  struct cellPair {
    size_t i;
    size_t j;
    size_t jn;
  };

  std::vector<AABB_2D> aabbs(cells.size());
  for (size_t ci = 0; ci < cells.size(); ci++) {
    aabbs[ci].set_single(cells[ci].nodes[0].pos);
    for (size_t in = 1; in < cells[ci].nodes.size(); in++) {
      aabbs[ci].add(cells[ci].nodes[in].pos);
    }
    aabbs[ci].enlarge(cells[ci].radius + 0.5 * distVerlet);
  }

  AABB_2D bigBox = aabbs[0];
  for (size_t ci = 1; ci < cells.size(); ci++) {
    bigBox.merge(aabbs[ci]);
  }

  vec2r subBoxSizes(linkCells_lx, linkCells_ly);
  linkCells2D lc(bigBox, subBoxSizes);
  for (size_t ci = 0; ci < cells.size(); ci++) {
    vec2r pos = aabbs[ci].getCenter();
    lc.add_body(ci, pos, aabbs[ci]);
  }

  std::vector<cellPair> cellPairs;

  AABB_2D_Cell *cc, *cv;
  for (size_t xc = 0; xc < lc.N.x; xc++) {
    for (size_t yc = 0; yc < lc.N.y; yc++) {

      cc = &(lc.cells[xc][yc]);

      for (size_t v = 0; v < lc.cells[xc][yc].pcells.size(); v++) {
        cv = cc->pcells[v];

        for (size_t bi = 0; bi < cc->bodies.size(); bi++) {
          size_t ci = cc->bodies[bi];
          for (size_t bj = 0; bj < cv->bodies.size(); bj++) {
            size_t cj = cv->bodies[bj];

            if (cj <= ci) {
              continue;
            }

            if (aabbs[ci].intersect(aabbs[cj])) {
              cellPair P;

              for (size_t jn = 0; jn < cells[cj].nodes.size(); jn++) {
                if (aabbs[ci].intersect(cells[cj].nodes[jn].pos)) {
                  P.i = ci;
                  P.j = cj;
                  P.jn = jn;
                  cellPairs.push_back(P);
                }
              } // for jn
            } // if intersect
          } // for bj
        } // for bi
      } // for v
    } // for yc
  } // for xc

  // Now we test if a too large bodies can collide another body in the
  cc = &(lc.oversized_bodies);
  for (size_t ix = 0; ix < lc.N.x; ++ix) {
    for (size_t iy = 0; iy < lc.N.y; ++iy) {
      cv = &(lc.cells[ix][iy]);

      for (size_t icc = 0; icc < cc->bodies.size(); ++icc) {
        size_t ci = cc->bodies[icc];
        for (size_t jcv = 0; jcv < cv->bodies.size(); ++jcv) {
          size_t cj = cv->bodies[jcv];

          if (cj <= ci)
            continue;

          if (aabbs[ci].intersect(aabbs[cj])) {
            cellPair P;

            for (size_t jn = 0; jn < cells[cj].nodes.size(); jn++) {
              if (aabbs[ci].intersect(cells[cj].nodes[jn].pos)) {
                P.i = ci;
                P.j = cj;
                P.jn = jn;
                cellPairs.push_back(P);
              }
            }
          }
        } // jcv
      } // icc
    } // iy
  } // ix

  // Now we test oversized vs oversized
  cc = &(lc.oversized_bodies);
  cv = cc;
  for (size_t icc = 0; icc < cc->bodies.size(); ++icc) {
    size_t ci = cc->bodies[icc];
    for (size_t jcv = 0; jcv < cv->bodies.size(); ++jcv) {
      size_t cj = cv->bodies[jcv];

      if (cj <= ci)
        continue;

      if (aabbs[ci].intersect(aabbs[cj])) {
        cellPair P;

        for (size_t jn = 0; jn < cells[cj].nodes.size(); jn++) {
          if (aabbs[ci].intersect(cells[cj].nodes[jn].pos)) {
            P.i = ci;
            P.j = cj;
            P.jn = jn;
            cellPairs.push_back(P);
          }
        }
      }
    } // jcv
  } // icc

  for (size_t ip = 0; ip < cellPairs.size(); ++ip) {

    size_t ci = cellPairs[ip].i;
    size_t cj = cellPairs[ip].j;
    size_t jn = cellPairs[ip].jn;

    for (size_t in = 0; in < cells[ci].nodes.size(); ++in) {
      addNodeToBarNeighbor(ci, cj, in, jn /*, cells[ci].radius */);
    }

    for (size_t in = 0; in < cells[ci].nodes.size(); in++) {
      addNodeToBarNeighbor(cj, ci, jn, in, cells[cj].radius);
    }
  }

  copy_neighbors_set_to_vec();
}

// ======================================================================================================
// Gluing the cells
// ======================================================================================================

///  Put a glued link where the distance is close enough
///
///  @param epsilonDist  maximum distance
///
void Lhyphen::glue(double epsilonDist, int modelGc) {
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
          Inter->glueState = modelGc;
          Inter->fn_coh = 0.0;
          Inter->ft_coh = 0.0;
        }

      } else {

        size_t jnext = cells[cj].nodes[jn].nextNode; // jnext est la fin de la barre (jn c'est le début)
        vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
        vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
        double u_length = u.normalize();
        double proj = b * u;

        if (proj < 0.0) {

          double sqrDist = norm2(b);
          double sumR = cells[ci].radius + cells[cj].radius + epsilonDist;
          if (sqrDist - sumR * sumR < 0.0) {
            Inter->glueState = modelGc;
            Inter->fn_coh = 0.0;
            Inter->ft_coh = 0.0;
          }

        } else if (proj > u_length) {

          b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
          double sqrDist = norm2(b);
          double sumR = cells[ci].radius + cells[cj].radius + epsilonDist;
          if (sqrDist - sumR * sumR < 0.0) {
            Inter->glueState = modelGc;
            Inter->fn_coh = 0.0;
            Inter->ft_coh = 0.0;
          }

        } else {

          vec2r T(-u.y, u.x);
          double sumR = cells[ci].radius + cells[cj].radius;
          double dist = b * T; // dot product
          double dn = fabs(dist) - sumR;
          if (dn < epsilonDist) {
            Inter->glueState = modelGc;
            Inter->fn_coh = 0.0;
            Inter->ft_coh = 0.0;
          }

        } // end if proj
      }
    }
  }

  associateGlue(modelGc);
}

///  Given a cell-cell interaction, compute the "position" of interaction.
///  This is a point such that the distance between this point and the
///  node in (ci,in) is equal to the radius of node in (ci,in).
///  This point is either on the bar (cj,jn) or one of the disks at its
///  extremities (cj,jn) or (cj,jnext).
///  This is used to compute the force in the modelGc.
///
///  @param ci      cell i
///  @param cj      cell j
///  @param in      node number in cell ci
///  @param jn      node number at the start of a bar in cell cj
///  @param pos     output position
///  @return        code indicating which case was used:
///                 0: no bar in cj
///                 1: interaction with disk j
///                 2: interaction with disk jnext
///                 3: interaction with bar
///
///
int Lhyphen::getPosition(size_t ci, size_t cj, size_t in, size_t jn, vec2r &pos) {
  if (cells[cj].nodes[jn].nextNode == null_size_t) { // This is the case where the bar in cj does not exist.

    vec2r b = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
    b *= cells[ci].radius / (cells[ci].radius + cells[cj].radius);
    pos = cells[ci].nodes[in].pos + b;
    return 0;

  } else { // a disk (ci,in) interacts with a bar (cj, jn -> jnext)

    size_t jnext = cells[cj].nodes[jn].nextNode;
    // jnext est le numéro du noeud à la fin de la barre dans la cellule cj (jn c'est le début)
    vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
    vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
    double u_length = u.normalize();
    double proj = b * u;
    double fact = cells[cj].radius / (cells[ci].radius + cells[cj].radius);

    if (proj <= 0.0) { // ====================== disque j du début

      vec2r ut = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
      pos = cells[cj].nodes[jn].pos + fact * ut;
      return 1;

    } else if (proj >= u_length) { // ====================== disque jnext (de fin)

      vec2r ut = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
      pos = cells[cj].nodes[jnext].pos + fact * ut;
      return 2;

    } else { // ====================== sur la barre

      pos = cells[cj].nodes[jn].pos + proj * u;
      vec2r ut = cells[ci].nodes[in].pos - pos;
      pos += fact * ut;
      return 3;
    }
  }
}

///  Associate a model of glue to all the neighbors between cells having a positive glueState.
///
///  @param modelGc the model of glue to associate
///
void Lhyphen::associateGlue(int modelGc, double activationLength) {

  struct Connect {
    size_t ic; // indice de la première cellule dans Sample::cells
    size_t jc; // indice de la seconde cellule dans Sample::cells (ic <= jc)
    size_t in; // indice du noeud de la première cellule dans Sample::cells[ic]
    size_t jn; // indice du noeud de la seconde cellule dans Sample::cells[jc]
    Neighbor *woami{nullptr};
    Connect(size_t t_ic, size_t t_jc, size_t t_in, size_t t_jn, Neighbor *t_woami)
        : ic(t_ic), jc(t_jc), in(t_in), jn(t_jn) {
      woami = t_woami;
    }
  };

  // on construit un vector de neighbors avec uniquement les points de colle
  std::vector<Connect> connects;
  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {
      Neighbor *Inter = const_cast<Neighbor *>(std::addressof(*InterIt));

      if (Inter->glueState > 0) {
        connects.push_back(Connect(ci, Inter->jc, Inter->in, Inter->jn, Inter));
      }
    }
  }

  // on réduit les plan collés à 2 points les plus éloignées et on calcul la distance (= surface collée)
  std::map<std::pair<size_t, size_t>, std::vector<size_t>> planMap;
  for (size_t i = 0; i < connects.size(); i++) {
    size_t cinf = std::min(connects[i].ic, connects[i].jc);
    size_t csup = std::max(connects[i].ic, connects[i].jc);
    std::pair<size_t, size_t> p = std::make_pair(cinf, csup);
    planMap[p].push_back(i);
  }

  for (const auto &elem : planMap) {

    double dmax = -1.0;
    size_t ii = 0, jj = 0;
    for (size_t ei = 0; ei < elem.second.size(); ei++) {
      for (size_t ej = ei + 1; ej < elem.second.size(); ej++) {

        vec2r p1;
        getPosition(connects[elem.second[ei]].ic, connects[elem.second[ei]].jc, connects[elem.second[ei]].in,
                    connects[elem.second[ei]].jn, p1);
        vec2r p2;
        getPosition(connects[elem.second[ej]].ic, connects[elem.second[ej]].jc, connects[elem.second[ej]].in,
                    connects[elem.second[ej]].jn, p2);

        double dst = (p2 - p1).length();
        if (dst > dmax) {
          ii = elem.second[ei];
          jj = elem.second[ej];
          dmax = dst;
        }
      }
    }

    for (size_t e = 0; e < elem.second.size(); e++) {
      connects[elem.second[e]].woami->glueState = 0;
    }

    if (dmax >= activationLength) {
      connects[ii].woami->glueState = modelGc;
      connects[jj].woami->glueState = modelGc;

      connects[ii].woami->length = dmax;
      connects[jj].woami->length = dmax;
    }
  }

  // on re-parcours les neighbors pour identifier les brothers (liens collés aux mêmes cellules)
  size_t cinfi, csupi;
  size_t cinfj, csupj;
  for (size_t i = 0; i < connects.size(); i++) {
    cinfi = std::min(connects[i].ic, connects[i].jc);
    csupi = std::max(connects[i].ic, connects[i].jc);
    for (size_t j = 0; j < connects.size(); j++) { // FIXME: commencer à i+1 ????
      if (i == j) {
        continue;
      }
      cinfj = std::min(connects[j].ic, connects[j].jc);
      csupj = std::max(connects[j].ic, connects[j].jc);
      if (cinfi == cinfj && csupi == csupj && connects[i].woami->glueState > 0 && connects[j].woami->glueState > 0) {
        connects[j].woami->brother = connects[i].woami;
        connects[i].woami->brother = connects[j].woami;
        break;
      }
    }
  }
}

// ======================================================================================================
// Computing
// ======================================================================================================

// Il faudra plus tard utiliser cette fonction pour éviter les duplications dans la fonction qui calcul les forces
// de contact
// fn_visc: instantaneous viscous normal force (does NOT enter breakage criterion, only force application)
void Lhyphen::glue_breakage(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jn, size_t jnext, double wbeg,
                            double wend, double fn_visc) {

  if (Inter->glueState == GLUE_FORCE_THRESHOLD) {

    // breakage criterion uses elastic spring forces only (not fn_visc)
    double zeta =
        -Inter->fn_coh / Inter->fn_coh_max + pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
    if (zeta > 0.0) {
      Inter->fn_coh = 0.0;
      Inter->ft_coh = 0.0;
      Inter->glueState = 0;
      if (Inter->brother != nullptr) {
        Inter->brother->fn_coh = 0.0;
        Inter->brother->ft_coh = 0.0;
        Inter->brother->glueState = 0;
      }
    } else {
      vec2r finc = (Inter->fn_coh + fn_visc) * Inter->n + Inter->ft_coh * Inter->T;
      cells[ci].nodes[in].force += finc;
      if (wbeg > 0.0) cells[cj].nodes[jn].force -= wbeg * finc;
      if (wend > 0.0) cells[cj].nodes[jnext].force -= wend * finc;
    }

  } else if (Inter->glueState == GLUE_GC) {

    // stored energy W uses elastic spring forces only (not fn_visc)
    double W = 0.0;
    if (Inter->fn_coh < 0.0) W += Inter->fn_coh * Inter->fn_coh / Inter->kn_coh; // heaviside(dn)
    W += Inter->ft_coh * Inter->ft_coh / Inter->kt_coh;

    if (Inter->brother != nullptr) {
      if (Inter->brother->fn < 0.0) W += Inter->brother->fn_coh * Inter->brother->fn_coh / Inter->brother->kn_coh;
      W += Inter->brother->ft_coh * Inter->brother->ft_coh / Inter->brother->kt_coh;

      double G = (Inter->length > 1.0e-12) ? W / (2.0 * Inter->length) : 0.0;

      if (G > Inter->Gc) {
        Inter->fn_coh = 0.0;
        Inter->ft_coh = 0.0;
        Inter->glueState = 0;
        Inter->brother->fn_coh = 0.0;
        Inter->brother->ft_coh = 0.0;
        Inter->brother->glueState = 0;
        cumulatedG += G;
        cumulatedL += Inter->length;
      } else {
        vec2r finc = (Inter->fn_coh + fn_visc) * Inter->n + Inter->ft_coh * Inter->T;
        cells[ci].nodes[in].force += finc;
        if (wbeg > 0.0) cells[cj].nodes[jn].force -= wbeg * finc;
        if (wend > 0.0) cells[cj].nodes[jnext].force -= wend * finc;
      }

    } else {
      // no brother: break immediately
      Inter->fn_coh = 0.0;
      Inter->ft_coh = 0.0;
      Inter->glueState = 0;
    }
  }
}

/// Calculation of interaction forces (between cell-elements)
///
/// !!!!!!!!!!!
/// Cette fonction est extrêmement longue et imbriquée
/// Il faudrait trouver une solution pour la rendre moins compliquée
/// !!!!!!!!!!!
///
void Lhyphen::computeInteractionForces_originalVersion() {
  START_TIMER("computeInteractionForces_originalVersion");

  vec2r vrel;

  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (size_t ii = 0; ii < cells[ci].vec_neighbors.size(); ++ii) {
      Neighbor *Inter = cells[ci].vec_neighbors[ii];

      size_t cj = Inter->jc;
      size_t in = Inter->in;
      size_t jn = Inter->jn;
      size_t jnext = cells[cj].nodes[jn].nextNode;

      if (jnext == null_size_t) {
        // =============================================================================
        // This is the case where the bar in cj does not exist.
        // In other words, it's an interaction between a disk of ci and a disc of cj
        // (which corresponds to the end of a bar)
        // =============================================================================

        // FIXME: Ce cas n'est pas correctement traité mais n'est pas utiliser en fait avec
        //        les cellules fermées

        vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
        double sqrDist = norm2(branch);
        double sumR = cells[ci].radius + cells[cj].radius;

        if (sqrDist - sumR * sumR < 0.0) { // overlap
          double b_length = sqrt(sqrDist);
          double dn = b_length - sumR;
          Inter->n = branch / b_length;
          Inter->fn = -kn * dn;
          if (adaptativeStiffness == 1) {
            Inter->fn *= sumR / (sumR + dn);
          }
          if (Inter->glueState == GLUE_NONE) {
            Inter->fn -= fadh;
          }
        } else {
          if (Inter->glueState > 0) {
            double b_length = sqrt(sqrDist);
            Inter->n = branch / b_length;
          }
          Inter->fn = 0.0;
          Inter->ft = 0.0;
        }

        if (Inter->glueState > 0) {
          // ...
          vrel = cells[cj].nodes[jn].vel - cells[ci].nodes[in].vel;
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
        // a disk (ci,in) interacts with a bar (cj, jn -> jnext)

        // jnext est le numéro du noeud à la fin de la barre dans la cellule cj (jn c'est le début)
        vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
        vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
        double u_length = u.normalize();
        double proj = b * u;

        // on doit gérer le calcul des forces (contact, frottement, cohesion...)
        // pour chacun des 3 cas suivants :
        //   cas 1 / disque (ci, in)   ---   disque (cj, jn) début de la barre
        //   cas 2 / disque (ci, in)   ---   disque (cj, jnext) fin de la barre
        //   cas 3 / disque (ci, in)   ---   barre (cj, jn--jnext)

        if (proj <= 0.0) { // ** cas 1 / ====================================================== disque j du début

          double sqrDist = norm2(b);
          double sumR = cells[ci].radius + cells[cj].radius;
          if (sqrDist - sumR * sumR < 0.0) {
            Inter->contactState = TOUCHING;
            double b_lenght = sqrt(sqrDist);
            double dn = b_lenght - sumR;
            Inter->n = b / b_lenght;
            Inter->fn = -kn * dn; // +- visc * vrel;
            if (adaptativeStiffness == LH_ENABLED) {
              Inter->fn *= sumR / (sumR + dn);
            }

            vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jn].vel;
            Inter->T.set(-Inter->n.y, Inter->n.x);
            Inter->ft -= kt * (vrel * Inter->T) * dt;
            double threshold = mu * Inter->fn;
            if (Inter->ft > threshold) {
              Inter->ft = threshold;
            } else if (Inter->ft < -threshold) {
              Inter->ft = -threshold;
            }

            if (Inter->glueState == GLUE_NONE) {
              Inter->fn -= fadh; // cette adhesion ne peut agir que si il n'y a pas de cohésion solide
            }

            // transfert des force vers les noeuds concernés
            vec2r finc = Inter->fn * Inter->n;
            cells[ci].nodes[in].force += finc;
            cells[cj].nodes[jn].force -= finc;
          } else {
            Inter->contactState = NOT_TOUCHING;
            Inter->fn = 0.0;
            Inter->ft = 0.0;
          }

          if (Inter->glueState >= 1) {
            if (Inter->contactState == NOT_TOUCHING) {
              double b_lenght = sqrt(sqrDist);
              Inter->n = b / b_lenght;
              vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jn].vel;
              Inter->T.set(-Inter->n.y, Inter->n.x);
            }

            Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
            Inter->ft_coh += -Inter->kt_coh * vrel * Inter->T * dt;

            if (Inter->fn_coh > 0.0) {
              Inter->fn_coh = 0.0;
            }

            // glue_breakage(Inter, ci, cj, in, jn, jnext, 1.0, 0.0);
            if (Inter->glueState == GLUE_FORCE_THRESHOLD) { // RUPTURE par forces
              // rupture
              double zeta = -Inter->fn_coh / Inter->fn_coh_max +
                            pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
              if (zeta > 0.0) {
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
                if (Inter->brother != nullptr) {
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;
                }
              } else {
                // transfert des force vers les noeuds concernés
                vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                cells[ci].nodes[in].force += finc;
                cells[cj].nodes[jn].force -= finc;
              }
            } else if (Inter->glueState == GLUE_GC) { // model rupture Gc
              double W = 0.0;
              if (Inter->fn_coh < 0.0) { // heaviside(dn)
                W += Inter->fn_coh * Inter->fn_coh / Inter->kn_coh;
              }
              W += Inter->ft_coh * Inter->ft_coh / Inter->kt_coh;

              if (Inter->brother != nullptr) {

                if (Inter->brother->fn < 0.0) { // heaviside(dn)
                  W += Inter->brother->fn_coh * Inter->brother->fn_coh / Inter->brother->kn_coh;
                }
                W += Inter->brother->ft_coh * Inter->brother->ft_coh / Inter->brother->kt_coh;

                double G = 0.0;
                if (Inter->length > 1.0e-12) {
                  G = W / (2.0 * Inter->length);
                }

                if (G > Inter->Gc) {
                  // breakage
                  Inter->fn_coh = 0.0;
                  Inter->ft_coh = 0.0;
                  Inter->glueState = 0;
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;

                  cumulatedG += G;
                  cumulatedL += Inter->length;
                } else {
                  // transfert des force vers les noeuds concernés
                  vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                  cells[ci].nodes[in].force += finc;

                  cells[cj].nodes[jn].force -= finc;
                }

              } else {
                // breakage if a 'brother' exists
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
              }
            } // fin glueState == 2
          } // if glueState == 1 ou 2

        } // fin de (disque j du début)
        else if (proj >= u_length) { // ** cas 2 / ============================================ disque jnext (de fin)

          b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
          double sqrDist = norm2(b);
          double sumR = cells[ci].radius + cells[cj].radius;

          if (sqrDist - sumR * sumR < 0.0) {
            Inter->contactState = TOUCHING;
            double b_lenght = sqrt(sqrDist);
            double dn = b_lenght - sumR;
            Inter->n = b / b_lenght;
            Inter->fn = -kn * dn;

            if (adaptativeStiffness == LH_ENABLED) {
              Inter->fn *= sumR / (sumR + dn);
            }

            vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jnext].vel;
            Inter->T.set(-Inter->n.y, Inter->n.x);
            Inter->ft -= kt * (vrel * Inter->T) * dt;
            double threshold = mu * Inter->fn;
            if (Inter->ft > threshold) {
              Inter->ft = threshold;
            } else if (Inter->ft < -threshold) {
              Inter->ft = -threshold;
            }

            if (Inter->glueState == GLUE_NONE) {
              Inter->fn -= fadh; // this adhesion can only act if there is no solid cohesion
            }

            // force transfer to the relevant nodes
            vec2r finc = Inter->fn * Inter->n;
            cells[ci].nodes[in].force += finc;
            cells[cj].nodes[jnext].force -= finc;

          } else {
            Inter->contactState = NOT_TOUCHING;
            Inter->fn = 0.0;
            Inter->ft = 0.0;
          }

          if (Inter->glueState >= 1) {
            if (Inter->contactState == NOT_TOUCHING) {
              double b_lenght = sqrt(sqrDist);
              Inter->n = b / b_lenght;
              vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jnext].vel;
              Inter->T.set(-Inter->n.y, Inter->n.x);
            }

            Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
            Inter->ft_coh += -Inter->kt_coh * vrel * Inter->T * dt;

            if (Inter->fn_coh > 0.0) {
              Inter->fn_coh = 0.0;
            }

            // glue_breakage(Inter, ci, cj, in, jnext, 0.0, 1.0);
            if (Inter->glueState == GLUE_FORCE_THRESHOLD) { // RUPTURE par forces

              double zeta = -Inter->fn_coh / Inter->fn_coh_max +
                            pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
              if (zeta > 0.0) {
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
                if (Inter->brother != nullptr) {
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;
                }
              } else {
                // transfert des force vers les noeuds concernés
                vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                cells[ci].nodes[in].force += finc;
                cells[cj].nodes[jnext].force -= finc;
              }
            } else if (Inter->glueState == GLUE_GC) { // model rupture Gc
              double W = 0.0;
              if (Inter->fn_coh < 0.0) { // heaviside(dn)
                W += Inter->fn_coh * Inter->fn_coh / Inter->kn_coh;
              }
              W += Inter->ft_coh * Inter->ft_coh / Inter->kt_coh;

              if (Inter->brother != nullptr) {

                if (Inter->brother->fn < 0.0) { // heaviside(dn)
                  W += Inter->brother->fn_coh * Inter->brother->fn_coh / Inter->brother->kn_coh;
                }
                W += Inter->brother->ft_coh * Inter->brother->ft_coh / Inter->brother->kt_coh;

                double G = 0.0;
                if (Inter->length > 1.0e-12) {
                  G = W / (2.0 * Inter->length);
                }

                if (G > Inter->Gc) {
                  // cassage
                  Inter->fn_coh = 0.0;
                  Inter->ft_coh = 0.0;
                  Inter->glueState = 0;
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;

                  cumulatedG += G;
                  cumulatedL += Inter->length;
                } else {
                  // transfert des force vers les noeuds concernés
                  vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                  cells[ci].nodes[in].force += finc;

                  cells[cj].nodes[jnext].force -= finc;
                }

              } else {
                // cassage si il n'y a pas de frère
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
              }
            } // fin glueState == 2
          } // if glueState == 1 ou 2

        } // fin "disque jnext (de fin)"
        else { // ** cas 3 / ====================================================== sur la barre

          vec2r urot(-u.y, u.x); // turn 90°
          double wend = 0.0;
          double wbeg = 0.0;
          double sumR = cells[ci].radius + cells[cj].radius;
          double dist = b * urot;
          double dn = fabs(dist) - sumR;

          if (dn < 0.0) { // overlap
            Inter->contactState = TOUCHING;
            Inter->n = urot;
            if (dist < 0.0) {
              Inter->n *= -1.0;
            }
            // here the normal vector is oriented from the bar to the disc

            Inter->fn = -kn * dn;
            if (adaptativeStiffness == LH_ENABLED) {
              Inter->fn *= sumR / (sumR + dn);
            }

            // vitesse relative noeud par rapport à barre
            wend = proj / u_length;
            wbeg = 1.0 - wend;
            // on neglige l'épaisseur des barres
            // à voir si c'est ok car on fait une approximation
            vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
            Inter->T.set(-Inter->n.y, Inter->n.x);
            Inter->ft -= kt * (vrel * Inter->T) * dt;
            double threshold = mu * Inter->fn;
            if (Inter->ft > threshold) {
              Inter->ft = threshold;
            } else if (Inter->ft < -threshold) {
              Inter->ft = -threshold;
            }

            if (Inter->glueState == GLUE_NONE) {
              Inter->fn -= fadh;
            }

            // transfert des forces vers les noeuds concernés
            // n est orienté de la barre vers le noeud
            // donc fn positif correspond à une force de la barre qui pousse le noeud

            vec2r finc = Inter->fn * Inter->n + Inter->ft * Inter->T;
            cells[ci].nodes[in].force += finc;

            cells[cj].nodes[jn].force -= wbeg * finc;
            cells[cj].nodes[jnext].force -= wend * finc;

          } else {
            Inter->contactState = NOT_TOUCHING;
            Inter->fn = 0.0;
            Inter->ft = 0.0;
          }

          // glue_breakage(Inter, ci, cj, in, jn, jnext, wbeg, wend);
          if (Inter->glueState >= 1) {

            // if there is no contact, some variables have not been calculated yet
            if (Inter->contactState == NOT_TOUCHING) {
              Inter->n = urot;
              if (dist < 0.0) {
                Inter->n *= -1.0;
              }

              wend = proj / u_length;
              wbeg = 1.0 - wend;
              vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
              Inter->T.set(-Inter->n.y, Inter->n.x);
            }

            Inter->fn_coh += -Inter->kn_coh * vrel * Inter->n * dt;
            Inter->ft_coh += -Inter->kt_coh * vrel * Inter->T * dt;

            if (Inter->fn_coh > 0.0) {
              Inter->fn_coh = 0.0;
            }

            // rupture
            // glue_breakage(Inter, ci, cj, in, jn, jnext, wbeg, wend);
            if (Inter->glueState == GLUE_FORCE_THRESHOLD) { // RUPTURE par forces
              double zeta = -Inter->fn_coh / Inter->fn_coh_max +
                            pow(fabs(Inter->ft_coh) / Inter->ft_coh_max, Inter->yieldPower) - 1.0;
              if (zeta > 0.0) {
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
                if (Inter->brother != nullptr) {
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;
                }
              } else {
                // transfert des force vers les noeuds concernés
                vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                cells[ci].nodes[in].force += finc;

                cells[cj].nodes[jn].force -= wbeg * finc;
                cells[cj].nodes[jnext].force -= wend * finc;
              }
            } else if (Inter->glueState == GLUE_GC) { // model rupture Gc
              double W = 0.0;
              if (Inter->fn_coh < 0.0) { // heaviside(dn)
                W += Inter->fn_coh * Inter->fn_coh / Inter->kn_coh;
              }
              W += Inter->ft_coh * Inter->ft_coh / Inter->kt_coh;

              if (Inter->brother != nullptr) {

                if (Inter->brother->fn < 0.0) { // heaviside(dn)
                  W += Inter->brother->fn_coh * Inter->brother->fn_coh / Inter->brother->kn_coh;
                }
                W += Inter->brother->ft_coh * Inter->brother->ft_coh / Inter->brother->kt_coh;

                double G = 0.0;
                if (Inter->length > 1.0e-12) {
                  G = W / (2.0 * Inter->length);
                }

                if (G > Inter->Gc) {
                  // cassage
                  Inter->fn_coh = 0.0;
                  Inter->ft_coh = 0.0;
                  Inter->glueState = 0;
                  Inter->brother->fn_coh = 0.0;
                  Inter->brother->ft_coh = 0.0;
                  Inter->brother->glueState = 0;

                  cumulatedG += G;
                  cumulatedL += Inter->length;
                } else {
                  // transfert des force vers les noeuds concernés
                  vec2r finc = Inter->fn_coh * Inter->n + Inter->ft_coh * Inter->T;
                  cells[ci].nodes[in].force += finc;

                  cells[cj].nodes[jn].force -= wbeg * finc;
                  cells[cj].nodes[jnext].force -= wend * finc;
                }

              } else {
                // cassage si il n'y a pas de frère
                Inter->fn_coh = 0.0;
                Inter->ft_coh = 0.0;
                Inter->glueState = 0;
              }
            } // fin glueState == 2
          }
        } // fin "sur la barre"
      }

    } // boucle sur les voisins dans ci
  } // boucle sur les cellules ci
}

// disk(ci,in) vs disk(cj,jn) — special case: jn has no next node
// FIXME: normal convention here is n = (cj-ci)/|...| (opposite to other handlers); not used with closed cells
void Lhyphen::handleDiskDiskContact(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jn) {
  vec2r branch = cells[cj].nodes[jn].pos - cells[ci].nodes[in].pos;
  double sqrDist = norm2(branch);
  double sumR = cells[ci].radius + cells[cj].radius;
  vec2r vrel;

  double cn = 0.0;
  if (viscnrate > 0.0) {
    double m_eff = 0.5 * (cells[ci].nodes[in].mass + cells[cj].nodes[jn].mass);
    cn = viscnrate * 2.0 * sqrt(m_eff * kn);
  }

  if (sqrDist < sumR * sumR) {
    Inter->contactState = TOUCHING;
    double b_length = sqrt(sqrDist);
    double dn = b_length - sumR;
    Inter->n = branch / b_length;
    Inter->fn = -kn * dn;
    if (adaptativeStiffness == LH_ENABLED)
      Inter->fn *= sumR / (sumR + dn);
    vrel = cells[cj].nodes[jn].vel - cells[ci].nodes[in].vel;
    Inter->fn -= cn * (vrel * Inter->n);
    if (Inter->glueState == GLUE_NONE)
      Inter->fn -= fadh;
    Inter->T.set(-Inter->n.y, Inter->n.x);
    Inter->ft -= kt * (vrel * Inter->T) * dt;
    double threshold = mu * Inter->fn;
    if (Inter->ft > threshold)
      Inter->ft = threshold;
    else if (Inter->ft < -threshold)
      Inter->ft = -threshold;
  } else {
    Inter->contactState = NOT_TOUCHING;
    if (Inter->glueState > 0)
      Inter->n = branch / sqrt(sqrDist);
    Inter->fn = 0.0;
    Inter->ft = 0.0;
  }

  vec2r finc = Inter->fn * Inter->n;
  cells[ci].nodes[in].force += finc;
  cells[cj].nodes[jn].force -= finc;

  if (Inter->glueState > 0) {
    if (Inter->contactState == NOT_TOUCHING) {
      vrel = cells[cj].nodes[jn].vel - cells[ci].nodes[in].vel;
      Inter->T.set(-Inter->n.y, Inter->n.x);
    }
    Inter->fn_coh += -Inter->kn_coh * (vrel * Inter->n) * dt;
    Inter->ft_coh += -Inter->kt_coh * (vrel * Inter->T) * dt;
    if (Inter->fn_coh > 0.0)
      Inter->fn_coh = 0.0;
    glue_breakage(Inter, ci, cj, in, jn, null_size_t, 1.0, 0.0, -cn * (vrel * Inter->n));
  }
}

// disk(ci,in) vs disk-end node jcontact (either jn or jnext of the bar)
// wbeg, wend: bar weights for distributing cohesive forces to jn and jnext
void Lhyphen::handleDiskEndContact(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jcontact, size_t jn,
                                   size_t jnext, double wbeg, double wend) {
  vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jcontact].pos;
  double sqrDist = norm2(b);
  double sumR = cells[ci].radius + cells[cj].radius;
  vec2r vrel;

  double cn = 0.0;
  if (viscnrate > 0.0) {
    double m_eff = 0.5 * (cells[ci].nodes[in].mass + cells[cj].nodes[jcontact].mass);
    cn = viscnrate * 2.0 * sqrt(m_eff * kn);
  }

  if (sqrDist < sumR * sumR) {
    Inter->contactState = TOUCHING;
    double b_length = sqrt(sqrDist);
    double dn = b_length - sumR;
    Inter->n = b / b_length;
    Inter->fn = -kn * dn;
    if (adaptativeStiffness == LH_ENABLED)
      Inter->fn *= sumR / (sumR + dn);
    vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jcontact].vel;
    Inter->fn -= cn * (vrel * Inter->n);
    Inter->T.set(-Inter->n.y, Inter->n.x);
    Inter->ft -= kt * (vrel * Inter->T) * dt;
    double threshold = mu * Inter->fn;
    if (Inter->ft > threshold)
      Inter->ft = threshold;
    else if (Inter->ft < -threshold)
      Inter->ft = -threshold;
    if (Inter->glueState == GLUE_NONE)
      Inter->fn -= fadh;
    vec2r finc = Inter->fn * Inter->n;
    cells[ci].nodes[in].force += finc;
    if (wbeg > 0.0) cells[cj].nodes[jn].force -= wbeg * finc;
    if (wend > 0.0) cells[cj].nodes[jnext].force -= wend * finc;
  } else {
    Inter->contactState = NOT_TOUCHING;
    Inter->fn = 0.0;
    Inter->ft = 0.0;
  }

  if (Inter->glueState >= 1) {
    if (Inter->contactState == NOT_TOUCHING) {
      Inter->n = b / sqrt(sqrDist);
      vrel = cells[ci].nodes[in].vel - cells[cj].nodes[jcontact].vel;
      Inter->T.set(-Inter->n.y, Inter->n.x);
    }
    Inter->fn_coh += -Inter->kn_coh * (vrel * Inter->n) * dt;
    Inter->ft_coh += -Inter->kt_coh * (vrel * Inter->T) * dt;
    if (Inter->fn_coh > 0.0)
      Inter->fn_coh = 0.0;
    glue_breakage(Inter, ci, cj, in, jn, jnext, wbeg, wend, -cn * (vrel * Inter->n));
  }
}

// disk(ci,in) vs bar segment (cj, jn->jnext), projection falls strictly between the two endpoints
// u, u_length, proj, b are pre-computed by the caller (handleDiskBarInteraction)
void Lhyphen::handleDiskOnBarContact(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jn, size_t jnext,
                                     const vec2r &u, double u_length, double proj, const vec2r &b) {
  vec2r urot(-u.y, u.x); // 90° rotation of bar direction: bar normal
  double sumR = cells[ci].radius + cells[cj].radius;
  double dist = b * urot; // signed distance from bar to disk center
  double dn = fabs(dist) - sumR;
  double wend = proj / u_length;
  double wbeg = 1.0 - wend;
  vec2r vrel;

  double cn = 0.0;
  if (viscnrate > 0.0) {
    double m_bar = wbeg * cells[cj].nodes[jn].mass + wend * cells[cj].nodes[jnext].mass;
    double m_eff = 0.5 * (cells[ci].nodes[in].mass + m_bar);
    cn = viscnrate * 2.0 * sqrt(m_eff * kn);
  }

  if (dn < 0.0) {
    Inter->contactState = TOUCHING;
    Inter->n = (dist >= 0.0) ? urot : -urot; // oriented from bar toward disk
    Inter->fn = -kn * dn;
    if (adaptativeStiffness == LH_ENABLED)
      Inter->fn *= sumR / (sumR + dn);
    vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
    Inter->fn -= cn * (vrel * Inter->n);
    Inter->T.set(-Inter->n.y, Inter->n.x);
    Inter->ft -= kt * (vrel * Inter->T) * dt;
    double threshold = mu * Inter->fn;
    if (Inter->ft > threshold)
      Inter->ft = threshold;
    else if (Inter->ft < -threshold)
      Inter->ft = -threshold;
    if (Inter->glueState == GLUE_NONE)
      Inter->fn -= fadh;
    vec2r finc = Inter->fn * Inter->n + Inter->ft * Inter->T;
    cells[ci].nodes[in].force += finc;
    cells[cj].nodes[jn].force -= wbeg * finc;
    cells[cj].nodes[jnext].force -= wend * finc;
  } else {
    Inter->contactState = NOT_TOUCHING;
    Inter->fn = 0.0;
    Inter->ft = 0.0;
  }

  if (Inter->glueState >= 1) {
    if (Inter->contactState == NOT_TOUCHING) {
      Inter->n = (dist >= 0.0) ? urot : -urot;
      vrel = cells[ci].nodes[in].vel - (wbeg * cells[cj].nodes[jn].vel + wend * cells[cj].nodes[jnext].vel);
      Inter->T.set(-Inter->n.y, Inter->n.x);
    }
    Inter->fn_coh += -Inter->kn_coh * (vrel * Inter->n) * dt;
    Inter->ft_coh += -Inter->kt_coh * (vrel * Inter->T) * dt;
    if (Inter->fn_coh > 0.0)
      Inter->fn_coh = 0.0;
    glue_breakage(Inter, ci, cj, in, jn, jnext, wbeg, wend, -cn * (vrel * Inter->n));
  }
}

// disk(ci,in) vs bar(cj, jn->jnext): projects the disk onto the bar and dispatches
// to the appropriate geometric sub-case (start cap, end cap, or bar body)
void Lhyphen::handleDiskBarInteraction(Neighbor *Inter, size_t ci, size_t cj, size_t in, size_t jn, size_t jnext) {
  // START_TIMER("handleDiskBarInteraction");
  vec2r b = cells[ci].nodes[in].pos - cells[cj].nodes[jn].pos;
  vec2r u = cells[cj].nodes[jnext].pos - cells[cj].nodes[jn].pos;
  double u_length = u.normalize();
  double proj = b * u;

  if (proj <= 0.0)
    handleDiskEndContact(Inter, ci, cj, in, jn, jn, jnext, 1.0, 0.0);
  else if (proj >= u_length)
    handleDiskEndContact(Inter, ci, cj, in, jnext, jn, jnext, 0.0, 1.0);
  else
    handleDiskOnBarContact(Inter, ci, cj, in, jn, jnext, u, u_length, proj, b);
}

void Lhyphen::computeInteractionForces() {
  START_TIMER("computeInteractionForces");
  for (size_t ci = 0; ci < cells.size(); ci++) {
    for (size_t ii = 0; ii < cells[ci].vec_neighbors.size(); ++ii) {
      Neighbor *Inter = cells[ci].vec_neighbors[ii];
      const size_t cj = Inter->jc, in = Inter->in, jn = Inter->jn;
      const size_t jnext = cells[cj].nodes[jn].nextNode;
      if (jnext == null_size_t)
        handleDiskDiskContact(Inter, ci, cj, in, jn);
      else
        handleDiskBarInteraction(Inter, ci, cj, in, jn, jnext);
    }
  }
}

///  Calculation of internal cell forces
///
void Lhyphen::computeNodeForces() {
  START_TIMER("computeNodeForces");

#ifdef _OPENMP
  int chunk_size = (int)cells.size() / nbThreads;
#endif

  // Forces normales dans les barres
  {
    START_TIMER("Bar_forces");

#pragma omp parallel for schedule(static, chunk_size)
    for (size_t c = 0; c < cells.size(); ++c) {
      for (size_t b = 0; b < cells[c].bars.size(); ++b) {
        size_t i = cells[c].bars[b].i;
        size_t j = cells[c].bars[b].j;
        vec2r posi = cells[c].nodes[i].pos;
        vec2r posj = cells[c].nodes[j].pos;

        vec2r branch = posj - posi;
        vec2r n = branch;
        double l = n.normalize();            // n est orienté de i vers j
        double dn = l - cells[c].bars[b].l0; // dn positif = allongement
        cells[c].bars[b].fn = -cells[c].bars[b].kn * dn;

        // viscosité linéaire
        if (cells[c].bars[b].visc > 0.0) {
          double vn = (cells[c].nodes[i].vel - cells[c].nodes[j].vel) * n;
          cells[c].bars[b].fn += cells[c].bars[b].visc * vn;
        }

        // TODO: ajouter plasticité dans les barres (?)

        vec2r finc = cells[c].bars[b].fn * n; // en cas d'allongement finc est orienté de j vers i
        cells[c].nodes[i].force -= finc;
        cells[c].nodes[j].force += finc;
      }
    }
  }

  // Moment au niveau des noeuds (entre les barres d'une même cellule).
  // Puisque les noeuds n'ont pas de ddl en rotation, imposer un moment est fait ici
  // en imposant une force à chaque extrémité de barre
  {
    START_TIMER("Node_moments");

#pragma omp parallel for schedule(static, chunk_size)
    for (size_t c = 0; c < cells.size(); ++c) {
      for (size_t n = 0; n < cells[c].nodes.size(); ++n) {

        size_t prev = cells[c].nodes[n].prevNode;
        size_t next = cells[c].nodes[n].nextNode;

        // si c'est une extrémité on ne calcul pas de moment
        if (prev == null_size_t || next == null_size_t) {
          continue;
        }

        vec2r nextDir = cells[c].nodes[next].pos - cells[c].nodes[n].pos;
        vec2r prevDir = cells[c].nodes[prev].pos - cells[c].nodes[n].pos;

        double lNextSqrInv = 1.0 / (nextDir * nextDir);
        vec2r TNext(-nextDir.y, nextDir.x);
        double omegaNext = ((cells[c].nodes[next].vel - cells[c].nodes[n].vel) * TNext * lNextSqrInv);

        double lPrevSqrInv = 1.0 / (prevDir * prevDir);
        vec2r TPrev(-prevDir.y, prevDir.x);
        double omegaPrev = ((cells[c].nodes[prev].vel - cells[c].nodes[n].vel) * TPrev * lPrevSqrInv);

        cells[c].nodes[n].mz += -cells[c].nodes[n].kr * (omegaNext - omegaPrev) * dt;

        // viscosité en rotation
        if (cells[c].nodes[n].visc > 0.0) {
          double vtheta = (omegaNext - omegaPrev);
          cells[c].nodes[n].mz += -cells[c].nodes[n].visc * vtheta;
        }

        // TODO: ajouter la possibilité de désactiver la plasticité (?)
        if (cells[c].nodes[n].mz > cells[c].nodes[n].mz_max) {
          cells[c].nodes[n].mz = cells[c].nodes[n].mz_max;
        } else if (cells[c].nodes[n].mz < -cells[c].nodes[n].mz_max) {
          cells[c].nodes[n].mz = -cells[c].nodes[n].mz_max;
        }

        // mz sur next
        double F = cells[c].nodes[n].mz * lNextSqrInv;
        vec2r finc = F * TNext;
        cells[c].nodes[next].force += finc;
        cells[c].nodes[n].force -= finc;

        // -mz sur prev
        F = cells[c].nodes[n].mz * lPrevSqrInv;
        finc = -F * TPrev;
        cells[c].nodes[prev].force += finc;
        cells[c].nodes[n].force -= finc;
      }
    }
  }

  {
    START_TIMER("globalViscosity");
    // viscosité globale
    if (globalViscosity > 0.0) {
#pragma omp parallel for schedule(static, chunk_size)
      for (size_t c = 0; c < cells.size(); c++) {
        for (size_t n = 0; n < cells[c].nodes.size(); n++) {
          cells[c].nodes[n].force -= globalViscosity * cells[c].nodes[n].vel;
        }
      }
    }
  }
}

/// Compute the node accelerations
///
///
void Lhyphen::nodeAccelerations() {
  START_TIMER("nodeAccelerations");

#ifdef _OPENMP
  int chunk_size = (int)cells.size() / nbThreads;
#endif

#pragma omp parallel for schedule(static, chunk_size)
  for (size_t c = 0; c < cells.size(); ++c) {
    for (size_t n = 0; n < cells[c].nodes.size(); ++n) {
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

  if (cellContent == CELL_CONSTANT_PV) {
    InternalGasPressureForce();
  } else if (cellContent == CELL_ELASTIC_PV) {
    InternalLiquidPressureForce();
  }

#pragma omp parallel for schedule(static, chunk_size)
  for (size_t c = 0; c < cells.size(); ++c) {
    for (size_t n = 0; n < cells[c].nodes.size(); ++n) {
      cells[c].nodes[n].acc += cells[c].nodes[n].force / cells[c].nodes[n].mass;
    }
  }
}

/// One step (velocity verlet scheme)
///
///
void Lhyphen::SingleStep() {
  START_TIMER("SingleStep");

#ifdef _OPENMP
  int chunk_size = (int)cells.size() / nbThreads;
#endif

#pragma omp parallel for schedule(static, chunk_size)
  for (size_t c = 0; c < cells.size(); ++c) {
    for (size_t n = 0; n < cells[c].nodes.size(); ++n) {
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
          cells[c].nodes[n].pos.x += cells[c].nodes[n].vel.x * dt;
        }

        if (ctrl->ymode == FORCE_CONTROL) {
          cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y * dt + cells[c].nodes[n].acc.y * dt2_2;
          cells[c].nodes[n].vel.y += cells[c].nodes[n].acc.y * dt_2;
        } else if (ctrl->ymode == VELOCITY_CONTROL) {
          double cycle = 1.0;
          if (cyclicVelPeriod > 0.0) {
            double d = t / cyclicVelPeriod;
            double r = d - floor(d);
            if (r >= 0.5)
              cycle = -1.0;
          }
          cells[c].nodes[n].vel.y = cycle * ctrl->yvalue;
          cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y * dt;
        }
      }
    }
  }

  nodeAccelerations();

#pragma omp parallel for schedule(static, chunk_size)
  for (size_t c = 0; c < cells.size(); ++c) {
    for (size_t n = 0; n < cells[c].nodes.size(); ++n) {
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
#pragma omp parallel for schedule(static, chunk_size)
    for (size_t c = 0; c < cells.size(); ++c) {
      for (size_t n = 0; n < cells[c].nodes.size(); ++n) {
        cells[c].nodes[n].vel *= (1.0 - numericalDissipation);
      }
    }
  }
}

///  Run the simulation!
///
void Lhyphen::integrate() {
  START_TIMER("integrate");

  std::vector<std::ofstream> cellFiles(followedCells.size());
  for (size_t c = 0; c < followedCells.size(); c++) {
    char name[256];
    snprintf(name, 256, "cell%ld.txt", followedCells[c]);
    cellFiles[c].open(name);
  }

  std::vector<std::ofstream> of(capturedNodes.size());
  for (size_t c = 0; c < capturedNodes.size(); c++) {
    of[c].open(capturedNodes[c].filename.c_str());
  }

  breakHistory.open("breakHistory.txt");
  breakHistory << "#ci cj xdeb ydeb xend yend released_NRJ" << std::endl;

  breakEvol.open("breakEvol.txt");
  breakEvol << "#time cumulated_breakage_NRJ cumulated_breaked_length" << std::endl;

  // Première construction de la liste de voisins, puis diagnostic avant de démarrer
  updateNeighbors();
  diagnostics();

  // === START THE LOOP ===
  for (int step = 0; step < nstep; step++) {

    if (step % nstepPeriodVerlet == 0) {
      updateNeighbors();
    }

    if (step % nstepPeriodRecord == 0) {
      for (size_t c = 0; c < followedCells.size(); c++) {
        size_t cid = followedCells[c];
        cells[cid].CellCenter();
        vec2r force;
        cells[cid].CellForce(force);
        cellFiles[c] << cells[cid].center << " " << force << '\n';
      }

      breakEvol << t << ' ' << cumulatedG << ' ' << cumulatedL << std::endl;
    }

    if (nstepPeriodSVG > 0 && step % nstepPeriodSVG == 0) {

      for (size_t c = 0; c < capturedNodes.size(); c++) {
        vec2r meanPos, sumForce, meanForce;
        for (size_t cn = 0; cn < capturedNodes[c].cellNodeIDs.size(); cn++) {
          size_t C = capturedNodes[c].cellNodeIDs[cn].c;
          size_t N = capturedNodes[c].cellNodeIDs[cn].n;
          meanPos += cells[C].nodes[N].pos;
          sumForce += cells[C].nodes[N].force;
        }
        if (capturedNodes[c].cellNodeIDs.size() > 0) {
          meanPos /= (double)capturedNodes[c].cellNodeIDs.size();
          meanForce = sumForce / (double)capturedNodes[c].cellNodeIDs.size();
        }
        of[c] << meanPos << ' ' << sumForce << ' ' << meanForce << '\n' << std::flush;
      }

      saveSVG(isvg);
      isvg++;
    }

    if (nstepPeriodConf > 0 && step % nstepPeriodConf == 0) {
      saveCONF(iconf);
      iconf++;
    }

    SingleStep();
    t += dt;
  }
}

// ======================================================================================================
// Volume internal forces
// ======================================================================================================

// loi des gaz parfaits à température constante
void Lhyphen::InternalGasPressureForce() {
  START_TIMER("InternalGasPressureForce");

  for (size_t c = 0; c < cells.size(); c++) {
    if (cells[c].close == true) {
      double PrevSurface = cells[c].surface; // plus utile
      double PrevPint = cells[c].p_int;
      cells[c].CellSurface();
      // P(k-1) S(k-1) = P(k)S(k) -> P(k) = P(k-1) S(k-1)/S(k)
      cells[c].p_int = PrevSurface * PrevPint / cells[c].surface;

      for (size_t b = 0; b < cells[c].bars.size(); b++) {
        // FIXME : pour l'instant  on fait l'hypothèse que les noeuds sont ordonnées (sens trigo ?)
        size_t i = cells[c].bars[b].i;
        size_t j = cells[c].bars[b].j;
        vec2r barVector = cells[c].nodes[i].pos - cells[c].nodes[j].pos;
        vec2r barDir = barVector;
        double barL = barDir.normalize();
        vec2r barDirRot(-barDir.y, barDir.x);
        vec2r Fp_bar = cells[c].p_int * barL * barDirRot;
        cells[c].nodes[i].force += Fp_bar * 0.5;
        cells[c].nodes[j].force += Fp_bar * 0.5;
      }
    } else {
      continue;
    }
  }
}

// remarque : COMPRESSIBLE plutot que liquide
void Lhyphen::InternalLiquidPressureForce() {
  START_TIMER("InternalLiquidPressureForce");

  for (size_t c = 0; c < cells.size(); c++) {
    if (cells[c].close == true) {

      cells[c].CellSurface();
      cells[c].p_int = -compressFactor * (cells[c].surface - cells[c].surface0) / cells[c].surface0;

      for (size_t b = 0; b < cells[c].bars.size(); b++) {
        // FIXME : pour l'instant  on fait l'hypothèse que les noeuds sont ordonnées (sens trigo ?)
        size_t i = cells[c].bars[b].i;
        size_t j = cells[c].bars[b].j;
        vec2r barVector = cells[c].nodes[i].pos - cells[c].nodes[j].pos;
        vec2r barDir = barVector;
        double barL = barDir.normalize();
        vec2r barDirRot(-barDir.y, barDir.x);
        vec2r Fp_bar = cells[c].p_int * barL * barDirRot;
        cells[c].nodes[i].force += Fp_bar * 0.5;
        cells[c].nodes[j].force += Fp_bar * 0.5;
      }
    } else {
      continue;
    }
  }
}

// ======================================================================================================
// Dumping (conf-files)
// ======================================================================================================

///  Save the current configuration in a file
///
///  @param fname  name of the file
///
void Lhyphen::saveCONF(const char *fname) {
  START_TIMER("saveCONF");

  std::ofstream file(fname);

  file << "kn " << kn << '\n';
  file << "kt " << kt << '\n';
  file << "adaptativeStiffness " << adaptativeStiffness << '\n';
  file << "mu " << mu << '\n';
  file << "gravity " << gravity << '\n';
  file << "numericalDissipation " << numericalDissipation << '\n';
  file << "globalViscosity " << globalViscosity << '\n';
  file << "distVerlet " << distVerlet << '\n';
  file << "t " << t << '\n';
  file << "dt " << dt << '\n';

  file << "nstep " << nstep << '\n';
  file << "nstepPeriodVerlet " << nstepPeriodVerlet << '\n';
  file << "nstepPeriodSVG " << nstepPeriodSVG << '\n';
  file << "nstepPeriodRecord " << nstepPeriodRecord << '\n';
  file << "nstepPeriodConf " << nstepPeriodConf << '\n';
  file << "isvg " << isvg << '\n';
  file << "iconf " << iconf << '\n';
  file << "limits " << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << '\n';

  file << "cells " << cells.size() << '\n';
  for (size_t i = 0; i < cells.size(); i++) {
    cells[i].CellSurface();
    file << cells[i].radius << ' ' << cells[i].nodes.size() << ' ' << cells[i].bars.size() << ' ' << cells[i].p_int
         << ' ' << cells[i].surface << ' ' << cells[i].surface0 << ' ' << (int)cells[i].close << '\n';
    for (size_t n = 0; n < cells[i].nodes.size(); n++) {
      file << cells[i].nodes[n].mass << ' ' << cells[i].nodes[n].pos << ' ' << cells[i].nodes[n].vel << ' '
           << cells[i].nodes[n].force << ' ';
      put_Ptr_size_t<'x'>(file, cells[i].nodes[n].ictrl);
      file << ' ';
      put_Ptr_size_t<'x'>(file, cells[i].nodes[n].prevNode);
      file << ' ';
      put_Ptr_size_t<'x'>(file, cells[i].nodes[n].nextNode);
      file << ' ';
      file << cells[i].nodes[n].kr << ' ' << cells[i].nodes[n].mz << ' ' << cells[i].nodes[n].mz_max << '\n';
    }
    for (size_t b = 0; b < cells[i].bars.size(); b++) {
      file << cells[i].bars[b].i << ' ' << cells[i].bars[b].j << ' ' << cells[i].bars[b].l0 << ' '
           << cells[i].bars[b].kn << ' ' << cells[i].bars[b].fn << '\n';
    }
  }

  file << "neighbors\n";
  for (size_t ci = 0; ci < cells.size(); ci++) {

    // count active neighbors
    size_t nbN = 0;
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {
      if (InterIt->glueState == 0 && InterIt->fn == 0.0) {
        continue;
      }
      nbN++;
    }

    file << nbN << '\n';
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {
      if (InterIt->glueState == 0 && InterIt->fn == 0.0) {
        continue;
      }

      file << InterIt->ic << ' ' << InterIt->jc << ' ' << InterIt->in << ' ' << InterIt->jn << ' ' << InterIt->n << ' '
           << InterIt->contactState << ' ' << InterIt->fn << ' ' << InterIt->ft << ' ' << InterIt->glueState << ' ';

      if (InterIt->glueState == 1) {
        file << InterIt->fn_coh << ' ' << InterIt->ft_coh << ' ' << InterIt->kn_coh << ' ' << InterIt->kt_coh << ' '
             << InterIt->fn_coh_max << ' ' << InterIt->ft_coh_max << ' ' << InterIt->yieldPower;
      } else if (InterIt->glueState == 2) {
        file << InterIt->fn_coh << ' ' << InterIt->ft_coh << ' ' << InterIt->kn_coh << ' ' << InterIt->kt_coh << ' '
             << InterIt->Gc;
      }
      file << '\n';
    }
  }

  if (!controlBoxAreas.empty()) {
    file << "controlBoxAreas " << controlBoxAreas.size() << std::endl;
    for (size_t i = 0; i < controlBoxAreas.size(); i++) {
      file << controlBoxAreas[i].xmin << ' ' << controlBoxAreas[i].xmax << ' ' << controlBoxAreas[i].ymin << ' '
           << controlBoxAreas[i].ymax << ' ' << controlBoxAreas[i].xmode << ' ' << controlBoxAreas[i].xvalue << ' '
           << controlBoxAreas[i].ymode << ' ' << controlBoxAreas[i].yvalue << std::endl;
    }
  }
}

///  Saves the CONF file.
///
///  @param ifile  the file number to save
///
void Lhyphen::saveCONF(int ifile) {
  char fname[256];
  snprintf(fname, 256, "conf%d", ifile);
  std::cout << "save CONF file: " << fname << '\n';
  saveCONF(fname);
}

///  Load a configuration that has been saved in a file
///
///  @param fname The name of the file
///
void Lhyphen::loadCONF(const char *fname) {
  static bool firstLoad = true;

  std::ifstream file(fname);

  if (exprParser == nullptr) {
    exprParser = new ExpressionParser();
  }

  int modelGc = 0;
  controlBoxAreas.clear();

  std::string token;
  file >> token;
  while (file.good()) {
    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(file, token); // ignore the rest of the current line
      file >> token;        // next token
      continue;
    } else if (token == "nbThreads") {
#ifdef _OPENMP
      nbThreads = 1;
      file >> nbThreads;
      omp_set_num_threads(nbThreads);
      std::cout << "OpenMP acceleration (Number of threads = " << nbThreads << ")\n";
#else
      std::cout << "No multithreading\n";
      file >> nbThreads;
      nbThreads = 1;
#endif
    } else if (token == "gravity") {
      exprParser->getValue(file, gravity.x);
      exprParser->getValue(file, gravity.y);
      if (firstLoad) {
        std::cout << "> gravity = " << gravity << std::endl;
      }
    } else if (token == "numericalDissipation") {
      exprParser->getValue(file, numericalDissipation);
      if (firstLoad)
        std::cout << "> numericalDissipation = " << numericalDissipation << std::endl;
    } else if (token == "globalViscosity") {
      exprParser->getValue(file, globalViscosity);
      if (firstLoad) {
        std::cout << "> globalViscosity = " << globalViscosity << std::endl;
      }
    } else if (token == "linkCells") {
      exprParser->getValue(file, linkCells_lx);
      exprParser->getValue(file, linkCells_ly);
      if (firstLoad) {
        std::cout << "* Update neighbors with linkCells: lx = " << linkCells_lx << ", ly = " << linkCells_ly
                  << std::endl;
      }
      updateNeighbors = [this]() { this->updateNeighbors_linkCells(); };
    } else if (token == "distVerlet") {
      exprParser->getValue(file, distVerlet);
      if (firstLoad) {
        std::cout << "> distVerlet = " << distVerlet << std::endl;
      }
    } else if (token == "t") {
      exprParser->getValue(file, t);
      if (firstLoad) {
        std::cout << "> t = " << t << std::endl;
      }
    } else if (token == "cyclicVelPeriod") {
      exprParser->getValue(file, cyclicVelPeriod);
      if (firstLoad) {
        std::cout << "> cyclicVelPeriod = " << cyclicVelPeriod << std::endl;
      }
    } else if (token == "dt") {
      double timestep;
      exprParser->getValue(file, timestep);
      if (firstLoad) {
        std::cout << "> dt = " << timestep << std::endl;
      }
      setTimeStep(timestep);
    } else if (token == "nstep") {
      exprParser->getValue(file, nstep);
      if (firstLoad) {
        std::cout << "> nstep = " << nstep << std::endl;
      }
    } else if (token == "nstepPeriodVerlet") {
      exprParser->getValue(file, nstepPeriodVerlet);
      if (firstLoad) {
        std::cout << "> nstepPeriodVerlet = " << nstepPeriodVerlet << std::endl;
      }
    } else if (token == "nstepPeriodSVG") {
      exprParser->getValue(file, nstepPeriodSVG);
      if (firstLoad) {
        std::cout << "> nstepPeriodSVG = " << nstepPeriodSVG << std::endl;
      }
    } else if (token == "nstepPeriodRecord") {
      exprParser->getValue(file, nstepPeriodRecord);
      if (firstLoad) {
        std::cout << "> nstepPeriodRecord = " << nstepPeriodRecord << std::endl;
      }
    } else if (token == "nstepPeriodConf") {
      exprParser->getValue(file, nstepPeriodConf);
      if (firstLoad) {
        std::cout << "> nstepPeriodConf = " << nstepPeriodConf << std::endl;
      }
    } else if (token == "isvg") {
      exprParser->getValue(file, isvg);
      if (firstLoad) {
        std::cout << "> isvg = " << isvg << std::endl;
      }
    } else if (token == "iconf") {
      exprParser->getValue(file, iconf);
      if (firstLoad) {
        std::cout << "> iconf = " << iconf << std::endl;
      }
    } else if (token == "kn") {
      exprParser->getValue(file, kn);
      if (firstLoad) {
        std::cout << "> kn = " << kn << std::endl;
      }
    } else if (token == "kt") {
      exprParser->getValue(file, kt);
      if (firstLoad) {
        std::cout << "> kt = " << kt << std::endl;
      }
    } else if (token == "viscnrate") {
      exprParser->getValue(file, viscnrate);
      if (firstLoad) {
        std::cout << "> viscnrate = " << viscnrate << std::endl;
      }
    } else if (token == "define") {
      std::string name;
      file >> name;
      double value = 0.0;
      exprParser->getValue(file, value);
      exprParser->addConstant(name, value);
      if (firstLoad)
        std::cout << "* Parameter defined as a constant: " << name << " " << value << std::endl;
    } else if (token == "adaptativeStiffness") {
      exprParser->getValue(file, adaptativeStiffness);
      if (firstLoad)
        std::cout << "> adaptativeStiffness = " << adaptativeStiffness << std::endl;
    } else if (token == "compressFactor") {
      exprParser->getValue(file, compressFactor);
      if (firstLoad)
        std::cout << "> compressFactor = " << compressFactor << std::endl;
    } else if (token == "setGlueSameProperties") {
      double kn_coh;
      double kt_coh;
      double fn_coh_max;
      double ft_coh_max;
      double yieldPower;
      exprParser->getValue(file, kn_coh);
      exprParser->getValue(file, kt_coh);
      exprParser->getValue(file, fn_coh_max);
      exprParser->getValue(file, ft_coh_max);
      exprParser->getValue(file, yieldPower);
      if (firstLoad) {
        std::cout << "* setGlueSameProperties: " << std::endl;
        std::cout << "    |     kn_coh = " << kn_coh << std::endl;
        std::cout << "    |     kt_coh = " << kt_coh << std::endl;
        std::cout << "    | fn_coh_max = " << fn_coh_max << std::endl;
        std::cout << "    | ft_coh_max = " << ft_coh_max << std::endl;
        std::cout << "    | yieldPower = " << yieldPower << std::endl;
      }
      setGlueSameProperties(kn_coh, kt_coh, fn_coh_max, ft_coh_max, yieldPower);
    } else if (token == "setGcGlueSameProperties") {
      double kn_coh;
      double kt_coh;
      double Gc;
      exprParser->getValue(file, kn_coh);
      exprParser->getValue(file, kt_coh);
      exprParser->getValue(file, Gc);
      if (firstLoad) {
        std::cout << "* setGcGlueSameProperties: " << std::endl;
        std::cout << "    | kn_coh = " << kn_coh << std::endl;
        std::cout << "    | kt_coh = " << kt_coh << std::endl;
        std::cout << "    |     Gc = " << Gc << std::endl;
      }
      setGcGlueSameProperties(kn_coh, kt_coh, Gc);
    } else if (token == "mu") {
      // file >> mu;
      exprParser->getValue(file, mu);
      if (firstLoad)
        std::cout << "> mu = " << mu << std::endl;
    } else if (token == "fadh") {
      exprParser->getValue(file, fadh);
      if (firstLoad)
        std::cout << "> fadh = " << fadh << std::endl;
    } else if (token == "limits") {
      file >> xmin >> xmax >> ymin >> ymax;
    } else if (token == "cells") {
      cells.clear();
      size_t nbCells;
      file >> nbCells;
      for (size_t i = 0; i < nbCells; i++) {
        Cell C;
        size_t nbNodes, nbBars;
        file >> C.radius >> nbNodes >> nbBars >> C.p_int >> C.surface >> C.surface0 >> C.close;
        Node N(0.0, 0.0);
        for (size_t n = 0; n < nbNodes; n++) {
          file >> N.mass >> N.pos >> N.vel >> N.force;
          N.ictrl = get_Ptr_size_t<'x'>(file);
          N.prevNode = get_Ptr_size_t<'x'>(file);
          N.nextNode = get_Ptr_size_t<'x'>(file);
          file >> N.kr >> N.mz >> N.mz_max;
          C.nodes.push_back(N);
        }
        Bar B(0, 0);
        for (size_t b = 0; b < nbBars; b++) {
          file >> B.i >> B.j >> B.l0 >> B.kn >> B.fn;

          C.bars.push_back(B);
        }
        cells.push_back(C);
      }
    } else if (token == "neighbors") {

      for (size_t ci = 0; ci < cells.size(); ci++) {
        cells[ci].neighbors.clear();
        size_t nbNeighbors;
        file >> nbNeighbors;
        for (size_t k = 0; k < nbNeighbors; k++) {
          Neighbor K(0, 0, 0, 0);
          file >> K.ic >> K.jc >> K.in >> K.jn >> K.n >> K.contactState >> K.fn >> K.ft >> K.glueState;
          if (K.glueState == 1) {
            modelGc = 1;
            file >> K.fn_coh >> K.ft_coh >> K.kn_coh >> K.kt_coh >> K.fn_coh_max >> K.ft_coh_max >> K.yieldPower;
          } else if (K.glueState == 2) {
            modelGc = 2;
            file >> K.fn_coh >> K.ft_coh >> K.kn_coh >> K.kt_coh >> K.Gc;
          }
          cells[ci].neighbors.insert(K);
        }
      }
    } else if (token == "addMultiLine") {
      named_arg_TwoPointsDataset h;
      named_arg_CellProperties p;
      int n;
      file >> h.xo >> h.yo >> h.xe >> h.ye >> h.barWidth >> n;
      file >> p.kn >> p.kr >> p.mz_max;
      // TODO: add exprParser ??
      p.p_int = 0.0;
      addMultiLine(h, p, n);
    } else if (token == "addRegularPolygonalCell") {
      named_arg_RegularCellDataset h;
      named_arg_CellProperties p;
      file >> h.nbFaces >> h.x >> h.y >> h.rot >> h.Rext >> h.barWidth;
      file >> p.kn >> p.kr >> p.mz_max;
      // TODO: add exprParser ??
      addRegularPolygonalCell(h, p);
    } else if (token == "addSquareBrickWallCells") {
      named_arg_CellProperties p;
      int nx, ny;
      double hdist;
      double xleft, ybottom;
      double barWidth;
      file >> nx >> ny >> hdist >> xleft >> ybottom >> barWidth;
      file >> p.kn >> p.kr >> p.mz_max;
      // TODO: add exprParser ??
      addSquareBrickWallCells(nx, ny, hdist, xleft, ybottom, barWidth, p);
    } else if (token == "findDisplayArea") {
      double d;
      file >> d;
      findDisplayArea(d);
    } else if (token == "glue") {
      double dist;
      file >> dist;
      modelGc = 1;
      glue(dist, modelGc);
    } else if (token == "GcGlue" || token == "distGcGlue") {
      double dist;
      file >> dist;
      modelGc = 2;
      glue(dist, modelGc);
    } else if (token == "setCellWallDampingRates") {
      double alpha_s, alpha_b;
      exprParser->getValue(file, alpha_s);
      exprParser->getValue(file, alpha_b);
      if (firstLoad) {
        std::cout << "* setCellWallDampingRates: " << std::endl;
        std::cout << "    | alpha_s = " << alpha_s << std::endl;
        std::cout << "    | alpha_b = " << alpha_b << std::endl;
      }
      setCellWallDampingRates(alpha_s, alpha_b);
    } else if (token == "setCellWallDampings") {
      double nu_s, nu_b;
      exprParser->getValue(file, nu_s);
      exprParser->getValue(file, nu_b);
      if (firstLoad) {
        std::cout << "* setCellWallDampingRates: " << std::endl;
        std::cout << "    | nu_s = " << nu_s << std::endl;
        std::cout << "    | nu_b = " << nu_b << std::endl;
      }
      setCellWallDampings(nu_s, nu_b);
    } else if (token == "setCellMasses") {
      double mass;
      exprParser->getValue(file, mass);
      if (firstLoad) {
        std::cout << "* setCellMasses: " << std::endl;
        std::cout << "    | mass = " << mass << std::endl;
      }
      setCellMasses(mass);
    } else if (token == "setNodeMasses") {
      double mass;
      exprParser->getValue(file, mass);
      if (firstLoad) {
        std::cout << "* setNodeMasses: " << std::endl;
        std::cout << "    | mass = " << mass << std::endl;
      }
      setNodeMasses(mass);
    } else if (token == "setCellWallDensities") {
      double rho, thickness;
      exprParser->getValue(file, rho);
      exprParser->getValue(file, thickness);
      if (firstLoad) {
        std::cout << "* setCellWallDensities: " << std::endl;
        std::cout << "    |       rho = " << rho << std::endl;
        std::cout << "    | thickness = " << thickness << std::endl;
      }
      setCellWallDensities(rho, thickness);
    } else if (token == "setCellDensities") {
      double rho, thickness;
      exprParser->getValue(file, rho);
      exprParser->getValue(file, thickness);
      if (firstLoad) {
        std::cout << "* setCellDensities: " << std::endl;
        std::cout << "    |       rho = " << rho << std::endl;
        std::cout << "    | thickness = " << thickness << std::endl;
      }
      setCellDensities(rho, thickness);
    } else if (token == "setNodeControl") {
      double xvalue, yvalue;
      int xmode, ymode;
      size_t c, n;
      file >> c >> n >> xmode >> xvalue >> ymode >> yvalue;
      // TODO: add exprParser ??
      setNodeControl(c, n, xmode, xvalue, ymode, yvalue);
    } else if (token == "setNodeControlInBox") {
      double xmin_, xmax_, ymin_, ymax_, xvalue, yvalue;
      int xmode, ymode;
      file >> xmin_ >> xmax_ >> ymin_ >> ymax_ >> xmode >> xvalue >> ymode >> yvalue;
      // TODO: add exprParser ??
      setNodeControlInBox(xmin_, xmax_, ymin_, ymax_, xmode, xvalue, ymode, yvalue);
      controlBoxAreas.push_back(ControlBoxArea(xmin_, xmax_, ymin_, ymax_, xmode, xvalue, ymode, yvalue));
    } else if (token == "controlBoxAreas") {
      double xmin_, xmax_, ymin_, ymax_, xvalue, yvalue;
      int xmode, ymode;
      size_t nb = 0;
      file >> nb;
      for (size_t i = 0; i < nb; i++) {
        file >> xmin_ >> xmax_ >> ymin_ >> ymax_ >> xmode >> xvalue >> ymode >> yvalue;
        controlBoxAreas.push_back(ControlBoxArea(xmin_, xmax_, ymin_, ymax_, xmode, xvalue, ymode, yvalue));
      }
    } else if (token == "captureNodes") {
      double xmin_, xmax_, ymin_, ymax_;
      std::string filename;
      file >> filename >> xmin_ >> xmax_ >> ymin_ >> ymax_;
      CapturedNodes CN(filename.c_str(), xmin_, xmax_, ymin_, ymax_);
      for (size_t c = 0; c < cells.size(); c++) {
        for (size_t n = 0; n < cells[c].nodes.size(); n++) {

          if (cells[c].nodes[n].pos.x >= xmin_ && cells[c].nodes[n].pos.x <= xmax_ &&
              cells[c].nodes[n].pos.y >= ymin_ && cells[c].nodes[n].pos.y <= ymax_) {
            CellNodeID CNID(c, n);
            CN.cellNodeIDs.push_back(CNID);
          }
        }
      }
      capturedNodes.push_back(CN);
    } else if (token == "followCell") {
      size_t cid;
      exprParser->getValue(file, cid);
      if (firstLoad) {
        std::cout << "* followCell cell.id = " << cid << std::endl;
      }
      followedCells.push_back(cid);
    } else if (token == "setCellControl") {
      size_t icell;
      int xmode, ymode;
      double xvalue, yvalue;
      file >> icell >> xmode >> xvalue >> ymode >> yvalue;
      // TODO: add exprParser ??
      setCellControl(icell, xmode, xvalue, ymode, yvalue);
    } else if (token == "reorder") {
      file >> reorder;
    } else if (token == "readNodeFile") {
      std::string fileName;
      double barWidth, Kn, Kr, Mz_max, p_int;
      file >> fileName;
      exprParser->getValue(file, barWidth);
      exprParser->getValue(file, Kn);
      exprParser->getValue(file, Kr);
      exprParser->getValue(file, Mz_max);
      exprParser->getValue(file, p_int);
      if (firstLoad) {
        std::cout << "* readNodeFile " << fileName << ":" << std::endl;
        std::cout << "    | barWidth = " << barWidth << std::endl;
        std::cout << "    |       Kn = " << Kn << std::endl;
        std::cout << "    |       Kr = " << Kr << std::endl;
        std::cout << "    |   Mz_max = " << Mz_max << std::endl;
        std::cout << "    |    p_int = " << p_int << std::endl;
      }
      readNodeFile(fileName.c_str(), barWidth, Kn, Kr, Mz_max, p_int);
    } else if (token == "setCellInternalPressure") {
      double p_int;
      size_t c;
      file >> c >> p_int;
      // TODO: add exprParser ??
      setCellInternalPressure(c, p_int);
    } else if (token == "enablePressures") {
      // enablePressures = 1;
      std::cout << "!!!! enablePressures has been replaced by cellContent" << std::endl;
    } else if (token == "disablePressures") {
      // enablePressures = 0;
      std::cout << "!!!! disablePressures has been replaced by cellContent" << std::endl;
    } else if (token == "cellContent") {
      file >> cellContent;
      if (firstLoad) {
        if (cellContent == CELL_EMPTY) {
          std::cout << "> cellContent = CELL_EMPTY" << std::endl;
        } else if (cellContent == CELL_CONSTANT_PV) {
          std::cout << "> cellContent = CELL_CONSTANT_PV" << std::endl;
        } else if (cellContent == CELL_ELASTIC_PV) {
          std::cout << "> cellContent = CELL_ELASTIC_PV" << std::endl;
        }
      }

    } else if (token == "setCellAsOpen") {
      size_t c;
      file >> c;
      cells[c].close = false;
    } else if (token == "setClose") {
      size_t c;
      file >> c;
      cells[c].close = true;
    } else {
      if (firstLoad)
        std::cout << "@Lhyphen::loadCONF, this token is not known: " << token << std::endl;
    }

    file >> token;
  }

  computeCellInitialSurfaces(); // surface des cellules
  copy_neighbors_set_to_vec();

  if (modelGc > 0) {
    associateGlue(modelGc);
  }

  firstLoad = false;
}

void Lhyphen::loadCONF(int ifile) {
  char fname[256];
  snprintf(fname, 256, "conf%d", ifile);
  loadCONF(fname);
}

void Lhyphen::computeCellInitialSurfaces() {
  for (size_t c = 0; c < cells.size(); c++) {
    if (cells[c].close == true && cells[c].surface0 == 0.0) {
      cells[c].CellSurface();
      cells[c].surface0 = cells[c].surface;
    }
  }
}

// ======================================================================================================
// Saving SVG
// ======================================================================================================

/// Find the boundaries of the area drawn in the svg files
///
/// @param factor multiplying size-factor
///
void Lhyphen::findDisplayArea(double factor) {
  if (cells.empty()) {
    std::cout << "@Lhyphen::findDisplayArea, cannot find the scene limits because no cell has been set!\n";
    return;
  }

  xmin = xmax = cells[0].nodes[0].pos.x;
  ymin = ymax = cells[0].nodes[0].pos.y;
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {

      if (cells[c].nodes[n].pos.x < xmin) {
        xmin = cells[c].nodes[n].pos.x;
      }
      if (cells[c].nodes[n].pos.x > xmax) {
        xmax = cells[c].nodes[n].pos.x;
      }
      if (cells[c].nodes[n].pos.y < ymin) {
        ymin = cells[c].nodes[n].pos.y;
      }
      if (cells[c].nodes[n].pos.y > ymax) {
        ymax = cells[c].nodes[n].pos.y;
      }
    }
  }

  std::cout << "limites des noeuds:\n";
  std::cout << "xmin = " << xmin << '\n';
  std::cout << "xmax = " << xmax << '\n';
  std::cout << "ymin = " << ymin << '\n';
  std::cout << "ymax = " << ymax << '\n';

  double halfL = 0.5 * (xmax - xmin);
  double halfH = 0.5 * (ymax - ymin);
  double x0 = 0.5 * (xmin + xmax);
  double y0 = 0.5 * (ymin + ymax);
  xmin = x0 - factor * halfL;
  xmax = x0 + factor * halfL;
  ymin = y0 - factor * halfH;
  ymax = y0 + factor * halfH;
}

///   Very basic function for drawing cells with blue lines,
///   Remarque: pour d'autre sorties graphiques, il vaut mieux le faire en post-traitement à partir des fichiers conf :)
///
///   @param num         numero du fichier
///   @param nameBase    nom de base (dans lequel sera inséré le numéro)
///   @param Canvaswidth largeur de canvas (par affichae rapide)
///
///   @see https://bestofcpp.com/repo/tomkwok-svgasm
///
void Lhyphen::saveSVG(int num, const char *nameBase, int CanvasWidth) {
  char name[256];
  snprintf(name, 256, nameBase, num);
  std::cout << "save file: " << name << '\n';
  std::ofstream ofs(name);
  SVGfile svg(ofs);

  int CanvasHeight = (int)(CanvasWidth * (ymax - ymin) / (xmax - xmin));
  svg.begin(CanvasWidth, CanvasHeight);
  viewZone vz(0, 0, CanvasWidth, CanvasHeight);
  double delta = (xmax - xmin) * 0.05; // c'est une bordure
  vz.adjustRange(xmin - delta, xmax + delta, ymin - delta, ymax + delta);

  char opt[256];

  colorRGBA col;

  switch (SVG_colorCells) {
  case 1: { // pression
    double pmin = 1e20;
    double pmax = -1e20;
    for (size_t c = 0; c < cells.size(); c++) {
      if (cells[c].p_int > pmax) {
        pmax = cells[c].p_int;
      }
      if (cells[c].p_int < pmin) {
        pmin = cells[c].p_int;
      }
    }

    double pp = std::max(fabs(pmax), fabs(pmin));
    ctNeg.setMinMax(-(float)pp, 0.0);
    ctPos.setMinMax(0.0, (float)pp);
  } break;

  case 2: { // NRJ elast

    double NRJmax = 0.0;
    for (size_t c = 0; c < cells.size(); c++) {
      double NRJ = cells[c].getElasticNRJ(compressFactor);
      if (NRJ > NRJmax) {
        NRJmax = NRJ;
      }
    }
    ctPos.setMinMax(0.0, (float)NRJmax);
  } break;
  }

  // draw the cells
  for (size_t c = 0; c < cells.size(); c++) {

    if (SVG_colorCells > 0 && cells[c].close == true) {

      switch (SVG_colorCells) {
      case 1: {
        if (cells[c].p_int >= 0.0) {
          ctPos.getRGB((float)cells[c].p_int, &col);
        } else {
          ctNeg.getRGB((float)cells[c].p_int, &col);
        }
      } break;

      case 2: {
        double NRJ = cells[c].getElasticNRJ(compressFactor);
        ctPos.getRGB((float)NRJ, &col);
      } break;
      }

      snprintf(opt, 256, "stroke:none;fill:rgb(%d,%d,%d)", col.r, col.g, col.b);

      std::vector<vec2r> contour;
      for (size_t s = 0; s < cells[c].nodes.size(); s++) {
        contour.push_back(cells[c].nodes[s].pos);
      }
      svg.polygon(vz, contour, opt);
    }

    // contour de la cellule
    for (size_t b = 0; b < cells[c].bars.size(); b++) {
      snprintf(opt, 256, "stroke:blue;fill:none;stroke-linecap:round;stroke-width:%g", 2 * cells[c].radius * vz.scalex);
      size_t i = cells[c].bars[b].i;
      size_t j = cells[c].bars[b].j;
      svg.line(vz, cells[c].nodes[i].pos.x, cells[c].nodes[i].pos.y, cells[c].nodes[j].pos.x, cells[c].nodes[j].pos.y,
               opt);
    }
  }

  // draw the forces
  if (SVG_cellForces > 0) {
    std::map<std::pair<size_t, size_t>, vec2r> force_map;

    std::pair<size_t, size_t> duo;
    double surfMin = 1.0e20;
    for (size_t c = 0; c < cells.size(); c++) {
      cells[c].CellCenter();
      if (cells[c].close == true && cells[c].surface0 < surfMin)
        surfMin = cells[c].surface0;
      for (auto &in : cells[c].neighbors) {
        duo.first = in.ic;
        duo.second = in.jc;
        if (cells[in.ic].close == false || cells[in.jc].close == false)
          continue;

        if (in.contactState == TOUCHING) {
          auto it = force_map.find(duo);
          if (it != force_map.end())
            it->second += in.fn * in.n;
          else
            force_map[duo] = in.fn * in.n;
        }
      }
    }

    double fnMax = 0.0;
    for (auto &fm : force_map) {
      double fn = norm(fm.second);
      if (fn > fnMax)
        fnMax = fn;
    }
    if (fabs(fnMax) > 1e-12) {
      double req = sqrt(surfMin / M_PI);

      for (auto &fm : force_map) {
        size_t ic = fm.first.first;
        size_t jc = fm.first.second;
        double width = norm(fm.second) / fnMax;
        snprintf(opt, 256, "stroke:black;fill:none;stroke-linecap:round;stroke-width:%g", req * width * vz.scalex);
        svg.line(vz, cells[ic].center.x, cells[ic].center.y, cells[jc].center.x, cells[jc].center.y, opt);
      }
    }
  }

  svg.end();
}
