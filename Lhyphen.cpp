#include "Lhyphen.hpp"
#include "null_size_t.hpp"

/**
 * @brief Construct a new Lhyphen object
 *
 */
Lhyphen::Lhyphen()
    : xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), dt(0.0), globalViscosity(0.0), numericalDissipation(0.0), gravity(),
      distVerlet(0.0), cellContent(CELL_CONTAIN_NOTHING), compressFactor(0.0), kn(1000.0), kt(1000.0), mu(0.0),
      fadh(0.0), nstep(1000), nstepPeriodVerlet(1), nstepPeriodSVG(10), nstepPeriodRecord(1), nstepPeriodConf(10),
      isvg(0), iconf(0) {

  SVG_colorCells = 2;
  SVG_colorTableRescale = 0;
  SVG_colorTableMin = -10000.0;
  SVG_colorTableMax = 10000.0;

  SVG_cellForces = 1;

  ctNeg.setSize(128);
  ctPos.setSize(128);
  ctNeg.rebuild_interp_rgba({0, 127}, {{204, 204, 230, 255}, {0, 0, 255, 255}});
  ctPos.rebuild_interp_rgba({0, 127}, {{204, 204, 230, 255}, {255, 0, 0, 255}});
}

/**
 * @brief affiche un petit entete sympatique
 *
 */
void Lhyphen::head() {
  std::cout << "    L-HYPHEN\n";
  std::cout << "    _    _\n";
  std::cout << "   / \\__/ \\\n";
  std::cout << "  |  o  <  |\n";
  std::cout << "   \\_\\ _/_/\n";
  std::cout << "     /_/\n";
  std::cout << "    /_/\n\n";
}

/**
 * @brief ajoute une cellule de forme polygonale régulière
 *
 * @param h les données pour définir un polyèdre régulier
 * @param p les propriétés mécaniques
 */
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
  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, true);
  // C.print();
  cells.push_back(C);
}

/**
 * @brief ajoute une boite à partir de deux points
 *
 * @param h les deux points qui definissent les dimensions extérieures de la boite
 * @param p les propriétés mécaniques
 */
void Lhyphen::addBoxCell(named_arg_TwoPointsDataset h, named_arg_CellProperties p) {
  Cell C;
  C.p_int = p.p_int;
  C.close = true;
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
void Lhyphen::addLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p) {
  Cell C;
  C.nodes.emplace_back(h.xo, h.yo);
  C.nodes.emplace_back(h.xe, h.ye);
  C.connectOrderedNodes(h.barWidth, p.kn, p.kr, p.mz_max, false);
  cells.push_back(C);
}

/**
 * @brief ajoute une ligne formée de plusieurs barres (cellule non fermée) entre deux points
 *
 * @param h les deux points
 * @param p les propriétés mécaniques
 * @param n nombre de barres
 */
void Lhyphen::addMultiLine(named_arg_TwoPointsDataset h, named_arg_CellProperties p, int n) {
  double dx = (h.xe - h.xo) / (double)n;
  double dy = (h.ye - h.yo) / (double)n;

  Cell C;
  C.nodes.emplace_back(h.xo, h.yo);

  for (int i = 1; i <= n; i++) {
    C.nodes.emplace_back(h.xo + i * dx, h.yo + i * dy);
  }
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
void Lhyphen::addSquareBrickWallCells(int nx, int ny, double horizontalDistance, double xleft, double ybottom,
                                      double barWidth, named_arg_CellProperties p) {
  addRegularPolygonalCellsOnTriangularGrid(nx, ny, horizontalDistance, horizontalDistance,
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
void Lhyphen::addHoneycombCells(int nx, int ny, double CellExternWidth, double xleft, double ybottom, double barWidth,
                                named_arg_CellProperties p) {
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
void Lhyphen::setTimeStep(double dt_) {
  dt = dt_;
  dt_2 = 0.5 * dt;
  dt2_2 = dt_2 * dt;
}

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

void Lhyphen::setCellMasses(double cellMass) {
  for (size_t c = 0; c < cells.size(); c++) {
    double nodeMass = cellMass / (double)cells[c].nodes.size();
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass = nodeMass;
    }
  }
}

void Lhyphen::setNodeMasses(double nodeMass) {
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      cells[c].nodes[n].mass = nodeMass;
    }
  }
}

void Lhyphen::setGlueSameProperties(double kn_coh, double kt_coh, double fn_coh_max, double ft_coh_max,
                                    double yieldPower) {

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
void Lhyphen::setNodeControl(size_t c, size_t n, int xmode, double xvalue, int ymode, double yvalue) {
  Control C(xmode, xvalue, ymode, yvalue);
  controls.push_back(C);
  cells[c].nodes[n].ictrl = controls.size() - 1;
}

/**
 * @brief ajoute un control à tous les noeuds de la cellule c
 *
 * @param c numéro de la cellule
 * @param xmode VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
 * @param xvalue valeur suivant x (une force ou une vitesse selon xmode)
 * @param ymode VELOCITY_CONTROL (0) ou FORCE_CONTROL (1)
 * @param yvalue valeur suivant y (une force ou une vitesse selon ymode)
 */
void Lhyphen::setCellControl(size_t c, int xmode, double xvalue, int ymode, double yvalue) {
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
void Lhyphen::setNodeControlInBox(double xmin_, double xmax_, double ymin_, double ymax_, int xmode, double xvalue,
                                  int ymode, double yvalue) {
  Control C(xmode, xvalue, ymode, yvalue);

  controls.push_back(C);
  size_t ictrl = controls.size() - 1;
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      vec2r pos = cells[c].nodes[n].pos;
      if (pos.x >= xmin_ && pos.x <= xmax_ && pos.y >= ymin_ && pos.y <= ymax_) {
        cells[c].nodes[n].ictrl = ictrl;
        std::cout << "@Lhyphen::setNodeControlInBox, node pos = " << pos << "\n";
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
 * @param epsilonEnds longueur à ignorer aux extremités de la barre
 */
void Lhyphen::addNodeToBarNeighbor(size_t ci, size_t cj, size_t in, size_t jn, double epsilonEnds) {
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
void Lhyphen::readNodeFile(const char *name, double barWidth, double Kn, double Kr, double Mz_max, double p_int) {
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

  for (size_t i = 0; i < cells.size(); i++) {
    cells[i].reorderNodes();
    // il faudrait un flag pour que cette fonction ne soit pas appelée si l'utilisateur le souhaite
    // Dans ce cas, il faudra que les noeuds du fichier soient positionnés dans le sens trigo (à vérifier)

    cells[i].connectOrderedNodes(barWidth, Kn, Kr, Mz_max, p_int, true);
  }
}

/**
 * @brief met à jour la liste des voisins
 *
 */
void Lhyphen::updateNeighbors() {
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
 * @param epsilonDist distance maximale
 */
void Lhyphen::glue(double epsilonDist) {
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
void Lhyphen::computeInteractionForces() {

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
            Inter->fn = -kn * dn;

            // ...vrel T ft (TODO)

            if (Inter->glueState == 0)
              Inter->fn -= fadh; // cette adhesion ne peut agir que si il n'y a pas de cohésion solide

            // transfert des force vers les noeuds concernés
            vec2r finc = Inter->fn * Inter->n;
            cells[ci].nodes[in].force += finc;
            cells[cj].nodes[jn].force -= finc;
          } else {
            Inter->contactState = 0;
            Inter->fn = 0.0;
            Inter->ft = 0.0;
          }

          if (Inter->glueState == 1) {
            // TODO
          }

        } else if (proj >= u_length) { // ====================== disque jnext (de fin)
                                       // ? doit on traiter ce cas ?
                                       // ! sinon il peut y avoir des forces d'interaction calculées deux fois !

          b = cells[ci].nodes[in].pos - cells[cj].nodes[jnext].pos;
          double sqrDist = norm2(b);
          double sumR = cells[ci].radius + cells[cj].radius;
          if (sqrDist - sumR * sumR < 0.0) {
            Inter->contactState = 1;
            double b_lenght = sqrt(sqrDist);
            double dn = b_lenght - sumR;
            Inter->n = b / b_lenght;
            Inter->fn = -kn * dn;

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
              Inter->fn -= fadh; // cette adhesion ne peut agir que si il n'y a pas de cohésion solide

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
void Lhyphen::computeNodeForces() {

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
      cells[c].bars[b].fn = -cells[c].bars[b].kn * dn;

      // TODO: ajouter plasticité dans les barres (?)

      vec2r finc = cells[c].bars[b].fn * n; // en cas d'allongement finc est orienté de j vers i
      cells[c].nodes[i].force -= finc;
      cells[c].nodes[j].force += finc;
    }
  }

  // Moment au niveau des noeuds (entre les barres d'une même cellule).
  // Puisque les noeuds n'ont pas de ddl en rotation, imposer un moment est fait ici
  // en imposant une force à chaque extrémité de barre
  for (size_t c = 0; c < cells.size(); c++) {
    for (size_t n = 0; n < cells[c].nodes.size(); n++) {
      size_t prev = cells[c].nodes[n].prevNode;
      size_t next = cells[c].nodes[n].nextNode;

      if (prev == null_size_t || next == null_size_t) { // si c'est une extrémité...
        continue;
      }

      double omegaNext = getRotationVelocityBar(cells[c].nodes[n].pos, cells[c].nodes[next].pos, cells[c].nodes[n].vel,
                                                cells[c].nodes[next].vel);
      double omegaPrev = getRotationVelocityBar(cells[c].nodes[n].pos, cells[c].nodes[prev].pos, cells[c].nodes[n].vel,
                                                cells[c].nodes[prev].vel);
      cells[c].nodes[n].mz += -cells[c].nodes[n].kr * (omegaNext - omegaPrev) * dt;

      // TODO: ajouter la possibilité de désactiver la plasticité
      if (cells[c].nodes[n].mz > cells[c].nodes[n].mz_max) {
        cells[c].nodes[n].mz = cells[c].nodes[n].mz_max;
      } else if (cells[c].nodes[n].mz < -cells[c].nodes[n].mz_max) {
        cells[c].nodes[n].mz = -cells[c].nodes[n].mz_max;
      }

      vec2r nextVector = cells[c].nodes[next].pos - cells[c].nodes[n].pos;
      vec2r prevVector = cells[c].nodes[n].pos - cells[c].nodes[prev].pos;

      // TODO: ici c'est optimisable...
      vec2r nextDir = nextVector;
      vec2r prevDir = prevVector;
      double nextL = nextDir.normalize();
      double prevL = prevDir.normalize();

      // mz sur next
      vec2r nextDirRot(-nextDir.y, nextDir.x);
      double F = cells[c].nodes[n].mz / nextL;
      vec2r finc = F * nextDirRot;
      cells[c].nodes[next].force += finc;
      cells[c].nodes[n].force -= finc;

      // -mz sur prev
      vec2r prevDirRot(-prevDir.y, prevDir.x);
      F = cells[c].nodes[n].mz / prevL;
      finc = F * prevDirRot;
      cells[c].nodes[prev].force += finc;
      cells[c].nodes[n].force -= finc;
    }
  }

  // viscosité globale
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
void Lhyphen::NodeAccelerations() {

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

  if (cellContent == CELL_CONTAIN_GAS)
    InternalGasPressureForce();
  else if (cellContent == CELL_CONTAIN_LIQUID)
    InternalLiquidPressureForce();

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
void Lhyphen::SingleStep() {

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
          // cells[c].nodes[n].force.x = 0.0;
          cells[c].nodes[n].pos.x += cells[c].nodes[n].vel.x * dt;
        }

        if (ctrl->ymode == FORCE_CONTROL) {
          cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y * dt + cells[c].nodes[n].acc.y * dt2_2;
          cells[c].nodes[n].vel.y += cells[c].nodes[n].acc.y * dt_2;
        } else if (ctrl->ymode == VELOCITY_CONTROL) {
          cells[c].nodes[n].vel.y = ctrl->yvalue;
          // cells[c].nodes[n].force.y = 0.0;
          cells[c].nodes[n].pos.y += cells[c].nodes[n].vel.y * dt;
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
void Lhyphen::integrate() {
  // std::ofstream of[capturedNodes.size()];
  std::vector<std::ofstream> of(capturedNodes.size());
  for (size_t c = 0; c < capturedNodes.size(); c++) {
    of[c].open(capturedNodes[c].filename.c_str());
  }

  for (int step = 0; step < nstep; step++) {
    if (step % nstepPeriodVerlet == 0) {
      updateNeighbors();
    }

    if (nstepPeriodSVG > 0 && step % nstepPeriodSVG == 0) {

      for (size_t c = 0; c < capturedNodes.size(); c++) {
        vec2r meanPos, meanForce;
        for (size_t cn = 0; cn < capturedNodes[c].cellNodeIDs.size(); cn++) {
          size_t C = capturedNodes[c].cellNodeIDs[cn].c;
          size_t N = capturedNodes[c].cellNodeIDs[cn].n;
          meanPos += cells[C].nodes[N].pos;
          meanForce += cells[C].nodes[N].force;
        }
        if (capturedNodes[c].cellNodeIDs.size() > 0) {
          meanPos /= (double)capturedNodes[c].cellNodeIDs.size();
          meanForce /= (double)capturedNodes[c].cellNodeIDs.size();
        }
        of[c] << meanPos << ' ' << meanForce << '\n' << std::flush;
      }

      saveSVG(isvg);
      isvg++;
    }

    if (nstepPeriodConf > 0 && step % nstepPeriodConf == 0) {
      saveCONF(iconf);
      iconf++;
    }

    SingleStep();
  }
}

/**
 * @brief Save the current configuration in a file
 *
 * @param fname Name of the file
 */
void Lhyphen::saveCONF(const char *fname) {
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
    file << cells[ci].neighbors.size() << '\n';
    for (std::set<Neighbor>::iterator InterIt = cells[ci].neighbors.begin(); InterIt != cells[ci].neighbors.end();
         ++InterIt) {

      file << InterIt->ic << ' ' << InterIt->jc << ' ' << InterIt->in << ' ' << InterIt->jn << ' ' << InterIt->n << ' '
           << InterIt->contactState << ' ' << InterIt->fn << ' ' << InterIt->ft << ' ' << InterIt->glueState << ' ';
      if (InterIt->glueState > 0) {
        file << InterIt->fn_coh << ' ' << InterIt->ft_coh << ' ' << InterIt->kn_coh << ' ' << InterIt->kt_coh << ' '
             << InterIt->fn_coh_max << ' ' << InterIt->ft_coh_max << ' ' << InterIt->yieldPower;
      }
      file << '\n';
    }
  }
}

void Lhyphen::saveCONF(int ifile) {
  char fname[256];
  snprintf(fname, 256, "conf%d", ifile);
  std::cout << "save CONF file: " << fname << '\n';
  saveCONF(fname);
}

/**
 * @brief Load a configuration that has been saved in a file
 *
 * @param fname The name of the file
 */
void Lhyphen::loadCONF(const char *fname) {
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
    } else if (token == "compressFactor") {
      file >> compressFactor;
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
          if (K.glueState > 0) {
            file >> K.fn_coh >> K.ft_coh >> K.kn_coh >> K.kt_coh >> K.fn_coh_max >> K.ft_coh_max >> K.yieldPower;
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
      addMultiLine(h, p, n);
    } else if (token == "addRegularPolygonalCell") {
      named_arg_RegularCellDataset h;
      named_arg_CellProperties p;
      file >> h.nbFaces >> h.x >> h.y >> h.rot >> h.Rext >> h.barWidth;
      file >> p.kn >> p.kr >> p.mz_max;
      addRegularPolygonalCell(h, p);
    } else if (token == "addSquareBrickWallCells") {
      named_arg_CellProperties p;
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
    } else if (token == "setCellMasses") {
      double mass;
      file >> mass;
      setCellMasses(mass);
    } else if (token == "setNodeMasses") {
      double mass;
      file >> mass;
      setNodeMasses(mass);
    } else if (token == "setCellWallDensities") {
      double rho, thickness;
      file >> rho >> thickness;
      setCellWallDensities(rho, thickness);
    } else if (token == "setCellDensities") {
      double rho, thickness;
      file >> rho >> thickness;
      setCellDensities(rho, thickness);
    } else if (token == "setNodeControl") {
      double xvalue, yvalue;
      int xmode, ymode;
      size_t c, n;
      file >> c >> n >> xmode >> xvalue >> ymode >> yvalue;
      setNodeControl(c, n, xmode, xvalue, ymode, yvalue);
    } else if (token == "setNodeControlInBox") {
      double xmin_, xmax_, ymin_, ymax_, xvalue, yvalue;
      int xmode, ymode;
      file >> xmin_ >> xmax_ >> ymin_ >> ymax_ >> xmode >> xvalue >> ymode >> yvalue;
      setNodeControlInBox(xmin_, xmax_, ymin_, ymax_, xmode, xvalue, ymode, yvalue);
    } else if (token == "captureNodes") {
      double xmin_, xmax_, ymin_, ymax_;
      std::string filename;
      file >> filename >> xmin_ >> xmax_ >> ymin_ >> ymax_;
      CapturedNodes CN(filename.c_str());
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
    } else if (token == "setCellControl") {
      size_t icell;
      int xmode, ymode;
      double xvalue, yvalue;
      file >> icell >> xmode >> xvalue >> ymode >> yvalue;
      setCellControl(icell, xmode, xvalue, ymode, yvalue);
    } else if (token == "readNodeFile") {
      std::string fileName;
      double barWidth, Kn, Kr, Mz_max, p_int;
      file >> fileName >> barWidth >> Kn >> Kr >> Mz_max >> p_int;
      readNodeFile(fileName.c_str(), barWidth, Kn, Kr, Mz_max, p_int);
    } else if (token == "setCellInternalPressure") {
      double p_int;
      size_t c;
      file >> c >> p_int;
      setCellInternalPressure(c, p_int);
    } else if (token == "enablePressures") {
      // enablePressures = 1;
      std::cout << "enablePressures has been replaced by cellContent\n";
    } else if (token == "disablePressures") {
      // enablePressures = 0;
      std::cout << "disablePressures has been replaced by cellContent\n";
    } else if (token == "cellContent") {
      file >> cellContent;
      if (cellContent == CELL_CONTAIN_NOTHING)
        std::cout << "cellContent = CELL_CONTAIN_NOTHING\n";
      else if (cellContent == CELL_CONTAIN_GAS)
        std::cout << "cellContent = CELL_CONTAIN_GAS\n";
      else if (cellContent == CELL_CONTAIN_LIQUID)
        std::cout << "cellContent = CELL_CONTAIN_LIQUID\n";
    } else if (token == "setCellAsOpen") {
      size_t c;
      file >> c;
      cells[c].close = false;
    } else if (token == "setClose") {
      size_t c;
      file >> c;
      cells[c].close = true;
    } else {
      std::cout << "@Lhyphen::loadCONF, this token is not known: " << token << '\n';
    }

    file >> token;
  }

  initCellInitialSurfaces();
}

void Lhyphen::loadCONF(int ifile) {
  char fname[256];
  snprintf(fname, 256, "conf%d", ifile);
  loadCONF(fname);
}

void Lhyphen::initCellInitialSurfaces() {
  for (size_t c = 0; c < cells.size(); c++) {

    if (cells[c].close == true && cells[c].surface0 == 0.0) {
      cells[c].CellSurface();
      cells[c].surface0 = cells[c].surface;
    }
  }
}

/**
 * @brief trouver les limites de la zone dessinée dans les fichiers svg
 *
 * @param factor facteur multiplicateur
 */
void Lhyphen::findDisplayArea(double factor) {
  if (cells.empty()) {
    std::cout << "@Lhyphen::findDisplayArea, cannot find the scene limits because no cell has been set!\n";
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
void Lhyphen::saveSVG(int num, const char *nameBase, int Canvaswidth) {
  char name[256];
  snprintf(name, 256, nameBase, num);
  std::cout << "save file: " << name << '\n';
  std::ofstream ofs(name);
  SVGfile svg(ofs);

  int CanvasHeight = Canvaswidth * (int)floor((ymax - ymin) / (xmax - xmin));
  svg.begin(Canvaswidth, CanvasHeight);
  viewZone vz(0, 0, Canvaswidth, CanvasHeight);
  double delta = (xmax - xmin) * 0.05; // c'est une bordure
  vz.adjustRange(xmin - delta, xmax + delta, ymin - delta, ymax + delta);

  char opt[256];

  colorRGBA col;

  switch (SVG_colorCells) {
  case 1: { // pression
    double pmin = 1e20;
    double pmax = -1e20;
    for (size_t c = 0; c < cells.size(); c++) {
      if (cells[c].p_int > pmax)
        pmax = cells[c].p_int;
      if (cells[c].p_int < pmin)
        pmin = cells[c].p_int;
    }

    double pp = std::max(fabs(pmax), fabs(pmin));
    ctNeg.setMinMax(-(float)pp, 0.0);
    ctPos.setMinMax(0.0, (float)pp);
  } break;

  case 2: { // NRJ elast

    double NRJmax = 0.0;
    for (size_t c = 0; c < cells.size(); c++) {
      double NRJ = cells[c].getElasticNRJ(compressFactor);
      // std::cout << "NRJ(" << c << ") = " << NRJ << '\n';
      if (NRJ > NRJmax)
        NRJmax = NRJ;
    }
    ctPos.setMinMax(0.0, (float)NRJmax);
    // std::cout << "NRJmax = " << NRJmax << '\n';
  } break;
  }

  // draw the cells
  for (size_t c = 0; c < cells.size(); c++) {

    if (SVG_colorCells > 0 && cells[c].close == true) {

      switch (SVG_colorCells) {
      case 1: {
        if (cells[c].p_int >= 0.0)
          ctPos.getRGB((float)cells[c].p_int, &col);
        else
          ctNeg.getRGB((float)cells[c].p_int, &col);
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

        if (in.contactState == 1) {
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
        // vec2r branch = cells[jc].center - cells[ic].center;
        double width = norm(fm.second) / fnMax;
        snprintf(opt, 256, "stroke:black;fill:none;stroke-linecap:round;stroke-width:%g", req * width * vz.scalex);
        svg.line(vz, cells[ic].center.x, cells[ic].center.y, cells[jc].center.x, cells[jc].center.y, opt);
      }
    }
  }

  svg.end();
}

/**
 * @brief Get the Rotation Velocity Bar object
 *
 * @param a Position point A
 * @param b Position point B
 * @param va Vitesse point A
 * @param vb Vitesse point B
 * @return double vitesse de rotation de solide rigide
 */
double Lhyphen::getRotationVelocityBar(vec2r &a, vec2r &b, vec2r &va, vec2r &vb) {
  vec2r u = b - a;
  double l = u.normalize();
  vec2r t(-u.y, u.x);
  return ((vb - va) * t / l);
}

// loi des gaz parfaits à température constante
void Lhyphen::InternalGasPressureForce() {
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
  for (size_t c = 0; c < cells.size(); c++) {
    if (cells[c].close == true) {

      cells[c].CellSurface();
      cells[c].p_int = -compressFactor * (cells[c].surface - cells[c].surface0);

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

void Lhyphen::setCellInternalPressure(size_t c, double p) {
  cells[c].p_int = p;
  // cells[c].close = false;
}
